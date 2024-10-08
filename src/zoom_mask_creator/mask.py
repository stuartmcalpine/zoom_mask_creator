import sys
import numpy as np
from mpi4py import MPI
from scipy import ndimage

from zoom_mask_creator.io import load_particles
from zoom_mask_creator.params import read_param_file
from zoom_mask_creator.peano import compute_ic_positions


class Mask:
    def __init__(self, param_file, comm=None):
        """
        Class to construct and store a mask.

        Upon instantiation, the parameter file is read and the mask is created.

        Parameters
        ----------
        param_file : str
            The name of the YAML parameter file defining the mask. If None
            (default), a parsed parameter structure must be provided as
            `params` instead.
        comm : mpi4py communicator, optional
        """

        # MPI info.
        self.comm = comm
        if comm is not None:
            self.comm_rank = comm.rank
            self.comm_size = comm.size
        else:
            self.comm_rank = 0
            self.comm_size = 1

        # Load the parameter file.
        self.params = self._load_parameters(param_file)

        # Load the particles.
        self._load_particles()

        # Compute the mask from the particles.
        self._compute_mask()

    def _load_parameters(self, param_file):
        """
        Rank 0 reads the parameter file and shares it.

        Parameters
        ----------
        param_file : str
        """

        if self.comm_rank == 0:
            params = read_param_file(param_file)
        else:
            params = None

        if self.comm_size > 1:
            params = self.comm.bcast(params)

        return params

    def build_basic_mask(self, r, max_width):
        """
        Build the basic mask for an input particle distribution.

        This is a cubic boolean array with an adaptively computed cell size and
        extent that stretches by at least `min_width` in each dimension.  The
        mask value is True for any cells that contain at least the specified
        threshold number of particles.

        The mask is based on particles on all MPI ranks.

        Parameters
        ----------
        r : ndarray(float) [N_p, 3]
            The coordinates of (local) particles for which to create the mask.
            They must be shifted such that they lie within +/- `min_width` from
            the origin in each dimension.
        max_width : float
            The maximum extent of the mask along all six directions from the
            origin. It may get shrunk if `max_width` is larger than the whole
            box, but the mask will always remain centred on the origin. Note
            that this value must be identical across MPI ranks.

        Returns
        -------
        mask : ndarray(Bool) [N_cell, 3]
            The constructed boolean mask.
        edges : ndarray(N_cell + 1)
            The cell edges (same along each dimension)

        """
        # Total number of particles.
        n_tot = len(r)

        # Find out how far from the origin we need to extend the mask
        width = min(max_width, self.params["snapshot"]["bs"])

        # Work out how many cells we need along each dimension so that the
        # cells remain below the specified threshold size
        num_bins = int(np.ceil(2 * width / self.params["mask"]["mask_cell_size"]))

        # Compute number of particles in each cell, across MPI ranks
        n_p, edges = np.histogramdd(r, bins=num_bins, range=[(-width, width)] * 3)
        if self.comm_size > 1:
            n_p = self.comm.allreduce(n_p, op=MPI.SUM)
            n_tot = self.comm.allreduce(n_tot, op=MPI.SUM)

        # Make sure we found all the particles.
        if int(np.sum(n_p)) != n_tot:
            raise ValueError(
                f"Did not find all particles, increase padding factor, "
                f"binned: {int(np.sum(n_p))} should have: {len(r)}"
            )

        # Convert particle counts to True/False mask
        mask = n_p >= self.params["mask"]["min_num_per_cell"]

        return mask, edges[0]  # edges is a 3-tuple

    def _load_particles(self):

        if self.params["ics"]["ic_type"] == "map_to_ics":
            if self.comm_size > 1:
                raise ValueError("No MPI with IC mapping")

            self.ic_coords = load_particles(
                self.params, self.comm, self.comm_rank, self.comm_size
            )

        elif self.params["ics"]["ic_type"] == "use_peano_ids":
            # Load IDs of particles within target high-res region from snapshot.
            # Note that only particles assigned to current MPI rank are loaded,
            # which may be none.
            ids = load_particles(self.params, self.comm, self.comm_rank, self.comm_size)

            # Find initial positions from particle IDs (recall that these are
            # really Peano-Hilbert indices). Coordinates are in the same units
            # as the box size, centred (and wrapped) on the high-res target region.
            if len(ids) > 0:
                self.ic_coords = compute_ic_positions(ids, self.params, self.comm_rank)
            else:
                self.ic_coords = np.empty((0, 3), dtype=np.float32)

        # How many total particles.
        ntot = len(self.ic_coords)
        if self.comm_size > 1:
            ntot = self.comm.reduce(ntot)
        if self.comm_rank == 0:
            print(f"Loaded {ntot} total particles")

    def _compute_mask(self):

        # Find the corners of a box enclosing all particles in the ICs.
        box, widths = self.compute_bounding_box(self.ic_coords)
        if self.comm_rank == 0:
            print(
                f"Determined bounding box edges in ICs (re-centred):\n"
                f"\t{box[0, 0]:.3f} / {box[0, 1]:.3f} / {box[0, 2]:.3f} --> "
                f"{box[1, 0]:.3f} / {box[1, 1]:.3f} / {box[1, 2]:.3f}"
            )

        # For simplicity, shift the coordinates relative to geometric box
        # center, so that particles extend equally far in each direction
        geo_centre = box[0, :] + widths / 2
        self.ic_coords -= geo_centre

        # Also need to keep track of the mask centre in the original frame.
        self.mask_centre = geo_centre  # self.params["region"]["coords"] + geo_centre

        # Build the basic mask. This is a cubic boolean array with an
        # adaptively computed cell size and extent that includes at least
        # twice the entire bounding box. It is True for any cells that contain
        # at least the specified threshold number of particles.
        #
        # `edges` holds the spatial coordinate of the lower cell edges. By
        # construction, this is the same along all three dimensions.
        mask, edges = self.build_basic_mask(self.ic_coords, np.max(widths))
        self.cell_size = edges[1] - edges[0]

        # We only need MPI rank 0 for the rest, since we are done working with
        # individual particles
        if self.comm_rank > 0:
            return None, None

        # Fill holes and extrude the mask. This has to be done separately
        # for each of the three projections.
        for idim, name in enumerate(["x-y", "y-z", "x-z"]):
            print(f"Topological extrision ({idim}/3, {name} plane)...")
            self.refine_mask(mask, idim)

        # Finally, we need to find the centre of all selected mask cells, and
        # the box enclosing all those cells
        ind_sel = np.where(mask)  # Note: 3-tuple of ndarrays!
        self.cell_coords = np.vstack(
            (edges[ind_sel[0]], edges[ind_sel[1]], edges[ind_sel[2]])
        ).T
        self.cell_coords += self.cell_size * 0.5

        # Find the box that (fully) encloses all selected cells, and the
        # side length of its surrounding cube
        self.mask_box, mask_widths = self.compute_bounding_box(
            self.cell_coords, serial_only=True
        )
        self.mask_box[0, :] -= self.cell_size * 0.5
        self.mask_box[1, :] += self.cell_size * 0.5
        mask_widths += self.cell_size

        # Special case of slab
        if self.params["region"]["shape"] == "slab":
            raise NotImplementedError()
            # self.params["region"]["slab_dim"] = np.argmin(mask_widths)
            # for ii in range(3):
            #    if ii == self.params["region"]["slab_dim"]:
            #        continue
            #    mask_widths[ii] = self.params["snapshot"]["bs"]

        self.mask_widths = mask_widths
        self.mask_extent = np.max(mask_widths)

        print(
            f"Encompassing dimensions of the mask cells:\n"
            f"\tx = {mask_widths[0]:.4f} Mpc/h\n"
            f"\ty = {mask_widths[1]:.4f} Mpc/h\n"
            f"\tz = {mask_widths[2]:.4f} Mpc/h\n"
            f"Bounding length: {self.mask_extent:.4f} Mpc/h"
        )

        box_volume = np.prod(mask_widths)
        n_sel = len(ind_sel[0])
        cell_fraction = n_sel * self.cell_size**3 / box_volume
        cell_fraction_cube = n_sel * self.cell_size**3 / self.mask_extent**3
        print(f"There are {n_sel:d} selected mask cells.")
        print(
            f"They fill {cell_fraction * 100:.3f} per cent of the bounding "
            f"box ({cell_fraction_cube * 100:.3f} per cent of bounding "
            f"cube)."
        )

        # Final sanity check: make sure that all particles are within cubic
        # bounding box. This is allowed however, as max_per_cell > 1 can put particles
        # outside the region, so we leave this as a warning.
        if not np.all(box[0, :] - geo_centre >= -self.mask_extent / 2) or not np.all(
            box[1, :] - geo_centre <= self.mask_extent / 2
        ):
            print(
                "***Warning*** "
                f"Cubic bounding box around final mask does not enclose all "
                f"input particles!"
                "This may just be the choice of max_per_cell, check output "
                "plot to make sure mask looks sensible."
            )

    def refine_mask(self, mask, idim):
        """
        Refine the mask by checking for holes and processing the morphology.

        The refinement is performed iteratively along all slices along the
        specified axis. It consists of the ndimage operations ...

        Parameters
        ----------
        mask : ndarray
            The raw mask before refinement (modified in place)
        idim : int
            The perpendicular to which slices of the mask are to be processed
            (0: x, 1: y, 2: z)

        Returns
        -------
        None

        """
        # Process each layer (slice) of the mask in turn
        for layer_id in range(mask.shape[idim]):

            # Since each dimension loops over a different axis, set up an
            # index for the current layer
            if idim == 0:
                index = np.s_[layer_id, :, :]
            elif idim == 1:
                index = np.s_[:, layer_id, :]
            elif idim == 2:
                index = np.s_[:, :, layer_id]
            else:
                raise ValueError(f"Invalid value idim={idim}!")

            # Step 1: fill holes in the mask
            if self.params["mask"]["topology_fill_holes"]:
                mask[index] = ndimage.binary_fill_holes(mask[index]).astype(bool)
            # Step 2: regularize the morphology
            if self.params["mask"]["topology_dilation_niter"] > 0:
                mask[index] = ndimage.binary_dilation(
                    mask[index],
                    iterations=self.params["mask"]["topology_dilation_niter"],
                ).astype(bool)
            if self.params["mask"]["topology_closing_niter"] > 0:
                mask[index] = ndimage.binary_closing(
                    mask[index],
                    iterations=self.params["mask"]["topology_closing_niter"],
                ).astype(bool)

    def compute_bounding_box(self, r, serial_only=False):
        """
        Find the corners of a box enclosing a set of points across MPI ranks.

        Parameters
        ----------
        r : ndarray(float) [N_part, 3]
            The coordinates of the `N_part` particles held on this MPI rank.
            The second array dimension holds the x/y/z components per point.
        serial_only : bool, optional
            Switch to disable cross-MPI comparison of box extent (default:
            False). If True, the return values will generally differ between
            different MPI ranks.

        Returns
        -------
        box : ndarray(float) [2, 3]
            The coordinates of the lower and upper vertices of the bounding
            box. These are stored in index 0 and 1 along the first dimension,
            respectively.

        widths : ndarray(float) [3]
            The width of the box along each dimension.

        """
        box = np.zeros((2, 3))

        # Find vertices of local particles (on this MPI rank). If there are
        # none, set lower (upper) vertices to very large (very negative)
        # numbers so that they will not influence the cross-MPI min/max.
        n_part = len(r)
        box[0, :] = np.min(r, axis=0) if n_part > 0 else sys.float_info.max
        box[1, :] = np.max(r, axis=0) if n_part > 0 else -sys.float_info.max

        # Now compare min/max values across all MPI ranks
        if (not serial_only) and (self.comm_size > 1):
            for idim in range(3):
                box[0, idim] = self.comm.allreduce(box[0, idim], op=MPI.MIN)
                box[1, idim] = self.comm.allreduce(box[1, idim], op=MPI.MAX)

        return box, box[1, :] - box[0, :]
