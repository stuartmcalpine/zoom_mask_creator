import os

import h5py
import numpy as np

HAVE_EAGLE = True
HAVE_SWIFT = True

# Do we have an EAGLE read routine?
try:
    from pyread_eagle import EagleSnapshot
except ImportError:
    HAVE_EAGLE = False

# Do we have a swift read routine?
try:
    from pyread_swift import SwiftSnapshot
except ImportError:
    HAVE_SWIFT = False

# Make sure we have at least one kind of read routine.
assert HAVE_EAGLE or HAVE_SWIFT, "No read libs found"


def periodic_wrapping(r, boxsize, return_copy=False):
    """
    Apply periodic wrapping to an input set of coordinates.

    Parameters
    ----------
    r : ndarray(float) [N, 3]
        The coordinates to wrap.
    boxsize : float
        The box size to wrap the coordinates to. The units must correspond to
        those used for `r`.
    return_copy : bool, optional
        Switch to return a (modified) copy of the input array, rather than
        modifying the input in place (which is the default).

    Returns
    -------
    r_wrapped : ndarray(float) [N, 3]
        The wrapped coordinates. Only returned if `return_copy` is True,
        otherwise the input array `r` is modified in-place.

    """
    if return_copy:
        r_wrapped = (r + 0.5 * boxsize) % boxsize - 0.5 * boxsize
        return r_wrapped

    # To perform the wrapping in-place, break it down into three steps
    r += 0.5 * boxsize
    r %= boxsize
    r -= 0.5 * boxsize


def _find_enclosing_frame(params, comm_rank):
    """
    Compute the bounding box that encloses the target zoom region at the target
    snapshot.

    This is only used to pre-select a region for particle reading from
    the snapshot, to make things more efficient.

    Parameters
    ----------
    params : dict
        The parameters of the run
    comm_rank : int
        Rank of this core

    Returns
    -------
    frame : np.ndarray(float)
        A 2x3 element array containing the lower and upper coordinates
        of the bounding region in the x, y, and z coordinates.

    """
    input_centre = params["input_centre"]
    frame = np.zeros((2, 3))

    # If the target region is a sphere, find the enclosing cube
    if params["shape"] == "sphere":
        frame[0, :] = input_centre - params["radius"]
        frame[1, :] = input_centre + params["radius"]

    # If the target region is a cuboid, simply transform from input_centre and
    # side length to lower and upper bounds along each coordinate
    elif params["shape"] in ["cuboid", "slab"]:
        frame[0, :] = input_centre - params["dim"] / 2.0
        frame[1, :] = input_centre + params["dim"] / 2.0

    if comm_rank == 0:
        print(
            f"Boundary frame in selection snapshot:\n"
            f"{frame[0, 0]:.2f} / {frame[0, 1]:.2f} / {frame[0, 2]:.2f}"
            " -- "
            f"{frame[1, 0]:.2f} / {frame[1, 1]:.2f} / {frame[1, 2]:.2f}"
        )

    return frame


def _convert_lengths_to_inverse_h(params):
    """
    Convert length parameters from h-free to 'h^-1' units (for IC-GEN).
    
    Modifies the length parameters in the 'params' dict directly.

    Parameters
    ----------
    params : dict
        The parameters of the run
    """
    h = params["h_factor"]

    for att in ["radius", "dim", "input_centre", "bs", "mask_cell_size"]:
        if att in params.keys():
            params[att] *= h

    params["length_unit"] = 'Mpc/h'

def load_particles(params, comm, comm_rank, comm_size):
    """
    Load the dark matter ParticleIDs from the target snapshot.

    In addition, relevant metadata are also loaded and stored in the `params`
    dict.

    Parameters
    ----------
    params : dict
        The parameters of the run
    comm : mpi4py comm
    comm_rank : int
    comm_size : int

    Returns
    -------
    ids : ndarray(int)
        The Peano-Hilbert keys of the particles.
    """

    # Find cuboidal frame enclosing the *target* high-resolution region
    load_region = _find_enclosing_frame(params, comm_rank)

    # To make life simpler, extractsome frequently used parameters
    cen = params["input_centre"]

    # First step: set up particle reader and load metadata.
    # This is different for SWIFT and GADGET simulations, but subsequent
    # loading
    if params["data_type"].lower() == "eagle":
        assert HAVE_EAGLE, "No EAGLE read routine found"
        snap = EagleSnapshot(params["snap_file"])
        params["bs"] = float(snap.boxsize)
        params["h_factor"] = 1.0
        params["length_unit"] = "Mph/h"
        params["redshift"] = -1.0
        snap.select_region(*load_region.T.flatten())
        if comm_size > 1:
            snap.split_selection(comm_rank, comm_size)
    elif params["data_type"].lower() == "swift":
        assert HAVE_SWIFT, "No SWIFT read routine found"
        snap = SwiftSnapshot(params["snap_file"], comm=comm)
        params["bs"] = float(snap.header["BoxSize"])
        params["h_factor"] = float(snap.header["h"])
        params["length_unit"] = "Mpc"
        params["redshift"] = snap.header["Redshift"]
        snap.select_region(1, *load_region.T.flatten())
        snap.split_selection()

    if comm_rank == 0:
        params["zred_snap"] = params["redshift"]
        print(f"Snapshot is at redshift z = {params['zred_snap']:.2f}.")

    # Load DM particle IDs and coordinates (uniform across GADGET/SWIFT)
    if comm_rank == 0:
        print("\nLoading particle data...")
    coords = snap.read_dataset(1, "Coordinates")

    # Shift coordinates relative to target centre, and wrap them to within
    # the periodic box (done by first shifting them up by half a box,
    # taking the modulus with the box size in each dimension, and then
    # shifting it back down by half a box)
    cen = params["input_centre"]
    coords -= cen
    periodic_wrapping(coords, params["bs"])

    # Select particles within target region
    if params['shape'] == "sphere":
        if comm_rank == 0:
            print(
                f"Clipping to sphere around {cen}, with radius "
                f"{params['radius']:.4f} {params['length_unit']}"
            )

        dists = np.linalg.norm(coords, axis=1)
        mask = np.where(dists <= params["radius"])[0]

    elif params["shape"] in ["cuboid", "slab"]:
        if comm_rank == 0:
            print(
                f"Clipping to {params['shape']} with "
                f"dx={params['dim'][0]:.2f} {params['length_unit']}, "
                f"dy={params['dim'][1]:.2f} {params['length_unit']}, "
                f"dz={params['dim'][2]:.2f} {params['length_unit']}\n"
                f"around {cen} {params['length_unit']}."
            )

        # To find particles within target cuboid, normalize each coordinate
        # offset by the maximum allowed extent in the corresponding
        # dimension, and find those where the result is between -1 and 1
        # for all three dimensions
        mask = np.where(np.max(np.abs(coords / (params["dim"] / 2)), axis=1) <= 1)[0]

    # Secondly, we need the IDs of particles lying in the mask region
    ids = snap.read_dataset(1, "ParticleIDs")[mask]

    # If IDs are Peano-Hilbert indices multiplied by two (as in e.g.
    # simulations with baryons), need to undo this multiplication here
    if params["divide_ids_by_two"]:
        ids = ids // 2

    print(f"[Rank {comm_rank}] Loaded {len(ids)} dark matter particles")

    # If the snapshot is from a user-friendly SWIFT simulation, all
    # lengths are in 'h-free' coordinates. Unfortunately, the ICs still
    # assume 'h^-1' units, so for consistency we now have to multiply
    # h factor back in...
    if params["data_type"].lower() == "swift":
        _convert_lengths_to_inverse_h(params)

    return ids


def save_mask(mask):
    """
    Save the generated mask.

    Note that, for consistency with IC GEN, all length dimensions must
    be in units of h^-1. This is already taken care of.

    Only rank 0 should do this, but we'll check just to be sure...

    Parameters
    -------
    mask : Mask object
        Has within it the cell positions from the computed mask
    """
    
    if mask.comm_rank != 0:
        return

    outloc = os.path.join(mask.params["output_dir"], mask.params["fname"]) + ".hdf5"

    with h5py.File(outloc, "w") as f:

        # Push parameter file data as file attributes
        g = f.create_group("Params")
        for param_attr in mask.params:
            g.attrs.create(param_attr, mask.params[param_attr])

        # Main output is the centres of selected mask cells
        ds = f.create_dataset(
            "Coordinates", data=np.array(mask.cell_coords, dtype="f8")
        )
        ds.attrs.create(
            "Description",
            "Coordinates of the centres of selected mask "
            "cells. The (uniform) cell width is stored as the "
            "attribute `grid_cell_width`.",
        )

        # Store attributes directly related to the mask as HDF5 attributes.
        ds.attrs.create("bounding_length", mask.mask_extent)
        ds.attrs.create("geo_centre", mask.mask_centre)
        ds.attrs.create("grid_cell_width", mask.cell_size)
        if mask.params["shape"] in ["cuboid", "slab"]:
            high_res_volume = np.prod(mask.params["dim"])

            if mask.params["shape"] == "slab":
                ds.attrs.create("slab_dim", mask.params["slab_dim"])
        else:
            high_res_volume = 4 / 3.0 * np.pi * mask.params["radius"] ** 3.0
        ds.attrs.create("high_res_volume", high_res_volume)
        ds.attrs.create("mask_widths", mask.mask_widths)
        ds.attrs.create("shape", mask.params["shape"])

    print(f"Saved mask data to file `{outloc}`.")
