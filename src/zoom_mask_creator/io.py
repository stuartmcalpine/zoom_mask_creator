import os

import h5py
import numpy as np

from .map import map_to_ics

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
    frame = np.zeros((2, 3))

    # If the target region is a sphere, find the enclosing cube
    if params["region"]["shape"] == "sphere":
        frame[0, :] = params["region"]["coords"] - params["region"]["radius"]
        frame[1, :] = params["region"]["coords"] + params["region"]["radius"]

    # If the target region is a cuboid, simply transform from input_centre and
    # side length to lower and upper bounds along each coordinate
    elif params["region"]["shape"] in ["cuboid", "slab"]:
        frame[0, :] = params["region"]["coords"] - params["region"]["dim"] / 2.0
        frame[1, :] = params["region"]["coords"] + params["region"]["dim"] / 2.0

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

    if params["snapshot"]["length_unit"] != "Mpc":
        raise ValueError(f"Should not be converting this type")

    h = params["snapshot"]["h_factor"]

    _to_change = {
        "region": ["radius", "dim", "coords"],
        "snapshot": ["bs"],
        "mask": ["mask_cell_size"],
    }

    for cat in _to_change.keys():
        for att in _to_change[cat]:
            if att in params[cat].keys():
                print(f"Converting {att} to h units")
                params[cat][att] *= h

    params["snapshot"]["length_unit"] = "Mpc/h"


def validate_snapshot_consistency(params, comm_rank):
    """
    Validate that all snapshot files have consistent properties
    """
    if len(params["snapshot"]["paths"]) <= 1:
        return True

    # Properties to check across files
    boxsizes = []
    h_factors = []
    redshifts = []

    for path in params["snapshot"]["paths"]:
        if params["snapshot"]["data_type"].lower() == "eagle":
            snap = EagleSnapshot(path)
            boxsizes.append(float(snap.boxsize))
            h_factors.append(1.0)  # EAGLE uses h=1
            redshifts.append(-1.0)  # EAGLE doesn't store this
        elif params["snapshot"]["data_type"].lower() == "swift":
            snap = SwiftSnapshot(path, comm=None)  # Just for metadata
            boxsizes.append(float(snap.header["BoxSize"]))
            h_factors.append(float(snap.header["h"]))
            redshifts.append(snap.header["Redshift"])

    # Check consistency
    if not np.allclose(boxsizes, boxsizes[0], rtol=1e-5):
        if comm_rank == 0:
            print("WARNING: Box sizes differ across snapshot files!")
            print(f"Box sizes: {boxsizes}")

    if not np.allclose(h_factors, h_factors[0], rtol=1e-5):
        if comm_rank == 0:
            print("WARNING: Hubble factors differ across snapshot files!")
            print(f"h factors: {h_factors}")

    if not np.allclose(redshifts, redshifts[0], rtol=1e-2):
        if comm_rank == 0:
            print("WARNING: Redshifts differ across snapshot files!")
            print(f"Redshifts: {redshifts}")

    return True


def load_particles(params, comm, comm_rank, comm_size):
    """
    Load the dark matter ParticleIDs from multiple target snapshots.

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
    ids : ndarray(int) or ndarray(float) [N, 3]
        The Peano-Hilbert keys of the particles or their coordinates in the ICs.
    """

    # Find cuboidal frame enclosing the *target* high-resolution region
    load_region = _find_enclosing_frame(params, comm_rank)

    # Check snapshot consistency if we have multiple files
    if len(params["snapshot"]["paths"]) > 1:
        validate_snapshot_consistency(params, comm_rank)

    # Handle different workflows depending on the IC type
    if params["ics"]["ic_type"] == "use_peano_ids":
        # For Peano-Hilbert IDs, we can just load and concatenate all IDs
        all_ids = []

        if comm_rank == 0:
            print(f"Processing {len(params['snapshot']['paths'])} snapshot files")

        # Process each snapshot file
        for file_idx, snapshot_path in enumerate(params["snapshot"]["paths"]):
            if comm_rank == 0:
                print(
                    f"\nProcessing file {file_idx+1}/{len(params['snapshot']['paths'])}: {snapshot_path}"
                )

            # Load particles from this file
            file_ids = _load_particles_from_file(
                snapshot_path,
                params,
                load_region,
                file_idx,
                len(params["snapshot"]["paths"]) - 1,
                comm,
                comm_rank,
                comm_size,
            )

            # Handle Peano-Hilbert IDs if needed
            if params["peano"].get("divide_ids_by_two", False):
                file_ids = file_ids // 2

            # Add this file's IDs to our collection
            all_ids.append(file_ids)

            if comm_rank == 0:
                print(f"  Found {len(file_ids)} particles in file {file_idx+1}")

        # Concatenate all particle IDs
        if len(all_ids) > 0:
            ids = np.concatenate(all_ids)
        else:
            ids = np.array([], dtype=np.int64)

        if comm_rank == 0:
            print(f"Total particles found across all files: {len(ids)}")

        print(f"[Rank {comm_rank}] Loaded {len(ids)} dark matter particles")
        return ids

    elif params["ics"]["ic_type"] == "map_to_ics":
        # For IC mapping, we need to process each snapshot-IC file pair
        # and combine the coordinates at the end
        all_ic_coords = []

        # Check if we have one IC file per snapshot or just one IC file
        if (
            len(params["map_to_ics"]["paths"]) == 1
            and len(params["snapshot"]["paths"]) > 1
        ):
            # One IC file for all snapshots - we'll load all IDs first, then map
            all_ids = []

            # Process each snapshot file to collect IDs
            for file_idx, snapshot_path in enumerate(params["snapshot"]["paths"]):
                if comm_rank == 0:
                    print(
                        f"\nProcessing file {file_idx+1}/{len(params['snapshot']['paths'])}: {snapshot_path}"
                    )

                # Load particles from this file
                file_ids = _load_particles_from_file(
                    snapshot_path,
                    params,
                    load_region,
                    file_idx,
                    len(params["snapshot"]["paths"]) - 1,
                    comm,
                    comm_rank,
                    comm_size,
                )

                # Add this file's IDs to our collection
                all_ids.append(file_ids)

                if comm_rank == 0:
                    print(f"  Found {len(file_ids)} particles in file {file_idx+1}")

            # Concatenate all particle IDs
            if len(all_ids) > 0:
                ids = np.concatenate(all_ids)
            else:
                ids = np.array([], dtype=np.int64)

            if comm_rank == 0:
                print(f"Total particles found across all files: {len(ids)}")

            # Map these IDs to coordinates using the single IC file
            ic_coords = map_to_ics(ids, params["map_to_ics"]["paths"][0])

            # Multiply by h as only swift snapshots can get this far
            ic_coords *= params["snapshot"]["h_factor"]

            return ic_coords

        else:
            # Multi-snapshot, multi-IC case
            # Determine which snapshot maps to which IC file
            if len(params["map_to_ics"]["paths"]) == len(params["snapshot"]["paths"]):
                # One-to-one mapping
                snapshot_to_ic_map = list(range(len(params["snapshot"]["paths"])))
            else:
                # Use modulo arithmetic to map snapshots to ICs
                snapshot_to_ic_map = [
                    i % len(params["map_to_ics"]["paths"])
                    for i in range(len(params["snapshot"]["paths"]))
                ]

            if comm_rank == 0:
                print(
                    f"Processing {len(params['snapshot']['paths'])} snapshot files "
                    f"with {len(params['map_to_ics']['paths'])} IC files"
                )

            # Process each snapshot file
            for file_idx, snapshot_path in enumerate(params["snapshot"]["paths"]):
                if comm_rank == 0:
                    print(
                        f"\nProcessing file {file_idx+1}/{len(params['snapshot']['paths'])}: {snapshot_path}"
                    )

                # Load particles from this file
                file_ids = _load_particles_from_file(
                    snapshot_path,
                    params,
                    load_region,
                    file_idx,
                    len(params["snapshot"]["paths"]) - 1,
                    comm,
                    comm_rank,
                    comm_size,
                )

                if len(file_ids) == 0:
                    if comm_rank == 0:
                        print(f"  No particles found in file {file_idx+1}, skipping")
                    continue

                # Determine which IC file to use for this snapshot
                ic_idx = snapshot_to_ic_map[file_idx]
                ic_path = params["map_to_ics"]["paths"][ic_idx]

                if comm_rank == 0:
                    print(
                        f"  Mapping {len(file_ids)} particles to IC file {ic_idx+1}: {ic_path}"
                    )

                # Map these IDs to coordinates
                file_ic_coords = map_to_ics(file_ids, ic_path)

                # Store these coordinates
                all_ic_coords.append(file_ic_coords)

                if comm_rank == 0:
                    print(
                        f"  Successfully mapped {len(file_ic_coords)} particles from file {file_idx+1}"
                    )

            # Combine all coordinates
            if len(all_ic_coords) > 0:
                ic_coords = np.vstack(all_ic_coords)
            else:
                ic_coords = np.zeros((0, 3), dtype=np.float32)

            if comm_rank == 0:
                print(f"Total particles mapped across all files: {len(ic_coords)}")

            # Multiply by h as only swift snapshots can get this far
            ic_coords *= params["snapshot"]["h_factor"]

            return ic_coords


def _load_particles_from_file(
    snapshot_path,
    params,
    load_region,
    file_idx,
    last_file_idx,
    comm,
    comm_rank,
    comm_size,
):
    """
    Helper function to load particles from a single snapshot file.

    Parameters
    ----------
    snapshot_path : str
        Path to the snapshot file
    params : dict
        The parameters of the run
    load_region : ndarray
        The region to load particles from
    file_idx : int
        Index of the file being processed
    last_file_idx : int
        Index of the last file to be preocessed
    comm : mpi4py comm
    comm_rank : int
    comm_size : int

    Returns
    -------
    file_ids : ndarray(int)
        The IDs of particles found in this file
    """
    # To make life simpler, extract frequently used parameters
    cen = params["region"]["coords"]

    # First step: set up particle reader and load metadata.
    # This is different for SWIFT and GADGET simulations
    if params["snapshot"]["data_type"].lower() == "eagle":
        assert HAVE_EAGLE, "No EAGLE read routine found"
        snap = EagleSnapshot(snapshot_path)

        # Only set these params from the first file if they're not already set
        if file_idx == 0:
            params["snapshot"]["bs"] = float(snap.boxsize)
            params["snapshot"]["h_factor"] = 1.0
            params["snapshot"]["length_unit"] = "Mph/h"
            params["snapshot"]["redshift"] = -1.0

        snap.select_region(*load_region.T.flatten())
        if comm_size > 1:
            snap.split_selection(comm_rank, comm_size)

    elif params["snapshot"]["data_type"].lower() == "swift":
        assert HAVE_SWIFT, "No SWIFT read routine found"
        snap = SwiftSnapshot(snapshot_path, comm=comm)

        # Only set these params from the first file if they're not already set
        if file_idx == 0:
            params["snapshot"]["bs"] = float(snap.header["BoxSize"])
            params["snapshot"]["h_factor"] = float(snap.header["h"])
            params["snapshot"]["length_unit"] = "Mpc"
            params["snapshot"]["redshift"] = snap.header["Redshift"]

        snap.select_region(1, *load_region.T.flatten())
        snap.split_selection()

    if file_idx == 0 and comm_rank == 0:
        params["snapshot"]["zred_snap"] = params["snapshot"]["redshift"]
        print(
            f"First snapshot is at redshift z = {params['snapshot']['zred_snap']:.2f}."
        )

    # Load DM particle coordinates and select within target region
    if comm_rank == 0:
        print(f"Loading particle data from {snapshot_path}...")

    coords = snap.read_dataset(1, "Coordinates")

    # Handle empty selection case
    if len(coords) == 0:
        return np.array([], dtype=np.int64)

    # Shift coordinates relative to target centre and wrap
    coords -= cen
    periodic_wrapping(coords, params["snapshot"]["bs"])

    # Select particles within target region
    if params["region"]["shape"] == "sphere":
        if comm_rank == 0:
            print(
                f"Clipping to sphere around {cen}, with radius "
                f"{params['region']['radius']:.4f} {params['snapshot']['length_unit']}"
            )

        dists = np.linalg.norm(coords, axis=1)
        mask = np.where(dists <= params["region"]["radius"])[0]

    elif params["region"]["shape"] in ["cuboid", "slab"]:
        if comm_rank == 0:
            print(
                f"Clipping to {params['region']['shape']} with "
                f"dx={params['region']['dim'][0]:.2f} {params['snapshot']['length_unit']}, "
                f"dy={params['region']['dim'][1]:.2f} {params['snapshot']['length_unit']}, "
                f"dz={params['region']['dim'][2]:.2f} {params['snapshot']['length_unit']}\n"
                f"around {cen} {params['snapshot']['length_unit']}."
            )

        # To find particles within target cuboid, normalize each coordinate
        # offset by the maximum allowed extent in the corresponding
        # dimension, and find those where the result is between -1 and 1
        # for all three dimensions
        mask = np.where(
            np.max(np.abs(coords / (params["region"]["dim"] / 2)), axis=1) <= 1
        )[0]

    # Get particle IDs for the selected particles
    del coords
    file_ids = snap.read_dataset(1, "ParticleIDs")[mask]

    # If the snapshot is from SWIFT in the first file, convert lengths to h^-1 units
    if params["snapshot"]["data_type"].lower() == "swift" and file_idx == last_file_idx:
        _convert_lengths_to_inverse_h(params)

    return file_ids


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

    outloc = os.path.join(mask.params["output"]["path"], "mask") + ".hdf5"

    with h5py.File(outloc, "w") as f:

        # Push parameter file data as file attributes
        g = f.create_group("Params")
        for cat in mask.params.keys():
            for param_attr in mask.params[cat]:
                g.attrs.create(f"{cat}:{param_attr}", mask.params[cat][param_attr])

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
        if mask.params["region"]["shape"] in ["cuboid", "slab"]:
            high_res_volume = np.prod(mask.params["region"]["dim"])

            if mask.params["region"]["shape"] == "slab":
                ds.attrs.create("slab_dim", mask.params["region"]["slab_dim"])
        else:
            high_res_volume = 4 / 3.0 * np.pi * mask.params["region"]["radius"] ** 3.0
        ds.attrs.create("high_res_volume", high_res_volume)
        ds.attrs.create("mask_widths", mask.mask_widths)
        ds.attrs.create("shape", mask.params["region"]["shape"])

    print(f"Saved mask data to file `{outloc}`.")
