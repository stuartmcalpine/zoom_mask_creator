import numpy as np
import h5py


def map_to_ics(ids, fname):
    """
    Find a given set of particles in an IC file, map, and find initial coordinates.

    Parameters
    ----------
    ids : ndarray
        List of ParticleIDs we want to match to
    fname : str
        Path to the IC file

    Returns
    -------
    coords : ndarray
        Coordinates of the particles in the ICs
    """
    # Load the ICs
    print(f"Loading IC file {fname}")
    with h5py.File(fname) as f:
        ics_ids = f["PartType1/ParticleIDs"][...]
        ics_coords = f["PartType1/Coordinates"][...]

    # Sort by ID for efficient searching
    print(f"Sorting {len(ics_ids)} IC particles by ID")
    idx = np.argsort(ics_ids)
    ics_ids = ics_ids[idx]
    ics_coords = ics_coords[idx]

    # Sort the input IDs for efficient matching
    print(f"Sorting {len(ids)} input particles")
    sort_idx = np.argsort(ids)
    sorted_ids = ids[sort_idx]

    # Prepare array for coordinates
    coords = np.zeros((len(ids), 3), dtype=ics_coords.dtype)

    # Match in batches to avoid memory issues with very large arrays
    print("Matching particles to ICs")
    batch_size = 1000000  # Adjust based on available memory
    found_count = 0

    for i in range(0, len(sorted_ids), batch_size):
        batch_end = min(i + batch_size, len(sorted_ids))
        batch_ids = sorted_ids[i:batch_end]

        # Find matching indices
        idx_in_ics = np.searchsorted(ics_ids, batch_ids)

        # Check which locations actually contain matching IDs
        valid_indices = (idx_in_ics < len(ics_ids)) & (ics_ids[idx_in_ics] == batch_ids)

        # Get the original indices in the unsorted ids array
        original_indices = sort_idx[i:batch_end][valid_indices]

        # Get the matched coordinates
        matched_coords = ics_coords[idx_in_ics[valid_indices]]

        # Assign coordinates
        coords[original_indices] = matched_coords
        found_count += np.sum(valid_indices)

        if i + batch_size >= len(sorted_ids) or (i > 0 and i % (10 * batch_size) == 0):
            print(f"  Matched {found_count}/{len(ids)} particles...")

    # Check if all particles were found
    if found_count < len(ids):
        missing_count = len(ids) - found_count
        missing_percent = missing_count / len(ids) * 100
        print(
            f"WARNING: {missing_count} particles ({missing_percent:.2f}%) not found in IC file!"
        )

        # Optionally, raise an error if too many particles are missing
        if missing_percent > 1.0:  # More than 1% missing
            raise ValueError(f"Too many particles missing: {missing_count}/{len(ids)}")

    return coords
