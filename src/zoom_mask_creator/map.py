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
    """

    # Load the ICs
    print(f"Loading IC file to map to {fname}")
    with h5py.File(fname) as f:
        ics_ids = f["PartType1/ParticleIDs"][...]
        ics_coords = f["PartType1/Coordinates"][...]

    # Sort by ID
    idx = np.argsort(ics_ids)
    ics_ids = ics_ids[idx]
    ics_coords = ics_coords[idx]

    ids = np.sort(ids)

    # Match
    idx = np.searchsorted(ics_ids, ids)

    assert np.array_equal(ics_ids[idx], ids)

    return ics_coords[idx]
