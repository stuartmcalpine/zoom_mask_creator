import os

import numpy as np
import yaml

import toml

# Default parameters
_DEFAULTS = {
    "region": {"shape": "sphere"},
    "mask": {
        "min_num_per_cell": 3,
        "mask_cell_size": 3.0,
        "topology_fill_holes": True,
        "topology_dilation_niter": 0,
        "topology_closing_niter": 0,
    },
}


def read_param_file(param_file):
    """
    Read parameters from a specified TOML file.

    See template file `param_files/template.yml` for a full listing and
    description of parameters. The values are return in the 'params` dict.

    If the parameter file specifies the target centre in terms of a
    particular Velociraptor halo, the centre and radius of the high-res
    region are determined internally.

    In the MPI version, the file is only read on one rank and the dict
    then broadcast to all other ranks.

    Parameters
    ----------
    param_file : str
        Path to TOML parameter file

    Returns
    -------
    params : dict
    """

    # Load parameters from TOML file
    params = toml.load(param_file)

    # Fill in the defaults
    for cat in _DEFAULTS.keys():
        if cat in params.keys():
            for att in _DEFAULTS[cat].keys():
                if att not in params[cat].keys():
                    params[cat][att] = _DEFAULTS[cat][att]
        else:
            params[cat] = _DEFAULTS[cat]

    # Make sure we have the minimum required params
    _required = {
        "snapshot": ["path", "data_type"],
        "ics": ["ic_type"],
        "output": ["path"],
    }

    for cat in _required:
        if cat not in params.keys():
            raise ValueError(f"Missing catagory {cat}")
        for att in _required[cat]:
            if att not in params[cat].keys():
                raise ValueError(f"Missing required param {cat}:{att}")

    if params["ics"]["ic_type"] == "use_peano_ids":
        _required = ["bits", "divide_ids_by_two"]

        if "peano" not in params.keys():
            raise ValueError(f"Need [peano] section with {_required}")
        for att in _required:
            if att not in params["peano"].keys():
                raise ValueError(f"Missing param peano:{att}")
    elif params["ics"]["ic_type"] == "map_to_ics":
        _required = ["path"]
        if "map_to_ics" not in params.keys():
            raise ValueError(f"Need [map_to_ics] section with {_required}")
        for att in _required:
            if att not in params["map_to_ics"].keys():
                raise ValueError(f"Missing params map_to_ics:{att}")
        # Needs to be swift snapshot for this mapping (assumed at least)
        if params["snapshot"]["data_type"].lower() != "swift":
            raise ValueError("Need to be swift snapshot to IC mapping")

    # Checks
    if params["region"]["shape"] not in ["sphere", "cuboid", "slab"]:
        raise ValueError(f"{params['region']['shape']} is a bad shape")

    if params["region"]["shape"] == "sphere":
        if "radius" not in params["region"].keys():
            raise ValueError("Need a radius for region shape of sphere")
    else:
        if "dim" not in params["region"].keys():
            raise ValueError("Need a dim for region shape of not sphere")

    _allowed_data_types = ["swift", "eagle"]
    if params["snapshot"]["data_type"] not in _allowed_data_types:
        raise ValueError(f"Allowed datatypes are {_allowed_data_types}")

    if not os.path.isfile(params["snapshot"]["path"]):
        raise FileNotFoundError(f"{params['snapshot']['path']} not found")

    _allowed_ic_types = ["use_peano_ids", "map_to_ics"]
    if params["ics"]["ic_type"] not in _allowed_ic_types:
        raise ValueError(f"Allowed ic tyoes are {_allowed_ic_types}")

    # Convert coordinates and cuboid/slab dimensions to ndarray
    params["region"]["coords"] = np.array(params["region"]["coords"], dtype="f8")
    if "dim" in params["region"].keys():
        if len(params["region"]["dim"]) != 3:
            raise ValueError("dim must be 3 dimentions [x,y,z]")
        params["region"]["dim"] = np.array(params["region"]["dim"], dtype="f8")

    # Create the output directory if it does not exist yet
    if not os.path.isdir(params["output"]["path"]):
        os.makedirs(params["output"]["path"])

    return params
