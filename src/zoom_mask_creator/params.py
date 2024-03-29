import os

import numpy as np
import yaml


def read_param_file(param_file):
    """
    Read parameters from a specified YAML file.

    See template file `param_files/template.yml` for a full listing and
    description of parameters. The values are return in the 'params` dict.

    If the parameter file specifies the target centre in terms of a
    particular Velociraptor halo, the centre and radius of the high-res
    region are determined internally.

    In the MPI version, the file is only read on one rank and the dict
    then broadcast to all other ranks.

    Parameters
    ----------
    param_file : string
        Path to YAML parameter file

    Returns
    -------
    params : dict
    """

    # Set default values for optional parameters
    params = {}
    params["min_num_per_cell"] = 3
    params["mask_cell_size"] = 3.0
    params["topology_fill_holes"] = True
    params["topology_dilation_niter"] = 0
    params["topology_closing_niter"] = 0

    # Read param file.
    read_params = yaml.safe_load(open(param_file))

    # Define a list of parameters that must be provided. An error
    # is raised if they are not found in the YAML file.
    required_params = [
        "fname",
        "snap_file",
        "bits",
        "data_type",
        "divide_ids_by_two",
        "output_dir",
        "input_centre",
        "shape",
    ]
    for att in required_params:
        if att not in read_params:
            raise KeyError(
                f"Need to provide a value for {att} in the parameter "
                f"file '{param_file}'!"
            )

    # Ingest.
    for att, v in read_params.items():
        params[att] = v

    # Consistency checks for manual target region selection
    if "input_centre" not in params.keys():
        raise KeyError(
            "Need to provide coordinates for the input_centre of the "
            "high-resolution region."
        )
    if "shape" not in params.keys():
        raise KeyError("Need to specify the shape of the target region!")
    if params["shape"] in ["cuboid", "slab"] and "dim" not in params.keys():
        raise KeyError(
            f"Need to provide dimensions of {params['shape']} "
            f"high-resolution region."
        )
    if params["shape"] == "sphere" and "radius" not in params.keys():
        raise KeyError(
            "Need to provide the radius of target high-resolution " "sphere!"
        )

    # Convert coordinates and cuboid/slab dimensions to ndarray
    params["input_centre"] = np.array(params["input_centre"], dtype="f8")
    if "dim" in params.keys():
        params["dim"] = np.array(params["dim"], dtype="f8")

    # Create the output directory if it does not exist yet
    if not os.path.isdir(params["output_dir"]):
        os.makedirs(params["output_dir"])

    return params
