[![python](https://img.shields.io/badge/Python-3.7-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org)
[![python](https://img.shields.io/badge/Python-3.8-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org)
[![python](https://img.shields.io/badge/Python-3.9-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org)
[![python](https://img.shields.io/badge/Python-3.10-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## Zoom-in simulations mask creator

A script to compute the Lagrangian positions (i.e., positions in the initial
conditions) of a set of particles at *z*=0 from a  [Virgo
Consortium](https://virgo.dur.ac.uk/)-like cosmological simulation. The
continuous region those particles occupy in the ICs is then masked, and stored,
which can then go on to be used as input to generate a particle load for a
"zoom-in" resimulation of the desired target region at *z*=0.

The code is compatible with simulation outputs that can be read by
[pyread_eagle](https://github.com/kyleaoman/pyread_eagle) or
[read_swift](https://github.com/stuartmcalpine/read_swift). The code supports two methods for finding the Lagrangian positions of the particles:

1. Using *Peano Hilbert* indexing in the parent simulation to map particle positions back to initial conditions
2. Matching particle IDs between the target snapshot and the corresponding IC file

This script partners with 
[zoom_particle_load_creator](https://github.com/stuartmcalpine/zoom_particle_load_creator),
which creates particle loads using masks generated from these scripts.

## Installation

### Requirements

* `OpenMPI` or other MPI library
* `python>=3.8,<3.11`

Recommended modules when working on COSMA7:

* `module load gnu_comp/14.1.0 openmpi/5.0.3 parallel_hdf5/1.14.4 fftw/3.3.10`
* `module load python/3.9.19`

### Installation from source

It is recommended you install the package within a virtual/conda environment.
Or alternatively, if you are installing on a shared access machine, to your
local profile by including the `--user` flag during the `pip` installation. You can ofcourse also install directly into your base Python environment if you prefer.

First make sure your `pip` is up-to-date:

* `python3 -m pip install --upgrade pip`

Then you can install the `zoom_mask_creator` package by typing the following in
the git directory: 

* `python3 -m pip install -e .`

which will install `zoom_mask_creator` and any dependencies (as an editable install).

### MPI installation for `read_swift`

If you are using `read_swift` and want to load large snapshots over MPI collectively
(i.e., multiple cores read in parallel from the same file), a bit of additional
setup is required.

Make sure you have MPI libraries installed on your machine (`OpenMPI` for example), and you have `hdf5` installed with **parallel** compatibility ([see here for details](https://docs.h5py.org/en/stable/mpi.html)).

First, uninstall any installed versions of `mpi4py` and `h5py`:

* `python3 -m pip uninstall mpi4py h5py`

Then reinstall `mpi4py` and `h5py` from source with MPI flags:

* `MPICC=mpicc CC=mpicc HDF5_MPI="ON" python3 -m pip install --no-binary=mpi4py,h5py mpi4py h5py`

If `pip` can't find your `HDF5` libraries automatically, e.g., `error: libhdf5.so: cannot open shared object file: No such file or directory`. You will have to specify the path to the HDF5 installation, i.e., `HDF5_DIR=/path/to/hdf5/lib` (see [here](https://docs.h5py.org/en/stable/build.html#building-against-parallel-hdf5) for more details).

For our COSMA7 setup, that will be:

`HDF5DIR="/cosma/local/parallel-hdf5//gnu_14.1.0_ompi_5.0.3/1.14.4/"`

## Usage

Once installed the `zoom-mask-creator` command will be available, which expects one argument, a parameter file, e.g.,

* `zoom-mask-creator ./examples/Eagle100_Group100.toml`

or in MPI:

* `mpirun -np XX zoom-mask-creator ./examples/Eagle100_Group100.toml`

When you run the `zoom-mask-creator` code, both a mask (in the form of a HDF5 file) and a plot of the mask get deposited into the `output_dir` directory. 

### Parameter file

All the parameters of the run are stored in a single TOML file. Here's an example of the parameter file structure:

```toml
[mask]
mask_cell_size = 5
min_num_per_cell = 3
topology_fill_holes = true
topology_dilation_niter = 0
topology_closing_niter = 0

[region]
coords = [250, 250, 250]
radius = 25.0
shape = "sphere"

[snapshot]
paths = [
    "/manticore_dir/2MPP_INNER_N128_DES_V2/R512/mcmc_0/swift_monofonic/snap_0001/snap_0001.0.hdf5",
]
data_type = "swift"

[ics]
ic_type = "map_to_ics"

[map_to_ics]
paths = [
    "/manticore_dir/2MPP_INNER_N128_DES_V2/R512/mcmc_0/monofonic/ics_swift.hdf5",
]

[output]
path = "masks/2MPP_INNER_N128_DES_V2/R512/RADIUS25MPC"
```

### Parameter descriptions

| Section | Parameter | Description | Allowed values | Default |
| --- | --- | --- | --- | --- |
| **region** | `shape` | Shape of the region to extract | "sphere", "cuboid", "slab" | "sphere" |
| | `coords` | 3D coordinates of the center of the region | [x, y, z] | |
| | `radius` | Radius of region (for sphere shape) | > 0 | |
| | `dim` | Dimensions (for cuboid or slab shape) | [x, y, z] | |
| **snapshot** | `paths` | Path(s) to target snapshot file(s) | List of strings | |
| | `data_type` | Type of snapshot | "swift" or "eagle" | |
| **ics** | `ic_type` | Method to map particles to initial conditions | "use_peano_ids" or "map_to_ics" | |
| **peano** | `bits` | Number of bits used in Peano Hilbert indexing | Integer | |
| | `divide_ids_by_two` | When particle IDs have been multiplied by 2 | Boolean | |
| **map_to_ics** | `paths` | Path(s) to IC file(s) for particle ID matching | List of strings | |
| **mask** | `mask_cell_size` | Mask cell size in simulation units | > 0 | 3.0 |
| | `min_num_per_cell` | Minimum number of particles per cell to include | ≥ 0 | 3 |
| | `topology_fill_holes` | Fill holes in the generated mask | Boolean | true |
| | `topology_dilation_niter` | Number of dilation iterations for buffering | ≥ 0 | 0 |
| | `topology_closing_niter` | Number of closing iterations | ≥ 0 | 0 |
| **output** | `path` | Output directory for mask and plot | String | |

### Method selection

The code now supports two methods for finding Lagrangian positions:

1. Using Peano Hilbert indexing (`ic_type = "use_peano_ids"`):
   - Requires the `bits` and `divide_ids_by_two` parameters in the `peano` section
   - Suitable for simulations where particle IDs are based on Peano-Hilbert space-filling curves

2. Matching particle IDs between snapshot and IC file (`ic_type = "map_to_ics"`):
   - Requires the `paths` parameter in the `map_to_ics` section pointing to the IC file(s)
   - Suitable when you have access to the original IC files and want to directly match particles

### An example: R=35 Mpc from Manticore-Mini

![image](https://github.com/user-attachments/assets/fc967109-83a6-4f64-8e69-9129b6f4e4d3)


The parameter file above (also in `./examples/` directory) generates a mask of the Lagrangian region created by the dark matter particles within a radius of `R <= 30 Mpc` of the centre of the Manticore-Mini parent volume.

When generating masks it is always a bit of trial and error, you may have to try a few different configurations and check the output mask looks reasonable (shown for this example above). The plot shows, in blue, (a subset of) the Lagrangian positions of the selected dark matter particles, the red squares are the cell positions of the generated mask, the white circle is the original target selection region at *z*=0 (remember everything shown is co-moving), and the red lines outline the minimum (and minimum symmetric) bounding boxes of the Lagrangian region.
