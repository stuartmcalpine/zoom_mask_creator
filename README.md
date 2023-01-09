[![python](https://img.shields.io/badge/Python-3.8-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org)
[![python](https://img.shields.io/badge/Python-3.9-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org)
[![python](https://img.shields.io/badge/Python-3.10-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## Zoom-in simulations mask creator

A simple script to find, and mask, the Lagrangian positions (i.e., positions in the
ICs) of a set of particles read from a [Virgo
Consortium](https://virgo.dur.ac.uk/)-like cosmological simulation. 

The code is compatible with simulation outputs that can be read by
[pyread_eagle](https://github.com/kyleaoman/pyread_eagle) or
[read_swift](https://github.com/stuartmcalpine/read_swift). Note we use the information from *Peano Hilbert* indexing in the ParticleIDs to find the Lagrangian positions, therefore the simulation must have used this indexing scheme.

The figure below shows an example. We want to know the Lagrangian positions of
all dark matter particles from the
*Sibelius-DARK* simulation that are within a radius of 5 Mpc from the Milky Way at redshift=0. The blue points are the computed
Lagrangian positions of the dark matter particles, the red symbols highlight
the constructed mask that outlines the Lagrangian region the particles form, and the white circle
is shows the size of the initial selection region at redshift=0 (note everything is in co-moving coordinates).

The mask can then be used to generate a "Particle Load" (using
[zoom_particle_load_creator](https://github.com/stuartmcalpine/zoom_particle_load_creator))
in order to perform a resimulation of the region.

<figure>
    <img src="/examples/Sibelius_5Mpc.png"
         alt="Sibelius 5Mpc region">
</figure>

## Installation

### Requirements

* `OpenMPI` or other MPI library
* `python>=3.8,<3.11`

Recommended modules when working on COSMA7:

* `module load gnu_comp/11.1.0 openmpi/4.1.4 parallel_hdf5/1.12.0 python/3.9.1-C7`

### Installation from source

It is recommended you install the package within a virtual/conda environment.
Or alternatively, if you are installing on a shared access machine, to your
local profile by including the `--user` flag during the `pip` installation.

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

* `MPICC=mpicc CC=mpicc HDF5_MPI="ON" python -m pip install --no-binary=mpi4py,h5py mpi4py h5py`

If `pip` can't find your `HDF5` libraries automatically, e.g., `error: libhdf5.so: cannot open shared object file: No such file or directory`. You will have to specify the path to the HDF5 installation, i.e.,

`HDF5_DIR=/path/to/hdf5/lib`

see [here](https://docs.h5py.org/en/stable/build.html#building-against-parallel-hdf5) for more details.

For our Cosma7 setup, that will be:

`HDF5DIR="/cosma/local/parallel-hdf5//gnu_11.1.0_ompi_4.1.4/1.12.0/"`

## Usage

Once installed the `zoom-mask-creator` command will be available. The command expects one argument, a parameter file, e.g.,

* `zoom-mask-creator ./examples/Sibelius_5Mpc.yml`

or in MPI:

* `mpirun -np XX zoom-mask-creator ./examples/Sibelius_5Mpc.yml`

### Expected entries in parameter file

All the parameters of the run are stored in a single YAML file, see `./examples/Sibelius_5Mpc.yml` as an example.

| Parameter | Description | Allowed values |
| --- | ----------- | ------------|
| `shape` | The shape of the region to extract from the target snapshot | sphere |
| `snap_file` | Full path to target snapshot (or snapshot part) file | |
| `bits` | Number of bits used in the Peano Hilbert indexing | |
| `fname` | Output filename for generated mask | |
| `input_centre` | 3D coordinates of the centre of the region we are extracting from the target snapshot | |
| `radius` | Used with shape=sphere. Radius of region to extract | |
| `data_type` | Type of snapshot | swift or eagle |
| `divide_ids_by_two` | When particles IDs have been multiplied for 2 (as was done for EAGLE gas simulations) ||
| `mask_cell_size` | Mask cell size in simulation units. Defines resolution of mask ||
| `min_num_per_cell` | Minimum number of dark matter particles in cell to be considered in mask ||
| `select_from_vr` | Ignore ||
| `output_dir` | Output directory to store mask ||
| `topology_dilation_niter` | Leave as 1 ||
