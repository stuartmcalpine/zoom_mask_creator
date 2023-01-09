[![python](https://img.shields.io/badge/Python-3.8-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org)
[![python](https://img.shields.io/badge/Python-3.9-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org)
[![python](https://img.shields.io/badge/Python-3.10-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## Zoom-in simulations mask creator

A simple script to find, and mask, the Lagrangian positions (positions in the
ICs) of a set of particles read from a [Virgo
Consortium](https://virgo.dur.ac.uk/)-like cosmological simulation. 

The code is compatible with simulation outputs that can be read by
[pyread_eagle](https://github.com/kyleaoman/pyread_eagle) or
[read_swift](https://github.com/stuartmcalpine/read_swift).  

The figure below shows an example. We want to know the Lagrangian positions of
all dark matter particles within a radius of 5 Mpc from the Milky Way from the
*Sibelius-DARK* simulation at redshift=0. The blue points are the computed
Lagrangian positions of the dark matter particles, the red symbols highlight
the constructed mask that outlines the Lagrangian region, and the white circle
is shows the initial selection region at redshift=0.

The mask can then be used to generate a "Particle Load" (using
[zoom_particle_load_creator](https://github.com/stuartmcalpine/zoom_particle_load_creator))
in order to perform a resimulation of the region.

<figure>
    <img src="/examples/Sibelius_5Mpc.png"
         alt="Sibelius 5Mpc region">
</figure>

### Requirements

* `OpenMPI` or other MPI library
* Python >=`3.8`,<`3.11`

Recommended modules when working on COSMA7:

* `module load gnu_comp/11.1.0 openmpi/4.1.4-romio-lustre parallel_hdf5/1.12.2 python/3.9.1-C7`

### Installation

It is recommended you install the package within a virtual/conda environment.
Or alternatively, if you are installing on a shared access machine, to your
local profile by including the `--user` flag during the `pip` installation.

First make sure your `pip` is up-to-date:

* `python3 -m pip install --upgrade pip`

Then you can install the `zoom-mask-creator` package by typing the following in
the git directory: 

* `python3 -m pip install .`

which will install `zoom-mask-creator` and any dependencies.

### MPI installation for `read_swift`

If you are using `read_swift` to load large snapshots over MPI collectively
(i.e., multiple cores read in parallel from the same file), a bit of additional
setup is required.

Make sure you have MPI libraries installed on your machine (`OpenMPI` for example), and you have `hdf5` installed with **parallel** compatibility ([see here for details](https://docs.h5py.org/en/stable/mpi.html)).

First, uninstall any installed versions of `h5py`:

* `python3 -m pip uninstall h5py`

Then reinstall `h5py` from source with MPI flags:

* `CC=mpicc HDF5_MPI="ON" python -m pip install --no-binary=h5py h5py`

If `pip` can't find your `HDF5` libraries automatically, e.g., `error: libhdf5.so: cannot open shared object file: No such file or directory`. You will have to specify the path to the HDF5 installation, i.e.,

`HDF5_DIR=/path/to/hdf5/lib`

see [here](https://docs.h5py.org/en/stable/build.html#building-against-parallel-hdf5) for more details.

## Usage

Once installed the `zoom-mask-creator` command will be available. The command expects one argument, a parameter file, e.g.,

* `zoom-mask-creator ./examples/Sibelius_5Mpc.yml`

or in MPI:

* `mpirun -np XX zoom-mask-creator ./examples/Sibelius_5Mpc.yml`
