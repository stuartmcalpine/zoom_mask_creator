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
[read_swift](https://github.com/stuartmcalpine/read_swift). Crucially, the
ParticleIDs in the parent simulation must have used *Peano Hilbert* indexing, which we
use to find the Lagrangian positions of the particles.

This script partners with 
[zoom_particle_load_creator](https://github.com/stuartmcalpine/zoom_particle_load_creator),
which creates particle loads using masks generated from these scripts.

## Installation

### Requirements

* `OpenMPI` or other MPI library
* `python>=3.8,<3.11`

Recommended modules when working on COSMA7:

* `module load gnu_comp/11.1.0 openmpi/4.1.4 parallel_hdf5/1.12.0 python/3.9.1-C7`

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

`HDF5DIR="/cosma/local/parallel-hdf5//gnu_11.1.0_ompi_4.1.4/1.12.0/"`

## Usage

Once installed the `zoom-mask-creator` command will be available, which expects one argument, a parameter file, e.g.,

* `zoom-mask-creator ./examples/Eagle100_Group100.yml`

or in MPI:

* `mpirun -np XX zoom-mask-creator ./examples/Eagle100_Group100.yml`

When you run the ``zoom-mask-creator` code, both a mask (in the form of a HDF5 file) and a plot of the mask get deposited into the `output_dir` directory. 

### Parameter file

All the parameters of the run are stored in a single YAML file, see `./examples/Eagle100_Group100.yml` as an example.

| Required parameters | Description | Allowed values |
| --- | ----------- | ------------|
| `shape` | The shape of the region to extract from the target snapshot | sphere |
| `snap_file` | Full path to target snapshot (or snapshot part) file | |
| `bits` | Number of bits used in the Peano Hilbert indexing | |
| `fname` | Output filename for generated mask | |
| `data_type` | Type of snapshot | swift or eagle |
| `divide_ids_by_two` | When particles IDs have been multiplied for 2 (as was done for EAGLE gas simulations) ||
| `output_dir` | Output directory to store mask ||
| `input_centre` | 3D coordinates of the center of the region we are extracting from the target snapshot | |

| Optional parameters | Description | Allowed values | Default |
| --- | ----------- | ------------| ------- |
| `radius` | Used with `shape`=sphere. Radius of region to extract. | >0 | |
| `mask_cell_size` | Mask cell size in simulation units (i.e., resolution of mask) | >0 | 3.0 |
| `min_num_per_cell` | Minimum number of dark matter particles in cell to be considered in mask | >= 0 | 3|
| `select_from_vr` | Not yet implemented | True/False | False |
| `topology_fill_holes` | Attempt to automatically fill in holes in the generated mask (see `scipy.ndimage.binary_fill_holes`) | True/False | True |
| `topology_dilation_niter` | Uses `scipy.ndimage.binary_dilation` to "buffer" the mask. I.e., `topology_dilation_niter=1` would generate a skin approximately one layer thick around the raw mask. Experiment with this if you want a bit of safety padding around your masked region.  | >= 0| 0 |
| `topology_closing_niter` | Basically the opposite of `topology_dialation_niter`, erodes away layers around the raw mask (see `scipy.ndimage.binary_closing`) | >= 0| 0 |


### An example: Group 100 from the Eagle 100 Mpc cosmological simulation

<figure>
    <img src="/docs/Eagle100_Group100.png"
         alt="Eagle100_Group100">
</figure>

The `Eagle100_Group100.yml` parameter file in the `./examples/` directory generates a mask of the Lagrangian region created by the dark matter particles within a radius of `R <= 2xR200_crit` of the 100th most massive halo in the Eagle 100 Mpc cosmological simulation. Note you will need access to the EAGLE data on COSMA to run this example.

When generating masks it is always a bit of trial and error, you may have to try a few different configurations and check the output mask looks reasonable (shown for this example above). The plot shows, in blue, (a subset of) the Lagrangian positions of the selected dark matter particles, the red squares are the cell positions of the generated mask, the white circle is the original target selection region at *z*=0 (remember everything shown is co-moving), and the red lines outline the minimum (and minimum symmetric) bounding boxes of the Lagrangian region.

What we considered generating this mask:

* `mask_cell_size=0.5` so we have a reasonable resolution for a mask covering this volume (a value too large and the mask is inefficient, too small and the mask is expensive and could be full of holes).

* `radius=0.7305606` is two times `R200_crit` of the desired halo (can find these details in the EAGLE database).

* `bits=14` was the choice of PH indexing when running EAGLE.

* `divide_ids_by_two=True` is needed because we are using the EAGLE hydro simulation (you would not need this if it is running from the DMO simulation).

* `topology_dilation_niter=0` as we wished a tight mask that only just selects the desired Lagrangian region and no more.

* `min_num_per_cell=100` to avoid including the spurious straggler particles separated from the main continuous region.
