## Installation

### Requirements

Python >=3.8,<3.11

### How to install

It is recommended you install the package within a virtual/conda enviromnent.
Or alternativley if you are installing on a shared access machine, to your
local profile by including the `--user` flag during the `pip` installation.

First make sure your `pip` is up-to-date:

`python3 -m pip install --upgrade pip`

Then you can install the `zoom-mask-creator` package by typing the following in
the git directory: 

`python3 -m pip install .`

which will install `zoom-mask-creator` and any dependencies.

### MPI installation for `read_swift`

If you are using `read_swift` to load large snapshots over MPI collectivley
(i.e., multiple cores read in parallel from the same file), a bit of additional
setup is required.

## Usage
