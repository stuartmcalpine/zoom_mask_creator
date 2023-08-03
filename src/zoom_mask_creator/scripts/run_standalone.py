from zoom_mask_creator.mask import Mask
from zoom_mask_creator.plot import plot_mask
from zoom_mask_creator.io import save_mask
from mpi4py import MPI

import sys

def main():
    # Generate mask using param file.
    comm = MPI.COMM_WORLD
    x = Mask(sys.argv[1], comm=comm)
    
    # Plot the mask.
    plot_mask(x)
    
    # Save the mask.
    save_mask(x)
