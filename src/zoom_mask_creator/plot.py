import os

import matplotlib

matplotlib.use("Agg")

import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np


def plot_mask(mask, max_npart_per_rank=int(1e5)):
    """
    Make an overview plot of the zoom-in region.

    Note that this function must be called on all MPI ranks, even though
    only rank 0 generates the actual plot. The others are still required
    to access (a subset of) the particles stored on them.

    Parameters
    ----------
    mask : Mask object
        Stores the coordinates from the simulation and the computed mask
    max_npart_per_rank : int (optional)
        The max number of DM particles to plot (per rank)
    """
    axis_labels = ["x", "y", "z"]

    # Select a random sub-sample of particle coordinates on each rank and
    # combine them all on rank 0
    np_ic = mask.ic_coords.shape[0]
    n_sample = int(min(np_ic, max_npart_per_rank))
    indices = np.random.choice(np_ic, n_sample, replace=False)
    plot_coords = mask.ic_coords[indices, :]
    if mask.comm_size > 1:
        plot_coords = mask.comm.gather(plot_coords)

    # Only need rank 0 from here on, combine all particles there.
    if mask.comm_rank != 0:
        return
    plot_coords = np.vstack(plot_coords)

    # Extract frequently needed attributes for easier structure
    bound = mask.mask_extent
    cell_size = mask.cell_size

    fig, axarr = plt.subplots(1, 3, figsize=(13, 4))

    # Plot each projection (xy, xz, yz) in a separate panel. `xx` and `yy`
    # denote the coordinate plotted on the x and y axis, respectively.
    for ii, (xx, yy) in enumerate(zip([0, 0, 1], [1, 2, 2])):
        ax = axarr[ii]
        ax.set_aspect("equal")

        # Draw the outline of the cubic bounding region
        rect = patches.Rectangle(
            [-bound / 2.0, -bound / 2.0],
            bound,
            bound,
            linewidth=1,
            edgecolor="maroon",
            facecolor="none",
        )
        ax.add_patch(rect)

        # Draw on outline of cuboidal bounding region
        box_corners = [mask.mask_box[:, xx], mask.mask_box[:, yy]]
        ax.plot(
            box_corners[0][[0, 1, 1, 0, 0]],
            box_corners[1][[0, 0, 1, 1, 0]],
            color="maroon",
            linestyle="--",
            linewidth=0.7,
        )

        # Plot particles.
        ax.scatter(
            plot_coords[:, xx],
            plot_coords[:, yy],
            s=0.5,
            c="blue",
            zorder=-100,
            alpha=0.3,
        )

        ax.set_xlim(-bound / 2.0 * 1.05, bound / 2.0 * 1.05)
        ax.set_ylim(-bound / 2.0 * 1.05, bound / 2.0 * 1.05)

        # Plot (the centres of) selected mask cells.
        ax.scatter(
            mask.cell_coords[:, xx],
            mask.cell_coords[:, yy],
            marker="x",
            color="red",
            s=5,
            alpha=0.2,
        )

        # Plot cell outlines if there are not too many of them.
        if mask.cell_coords.shape[0] < 10000:
            for e_x, e_y in zip(mask.cell_coords[:, xx], mask.cell_coords[:, yy]):
                rect = patches.Rectangle(
                    (e_x - cell_size / 2, e_y - cell_size / 2),
                    cell_size,
                    cell_size,
                    linewidth=0.5,
                    edgecolor="r",
                    facecolor="none",
                    alpha=0.2,
                )
                ax.add_patch(rect)

        ax.set_xlabel(f"${axis_labels[xx]}$ [$h^{{-1}}$ Mpc]")
        ax.set_ylabel(f"${axis_labels[yy]}$ [$h^{{-1}}$ Mpc]")

        # Plot target high-resolution sphere (if that is our shape).
        if mask.params["region"]["shape"] == "sphere":
            circle = patches.Circle(
                (0, 0),
                mask.params["region"]["radius"],
                linewidth=2,
                edgecolor="black",
                facecolor="none",
                linestyle="--",
            )
            ax.add_patch(circle)
            text_where = mask.params["region"]["radius"]
        # Plot target resolution cuboid/slab
        else:
            dim = mask.params["region"]["dim"]
            square = patches.Rectangle(
                (-dim[xx] / 2.0, -dim[yy] / 2.0),
                dim[xx],
                dim[yy],
                linewidth=2,
                edgecolor="black",
                facecolor="none",
                linestyle="--",
            )
            ax.add_patch(square)
            text_where = dim[yy] / 2.0

        ax.text(
            0,
            text_where,
            f"z = ${mask.params['snapshot']['zred_snap']:.2f}$",
            color="grey",
            fontsize=6,
            va="bottom",
            ha="center",
            bbox={
                "facecolor": "white",
                "edgecolor": "black",
                "pad": 0.25,
                "boxstyle": "round",
                "linewidth": 0.3,
            },
        )

    # Save the plot
    plt.subplots_adjust(left=0.05, right=0.99, bottom=0.15, top=0.99)
    plotloc = os.path.join(mask.params["output"]["path"], "mask_plot") + ".png"
    plt.savefig(plotloc, dpi=200)
    plt.close()
