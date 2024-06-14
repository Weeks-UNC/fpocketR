"""Plot distance histograms."""

import numpy as np
from rnavigate import plots


class DistHist(plots.Plot):
    """Create a distance histogram plot.

    Parameters
    ----------
    num_samples : int
        Number of samples to plot.
    **kwargs
        Keyword arguments passed to :class:`rnavigate.plots.Plot`.

    Attributes
    ----------
    fig : matplotlib.figure.Figure
        Figure object.
    axes : numpy.ndarray of matplotlib.axes.Axes
        Array of axes objects.
    axes2 : dict
        Dictionary of twin axes objects. Keys are axes objects from
        :attr:`rnavigate.plots.DistHist.axes`.
    i : int
        Index of current axes object. Increments with each call to
        :meth:`rnavigate.plots.DistHist.plot_data`.
    """

    def __init__(self, num_samples, **plot_kwargs):
        """Initialize DistHist object."""
        super().__init__(num_samples, sharey=True, **plot_kwargs)
        self.axes2 = {}
        base_ax2 = None
        for row in self.axes:
            for ax in row:
                ax2 = ax.twinx()
                self.axes2[ax] = ax2
                if base_ax2 is None:
                    base_ax2 = ax2
                else:
                    ax2.sharey(base_ax2)

    def set_figure_size(
        self,
        height_ax_rel=None,
        width_ax_rel=None,
        width_ax_in=2,
        height_ax_in=2,
        height_gap_in=1,
        width_gap_in=0.4,
        top_in=1,
        bottom_in=1,
        left_in=1,
        right_in=1,
    ):
        """Set figure size.

        Parameters
        ----------
        height_ax_rel : float, defaults to None
            Height of axes objects relative to y-axis limits.
        width_ax_rel : float, defaults to None
            Width of axes objects relative to x-axis limits.
        width_ax_in : float, defaults to 2
            Width of axes objects in inches.
        height_ax_in : float, defaults to 2
            Height of axes objects in inches.
        height_gap_in : float, defaults to 1
            Height of gap between axes objects in inches.
        width_gap_in : float, defaults to 0.4
            Width of gap between axes objects in inches.
        top_in : float, defaults to 1
            Height of top margin in inches.
        bottom_in : float, defaults to 1
            Height of bottom margin in inches.
        left_in : float, defaults to 1
            Width of left margin in inches.
        right_in : float, defaults to 1
            Width of right margin in inches.
        """
        super().set_figure_size(
            height_ax_rel=height_ax_rel,
            width_ax_rel=width_ax_rel,
            width_ax_in=width_ax_in,
            height_ax_in=height_ax_in,
            height_gap_in=height_gap_in,
            width_gap_in=width_gap_in,
            top_in=top_in,
            bottom_in=bottom_in,
            left_in=left_in,
            right_in=right_in,
        )

    def plot_data(
        self, structure, interactions, bg_interactions, label, atom="O2'", ax=None
    ):
        """Plot data on the current (or specified) axes object.

        Parameters
        ----------
        structure : rnavigate.data.SecondaryStructure or rnavigate.data.PDB
            Structure object to compute contact distances or 3D distances.
        interactions : rnavigate.data.Interactions
            Filtered Interactions object to to visualize pairwise distances.
        bg_interactions : rnavigate.data.Interactions
            Filtered Interactions object to to visualize pairwise background distances.
        label : str
            Label for the current axes object.
        atom : str, defaults to "O2'"
            Atom to compute distances from.
        ax : matplotlib.axes.Axes, defaults to None
            Axes object to plot on. If None, use the current axes object.
        """
        if ax is None:
            ax = self.get_ax()
        ax2 = self.axes2[ax]
        if bg_interactions is not None:
            self.plot_experimental_distances(
                ax=ax2,
                structure=structure,
                interactions=bg_interactions,
                atom=atom,
                histtype="step",
            )
        else:
            self.plot_structure_distances(ax=ax2, structure=structure, atom=atom)
        self.plot_experimental_distances(
            ax=ax, structure=structure, interactions=interactions, atom=atom
        )
        ax.set(title=label)
        self.i += 1
        if self.i == self.length:
            for leftmost_axis in self.axes[:, 0]:
                leftmost_axis.set(ylabel="Experimental")
            for rightmost_axis in self.axes[:, -1]:
                self.axes2[rightmost_axis].set(ylabel="Pairwise")
            for top_axis in self.axes[-1, :]:
                top_axis.set(xlabel="3D distance")
            for row in self.axes:
                for not_rightmost_axis in row[:-1]:
                    self.axes2[not_rightmost_axis].yaxis.set_tick_params(
                        labelright=False
                    )

    def plot_structure_distances(self, ax, structure, atom):
        """Plot all distances in the structure.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes object to plot on.
        structure : rnavigate.data.SecondaryStructure or rnavigate.data.PDB
            Structure object to compute contact distances or 3D distances.
        atom : str
            Atom to compute distances from.
        """
        matrix = structure.get_distance_matrix(atom=atom)
        dists = []
        for i in range(len(matrix) - 6):
            for j in range(i + 6, len(matrix)):
                if not np.isnan(matrix[i, j]):
                    dists.append(matrix.item(i, j))
        ax.hist(
            dists,
            bins=range(0, int(max(dists)) + 5, 5),
            histtype="step",
            color="0.5",
            label="All distances",
        )

    def plot_experimental_distances(
        self, ax, structure, interactions, atom, histtype="bar"
    ):
        """Plot pairwise distances from the interactions object.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes object to plot on.
        structure : rnavigate.data.SecondaryStructure or rnavigate.data.PDB
            Structure object to compute contact distances or 3D distances.
        interactions : rnavigate.data.Interactions
            Filtered Interactions object to to visualize pairwise distances.
        atom : str
            Atom to compute distances from.
        histtype : "bar" or "step", defaults to "bar"
            Type of histogram to plot.
        """
        interactions.set_3d_distances(structure, atom)
        ij_dists = interactions.data.loc[interactions.data["mask"], "Distance"]
        if (len(ij_dists) > 0) and (histtype == "bar"):
            ax.hist(
                ij_dists, bins=range(0, int(max(ij_dists)) + 5, 5), width=5, ec="none"
            )
        elif (len(ij_dists) > 0) and (histtype == "step"):
            ax.hist(
                ij_dists,
                bins=range(0, int(max(ij_dists)) + 5, 5),
                histtype="step",
                color="0.5",
                label="All distances",
            )
