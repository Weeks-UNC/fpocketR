"""Heatmap plot class."""

from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import numpy as np
from rnavigate import plots
from rnavigate import data


class Heatmap(plots.Plot):
    """Create a heatmap plot.

    Parameters
    ----------
    num_samples : int
        Number of samples to plot.
    structure : rnavigate.data.SecondaryStructure or rnavigate.data.PDB
        Structure object.
    **kwargs
        Keyword arguments passed to :class:`rnavigate.plots.Plot`.

    Attributes
    ----------
    structure : rnavigate.data.SecondaryStructure or rnavigate.data.PDB
        Structure object.
    fig : matplotlib.figure.Figure
        Figure object.
    axes : numpy.ndarray of matplotlib.axes.Axes
        Array of axes objects.
    i : int
        Index of current axes object. Increments with each call to
        :meth:`rnavigate.plots.Heatmap.plot_data`.
    """

    def __init__(self, num_samples, structure, **plot_kwargs):
        """Initialize Heatmap object."""
        super().__init__(num_samples, sharey=True, **plot_kwargs)
        self.structure = structure
        for i in range(num_samples):
            ax = self.get_ax(i)
            ax.set_aspect("equal")
            mn_mx = (0.5, structure.length + 0.5)
            ax.set(
                ylim=mn_mx, xlim=mn_mx, xticks=list(range(50, structure.length + 1, 50))
            )

    def set_figure_size(
        self,
        height_ax_rel=None,
        width_ax_rel=None,
        width_ax_in=2,
        height_ax_in=2,
        height_gap_in=1,
        width_gap_in=0.5,
        top_in=1,
        bottom_in=0.5,
        left_in=0.5,
        right_in=0.5,
    ):
        """Set figure size.

        Parameters
        ----------
        height_ax_rel : float, defaults to None
            Height of each axes object relative to y-axis limits.
        width_ax_rel : float, defaults to None
            Width of each axes object relative to x-axis limits.
        width_ax_in : float, defaults to 2
            Width of each axes object in inches.
        height_ax_in : float, defaults to 2
            Height of each axes object in inches.
        height_gap_in : float, defaults to 1
            Height of gap between axes objects in inches.
        width_gap_in : float, defaults to 0.5
            Width of gap between axes objects in inches.
        top_in : float, defaults to 1
            Top margin in inches.
        bottom_in : float, defaults to 0.5
            Bottom margin in inches.
        left_in : float, defaults to 0.5
            Left margin in inches.
        right_in : float, defaults to 0.5
            Right margin in inches.
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
        self,
        interactions,
        label,
        levels=None,
        regions=None,
        interpolation=None,
        atom="O2'",
        plot_type="heatmap",
        weights=None,
    ):
        """Plot heatmap and contour data on the current (or specified) axes.

        Parameters
        ----------
        interactions : rnavigate.data.Interactions
            Interactions object containing inter-nucleotide data to plot as heatmap.
        label : str
            Label for plot title.
        levels : list of float, defaults to None
            Levels for contour lines. If None, levels are set automatically.
            [5] for SecondaryStructure, [20] for PDB.
        regions : list of tuple of int, defaults to None
            Regions to plot as contour lines. If None, distances are plotted.
        interpolation : str, defaults to None
            Interpolation method for heatmap. If None, no interpolation is used.
        atom : str, defaults to "O2'"
            Atom to use for calculating distances.
        plot_type : "heatmap" or "kde", defaults to "heatmap"
            Type of plot to create.
        weights : array-like, defaults to None
            A column name of `interactions` to use as weights for KDE plot.
        """
        ax = self.get_ax()
        if regions is not None:
            self.plot_contour_regions(ax, interactions, regions)
        else:
            self.plot_contour_distances(ax, levels, atom)
        assert plot_type in ["heatmap", "kde"], "plot_type must be heatmap or kde"
        if plot_type == "heatmap":
            self.plot_heatmap_data(ax, interactions, interpolation)
        elif plot_type == "kde":
            self.plot_kde_data(ax, interactions, weights=weights)
        ax.set_title(label)
        self.i += 1

    def plot_contour_regions(self, ax, interactions, regions):
        """Plot contour regions on the current (or specified) axes.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes object to plot on.
        interactions : rnavigate.data.Interactions
            Interactions object containing inter-nucleotide data to plot as heatmap.
        regions : list of tuple of int
            Regions to plot as contour lines.
        """
        matrix = np.full([interactions.length, interactions.length], 0)
        for (mn1, mx1), (mn2, mx2) in regions:
            matrix[mn1 - 1 : mx1 - 1, mn2 - 1 : mx2 - 1] += 1
        levels = [0.5]
        cmap = LinearSegmentedColormap.from_list("contours", ["black", "gray"])
        x_y = list(range(1, interactions.length + 1))
        ax.contour(x_y, x_y, matrix, levels=levels, cmap=cmap, linewidths=0.3)

    def plot_contour_distances(self, ax, levels, atom):
        """Plot contour distances on the current (or specified) axes.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes object to plot on.
        levels : list of float, defaults to None
            Levels for contour lines. If None, levels are set automatically.
            [5] for SecondaryStructure, [20] for PDB.
        atom : str, defaults to "O2'"
            Atom to use for calculating distances.
        """
        structure = self.structure
        distances = structure.get_distance_matrix(atom)
        for i in range(structure.length):
            for j in range(i, structure.length):
                distances[i, j] = 0
        if (levels is None) and isinstance(structure, data.SecondaryStructure):
            levels = [5]
        elif (levels is None) and isinstance(structure, data.PDB):
            levels = [20]
        cmap = LinearSegmentedColormap.from_list("contours", ["black", "gray"])
        x_y = list(range(1, structure.length + 1))
        ax.contour(x_y, x_y, distances, levels=levels, cmap=cmap, linewidths=0.3)

    def plot_heatmap_data(self, ax, interactions, interpolation):
        """Plot heatmap data on the specified axes.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes object to plot on.
        interactions : rnavigate.data.Interactions
            Interactions object containing inter-nucleotide data to plot as heatmap.
        interpolation : str, defaults to None
            Interpolation method for heatmap. If None, no interpolation is used.
        """
        data = interactions.get_sorted_data()
        metric = interactions.metric
        data_im = np.full([interactions.length] * 2, np.nan)
        window = interactions.window
        for _, row in data.iterrows():
            i = int(row["i"] - 1)
            j = int(row["j"] - 1)
            data_im[j : j + window, i : i + window] = row[metric]
        ax.imshow(
            data_im,
            cmap=interactions.cmap.cmap,
            norm=interactions.cmap.norm,
            interpolation=interpolation,
        )

    def plot_kde_data(self, ax, interactions, weights=None, **kwargs):
        """Plot KDE data on the specified axes.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes object to plot on.
        interactions : rnavigate.data.Interactions
            Interactions object containing inter-nucleotide data to plot as KDE.
        weights : array-like, defaults to None
            A column name of `interactions` to use as weights for KDE plot.
        **kwargs
            Keyword arguments passed to :func:`seaborn.kdeplot`.
        """
        data = interactions.get_sorted_data()
        sns.kdeplot(
            ax=ax,
            data=data,
            x="i",
            y="j",
            fill=True,
            levels=5,
            bw_adjust=0.2,
            cmap=interactions.cmap.cmap,
            common_norm=True,
            weights=weights,
            **kwargs
        )
        ax.set(xlabel="Nucleotide position", ylabel="Nucleotide position")
