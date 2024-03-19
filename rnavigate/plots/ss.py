from rnavigate import plots
import numpy as np


class SS(plots.Plot):
    """Plot secondary structure diagrams.

    Parameters
    ----------
    num_samples : int
        Number of samples to plot.
    **kwargs
        Keyword arguments passed to `rnavigate.plots.Plot`.

    Attributes
    ----------
    fig : matplotlib.figure.Figure
        Figure object.
    axes : numpy.ndarray of matplotlib.axes.Axes
        Array of axes objects.
    xlims : list of float
        x limits of the plot.
    ylims : list of float
        y limits of the plot.
    i : int
        Index of the current plot.
    """

    def __init__(self, num_samples, **kwargs):
        """Initialize the plot."""
        super().__init__(num_samples, sharey=True, sharex=True, **kwargs)
        for i in range(self.length):
            ax = self.get_ax(i)
            ax.axis("off")
        self.xlims = [-2, 2]
        self.ylims = [-2, 2]

    def set_figure_size(
        self,
        height_ax_rel=0.2,
        width_ax_rel=0.2,
        width_ax_in=None,
        height_ax_in=None,
        height_gap_in=0.5,
        width_gap_in=0.2,
        top_in=1,
        bottom_in=0.5,
        left_in=0.5,
        right_in=0.5,
    ):
        """Set the figure size.

        Parameters
        ----------
        height_ax_rel : float, defaults to 0.2
            Height of the axes relative to the y-axis limits.
        width_ax_rel : float, defaults to 0.2
            Width of the axes relative to the x-axis limits.
        width_ax_in : float
            Width of the axes in inches.
        height_ax_in : float
            Height of the axes in inches.
        height_gap_in : float, defaults to 0.5
            Height of the gap between axes in inches.
        width_gap_in : float, defaults to 0.2
            Width of the gap between axes in inches.
        top_in : float, defaults to 1
            Height of the top margin in inches.
        bottom_in : float, defaults to 0.5
            Height of the bottom margin in inches.
        left_in : float, defaults to 0.5
            Width of the left margin in inches.
        right_in : float, defaults to 0.5
            Width of the right margin in inches.
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
        structure,
        interactions=None,
        interactions2=None,
        profile=None,
        annotations=None,
        label="",
        colors=None,
        nt_ticks=None,
        bp_style="dotted",
    ):
        """Plot the data on the current axes.

        Parameters
        ----------
        structure : rnavigate.data.SecondaryStructure
            Structure object with diagram drawing coordinates.
        interactions : rnavigate.data.Interactions, optional
            Interactions object to plot as lines between nucleotides.
        interactions2 : rnavigate.data.Interactions, optional
            Interactions object to plot as lines between nucleotides.
        profile : rnavigate.data.Profile, optional
            Profile object used to color nucleotides.
        annotations : list of rnavigate.data.Annotation, optional
            Annotation objects to highlight regions or nucleotides of interest.
        label : str, optional
            Label for the plot title.
        colors : dict, optional
            Dictionary of colors for each plot element. Keys are "sequence",
            "nucleotides", "structure", and "basepairs". Values are either
            matplotlib colors or strings specifying the color scheme.
        nt_ticks : int, optional
            Number of nucleotides between tick marks.
        bp_style : "dotted", "solid", or "conventional", defaults to "dotted"
            Style of base pair lines.
        """
        if annotations is None:
            annotations = []
        ax = self.get_ax()
        if colors is None:
            colors = {}
        colors = {
            "sequence": None,
            "nucleotides": "sequence",
            "structure": "grey",
            "basepairs": "grey",
        } | colors
        for key in ["nucleotides", "structure", "basepairs", "sequence"]:
            if isinstance(colors[key], str) and colors[key] == "contrast":
                continue
            elif colors[key] is None:
                continue
            colors[key], colormap = structure.get_colors(
                colors[key],
                profile=profile,
                structure=structure,
                annotations=annotations,
            )
            self.add_colorbar_args(colormap)
        if isinstance(colors["sequence"], np.ndarray):
            colors["nucleotides"], _ = structure.get_colors("white")
        elif colors["sequence"] == "contrast":
            colors["sequence"] = plots.get_contrasting_colors(colors["nucleotides"])
        if colors["sequence"] is not None:
            plots.plot_sequence_ss(ax, structure, colors["sequence"])
        if colors["nucleotides"] is not None:
            plots.plot_nucleotides_ss(ax, structure, colors["nucleotides"])
        if colors["structure"] is not None:
            plots.plot_structure_ss(ax, structure, colors["structure"])
        if colors["basepairs"] is not None:
            plots.plot_basepairs_ss(ax, structure, bp_style)
        if nt_ticks is not None:
            plots.plot_positions_ss(ax, structure, nt_ticks)
        if interactions is not None:
            plots.plot_interactions_ss(ax, structure, interactions)
            self.add_colorbar_args(interactions.cmap)
        if interactions2 is not None:
            plots.plot_interactions_ss(ax, structure, interactions2)
            self.add_colorbar_args(interactions2.cmap)
        for annotation in annotations:
            plots.plot_annotation_ss(ax, structure, annotation)
        ax.set_title(label)
        x = structure.xcoordinates
        y = structure.ycoordinates
        xlims = ax.get_xlim()
        ylims = ax.get_ylim()
        xlims = [min([xlims[0], min(x) - 2]), max([xlims[1], max(x) + 2])]
        ylims = [min([ylims[0], min(y) - 2]), max([ylims[1], max(y) + 2])]
        ax.set(ylim=ylims, xlim=xlims)
        self.i += 1
