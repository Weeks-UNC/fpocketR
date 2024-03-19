"""Plotting class for sequence alignments."""

from rnavigate import plots


class Alignment(plots.Plot):
    """Class for plotting sequence alignments.

    Parameters
    ----------
    num_samples : int
        The number of samples to plot. Always 2.
    rows : int, optional
        The number of rows of plots. Always 1.
    cols : int, optional
        The number of columns of plots. Always 1.
    **kwargs : dict
        Additional keyword arguments to pass to plots.Plot.

    Attributes
    ----------
    region : tuple
        The region of the alignment to plot. Defaults to (1, len(alignment1))
    fig : matplotlib.figure.Figure
        The figure containing the plot.
    axes : numpy.array
        A 1x1 array of the axes containing the plot.
    """

    def __init__(self, num_samples, rows=None, cols=1, **kwargs):
        """Create a new Alignment object."""
        super().__init__(num_samples, rows, cols, **kwargs)

    def plot_data(self, alignment, label, ax=None):
        """Add the alignment to the next (or given) axes.

        Parameters
        ----------
        alignment : rnavigate.Alignment
            The alignment to plot.
        label : str
            The label for the alignment.
        ax : matplotlib.axes.Axes, optional
            The axes containing the plot.
            Defaults to None.
        """
        if ax is None:
            ax = self.get_ax()
        alignment1 = alignment.alignment1
        self.region = (1, len(alignment1))
        ax.set(xlim=(0, len(alignment1) + 1), ylim=(-5, 5), yticks=[])
        plots.plot_sequence_alignment(ax, alignment, label)
        for spine in ["top", "bottom", "left", "right"]:
            ax.spines[spine].set_color(None)

    def set_figure_size(
        self,
        fig=None,
        ax=None,
        rows=None,
        cols=None,
        height_ax_rel=0.03,
        width_ax_rel=0.03,
        width_ax_in=None,
        height_ax_in=None,
        height_gap_in=1,
        width_gap_in=0.5,
        top_in=1,
        bottom_in=0.5,
        left_in=0.5,
        right_in=0.5,
    ):
        """Set the figure size for the plot.

        Parameters
        ----------
        fig : matplotlib.figure.Figure, optional
            The figure containing the plot.
            Defaults to None.
        ax : matplotlib.axes.Axes, optional
            The axes containing the plot.
            Defaults to None.
        rows : int, optional
            The number of rows of plots.
            Defaults to None.
        cols : int, optional
            The number of columns of plots.
            Defaults to None.
        height_ax_rel : float, optional
            The relative height of each axes.
            Defaults to 0.03.
        width_ax_rel : float, optional
            The relative width of each axes.
            Defaults to 0.03.
        width_ax_in : float, optional
            The width of each axes in inches.
            Defaults to None.
        height_ax_in : float, optional
            The height of each axes in inches.
            Defaults to None.
        height_gap_in : float, optional
            The height of the gap between axes in inches.
            Defaults to 1.
        width_gap_in : float, optional
            The width of the gap between axes in inches.
            Defaults to 0.5.
        top_in : float, optional
            The top margin in inches.
            Defaults to 1.
        bottom_in : float, optional
            The bottom margin in inches.
            Defaults to 0.5.
        left_in : float, optional
            The left margin in inches.
            Defaults to 0.5.
        right_in : float, optional
            The right margin in inches.
            Defaults to 0.5.
        """
        super().set_figure_size(
            fig=fig,
            ax=ax,
            rows=rows,
            cols=cols,
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
