"""Plot nucleotide distribution of a profile."""

import numpy as np
import seaborn as sns
from rnavigate import plots, styles


class NucleotideDistribution(plots.Plot):
    """Plot nucleotide distribution of a profile.

    Parameters
    ----------
    num_samples : int
        Number of samples to plot.
    sharex : bool, optional
        Whether to share the x-axis between plots.
    cols : int, optional
        Number of columns in the plot.
    **plot_kwargs
        Keyword arguments passed to the plot function.

    Attributes
    ----------
    fig : matplotlib.figure.Figure
        Figure object.
    axes : numpy.ndarray of matplotlib.axes.Axes
        Array of axes objects.
    i : int
        Index of the current plot.
    """

    def __init__(self, num_samples, **plot_kwargs):
        """Initialize the plot."""
        super().__init__(num_samples, sharex=True, cols=1, **plot_kwargs)

    def set_figure_size(
        self,
        height_ax_rel=None,
        width_ax_rel=None,
        width_ax_in=2,
        height_ax_in=2,
        height_gap_in=0.2,
        width_gap_in=0.4,
        top_in=1,
        bottom_in=1,
        left_in=1,
        right_in=1,
    ):
        """Set the figure size.

        Parameters
        ----------
        height_ax_rel : float, optional
            Height of the axes relative to the y-axis limits.
        width_ax_rel : float, optional
            Width of the axes relative to the x-axis limits.
        width_ax_in : float, defaults to 2
            Width of the axes in inches.
        height_ax_in : float, defaults to 2
            Height of the axes in inches.
        height_gap_in : float, defaults to 0.2
            Height of the gap between axes in inches.
        width_gap_in : float, defaults to 0.4
            Width of the gap between axes in inches.
        top_in : float, defaults to 1
            Top margin in inches.
        bottom_in : float, defaults to 1
            Bottom margin in inches.
        left_in : float, defaults to 1
            Left margin in inches.
        right_in : float, defaults to 1
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

    def plot_data(self, profile, label, column=None, normalize=None, ax=None):
        """Plot data to the current (or specified) axes.

        Parameters
        ----------
        profile : rnavigate.profile.Profile
            Profile object.
        label : str
            label for the y-axis.
        column : str, optional
            Column of `profile.data` to plot.
        normalize : dict, optional
            Keyword arguments passed to `profile.normalize`.
        ax : matplotlib.axes.Axes, optional
            Axes object to plot to. If not specified, the current axes is used.
        """
        profile = profile.copy()
        if ax is None:
            ax = self.get_ax()
            self.i += 1
        if column is None:
            column = profile.metric
        if normalize is not None:
            profile.normalize(**normalize)
        data = profile.data
        nt_idx = data[column] > 0
        sns.kdeplot(
            ax=ax,
            data=data.loc[nt_idx],
            x=np.log10(data.loc[nt_idx, column]),
            hue="Sequence",
            hue_order=["A", "U", "C", "G"],
            palette={nt: styles.get_nt_color(nt) for nt in "AUGC"},
            common_norm=False,
        )
        ax.set(
            xlim=(-2.5, 1.5),
            xticks=[-2, -1, 0, 1],
            xticklabels=[0.01, 0.1, 1, 10],
            xlabel="Normalized profile (log scale)",
            yticks=[],
            ylabel=label,
        )
        ax.axvspan(-1, 0, color="lightgrey", alpha=0.4, ec="none")
