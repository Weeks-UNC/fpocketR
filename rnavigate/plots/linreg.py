"""Plot linear regression scatter plots."""

from scipy import stats
import seaborn as sns
import numpy as np
from rnavigate import plots


class LinReg(plots.Plot):
    """Plot a linear regression scatter plot, pairwise, between profiles.

    Parameters
    ----------
    num_samples : int
        Number of samples to plot.
    scale : str, optional
        Scale of the plot, either 'linear' or 'log'. The default is 'linear'.
    regression : 'pearson' or 'spearman', Defaults to 'pearson'
        Type of regression to perform.
    kde : bool, optional
        Whether to plot a kernel density estimate instead of a scatter plot.
        The default is False.
    region : 'all' or tuple of 2 integers, defaults to 'all'
        Start and end positions of the region to plot. If 'all', plot the
        entire profile.

    Attributes
    ----------
    fig : matplotlib.figure.Figure
        Figure of the plot.
    axes : numpy.ndarray of matplotlib.axes.Axes
        Axes of the plot.
    length : int
        Number of samples to plot.
    lims : list of 2 floats
        Limits of the plot applied to all x and y axes.
    profiles : list of numpy.ndarray
        Each sample's per-nucleotide values over the region and column of interest.
    colors : list of numpy.ndarray
        A color for each nucleotide in the region applied to the scatter plot.
    labels : list of str
        A label for each sample.
    scale : 'linear' or 'log'
        Scale of the plot axes.
    kde : bool
        Whether to plot a kernel density estimate instead of a scatter plot.
    regression : function
        Regression function to use.
    region : tuple of 2 integers
        Start and end positions of the region to plot.
    """

    def __init__(
        self, num_samples, scale="linear", regression="pearson", kde=False, region="all"
    ):
        """Initialize the linear regression plot."""
        super().__init__(num_samples, rows=num_samples - 1, cols=num_samples - 1)
        self.region = region
        linreg_axes = []
        for row in range(num_samples - 1):
            for col in range(row, num_samples - 1):
                linreg_axes.append(self.axes[row, col])
        self.axes[0, 0].get_shared_x_axes().join(*linreg_axes)
        self.axes[0, 0].get_shared_y_axes().join(*linreg_axes)
        self.lims = [0, 0]
        self.scale = scale
        self.kde = kde
        if regression == "pearson":
            self.regression = stats.pearsonr
        elif regression == "spearman":
            self.regression = stats.spearmanr
        self.profiles = []
        self.colors = []
        self.labels = []

    def set_figure_size(
        self,
        height_ax_rel=None,
        width_ax_rel=None,
        width_ax_in=2,
        height_ax_in=2,
        height_gap_in=0.3,
        width_gap_in=0.3,
        top_in=1,
        bottom_in=0.5,
        left_in=0.5,
        right_in=0.5,
    ):
        """Set the figure size.

        Parameters
        ----------
        height_ax_rel : float, defaults to None
            Height of the axes relative to the y-axis limits.
        width_ax_rel : float, defaults to None
            Width of the axes relative to the x-axis limits.
        width_ax_in : float, defaults to 2
            Width of the axes in inches.
        height_ax_in : float, defaults to 2
            Height of the axes in inches.
        height_gap_in : float, defaults to 0.3
            Height of the gap between axes in inches.
        width_gap_in : float, defaults to 0.3
            Width of the gap between axes in inches.
        top_in : float, defaults to 1
            Top margin of the figure in inches.
        bottom_in : float, defaults to 0.5
            Bottom margin of the figure in inches.
        left_in : float, defaults to 0.5
            Left margin of the figure in inches.
        right_in : float, defaults to 0.5
            Right margin of the figure in inches.
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
        self, structure, profile, annotations, label, column=None, colors="sequence"
    ):
        """Add profile data. If all samples have been added, plot the regression.

        Parameters
        ----------
        structure : Structure
            Structure object.
        profile : Profile
            Profile object.
        annotations : Annotations
            Annotations object.
        label : str
            Label for the sample.
        column : str, optional
            Column of the profile to plot. The default is None.
        colors : str, optional
            Color scheme to use. The default is 'sequence'.
        """
        if self.region == "all":
            start, end = 1, profile.length
        else:
            start, end = self.region
        if column is None:
            column = profile.metric
        self.labels.append(label)
        colors, colormap = profile.get_colors(
            colors,
            profile=profile,
            structure=structure,
            annotations=annotations,
        )
        self.colors.append(colors[start - 1 : end])
        self.add_colorbar_args(colormap)
        self.profiles.append(profile.data[column].to_numpy(copy=True)[start - 1 : end])
        if len(self.profiles) == self.length:
            for row in range(self.length - 1):
                for col in range(self.length - 1):
                    if row <= col:
                        self.plot_regression(i=row, j=col)
            self.finalize()

    def finalize(self):
        """Finalize the plot by formatting the axes and adding labels."""
        # change the plotting buffer for linear scale
        if self.scale == "log" and self.kde:
            # change the tick labels to 'fake' a log scale plot
            ticks = []
            for power in range(-5, 5):
                power = 10**power
                if self.lims[0] > power or power > self.lims[1]:
                    continue
                ticks.append(power)
                # 10^n formatted tick labels
            for i in range(len(self.axes)):
                self.axes[i, i].set(xticks=ticks, yticks=ticks)
        if self.scale == "linear":
            buffer = 0.05 * (self.lims[1] - self.lims[0])
            self.lims[0] -= buffer
            self.lims[1] += buffer
            self.axes[0, 0].set(ylim=self.lims, xlim=self.lims)
        # format axis spines and labels
        ticks = self.axes[0, 0].get_yticks()[1:-1]
        for row in range(self.length - 1):
            for col in range(self.length - 1):
                ax = self.axes[row, col]
                if row <= col:
                    sns.despine(ax=ax)
                    ax.set_xticks(ticks)
                    ax.set_yticks(ticks)
                else:
                    self.axes[row, col].set_axis_off()
                if row < col:
                    ax.set_xticklabels([" "] * len(ticks))
                    ax.set_yticklabels([" "] * len(ticks))
                if row == 0:
                    ax.annotate(
                        self.labels[col + 1],
                        xy=(0.5, 1),
                        xytext=(0, 5),
                        xycoords="axes fraction",
                        textcoords="offset points",
                        size="large",
                        ha="center",
                        va="baseline",
                    )
                if col == (self.length - 2):
                    ax.annotate(
                        self.labels[row],
                        xy=(1, 0.5),
                        xytext=(5, 0),
                        rotation=-90,
                        xycoords="axes fraction",
                        textcoords="offset points",
                        size="large",
                        ha="left",
                        va="center",
                    )

    def plot_regression(self, i, j):
        """Plot a linear regression between two profiles.

        Profiles must already have been added using `plot_data`.

        Parameters
        ----------
        i : int
            Index of the first profile.
        j : int
            Index of the second profile.
        """
        ax = self.axes[i, j]
        p1 = self.profiles[j + 1]
        p2 = self.profiles[i]
        colors = self.colors[i]

        # ax.plot([0, 1], [0, 1], color='black')
        if self.scale == "log":
            notNans = (p1 > 0) & (p2 > 0)
            p1 = p1[notNans]
            p2 = p2[notNans]
            colors = colors[notNans]
            r_value, _ = self.regression(np.log10(p2), np.log10(p1))
        if self.scale == "linear":
            notNans = ~np.isnan(p1) & ~np.isnan(p2)
            p1 = p1[notNans]
            p2 = p2[notNans]
            colors = colors[notNans]
            r_value, _ = self.regression(p2, p1)
        self.lims[0] = min([self.lims[0], min(p1), min(p2)])
        self.lims[1] = max([self.lims[1], max(p1), max(p2)])
        ax.text(
            0.1,
            0.95,
            f"r = {r_value:.2f}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            bbox=dict(fc="white", alpha=0.5, ec="black"),
        )
        if self.kde:
            if self.scale == "log":
                sns.kdeplot(
                    ax=ax,
                    x=p1,
                    y=p2,
                    fill=True,
                    log_scale=True,
                    levels=np.arange(1, 11) / 10,
                )
            elif self.scale == "linear":
                sns.kdeplot(ax=ax, x=p1, y=p2, fill=True, levels=np.arange(1, 11) / 10)
        else:
            ax.scatter(p2, p1, c=colors, marker=".")
            ax.set(xscale=self.scale, yscale=self.scale)
