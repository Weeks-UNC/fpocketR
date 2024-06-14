from rnavigate import plots, styles
import seaborn as sns


class Skyline(plots.Plot):
    """Plot per-nucleotide measurements as stepped line graphs (skyline plots).

    Parameters
    ----------
    num_samples : int
        Number of samples to plot.
    nt_length : int
        Length of the nucleotide sequence.
    region : tuple of int, defaults to "all" (entire sequence)
        start and end position of the region to plot. If "all", plot the entire
        sequence.
    **kwargs
        Keyword arguments passed to `rnavigate.plots.Plot`.

    Attributes
    ----------
    nt_length : int
        Length of the nucleotide sequence.
    region : tuple of int
        start and end position of the region to plot.
    track_height : float
        Height of the tracks in the plot.
    fig : matplotlib.figure.Figure
        Figure object.
    ax : matplotlib.axes.Axes
        Axes object.
    axes : numpy.ndarray of matplotlib.axes.Axes
        Array of axes objects.
    i : int
        Index of the current plot.
    """

    def __init__(self, num_samples, nt_length, region="all", **kwargs):
        """Initialize the plot."""
        if region == "all":
            self.nt_length = nt_length
            self.region = (1, nt_length)
        else:
            self.nt_length = region[1] - region[0] + 1
            self.region = region
        super().__init__(num_samples=num_samples, **kwargs)
        self.ax = self.axes[0, 0]
        self.track_height = 0

    def set_figure_size(
        self,
        height_ax_rel=None,
        width_ax_rel=0.03,
        width_ax_in=None,
        height_ax_in=2,
        height_gap_in=1,
        width_gap_in=0.5,
        top_in=1,
        bottom_in=1,
        left_in=1,
        right_in=1,
    ):
        """Set the figure size.

        Parameters
        ----------
        height_ax_rel : float
            Height of the axes relative to the y-axis limits.
        width_ax_rel : float
            Width of the axes relative to the x-axis limits.
        width_ax_in : float
            Width of the axes in inches.
        height_ax_in : float
            Height of the axes in inches.
        height_gap_in : float
            Height of the gap between axes in inches.
        width_gap_in : float
            Width of the gap between axes in inches.
        top_in : float
            Height of the top margin in inches.
        bottom_in : float
            Height of the bottom margin in inches.
        left_in : float
            Width of the left margin in inches.
        right_in : float
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
            bottom_in=bottom_in + 2 * self.track_height,
            left_in=left_in,
            right_in=right_in,
        )

    def get_rows_columns(self, rows=None, cols=None):
        """Get the number of rows and columns.

        Parameters
        ----------
        rows : int
            Number of rows. This is ignored.
        cols : int
            Number of columns. This is ignored.

        Returns
        -------
        rows : int
            Number of rows. This is always 1.
        cols : int
            Number of columns. This is always 1.
        """
        return (1, 1)

    def get_ax(self, i=None):
        """Get the current axes object.

        Parameters
        ----------
        i : int, optional
            Index of the axes object. This is ignored.

        Returns
        -------
        ax : matplotlib.axes.Axes
            Axes object.
        """
        return self.ax

    def plot_data(
        self,
        profile,
        annotations=None,
        domains=None,
        label=None,
        columns=None,
        seqbar=True,
        errors=None,
        annotations_mode="track",
        nt_ticks=(20, 5),
    ):
        """Add data to the axes.

        Parameters
        ----------
        profile : rnavigate.profile.Profile
            Profile object.
        annotations : list of rnavigate.annotation.Annotation, optional
            List of annotation objects.
        domains : list of rnavigate.domain.Domain, optional
            List of domain objects.
        label : str, optional
            Sample name.
        columns : str or list of str, optional
            Which columns to plot. If None, plot the metric column.
        seqbar : bool, defaults to True
            Whether to plot a sequence bar.
        errors : str, optional
            Which error columns to plot. If None, do not plot errors.
        annotations_mode : "track" or "bar", defaults to "track"
            Whether to plot annotations as a track or as vertical bars.
        nt_ticks : tuple of int, defaults to (20, 5)
            Major and minor tick frequency for nucleotide positions.
        """
        if columns is None:
            columns = profile.metric
        track_unit = 0.03
        annotations_track = 2 * track_unit * len(annotations)
        annotations_track *= annotations_mode == "track"
        domains_track = 3 * track_unit * (domains is not None)
        sequence_track = track_unit * seqbar
        self.track_height = sequence_track + annotations_track + domains_track

        ax = self.get_ax()
        plots.plot_profile_skyline(ax, profile, label, columns, errors)
        if self.i == 0:
            if seqbar:
                plots.plot_sequence_track(
                    ax=ax,
                    sequence=profile.sequence,
                    yvalue=self.track_height * -1,
                    height=sequence_track,
                    region=self.region,
                    ytrans="axes",
                )
                self.add_colorbar_args(styles.get_nt_cmap())
            if domains is not None:
                for domain in domains:
                    plots.plot_domain_track(
                        ax=ax,
                        spans=domain,
                        yvalue=sequence_track - self.track_height,
                        height=domains_track,
                        region=self.region,
                        ytrans="axes",
                    )
            for i, annotation in enumerate(annotations):
                plots.plot_annotation_track(
                    ax=ax,
                    annotation=annotation,
                    mode=annotations_mode,
                    yvalue=(
                        sequence_track
                        + domains_track
                        - self.track_height
                        + (i + 0.5) * (2 * track_unit)
                    ),
                    height=2 * track_unit,
                    ytrans="axes",
                    region=self.region,
                )

        self.i += 1
        if self.i == self.length:
            if isinstance(columns, list):
                ylabel = [column.replace("_", " ") for column in columns]
                ylabel = ", ".join(ylabel)
            else:
                ylabel = columns.replace("_", " ")
            self.set_labels(ax=ax, ylabel=ylabel, axis_title=None, legend_title=None)
            self.set_axis(ax=ax, sequence=profile, nt_ticks=nt_ticks)

    def set_axis(self, ax, sequence, nt_ticks):
        """Set the axis limits and ticks.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes object.
        sequence : rnavigate.data.Sequence
            The sequence on which position labels are based. Dashes are ignored.
        nt_ticks : tuple of int
            Major and minor tick frequency for nucleotide positions.
        """
        xlim = self.region
        ax.set_xlim([xlim[0] - 0.5, xlim[1] + 0.5])
        plots.set_nt_ticks(
            ax=ax,
            sequence=sequence,
            region=self.region,
            major=nt_ticks[0],
            minor=nt_ticks[1],
        )
        ax.spines["bottom"].set(
            position=("axes", self.track_height * -1), visible=False
        )
        ax.spines["left"].set_visible(False)
        ax.grid(axis="y", visible=True)
        ax.tick_params(axis="y", length=0, grid_alpha=0.4)
        sns.despine(ax=ax, bottom=True, left=True)

    def set_labels(
        self,
        ax,
        axis_title="Raw Reactivity Profile",
        legend_title="Samples",
        xlabel="Nucleotide Position",
        ylabel="Profile",
    ):
        """Set the axis labels and legend.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes object.
        axis_title : str, optional
            Title of the axes.
        legend_title : str, optional
            Title of the legend.
        xlabel : str, optional
            Label of the x-axis.
        ylabel : str, optional
            Label of the y-axis.
        """
        ax.set_title(axis_title, loc="left")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(
            title=legend_title, labelcolor="linecolor", frameon=False, handlelength=0
        )
