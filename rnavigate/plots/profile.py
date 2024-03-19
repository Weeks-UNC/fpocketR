import seaborn as sns
from rnavigate import plots, styles


class Profile(plots.Plot):
    """Plot per-nucleotide measurements as colored bars.

    Parameters
    ----------
    num_samples : int
        Number of samples to plot.
    nt_length : int
        Length of the nucleotide sequence.
    region : tuple of int, optional
        Region of the nucleotide sequence to plot.
    **kwargs
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

    def __init__(self, num_samples, nt_length, region="all", **kwargs):
        """Initialize the plot."""
        if region == "all":
            self.nt_length = nt_length
            self.region = (1, nt_length)
        else:
            self.nt_length = region[1] - region[0] + 1
            self.region = region
        super().__init__(num_samples, sharey=True, **kwargs)
        self.track_height = 0

    def get_rows_columns(self, rows=None, cols=None):
        """Get the number of rows and columns in the plot.

        Parameters
        ----------
        rows : int, optional
            Number of rows in the plot.
        cols : int, optional
            Number of columns in the plot. This is ignored and set to 1.

        Returns
        -------
        rows : int
            Number of rows in the plot.
        cols : int
            Number of columns in the plot.
        """
        return super().get_rows_columns(rows=rows, cols=1)

    def set_labels(
        self,
        ax,
        axis_title="Reactivity Profile",
        xlabel="Nucleotide Position",
        ylabel="Reactivity",
    ):
        """Set the labels of the plot.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes object.
        axis_title : str, optional
            Title of the axis.
        xlabel : str, optional
            Label of the x-axis.
        ylabel : str, optional
            Label of the y-axis.
        """
        ax.set_title(axis_title, loc="left")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

    def plot_data(
        self,
        profile,
        annotations,
        domains,
        label,
        plot_error=True,
        column=None,
        seqbar=True,
        annotations_mode="track",
        nt_ticks=(20, 5),
    ):
        """Plot data to the current (or specified) axes.

        Parameters
        ----------
        profile : rnavigate.profile.Profile
            Profile object.
        annotations : list of rnavigate.data.Annotation
            List of annotation objects to plot along the sequence.
        domains : list of rnavigate.data.Annotation, optional
            List of domains to plot along the sequence.
        label : str
            label for the y-axis.
        plot_error : bool, optional
            Whether to plot the error bars.
        column : str, optional
            Column of `profile` to plot.
        seqbar : bool, optional
            Whether to plot the sequence track.
        annotations_mode : "track" or "bar", defaults to "track"
            Mode of the annotations track.
        nt_ticks : tuple of int, optional
            Major and minor tick interval for the nucleotide axis.
        """
        ax = self.get_ax()
        if column is not None:
            profile.metric = column
        column = profile.metric
        plots.plot_profile_bars(
            ax=ax, profile=profile, plot_error=plot_error, region=self.region
        )
        self.add_colorbar_args(profile.cmap)

        track_unit = 0.03
        annotations_track = 2 * track_unit * len(annotations)
        annotations_track *= annotations_mode == "track"
        domains_track = 3 * track_unit * (domains is not None)
        sequence_track = track_unit * seqbar
        self.track_height = sequence_track + annotations_track + domains_track
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
                    + (i + 0.5) * (4 * track_unit)
                ),
                height=2 * track_unit,
                ytrans="axes",
                region=self.region,
            )

        self.i += 1
        ylabel = column.replace("_", " ")
        self.set_labels(ax=ax, ylabel=ylabel, axis_title=label)
        self.set_axis(ax=ax, sequence=profile.sequence, nt_ticks=nt_ticks)

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
        """Set the size of the figure.

        Parameters
        ----------
        height_ax_rel : float, optional
            Height of the axes relative to the y-axis limits.
        width_ax_rel : float, defaults to 0.03
            Width of the axes relative to the x-axis limits.
        width_ax_in : float, optional
            Width of the axes in inches.
        height_ax_in : float, defaults to 2
            Height of the axes in inches.
        height_gap_in : float, defaults to 1
            Height of the gap between axes in inches.
        width_gap_in : float, defaults to 0.5
            Width of the gap between axes in inches.
        top_in : float, defaults to 1
            Height of the top margin in inches.
        bottom_in : float, defaults to 1
            Height of the bottom margin in inches.
        left_in : float, defaults to 1
            Width of the left margin in inches.
        right_in : float, defaults to 1
            Width of the right margin in inches.
        """
        super().set_figure_size(
            height_ax_rel=height_ax_rel,
            width_ax_rel=width_ax_rel,
            width_ax_in=width_ax_in,
            height_ax_in=height_ax_in,
            height_gap_in=height_gap_in + 2 * self.track_height,
            width_gap_in=width_gap_in,
            top_in=top_in,
            bottom_in=bottom_in + 2 * self.track_height,
            left_in=left_in,
            right_in=right_in,
        )

    def set_axis(self, ax, sequence, nt_ticks=(20, 5)):
        """Set up axis properties for aesthetics.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes object.
        sequence : str
            Nucleotide sequence.
        nt_ticks : tuple of int, optional
            Major and minor tick interval for the nucleotide axis.
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
