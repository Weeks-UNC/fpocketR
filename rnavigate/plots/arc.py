"""Module for plotting arc plots."""

from rnavigate import plots, styles


class AP(plots.Plot):
    """Class for plotting arc plots

    Parameters
    ----------
    num_samples : int
        Number of samples to plot
    nt_length : int
        Length of the sequence
    region : tuple of 2 integers, optional
        starting and ending positions of the region to plot.
        Default is "all", which plots the entire sequence.
    track_labels : bool, optional
        Whether to plot track labels. Default is True.
    **kwargs
        Additional keyword arguments are passed to plots.Plot

    Attributes
    ----------
    nt_length : int
        Length of the sequence
    region : tuple of 2 integers
        starting and ending positions of the region to plot.
    track_labels : bool
        Whether to plot track labels.
    fig : matplotlib.figure.Figure
        Figure object containing the plot
    axes : numpy.ndarray of matplotlib.axes.Axes
        Array of axes objects containing the plots
    i : int
        Index of the current plot
    """

    def __init__(
        self, num_samples, nt_length, region="all", track_labels=True, **kwargs
    ):
        """Initialize AP object."""
        self.track_labels = track_labels
        if region == "all":
            self.nt_length = nt_length
            self.region = (1, nt_length)
        else:
            self.nt_length = region[1] - region[0] + 1
            self.region = region
        super().__init__(num_samples, **kwargs)

    def plot_data(
        self,
        sequence,
        structure=None,
        structure2=None,
        interactions=None,
        interactions2=None,
        profile=None,
        annotations=None,
        domains=None,
        label="",
        ax=None,
        seqbar=True,
        title=True,
        panels=None,
        annotation_mode="track",
        track_height=None,
        profile_scale_factor=1,
        plot_error=False,
        nt_ticks=(20, 5),
    ):
        """Add data to the next (or specified) plot axes.

        This function assumes data has already been aligned to a common sequence.
        rnavigate.plot_ functions can be used to automatically align data.

        Parameters
        ----------
        sequence : rnavigate.data.Sequence
            Sequence object containing the sequence to plot
        structure : rnavigate.data.Structure, optional
            Structure object containing a structure to plot
        structure2 : rnavigate.data.Structure, optional
            Structure object containing a structure to compare to the first
        interactions : rnavigate.data.Interactions, optional
            Interactions object containing inter-nucleotide data to plot
        interactions2 : rnavigate.data.Interactions, optional
            Interactions object containing other inter-nucleotide data to plot
        profile : rnavigate.data.Profile, optional
            Profile object containing per-nucleotide data to plot
        annotations : list of rnavigate.data.Annotation, optional
            List of Annotation objects containing annotations to plot
        domains : list of rnavigate.data.Spans, optional
            List of Spans objects containing domains to plot
        label : str, defaults to ""
            Label for the title of the plot.
        ax : matplotlib.axes.Axes, optional
            Axes object to plot on. If None, the next axes in the figure will be used.
        seqbar : bool, Defaults to True
            Whether to plot the sequence track.
        title : bool, defaults to True
            Whether to show the title.
        panels : dict, optional
            Dictionary of panels to plot, with keys being the panel name and values
            being the panel location. Default is {"interactions": "bottom",
            "interactions2": "bottom", "structure": "top", "profile": "top"}
        annotation_mode : "track" or "bar", defaults to "track"
            Mode for plotting annotations.
        track_height : int, optional
            Height of the track. If None, the height is automatically determined.
        profile_scale_factor : float, defaults to 1
            Scale factor for the profile track.
        plot_error : bool, defaults to False
            Whether to plot the error bars for the profile track.
        nt_ticks : tuple of 2 ints, optional
            Major and minor tick spacing for the nucleotide axis. Default is (20, 5).
        """
        ax = self.get_ax(ax)
        if panels is None:
            panels = {}
        panels = {
            "interactions": "bottom",
            "interactions2": "bottom",
            "structure": "top",
            "profile": "top",
        } | panels
        if annotations is None:
            annotations = []

        annotation_gap = 4 * len(annotations) * (annotation_mode == "track")
        seqbar_height = 2 * seqbar
        domains_height = 6 * (domains is not None)
        if track_height is None:
            track_height = annotation_gap + seqbar_height + domains_height

        yvalues = {"bottom": 0, "top": track_height}
        if structure is not None:
            structure = structure.as_interactions(structure2)
            plots.plot_interactions_arcs(
                ax=ax,
                interactions=structure,
                panel=panels["structure"],
                yvalue=yvalues[panels["structure"]],
                region=self.region,
            )
            self.add_colorbar_args(structure.cmap)
        if interactions is not None:
            plots.plot_interactions_arcs(
                ax=ax,
                interactions=interactions,
                panel=panels["interactions"],
                yvalue=yvalues[panels["interactions"]],
                region=self.region,
            )
            self.add_colorbar_args(interactions.cmap)
        if interactions2 is not None:
            plots.plot_interactions_arcs(
                ax=ax,
                interactions=interactions2,
                panel=panels["interactions2"],
                yvalue=yvalues[panels["interactions2"]],
                region=self.region,
            )
            self.add_colorbar_args(interactions2.cmap)
        if profile is not None:
            if panels["profile"] == "bottom":
                scale_factor = profile_scale_factor * -1
                bottom = 0
            else:
                scale_factor = profile_scale_factor
                bottom = track_height
            plots.plot_profile_bars(
                ax=ax,
                profile=profile,
                bottom=bottom,
                scale_factor=scale_factor,
                plot_error=plot_error,
                region=self.region,
            )
            self.add_colorbar_args(profile.cmap)

        yticks, ylabels = [], []
        if seqbar:
            plots.plot_sequence_track(
                ax,
                sequence.sequence,
                yvalue=0,
                height=2,
                ytrans="data",
                region=self.region,
            )
            yticks.append(1)
            ylabels.append("sequence")
            self.add_colorbar_args(styles.get_nt_cmap())
        if domains is not None:
            for domain in domains:
                plots.plot_domain_track(
                    ax=ax,
                    spans=domain,
                    yvalue=seqbar_height,
                    height=6,
                    region=self.region,
                )
            yticks.append(seqbar_height + 3)
            ylabels.append("domains")
        for i, annotation in enumerate(annotations):
            yvalue = seqbar_height + domains_height + 4 * (i + 1) - 2
            plots.plot_annotation_track(
                ax,
                annotation=annotation,
                yvalue=yvalue,
                height=4,
                mode=annotation_mode,
                region=self.region,
            )
            yticks.append(yvalue)
            ylabels.append(annotation.name)
        if not self.track_labels:
            yticks, ylabels = [], []
        self.set_axis(
            ax=ax,
            sequence=sequence.sequence,
            track_height=track_height,
            nt_ticks=nt_ticks,
            yticks=yticks,
            ylabels=ylabels,
        )
        if title:
            ax.set_title(label)
        self.i += 1

    def set_axis(
        self,
        ax,
        sequence,
        track_height=0,
        nt_ticks=(20, 5),
        max_height=300,
        yticks=None,
        ylabels=None,
    ):
        """Set up the plotting axis settings for an aesthetic arc plot.

        Sets the following properties of the given axis:
        1. spine positions
        2. x-axis and y-axis limits
        3. x-axis tick labels and positions according to `sequence`
        4. background boxes for x-axis tick labels
        5. y-axis tick labels and positions according to `track_height` and `ylabels`

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes object to set up.
        sequence : str
            Sequence to plot.
        track_height : int, optional
            Height of the track. If None, the height is automatically determined.
        nt_ticks : tuple of 2 ints, optional
            Major and minor tick spacing for the nucleotide axis. Default is (20, 5).
        max_height : int, optional
            Maximum height of the plot. Default is 300.
        yticks : list of ints, optional
            List of ytick positions. If None, the yticks are automatically determined.
        ylabels : list of str, optional
            List of ytick labels. If None, the ylabels are automatically determined.
        """
        ax.spines["left"].set_color("none")
        ax.spines["right"].set_color("none")
        ax.spines["bottom"].set(position=("data", 0), visible=False)
        ax.spines["top"].set_color("none")
        height = min(max_height, self.nt_length / 2)
        mn, mx = self.region
        ax.set(
            xlim=(mn - 0.5, mx + 0.5),
            ylim=(-height - 1, height + 1 + track_height),
            yticks=yticks,
            yticklabels=ylabels,
            axisbelow=False,
        )
        plots.set_nt_ticks(
            ax=ax,
            sequence=sequence,
            region=self.region,
            major=nt_ticks[0],
            minor=nt_ticks[1],
        )
        for label in ax.get_xticklabels():
            label.set_bbox(
                {
                    "facecolor": "white",
                    "edgecolor": "None",
                    "alpha": 0.5,
                    "boxstyle": "round,pad=0.1,rounding_size=0.2",
                }
            )

    def set_figure_size(
        self,
        height_ax_rel=0.03,
        width_ax_rel=0.03,
        width_ax_in=None,
        height_ax_in=None,
        height_gap_in=0.5,
        width_gap_in=0.5,
        top_in=1,
        bottom_in=1,
        left_in=1,
        right_in=1,
    ):
        """Set the figure size for an arc plot.

        Parameters
        ----------
        height_ax_rel : float, Default is 0.03.
            Relative height of each axes.
        width_ax_rel : float, Default is 0.03.
            Relative width of each axes.
        width_ax_in : float, optional
            Absolute width of each axes in inches. If None, the width is automatically
            determined.
        height_ax_in : float, optional
            Absolute height of each axes in inches. If None, the height is
            automatically determined.
        height_gap_in : float, Default is 0.5.
            Absolute height of the gap between axes in inches.
        width_gap_in : float, Default is 0.5.
            Absolute width of the gap between axes in inches.
        top_in : float, Default is 1.
            Absolute top margin in inches.
        bottom_in : float, Default is 1.
            Absolute bottom margin in inches.
        left_in : float, Default is 1.
            Absolute left margin in inches.
        right_in : float, Default is 1.
            Absolute right margin in inches.
        """
        super().set_figure_size(
            height_ax_rel=height_ax_rel,
            width_ax_rel=width_ax_rel,
            width_ax_in=width_ax_in,
            height_ax_in=height_ax_in,
            height_gap_in=height_gap_in,
            width_gap_in=width_gap_in + self.track_labels * 0.5,
            top_in=top_in,
            bottom_in=bottom_in,
            left_in=left_in,
            right_in=right_in,
        )
