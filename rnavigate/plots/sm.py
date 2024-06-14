import numpy as np
from rnavigate import plots
from matplotlib.patches import Rectangle
from rnavigate.styles import rx_color, bg_color, dc_color


class SM(plots.Plot):
    """Plot classic ShapeMapper plots.

    Parameters
    ----------
    nt_length : int
        Length of the nucleotide sequence.
    region : tuple of int, defaults to "all" (entire sequence)
        start and end position of the region to plot. If "all", plot the entire
        sequence.
    panels : list of str, defaults to ["profile", "rates", "depth"]
        Which panels to plot. Options are "profile", "rates", and "depth".

    Attributes
    ----------
    nt_length : int
        Length of the nucleotide sequence.
    region : tuple of int
        start and end position of the region to plot.
    fig : matplotlib.figure.Figure
        Figure object.
    axes : numpy.ndarray of matplotlib.axes.Axes
        Array of axes objects.
    i : int
        Index of the current plot.
    """

    def __init__(self, nt_length, region=None, panels=["profile", "rates", "depth"]):
        """Initialize the plot."""
        if region is None:
            self.nt_length = nt_length
            self.region = [1, nt_length]
        super().__init__(len(panels), len(panels), cols=1)
        self.panels = panels

    def set_figure_size(
        self,
        height_ax_rel=None,
        width_ax_rel=0.03,
        width_ax_in=None,
        height_ax_in=2,
        height_gap_in=0.5,
        width_gap_in=0.5,
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
        width_ax_rel : float, defaults to 0.03
            Width of the axes relative to the x-axis limits.
        width_ax_in : float, defaults to None
            Width of the axes in inches.
        height_ax_in : float, defaults to 2
            Height of the axes in inches.
        height_gap_in : float, defaults to 0.5
            Height of the gap between axes in inches.
        width_gap_in : float, defaults to 0.5
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

    def plot_data(self, profile, label):
        """Plot the data.

        Parameters
        ----------
        profile : rnavigate.profile.Profile
            Profile object.
        label : str
            Label for the plot.
        """
        self.fig.suptitle(label)
        for i, plot in enumerate(self.panels):
            ax = self.get_ax(i)
            plot_func = {
                "profile": self.plot_sm_profile,
                "rates": self.plot_sm_rates,
                "depth": self.plot_sm_depth,
            }[plot]
            plot_func(ax, profile)
        self.axes[-1, 0].set_xlabel("Nucleotide position")

    def plot_sm_profile(self, ax, profile):
        """Plots classic ShapeMapper profile on the given axis.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes object.
        profile : rnavigate.profile.Profile
            Profile object.
        """
        ymin, ymax = (-0.5, 4)
        near_black = (0, 0, 1 / 255.0)
        orange_thresh = 0.4
        red_thresh = 0.85
        colors, _ = profile.get_colors("profile", profile=profile)
        sample = profile.data["Norm_profile"].copy()
        sample[np.isnan(sample)] = -1
        ax.bar(
            profile.data["Nucleotide"],
            sample,
            align="center",
            width=1.05,
            color=colors,
            edgecolor=colors,
            linewidth=0.0,
            yerr=profile.data["Norm_stderr"],
            ecolor=near_black,
            capsize=1,
        )
        ax.set_title("Normalized Profile")
        ax.set_ylim(ymin, ymax)
        ax.set_xlim(1, profile.length)
        ax.yaxis.grid(True)
        ax.set_axisbelow(True)
        ax.set_ylabel("Shape Reactivity")
        # add a SHAPE colorbar to the vertical ax
        # uses a little transformation magic to place correctly
        inv = ax.transData.inverted()
        for loc, spine in list(ax.spines.items()):
            if loc == "left":
                trans = spine.get_transform()
        tp = trans.transform_point([0, 0])
        tp2 = inv.transform_point(tp)
        rectX = tp2[0]
        tpA = (0, 0)
        tpB = (6, 0)
        tpA2 = inv.transform_point(tpA)
        tpB2 = inv.transform_point(tpB)
        rectW = tpB2[0] - tpA2[0]
        rect = Rectangle(
            (rectX, -0.5),
            rectW,
            orange_thresh + 0.5,
            facecolor="black",
            edgecolor="none",
        )
        ax.add_patch(rect)
        rect.set_clip_on(False)
        rect = Rectangle(
            (rectX, orange_thresh),
            rectW,
            red_thresh - orange_thresh,
            facecolor="orange",
            edgecolor="none",
        )
        ax.add_patch(rect)
        rect.set_clip_on(False)
        rect = Rectangle(
            (rectX, red_thresh),
            rectW,
            4 - red_thresh,
            facecolor="red",
            edgecolor="none",
        )
        ax.add_patch(rect)
        rect.set_clip_on(False)
        ax.get_xaxis().tick_bottom()  # remove unneeded ticks
        ax.get_yaxis().tick_left()
        ax.tick_params(axis="y", which="minor", left=False)
        ax.minorticks_on()
        ax.set(yticks=[0, 1, 2, 3, 4])

        for line in ax.get_yticklines():
            line.set_markersize(6)
            line.set_markeredgewidth(1)

        for line in ax.get_xticklines():
            line.set_markersize(7)
            line.set_markeredgewidth(2)

        for line in ax.xaxis.get_ticklines(minor=True):
            line.set_markersize(5)
            line.set_markeredgewidth(1)

        # put nuc sequence below ax
        plots.plot_sequence_track(ax, profile.sequence, yvalue=0, ytrans="axes")

    def plot_sm_depth(self, ax, profile):
        """Plots classic ShapeMapper read depth on the given axis.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes object.
        profile : rnavigate.profile.Profile
            Profile object.
        """
        sample = profile.data
        ax.plot(
            sample["Nucleotide"],
            sample["Modified_read_depth"],
            linewidth=1.5,
            color=rx_color,
            alpha=1.0,
            label="Modified",
        )
        ax.plot(
            sample["Nucleotide"],
            sample["Untreated_read_depth"],
            linewidth=1.5,
            color=bg_color,
            alpha=1.0,
            label="Untreated",
        )
        ax.plot(
            sample["Nucleotide"],
            sample["Denatured_read_depth"],
            linewidth=1.5,
            color=dc_color,
            alpha=1.0,
            label="Denatured",
        )
        ax.set_xlim(1, profile.length)
        ax.legend(
            title="Effective depth\nin lighter colors",
            loc=2,
            borderpad=0.8,
            handletextpad=0.2,
            framealpha=0.75,
        )
        ax.plot(
            sample["Nucleotide"],
            sample["Modified_effective_depth"],
            linewidth=1.0,
            color=rx_color,
            alpha=0.3,
        )
        ax.plot(
            sample["Nucleotide"],
            sample["Untreated_effective_depth"],
            linewidth=1.0,
            color=bg_color,
            alpha=0.3,
        )
        ax.plot(
            sample["Nucleotide"],
            sample["Denatured_effective_depth"],
            linewidth=1.0,
            color=dc_color,
            alpha=0.3,
        )
        xmin, xmax, ymin, ymax = ax.axis()
        ax.set_ylim(0, ymax)
        ax.minorticks_on()
        ax.tick_params(axis="y", which="minor", left=False)
        yticks = [int(y) for y in ax.get_yticks()]
        formatted_ticks = [self.metric_abbreviate(val) for val in yticks]
        ax.set(yticks=yticks, yticklabels=formatted_ticks)
        for line in ax.get_yticklines():
            line.set_markersize(6)
            line.set_markeredgewidth(1)
        for line in ax.get_xticklines():
            line.set_markersize(7)
            line.set_markeredgewidth(2)
        for line in ax.xaxis.get_ticklines(minor=True):
            line.set_markersize(5)
            line.set_markeredgewidth(1)
        ax.yaxis.grid(True)
        ax.set_axisbelow(True)
        ax.set_ylabel("Read depth")
        # tried to make an offset, smaller font note about effective depths,
        # but couldn't get positioning/transforms to work properly.
        # For now just putting in xaxis label

    def plot_sm_rates(self, ax, profile):
        """Plots classic ShapeMapper mutation rates on the given axis.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes object.
        profile : rnavigate.profile.Profile
            Profile object.
        """
        sample = profile.data
        # choose a decent range for ax, excluding high-background positions
        temp_rates = sample.loc[sample["Untreated_rate"] <= 0.05, "Modified_rate"]
        near_top_rate = np.nanpercentile(temp_rates, 98.0)
        maxes = np.array([0.32, 0.16, 0.08, 0.04, 0.02, 0.01])
        ymax = np.amin(maxes[maxes > near_top_rate])
        rx_err = sample["Modified_rate"] / sample["Modified_effective_depth"]
        rx_upper = sample["Modified_rate"] + rx_err
        rx_lower = sample["Modified_rate"] - rx_err
        bg_err = sample["Untreated_rate"] / sample["Untreated_effective_depth"]
        bg_upper = sample["Untreated_rate"] + bg_err
        bg_lower = sample["Untreated_rate"] - bg_err
        dc_err = sample["Denatured_rate"] / sample["Denatured_effective_depth"]
        dc_upper = sample["Denatured_rate"] + dc_err
        dc_lower = sample["Denatured_rate"] - dc_err
        ax.set_ylabel("Mutation rate (%)")
        ax.plot(
            sample["Nucleotide"],
            sample["Modified_rate"],
            zorder=3,
            color=rx_color,
            linewidth=1.5,
            label="Modified",
        )
        ax.plot(
            sample["Nucleotide"],
            sample["Untreated_rate"],
            zorder=2,
            color=bg_color,
            linewidth=1.5,
            label="Untreated",
        )
        ax.plot(
            sample["Nucleotide"],
            sample["Denatured_rate"],
            zorder=2,
            color=dc_color,
            linewidth=1.5,
        )
        ax.fill_between(
            sample["Nucleotide"],
            rx_lower,
            rx_upper,
            edgecolor="none",
            alpha=0.5,
            facecolor=rx_color,
        )
        ax.fill_between(
            sample["Nucleotide"],
            bg_lower,
            bg_upper,
            edgecolor="none",
            alpha=0.5,
            facecolor=bg_color,
        )
        ax.fill_between(
            sample["Nucleotide"],
            dc_lower,
            dc_upper,
            edgecolor="none",
            alpha=0.5,
            facecolor=dc_color,
        )
        ax.legend(loc=2, borderpad=0.8, handletextpad=0.2, framealpha=0.75)
        ax.set_xlim((1, len(sample["Modified_rate"])))
        ax.set_ylim((0, ymax))
        ax.minorticks_on()
        ax.tick_params(axis="y", which="minor", left=False)
        yticks = ax.get_yticks()
        yticklabels = [str(int(x * 100)) for x in yticks]
        ax.set(yticks=yticks, yticklabels=yticklabels)
        for line in ax.get_yticklines():
            line.set_markersize(6)
            line.set_markeredgewidth(1)
        for line in ax.get_xticklines():
            line.set_markersize(7)
            line.set_markeredgewidth(2)
        for line in ax.xaxis.get_ticklines(minor=True):
            line.set_markersize(5)
            line.set_markeredgewidth(1)
        ax.yaxis.grid(True)
        ax.set_axisbelow(True)

    def metric_abbreviate(self, num):
        """takes a large number and applies an appropriate abbreviation

        Parameters
        ----------
            num : int
                number to be abbreviated

        Returns
        -------
            str
                abbreviated number
        """
        suffixes = {3: "k", 6: "M", 9: "G"}
        s = str(num)
        # replace trailing zeros with metric abbreviation
        zero_count = len(s) - len(s.rstrip("0"))
        suffix = ""
        new_string = str(s)
        for num_zeros in sorted(suffixes.keys()):
            if num_zeros <= zero_count:
                suffix = suffixes[num_zeros]
                new_string = s[:-num_zeros]
                new_string = new_string + suffix
        return new_string
