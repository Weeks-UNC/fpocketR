import numpy as np
import matplotlib as mp
import matplotlib.patches as mp_patches
from rnavigate import styles


# 1-dimensional x-axis tracks
def plot_sequence_track(
    ax,
    sequence,
    yvalue=-0.05,
    height=0.05,
    ytrans="data",
    verticalalignment="bottom",
    region="all",
):
    """Plot a sequence track along the x-axis of a plot.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot the sequence track on.
    sequence : str
        The sequence to plot.
    yvalue : float, defaults to -0.05
        The y-value of the sequence track.
    height : float, defaults to 0.05
        The height of the sequence track.
    ytrans : "data" or "axes", defaults to "data"
        The y-axis coordinate system.
    verticalalignment : "top" or "bottom", defaults to "bottom"
        The vertical alignment of the sequence track.
    region : list of 2 int, defaults to "all"
        Start and end positions of the region to plot. If "all", plot the entire
        sequence.
    """
    style = styles.settings["sequence_bar"]
    ymin, ymax = ax.get_ylim()
    if ytrans == "axes":
        yvalue = (ymax - ymin) * yvalue + ymin
        height = (ymax - ymin) * height
    if region == "all":
        region = [1, len(sequence)]
    sequence = sequence[region[0] - 1 : region[1]]
    if style == "alphabet":
        for i, nt in enumerate(sequence):
            # set font style and colors for each nucleotide
            font_prop = mp.font_manager.FontProperties(
                family="monospace", style="normal", weight="bold", size="4"
            )
            col = styles.get_nt_color(nt)
            ax.text(
                i + region[0],
                yvalue,
                nt,
                fontproperties=font_prop,
                color=col,
                horizontalalignment="center",
                verticalalignment=verticalalignment,
            )
    if style == "bars":
        sequence = list(sequence)
        colors = [styles.get_nt_color(nt) for nt in sequence]
        is_uppercase = np.char.isupper(sequence)
        heights = np.full(len(sequence), height / 2)
        heights *= is_uppercase + 1
        heights *= np.char.not_equal(sequence, "-")
        if verticalalignment == "top":
            heights *= -1
        x = np.arange(region[0], region[1] + 1)
        ax.bar(
            x, heights, bottom=yvalue, width=1, color=colors, ec="none", clip_on=False
        )
    ax.set_ylim(ymin, ymax)


def plot_annotation_track(
    ax,
    annotation,
    yvalue,
    height,
    mode,
    region="all",
    ytrans="data",
):
    """Plot an annotation track along the x-axis of a plot.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot the annotation track on.
    annotation : rnavigate.data.Annotation
        The annotation to plot.
    yvalue : float
        The y-value of the annotation track.
    height : float
        The height of the annotation track.
    mode : "track" or "bar"
        The annotation mode.
    region : list of 2 int, defaults to "all"
        Start and end positions of the region to plot. If "all", plot the entire
        sequence.
    ytrans : "data" or "axes", defaults to "data"
        The y-axis coordinate system.
    """
    if region == "all":
        region = [1, annotation.length]
    mn, mx = region
    if ytrans == "axes":
        ymin, ymax = ax.get_ylim()
        yvalue = (ymax - ymin) * yvalue + ymin
        height = (ymax - ymin) * height
    color = annotation.color
    modes = ["track", "bar"]

    def plot_track(*s, alpha=1):
        start, end = min(s) - 0.5, max(s) + 0.5
        if mode == "track":
            ax.plot(
                [start, end],
                [yvalue] * 2,
                color=color,
                alpha=alpha,
                lw=5,
                solid_capstyle="butt",
                clip_on=False,
            )
        elif mode == "bar" and len(s) == 2:
            ax.axvspan(start - 0.5, end + 0.5, fc=color, ec="none", alpha=0.1)
        elif mode == "bar" and len(s) == 1:
            ax.axvline(start + 0.5, color=color, ls=":")

    assert mode in modes, f"annotation mode must be one of: {modes}"
    a_type = annotation.annotation_type
    if a_type == "spans" or (a_type == "primers" and mode == "bar"):
        for start, end in annotation:
            plot_track(start, end)
    elif a_type == "sites":
        for site in annotation.data["site"]:
            plot_track(site)
    elif a_type == "group":
        sites = annotation.data["site"]
        plot_track(min(sites), max(sites), alpha=0.5)
        for site in sites:
            plot_track(site)
    elif a_type == "primers":
        for start, end in annotation:
            if start > end:
                start, end = start + 0.5, end - 0.5
            else:
                start, end = start - 0.5, end + 0.5
            ax.arrow(
                x=start,
                y=yvalue,
                dx=end - start + 1,
                dy=0,
                color=color,
                shape="right",
                clip_on=False,
                head_width=height * 0.8,
                head_length=3,
                length_includes_head=True,
            )


def plot_domain_track(ax, spans, yvalue, height, region="all", ytrans="data"):
    """Plot a domain track along the x-axis of a plot.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot the domain track on.
    spans : rnavigate.data.Spans
        The spans to plot.
    yvalue : float
        The y-value of the domain track.
    height : float
        The height of the domain track.
    region : list of 2 int, defaults to "all"
        Start and end positions of the region to plot. If "all", plot the entire
        sequence.
    ytrans : "data" or "axes", defaults to "data"
        The y-axis coordinate system.
    """
    if region == "all":
        region = [1, spans.length]
    mn, mx = region
    text_height = height
    if ytrans == "axes":
        ymin, ymax = ax.get_ylim()
        yvalue = (ymax - ymin) * yvalue + ymin
        height = (ymax - ymin) * height
        text_height = 6
    name = spans.name
    color = spans.color
    for start, end in spans:
        start = max([mn, start])
        end = min([mx, end])
        rect = mp_patches.Rectangle(
            xy=(start - 0.5, yvalue),
            width=end - start + 1,
            height=height,
            linewidth=1,
            edgecolor="black",
            facecolor=color,
            clip_on=False,
        )
        ax.add_patch(rect)
        ax.text(
            x=(end - start) / 2 + start,
            y=yvalue + height / 2,
            s=name,
            zorder=25,
            horizontalalignment="center",
            verticalalignment="center",
            fontdict={"size": text_height * 1.5},
            clip_on=False,
        )
