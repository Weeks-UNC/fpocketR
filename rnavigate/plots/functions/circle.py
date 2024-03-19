from matplotlib.collections import PatchCollection
from matplotlib.patches import Path, PathPatch
import numpy as np
from rnavigate import styles


def plot_interactions_circle(ax, seq_circle, interactions):
    """Plot interactions on a circle plot.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on.
    seq_circle : rnavigate.data.SequenceCircle
        The sequence circle object containing nucleotide positions.
    interactions : rnavigate.data.Interactions
        The interactions to be plotted as arcs between nucleotides.
    """
    radius = seq_circle.radius
    theta = seq_circle.data["Theta"]
    ij_colors = interactions.get_ij_colors()
    patches = []
    for i, j, color in zip(*ij_colors):
        if j < i:
            i, j = j, i
        theta_i = theta[i - 1]
        theta_j = theta[j - 1]
        # scaling the center point towards zero depending on angle(i,j)
        theta_center = (theta_i + theta_j) / 2
        theta_diff = theta_j - theta_i
        if theta_diff > np.pi:
            theta_diff = 2 * np.pi - theta_diff
            theta_center -= np.pi
        r_center = radius * (1 - (theta_diff / 2 / np.pi)) ** 4
        verts = [[theta_i, radius], [theta_center, r_center], [theta_j, radius]]
        codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
        patches.append(PathPatch(Path(verts, codes), fc="none", ec=color))
    ax.add_collection(PatchCollection(patches, match_original=True))


def plot_annotation_circle(ax, seq_circle, annotation, offset=1):
    """Plot annotations on a circle plot.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on.
    seq_circle : rnavigate.data.SequenceCircle
        The sequence circle object containing nucleotide positions.
    annotation : rnavigate.data.Annotation
        The annotation to be plotted.
    offset : float, optional
        The offset from the circle circumference to plot the annotation.
    """
    radius = seq_circle.radius + offset
    color = annotation.color
    if annotation.annotation_type in ["spans", "primers"]:
        for start, end in annotation:
            if start > end:
                start, end = end, start
            theta = seq_circle.data.loc[start - 1 : end, "Theta"]
            if start == end:
                ax.scatter(theta, radius, color=color, **styles.settings["ss"]["sites"])
            else:
                ax.plot(theta, radius, color=color, **styles.settings["ss"]["spans"])
    elif annotation.annotation_type in ["sites", "group"]:
        sites = annotation.data["Site"] - 1
        theta = seq_circle.data.loc[sites, "Theta"]
        if annotation.annotation_type == "group":
            start, end = min(sites), max(sites)
            thetas = seq_circle.data.loc[start : end - 1, "Thetas"]
            ax.plot(thetas, radius, color=color, **styles.settings["ss"]["spans"])
        ax.scatter(theta, radius, color=color, **styles.settings["ss"]["sites"])
