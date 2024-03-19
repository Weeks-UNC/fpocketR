import numpy as np
import matplotlib.patches as mp_patches
import matplotlib.collections as mp_collections
from scipy.spatial.distance import cdist
from rnavigate import styles, data, plots


# Secondary structure diagram ploting functions
def plot_structure_ss(ax, structure, colors):
    """Plot the structure of a secondary structure diagram.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on.
    structure : rnavigate.data.SecondaryStructure
        Secondary structure with diagram drawing coordinates to plot.
    colors : list
        List of colors to use for each nucleotide in the structure.
    """
    x = structure.xcoordinates
    y = structure.ycoordinates

    # ax.scatter(x, y, marker=".", c=colors, **styles.settings['ss']['points'])

    path = mp_patches.Path(np.column_stack([x, y]))
    verts = path.interpolated(steps=2).vertices
    x, y = verts[:, 0], verts[:, 1]
    colors = [x for c in colors for x in [c, c]][1:-1]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    current_color = None
    new_segments = []
    new_colors = []
    for segment, color in zip(segments, colors):
        if current_color is None:
            current_color = color
            current_segment = list(segment)
        elif current_color == color:
            current_segment.append(segment[-1])
        else:
            new_segments.append(current_segment)
            new_colors.append(current_color)
            current_color = color
            current_segment = list(segment)
    new_segments.append(current_segment)
    new_colors.append(current_color)

    lc = mp_collections.LineCollection(
        new_segments, colors=new_colors, **styles.settings["ss"]["structure"]
    )
    ax.add_collection(lc)


def plot_basepairs_ss(ax, structure, bp_style):
    """Plot the basepairs of a secondary structure diagram.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on.
    structure : rnavigate.data.SecondaryStructure
        Secondary structure with diagram drawing coordinates to plot.
    bp_style : "conventional", "dotted" or "line"
        Style of basepairs to plot.
    """
    if bp_style not in ["conventional", "dotted", "line"]:
        raise ValueError('bp_style must be one of "conventional", "dotted" or "line".')
    x = structure.xcoordinates
    y = structure.ycoordinates
    zorder = styles.settings["ss"]["basepairs"]["zorder"]
    for pair in structure.get_pairs():
        x = structure.xcoordinates[[p - 1 for p in pair]]
        y = structure.ycoordinates[[p - 1 for p in pair]]
        xdist = x[1] - x[0]
        ydist = y[1] - y[0]
        if xdist != 0:
            angle_xy = np.arctan(ydist / xdist)
        elif ydist > 0:
            angle_xy = 1 / 2 * np.pi
        elif ydist < 0:
            angle_xy = 3 / 2 * np.pi
        if xdist < 0:
            angle_xy += np.pi
        x_offset = np.cos(angle_xy) / 3
        y_offset = np.sin(angle_xy) / 3
        x_caps = [x[0] + x_offset, x[1] - x_offset]
        y_caps = [y[0] + y_offset, y[1] - y_offset]
        if bp_style == "dotted":
            x_caps_dist = x_caps[0] - x_caps[1]
            y_caps_dist = y_caps[0] - y_caps[1]
            caps_dist = (x_caps_dist**2 + y_caps_dist**2) ** 0.5
            segs = int(max([1 / 3, caps_dist]) * 3) * 2
            x_dots = [x_caps[0] - x_caps_dist * i / segs for i in range(segs + 1)]
            y_dots = [y_caps[0] - y_caps_dist * i / segs for i in range(segs + 1)]
            ax.scatter(x_dots, y_dots, c="grey", marker=".", s=4, zorder=zorder)
        if bp_style == "line":
            ax.plot(x_caps, y_caps, color="grey", zorder=zorder)
        if bp_style == "conventional":
            nts = "".join([structure.sequence[p - 1] for p in pair]).upper()
            # x_caps = [x[0] + i * xdist / 7 for i in [2, 5]]
            # y_caps = [y[0] + i * ydist / 7 for i in [2, 5]]
            if nts in ["UA", "AU"]:
                ax.plot(
                    x_caps, y_caps, color="grey", zorder=zorder, solid_capstyle="butt"
                )
            if nts in ["GU", "UG", "AG", "GA"]:
                x_center = x[0] + xdist / 2
                y_center = y[0] + ydist / 2
                if nts in ["GU", "UG"]:
                    fc = "white"
                    ec = "grey"
                elif nts in ["AG", "GA"]:
                    fc = "grey"
                    ec = "none"
                ax.scatter(
                    x_center,
                    y_center,
                    zorder=zorder + 1,
                    s=3**2,
                    color="grey",
                    marker="o",
                    fc=fc,
                    ec=ec,
                )
            if nts in ["GC", "CG"]:
                ax.plot(
                    x_caps,
                    y_caps,
                    color="grey",
                    zorder=zorder,
                    linewidth=2,
                    solid_capstyle="butt",
                )
                ax.plot(x_caps, y_caps, color="white", zorder=zorder, linewidth=0.7)


def plot_sequence_ss(ax, structure, colors):
    """Plot the sequence of a secondary structure diagram.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on.
    structure : rnavigate.data.SecondaryStructure
        Secondary structure with diagram drawing coordinates to plot.
    colors : list
        List of colors to use for each nucleotide in the structure.
    """
    if isinstance(structure, data.SequenceCircle):
        x = structure.data["Theta"]
        y = np.full(x.shape, structure.radius)
    else:
        x = structure.xcoordinates
        y = structure.ycoordinates
    if colors is None:
        colors = structure.get_colors("gray")
    for nuc in "GUTACgutac":
        mask = [nt == nuc for nt in structure.sequence]
        xcoords = x[mask]
        ycoords = y[mask]
        marker = "$\mathsf{" + nuc + "}$"
        ax.scatter(
            xcoords,
            ycoords,
            marker=marker,
            c=colors[mask],
            **styles.settings["ss"]["sequence"]
        )


def plot_nucleotides_ss(ax, structure, colors):
    """Plot the nucleotides of a secondary structure diagram.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on.
    structure : rnavigate.data.SecondaryStructure
        Secondary structure with diagram drawing coordinates to plot.
    colors : list
        List of colors to use for each nucleotide in the structure.
    """
    if isinstance(structure, data.SequenceCircle):
        x = structure.data["Theta"]
        y = np.full(x.shape, structure.radius)
    else:
        x = structure.xcoordinates
        y = structure.ycoordinates
    if colors is None:
        colors = structure.get_colors("grey")
    ax.scatter(x, y, c=colors, **styles.settings["ss"]["nucleotides"])


def plot_positions_ss(ax, structure, xticks=20):
    """Plot the positions of a secondary structure diagram.

    Label locations are chosen from a point on a circle around each position that is
    the furthest from any other nucleotides. This sometimes causes tick marks and
    labels to overlap with other plot elements.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on.
    structure : rnavigate.data.SecondaryStructure
        Secondary structure with diagram drawing coordinates to plot.
    xticks : int
        Spacing between position labels.
    """
    xs = structure.xcoordinates
    ys = structure.ycoordinates
    thetas = np.pi / 32 * np.arange(64)
    x_shift = np.sin(thetas)
    y_shift = np.cos(thetas)
    ticks, labels = plots.get_nt_ticks(
        sequence=structure.sequence, region=(1, structure.length), gap=xticks
    )
    for nt, label in zip([t - 1 for t in ticks], labels):
        x_nt = xs[nt]
        y_nt = ys[nt]
        x_pos = xs[nt] + (x_shift * 1.5)
        y_pos = ys[nt] + (y_shift * 1.5)
        not_nt = np.arange(structure.length) != nt
        nt_box = (np.abs(xs - xs[nt]) < 4) & (np.abs(ys - ys[nt]) < 4)
        # print(np.vstack((x_pos, y_pos)).T)
        dists = cdist(
            np.vstack((x_pos, y_pos)).T,
            np.vstack((xs[not_nt & nt_box], ys[not_nt & nt_box])).T,
            "euclidean",
        )
        min_dists = np.min(dists, axis=1)
        which_pos = np.where(min_dists == np.max(min_dists))[0][0]
        x_label = x_pos[which_pos]
        y_label = y_pos[which_pos]
        bbox = dict(boxstyle="round", fc="w", ec="w", pad=0.1)
        ax.plot([x_nt, x_label], [y_nt, y_label], color="k", lw=1)
        ax.text(
            x_label,
            y_label,
            label,
            ha="center",
            va="center",
            bbox=bbox,
            **styles.settings["ss"]["positions"]
        )


def plot_interactions_ss(ax, structure, interactions):
    """Plot the interactions as lines over a secondary structure diagram.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on.
    structure : rnavigate.data.SecondaryStructure
        Secondary structure with diagram drawing coordinates.
    interactions : rnavigate.data.Interactions
        Interactions to plot.
    """
    ij_colors = interactions.get_ij_colors()
    segments = []
    colors = []
    for i, j, color in zip(*ij_colors):
        segments.append(
            [
                [structure.xcoordinates[i - 1], structure.ycoordinates[i - 1]],
                [structure.xcoordinates[j - 1], structure.ycoordinates[j - 1]],
            ]
        )
        colors.append(color)
    ax.add_collection(
        mp_collections.LineCollection(
            segments=segments, colors=colors, **styles.settings["ss"]["interactions"]
        )
    )


def plot_annotation_ss(ax, structure, annotation):
    """Highlight regions or nucleotides of interest on a secondary structure diagram.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on.
    structure : rnavigate.data.SecondaryStructure
        Secondary structure with diagram drawing coordinates.
    annotation : rnavigate.data.Annotation
        Annotation to plot.
    """
    color = annotation.color
    if annotation.annotation_type in ["spans", "primers"]:
        for start, end in annotation:
            if start > end:
                start, end = end, start
            x = structure.xcoordinates[start - 1 : end]
            y = structure.ycoordinates[start - 1 : end]
            if start == end:
                ax.scatter(x, y, color=color, **styles.settings["ss"]["sites"])
            else:
                ax.plot(x, y, color=color, **styles.settings["ss"]["spans"])
    elif annotation.annotation_type in ["sites", "group"]:
        sites = np.array(annotation[:]) - 1
        x = structure.xcoordinates[sites]
        y = structure.ycoordinates[sites]
        if annotation.annotation_type == "group":
            ax.plot(x, y, color=color, **styles.settings["ss"]["spans"])
        ax.scatter(x, y, color=color, **styles.settings["ss"]["sites"])
