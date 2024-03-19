"""Contains global plot display settings."""
from functools import wraps
import matplotlib as mpl
import seaborn as sns
from rnavigate import data

# STYLE SHEET
###############################################################################

__all__ = [
    "update_copy",
    "Settings",
    "set_defaults",
    "get_nt_color",
    "get_nt_cmap",
    "apply_style",
]

settings = {
    "sequence_bar": "alphabet",  # "bars"
    "sequence_colors": "rnavigate",  # "old"
    "ss": {
        "structure": {
            "linewidth": 1,
            "zorder": 2,
        },
        "interactions": {
            "linewidth": 3,
            "alpha": None,
            "zorder": 15,
        },
        "spans": {
            "linewidth": 10,
            "alpha": 0.4,
            "zorder": 5,
        },
        "sites": {
            "marker": "o",
            "edgecolor": "none",
            "s": 10**2,
            "alpha": 0.4,
            "zorder": 5,
        },
        "basepairs": {"zorder": 0},
        "nucleotides": {
            "marker": "o",
            "edgecolor": "none",
            "s": 5**2,
            "zorder": 10,
        },
        "sequence": {
            "linewidth": 0.3,
            "s": 3**2,
            "zorder": 20,
        },
        "positions": {"zorder": 25},
    },
}


def update_copy(original_settings, user_settings):
    """Recursively updates and returns a copy of og settings with new settings applied.

    Parameters
    ----------
        original_settings (dict)
            a default settings dictionary, usually rnav.settings
        user_settings (dict)
            a dictionary with only the fields from the original_settings that
            are to be changed

    Returns
    -------
        settings dict
            the original_settings dictionary with the new_settings dictionary
            values recursively applied
    """
    new_settings = dict()
    for k, v in original_settings.items():
        if isinstance(v, dict):
            try:
                new_settings[k] = update_copy(v, user_settings[k])
            except KeyError:
                new_settings[k] = v | {}
        else:
            try:
                new_settings[k] = user_settings[k]
            except KeyError:
                new_settings[k] = v
    return new_settings


class Settings(dict):
    """Context manager for temporarily changing global settings.

    Parameters
    ----------
        user_settings : dict
            a dictionary with only the fields from the original_settings that
            are to be changed
    """

    def __init__(self, user_settings):
        self.original_settings = update_copy(settings, {})
        self.user_settings = update_copy(settings, user_settings)

    def __enter__(self):
        settings.update(self.user_settings)

    def __exit__(self, *args, **kwargs):
        settings.update(self.original_settings)


def set_defaults(context="paper", style="ticks", colors="default", dpi=140):
    """Set or reset the major global style settings.

    Parameters
    ----------
        context : str, defaults to "paper"
            Passed to seaborn.set_context
            Defaults to "paper"
        style : str, defaults to "ticks"
            Passed to seaborn.set_style
        colors : str, defaults to "default"
            Passed to seaborn.set_palette
        dpi : int, defaults to 140
            Sets the dots-per-inch for inline and exported images
    """
    if colors == "default":
        colors = [
            "#0092edff",  # Blue
            "#ff8300ff",  # Orange
            "#a100ffff",  # Purple
            "#edc600ff",  # Yellow
            "#ff48e9ff",  # Pink
            "#3fd125ff",  # Green
        ]
    sns.set_context(context)
    sns.set_style(style)
    sns.set_palette(colors)
    mpl.rcParams["font.sans-serif"].insert(0, "Arial")
    mpl.rcParams["figure.dpi"] = dpi
    mpl.rcParams["svg.fonttype"] = "none"


def get_nt_color(nt, colors=None):
    """Get the RNAvigate color for a given nucleotide.

    Invalid nucleotides are set to gray

    Parameters
    ----------
        nt : str
            a nucleotide letter
        colors "rnavigate" or "old", defaults to settings["sequence_colors"]
            "rnavigate" uses blue, light blue, red, light red for "AUCG"
            "old" uses traditional red, yellow, blue, green for "AUCG"

    Returns:
        color : str
            a hex color string
    """
    if colors is None:
        colors = settings["sequence_colors"]
    try:
        return {
            "old": {
                "A": "#f20000",  # red
                "U": "#f28f00",  # yellow
                "G": "#00509d",  # blue
                "C": "#00c200",  # green
            },
            "rnavigate": {
                "A": "#366ef0",  # blue
                "U": "#9bb9ff",  # light blue
                "G": "#f04c4c",  # red
                "C": "#ffa77c",  # light red
            },
        }[colors][nt.upper().replace("T", "U")]
    except KeyError:
        return "#aaaaaa"


def get_nt_cmap(colors=None):
    """Get an rnavigate color map for nucleotides.

    Parameters
    ----------
        colors "rnavigate" or "old", defaults to settings["sequence_colors"]
            "rnavigate" uses blue, light blue, red, light red for "AUCG"
            "old" uses traditional red, yellow, blue, green for "AUCG"

    Returns
    -------
        cmap : rnavigate.data.colors.ScalarMappable (matplotlib.cm.ScalarMappable)
            a color map for nucleotides
    """
    return data.ScalarMappable(
        cmap=[get_nt_color(nt, colors=colors) for nt in "AUGC"],
        normalization="none",
        values=None,
        title="sequence",
        tick_labels=["A", "U", "G", "C"],
    )


def apply_style(style_dict):
    """Decorator for applying matplotlib style settings to a function."""

    def decorator(function):
        @wraps(function)
        def wrapper(*args, **kwargs):
            with mpl.rc_context(style_dict):
                return function(*args, **kwargs)

        return wrapper

    return decorator


# TODO: these belong in settings dict
# ShapeMapper Plot Styles
sm = {
    "font.family": "sans-serif",
    "pdf.fonttype": 42,
    # use TrueType fonts when exporting PDFs
    # (embeds most fonts - this is especially
    #  useful when opening in Adobe Illustrator)
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.fontsize": 14,
    "grid.color": ".8",
    "grid.linestyle": "-",
    "grid.linewidth": 1,
}

rx_color = "red"
bg_color = "blue"
dc_color = "darkgoldenrod"


# set default styles for plotting
set_defaults()
