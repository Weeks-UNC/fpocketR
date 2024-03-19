import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mpc


class ScalarMappable(cm.ScalarMappable):
    """Used to map scalar values to a color and to create a colorbar plot.

    Parameters
    ----------
    cmap : str, tuple, float, or list
        A valid mpl color, list of valid colors or a valid colormap name
    normalization : "min_max", "0_1", "none", or "bins"
        The type of normalization to use when mapping values to colors
    values : list
        The values to use when normalizing the data
    title : str, defaults to ""
        The title of the colorbar.
    tick_labels : list, defaults to None
        The labels to use for the colorbar ticks. If None, values are
        determined automatically.
    **cbar_args : dict
        Additional arguments to pass to the colorbar function

    Attributes
    ----------
    rnav_norm : str
        The type of normalization to use when mapping values to colors
    rnav_vals : list
        The values to use when normalizing the data
    rnav_cmap : list
        The colors to use when mapping values to colors
    cbar_args : dict
        Additional arguments to pass to the colorbar function
    tick_labels : list
        The labels to use for the colorbar ticks. If None, values are
        determined automatically.
    title : str
        The title of the colorbar.
    """

    def __init__(
        self, cmap, normalization, values, title="", tick_labels=None, **cbar_args
    ):
        """Initialize the ScalarMappable object"""
        cmap = self.get_cmap(cmap)
        norm = self.get_norm(normalization, values, cmap)
        super().__init__(norm, cmap)
        self.rnav_norm = normalization
        self.rnav_vals = values
        self.rnav_cmap = cmap(np.arange(cmap.N))
        self.cbar_args = cbar_args
        self.tick_labels = tick_labels
        self.title = title

    def is_equivalent_to(self, cmap2):
        """Check if two ScalarMappable objects are equivalent.

        Parameters
        ----------
        cmap2 : ScalarMappable
            The ScalarMappable object to compare to

        Returns
        -------
        bool
            True if the two ScalarMappable objects are equivalent, False
            otherwise
        """
        return (
            (self.rnav_norm == cmap2.rnav_norm)
            and (self.rnav_vals == cmap2.rnav_vals)
            and (len(self.rnav_cmap) == len(cmap2.rnav_cmap))
            and mpc.same_color(self.rnav_cmap, cmap2.rnav_cmap)
            and (self.tick_labels == cmap2.tick_labels)
            and (self.cbar_args == cmap2.cbar_args)
            and (self.title == cmap2.title)
        )

    def values_to_hexcolors(self, values, alpha=1.0):
        """Map values to colors and return a list of hex colors.

        Parameters
        ----------
        values : list
            The values to map to colors
        alpha : float, defaults to 1.0
            The alpha value to use for the colors

        Returns
        -------
        list of strings
            A list of hex colors
        """
        colors = super().to_rgba(x=values, alpha=alpha)
        return np.array([mpc.to_hex(c, keep_alpha=True) for c in colors])

    def get_norm(self, normalization, values, cmap):
        """Given a normalization type and values, return a normalization object.

        Parameters
        ----------
        normalization : "min_max", "0_1", "none", or "bins"
            The type of normalization to use when mapping values to colors
        values : list
            The values to use when normalizing the data
        cmap : matplotlib colormap
            The colormap to use when normalizing the data

        Returns
        -------
        matplotlib.colors normalization object
            Used to normalize data before mapping to colors
        """
        if normalization == "min_max":
            return mpc.Normalize(values[0], values[1])
        elif normalization == "0_1":
            return mpc.Normalize()
        elif normalization == "none":
            return mpc.NoNorm()
        elif normalization == "bins":
            if len(values) == (cmap.N - 1):
                return mpc.BoundaryNorm(values, cmap.N, extend="both")
            if len(values) == (cmap.N + 1):
                return mpc.BoundaryNorm(values, cmap.N)

    def get_cmap(self, cmap):
        """Converts a cmap specification to a matplotlib colormap object.

        Parameters
        ----------
            cmap : string, tuple, float, or list
                A valid mpl color, list of valid colors or a valid colormap name

        Returns
        -------
            matplotlib colormap
                a colormap matching the input
        """
        if mpc.is_color_like(cmap):
            cmap = mpc.ListedColormap([cmap])
            cmap.set_bad("grey")
            return cmap
        elif isinstance(cmap, list) and all(mpc.is_color_like(c) for c in cmap):
            cmap = mpc.ListedColormap(cmap)
            cmap.set_bad("grey")
            return cmap
        try:
            return cm.get_cmap(cmap)
        except ValueError as exception:
            print(
                "cmap must be one of: valid mpl color, list of mpl colors, or "
                f"mpl colormap:\n{str(cmap)}"
            )
            raise exception
