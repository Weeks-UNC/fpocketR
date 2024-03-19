"""Classes for storing and manipulating data for RNAvigate.

This module contains the base classes for RNAvigate data classes:
    Sequence: represents a nucleotide sequence
    Data: represents a data table with a sequence
"""

from os.path import isfile
import Bio.SeqIO
import numpy as np
import pandas as pd
import matplotlib.colors as mpc
from rnavigate import styles
from rnavigate import data


def normalize_sequence(sequence, t_or_u="U", uppercase=True):
    """Returns sequence as all uppercase nucleotides and/or corrects T or U.

    Parameters
    ----------
    sequence : string or RNAvigate Sequence)
        The sequence
        If given an RNAvigate Sequence, the sequence string is retrieved
    t_or_u : "T", "U", or False, defaults to "U"
        "T" converts "U"s to "T"s
        "U" converts "T"s to "U"s
        False does nothing
    uppercase bool, defaults to True
        Whether to make sequence all uppercase

    Returns
    -------
        string
            the cleaned-up sequence string
    """
    if isinstance(sequence, Sequence):
        sequence = sequence.sequence
    if uppercase:
        sequence = sequence.upper()
    if not t_or_u:
        pass
    elif t_or_u.upper() == "U":
        sequence = sequence.replace("t", "u").replace("T", "U")
    elif t_or_u.upper() == "T":
        sequence = sequence.replace("u", "t").replace("U", "T")
    return sequence


class Sequence:
    """A class for storing and manipulating RNA sequences.

    Parameters
    ----------
    sequence : string or pandas.DataFrame
        sequence string, fasta file, or a Pandas dataframe containing a
        "Sequence" column
    name : string, optional
        The name of the sequence, defaults to None

    Attributes
    ----------
    sequence : string
        The sequence string
    name : string
        The name of the sequence
    other_info : dict
        A dictionary of additional information about the sequence
    null_alignment : SequenceAlignment
        An alignment of the sequence to itself
    """

    def __init__(self, input_data, name=None):
        """Initialize the Sequence object."""
        self.name = name
        self.other_info = dict()
        if isinstance(input_data, str):
            if isfile(input_data):
                self.sequence = self.read_fasta(input_data)
            else:
                self.sequence = input_data
        elif isinstance(input_data, pd.DataFrame):
            self.sequence = self.get_seq_from_dataframe(input_data)
        elif isinstance(input_data, Sequence):
            self.sequence = input_data.sequence
        self.null_alignment = data.SequenceAlignment(self, self)

    def __str__(self):
        """Return the name of the sequence."""
        if self.name is None:
            return "seq-object"
        return self.name

    def read_fasta(self, fasta):
        """Parse a fasta file for the first sequence.

        Parameters
        ----------
        fasta : string
            path to fasta file
        """
        with open(fasta, "r") as file:
            fasta = list(Bio.SeqIO.parse(file, "fasta"))
        return str(fasta[0].seq).replace("T", "U")

    def get_seq_from_dataframe(self, dataframe):
        """Parse a dataframe for the sequence string, store as self.sequence.

        Parameters
        ----------
        dataframe : pandas.DataFrame
            must contain a "Sequence" column
        """
        sequence = "".join(dataframe["Sequence"].values)
        return sequence.replace("T", "U").replace("t", "u")

    @property
    def length(self):
        """Get the length of the sequence

        Returns
        -------
        length : int
            the length of self.sequence
        """
        return len(self.sequence)

    def normalize_sequence(self, t_or_u="U", uppercase=True):
        """Converts sequence to all uppercase nucleotides and corrects T or U.

        Parameters
        ----------
        t_or_u : "T", "U", or False, defaults to "U"
            "T" converts "U"s to "T"s
            "U" converts "T"s to "U"s
            False does nothing.
        uppercase : bool, defaults to True
            Whether to make sequence all uppercase
        """
        self.sequence = normalize_sequence(self, t_or_u=t_or_u, uppercase=uppercase)

    def get_aligned_data(self, alignment):
        """Get a copy of the sequence positionally aligned to another sequence.

        Parameters
        ----------
        alignment : rnavigate.data.Alignment
            the alignment to use

        Returns
        -------
        aligned_sequence : rnavigate.data.Sequence
            the aligned sequence
        """
        return Sequence(alignment.target_sequence)

    def get_colors_from_sequence(self):
        """Get a colors and colormap representing the nucleotide sequence.

        Returns
        -------
        colors : numpy array
            one matplotlib color-like value for each nucleotide in self.sequence
        colormap : rnavigate.data.ScalarMappable
            a colormap used for creating a colorbar
        """
        colors = np.array([styles.get_nt_color(nt) for nt in self.sequence])
        colormap = data.ScalarMappable(
            cmap=[styles.get_nt_color(nt) for nt in "AUGC"],
            normalization="none",
            values=None,
            extend="neither",
            title="Nucleotide identity",
            alpha=1,
            ticks=[0, 1, 2, 3],
            tick_labels=["A", "U", "G", "C"],
        )
        return colors, colormap

    def get_colors_from_positions(self, pos_cmap="rainbow"):
        """Get colors and colormap representing the nucleotide position.

        Parameters
        ----------
        pos_cmap : str, defaults to "rainbow"
            cmap used for position colors

        Returns
        -------
        colors : numpy array
            one matplotlib color-like value for each nucleotide in self.sequence
        colormap : rnavigate.data.ScalarMappable
            a colormap used for creating a colorbar
        """
        colormap = data.ScalarMappable(
            cmap=pos_cmap,
            normalization="min_max",
            values=[1, self.length],
            extend="neither",
            title="Nucleotide position",
            alpha=1,
        )
        colors = colormap.values_to_hexcolors(np.arange(self.length))
        return colors, colormap

    def get_colors_from_profile(self, profile):
        """Get colors and colormap representing per-nucleotide data.

        Parameters
        ----------
        profile : rnavigate.data.Profile
            the profile to use to get colors.

        Returns
        -------
        colors : numpy array
            one matplotlib color-like value for each nucleotide in self.sequence
        colormap : rnavigate.data.ScalarMappable
            a colormap used for creating a colorbar
        """
        alignment = data.SequenceAlignment(profile, self)
        colors = alignment.map_values(profile.colors, fill="#808080")
        colormap = profile.cmap
        return colors, colormap

    def get_colors_from_annotations(self, annotations, default_color="gray"):
        """Get colors and colormap representing sequence annotations.

        Parameters
        ----------
        annotations : list of rnavigate.data.Annotations
            the annotations to use to get colors.
        default_color : matplotlib color-like, defaults to "gray"
            the color to use for nucleotides not in any annotation

        Returns
        -------
        colors : numpy array
            one matplotlib color-like value for each nucleotide in self.sequence
        colormap : rnavigate.data.ScalarMappable
            a colormap used for creating a colorbar
        """
        colors = np.full(self.length, default_color, dtype="<U16")
        cmap = [default_color]
        tick_labels = ["other"]
        for annotation in annotations:
            cmap.append(annotation.color)
            tick_labels.append(annotation.name)
            annotation = annotation.get_aligned_data(
                data.SequenceAlignment(annotation, self)
            )
            for site in annotation.get_sites():
                colors[site - 1] = annotation.color
        colormap = data.ScalarMappable(
            cmap=cmap,
            normalization="none",
            values=None,
            extend="neither",
            title="Annotations",
            alpha=1,
            ticks=list(range(len(annotations) + 1)),
            tick_labels=tick_labels,
        )
        return colors, colormap

    def get_colors_from_structure(self, structure):
        """Get colors and colormap representing base-pairing status.

        Parameters
        ----------
        structure : rnavigate.data.SecondaryStructure
            the structure to use to get colors.

        Returns
        -------
        colors : numpy array
            one matplotlib color-like value for each nucleotide in self.sequence
        colormap : rnavigate.data.ScalarMappable
            a colormap used for creating a colorbar
        """
        cmap = ["darkOrange", "darkOrchid", "gray"]
        ct_colors = [cmap[int(pair == 0)] for pair in structure.pair_nts]
        alignment = data.SequenceAlignment(structure, self)
        colors = alignment.map_values(ct_colors, fill="gray")
        colormap = data.ScalarMappable(
            cmap=cmap,
            normalization="none",
            values=None,
            extend="neither",
            title="Base-pairing status",
            alpha=1,
            ticks=[0, 1, 2],
            tick_labels=["unpaired", "paired", "unaligned"],
        )
        return colors, colormap

    def get_colors(
        self, source, pos_cmap="rainbow", profile=None, structure=None, annotations=None
    ):
        """Get colors and colormap representing information about the sequence.

        Parameters
        ----------
        source : str, list, or matplotlib color-like
            the source of the color information
            if a string, must be one of:
                "sequence", "position", "profile", "structure", "annotations"
            if a list, must be a list of matplotlib color-like values, colormap
                will be None.
            if a matplotlib color-like value, all nucleotides will be colored
                that color, colormap will be None.
        pos_cmap : str, defaults to "rainbow"
            cmap used for position colors if source is "position"
        profile : rnavigate.data.Profile, optional
            the profile to use to get colors if source is "profile"
        structure : rnavigate.data.SecondaryStructure, optional
            the structure to use to get colors if source is "structure"
        annotations : list of rnavigate.data.Annotations, optional
            the annotations to use to get colors if source is "annotations"

        Returns
        -------
        colors : numpy array
            one matplotlib color-like value for each nucleotide in self.sequence
        colormap : rnavigate.data.ScalarMappable
            a colormap used for creating a colorbar
        """
        if isinstance(source, str) and (source == "sequence"):
            return self.get_colors_from_sequence()
        elif isinstance(source, str) and (source == "position"):
            return self.get_colors_from_positions(pos_cmap=pos_cmap)
        elif isinstance(source, str) and (source == "profile"):
            return self.get_colors_from_profile(profile=profile)
        elif isinstance(source, str) and (source == "annotations"):
            return self.get_colors_from_annotations(annotations=annotations)
        elif isinstance(source, str) and (source == "structure"):
            return self.get_colors_from_structure(structure)
        elif mpc.is_color_like(source):
            return np.full(self.length, source, dtype="<U16"), None
        elif (len(source) == self.length) and all(mpc.is_color_like(c) for c in source):
            return np.array(list(source)), None
        else:
            print("Invalid colors:")
            print("\tchoices are profile, sequence, position, structure,")
            print("annotations, a list of mpl colors, or a single mpl color.")
            print("Defaulting to sequence.")
            return self.get_colors_from_sequence()


class Data(Sequence):
    """The base class for RNAvigate Profile and Interactions classes.

    Parameters
    ----------
    input_data : pandas.DataFrame or str
        a pandas dataframe or path to a data file
    sequence : string or rnavigate.data.Sequence
        the sequence to use for the data
    metric : string or dict
        the column of the dataframe to use as the default metric to visualize
    metric_defaults : dict
        a dictionary of metric defaults
    read_table_kw : dict, optional
        kwargs dictionary passed to pd.read_table
    name : string, optional
        the name of the data, defaults to None

    Attributes
    ----------
    data : pandas.DataFrame
        the data table
    filepath : string
        the path to the data file
    sequence : string or rnavigate.data.Sequence
        the sequence to use for the data
    metric : string or dict
        the column of the dataframe to use as the metric to visualize
    metric_defaults : dict
        A dictionary of metric values and default settings for visualization
    default_metric : string
        the default metric to use for visualization
    """

    def __init__(
        self,
        input_data,
        sequence,
        metric,
        metric_defaults,
        read_table_kw=None,
        name=None,
    ):
        """Initialize the Data object."""
        if read_table_kw is None:
            read_table_kw = {}
        # assign data
        if isinstance(input_data, pd.DataFrame):
            self.data = input_data
            self.filepath = "dataframe"
        elif isfile(input_data):
            self.data = self.read_file(input_data, read_table_kw)
            self.filepath = input_data
        else:
            print(f"{self} initialized without data.")
        # assign sequence
        if sequence is None:
            sequence = self.data
        super().__init__(sequence, name=name)
        # assign metrics
        self.metric_defaults = {
            "default": {
                "metric_column": "Profile",
                "error_column": None,
                "color_column": None,
                "cmap": "viridis",
                "normalization": "0_1",
                "values": None,
                "extend": "neither",
                "title": "Type: metric",
                "alpha": 0.7,
            },
            "Distance": {
                "metric_column": "Distance",
                "error_column": None,
                "color_column": None,
                "cmap": "cool",
                "normalization": "min_max",
                "values": [5, 50],
                "extend": "both",
                "title": "3D distance",
                "alpha": 0.7,
            },
        }
        self.add_metric_defaults(metric_defaults)
        self.default_metric = metric
        self._metric = None
        self.metric = metric

    def add_metric_defaults(self, metric_defaults):
        """Add metric defaults to self.metric_defaults"""
        default_defaults = self.metric_defaults["default"]
        for metric, defaults in metric_defaults.items():
            self.metric_defaults[metric] = default_defaults | defaults

    @property
    def metric(self):
        """Get the column of the dataframe to use as the metric for visualization."""
        return self._metric["metric_column"]

    @metric.setter
    def metric(self, value):
        """Set the column of the dataframe and other settings to use for visualization.

        Parameters
        ----------
        value : string or dict
            the column of the dataframe to use as the metric for visualization
            if a string, must be one of:
                "default", "Distance", or a column name in self.data
            if a dict, must contain at least one of:
                "metric_column", "error_column", "color_column", "cmap",
                "normalization", "values", "extend", "title", or "alpha"
        """
        if isinstance(value, str):
            try:
                self._metric = self.metric_defaults[value]
                return
            except KeyError as exception:
                if value not in self.data.columns:
                    print(
                        f"metric ({value}) not found in data:\n"
                        + str(list(self.data.columns))
                    )
                    raise exception
                self._metric = self.metric_defaults["default"]
                self._metric["metric_column"] = value
                return
        try:
            defaults = self.metric_defaults[value["metric_column"]]
        except KeyError:
            defaults = self.metric_defaults["default"]
        for key in value:
            if key not in defaults:
                raise ValueError(
                    f"{key} is not an expected value for metric "
                    f"setting:\n{list(defaults.keys())}"
                )
        self._metric = defaults | value

    @property
    def error_column(self):
        """Get the column of the dataframe to use as the error for visualization."""
        if self._metric["error_column"] in self.data.columns:
            return self._metric["error_column"]
        print(f"Warning: {self} missing expected error column")
        return None

    @property
    def color_column(self):
        """Get the column of the dataframe to use as the color for visualization."""
        if self._metric["color_column"] in self.data.columns:
            return self._metric["color_column"]
        return self._metric["metric_column"]

    @property
    def cmap(self):
        """Get the colormap to use for colorbars and to retrieve colors."""
        cmap_kwargs = {}
        for k, v in self._metric.items():
            if not k.endswith("column"):
                cmap_kwargs[k] = v
        return data.ScalarMappable(**cmap_kwargs)

    @property
    def colors(self):
        """Get one matplotlib color-like value for each nucleotide in self.sequence."""
        values = self.data[self.color_column]
        if values.dtype.name == "bool":
            values = values.astype("int")
        return self.cmap.values_to_hexcolors(np.ma.masked_invalid(values))

    def read_file(self, filepath, read_table_kw):
        """Convert data file to pandas dataframe and store as self.data

        Parameters
        ----------
        filepath : string
            path to data file containing interactions
        read_table_kw : dict
            kwargs dictionary passed to pd.read_table

        Returns
        -------
        dataframe : pandas.DataFrame
            the data table
        """
        return pd.read_table(filepath, **read_table_kw)
