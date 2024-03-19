"""The RNAvigate Sample object for parsing, assigning, and organizing data."""

from rnavigate import data
from rnavigate.data_loading import create_data


class Sample:
    """Loads and organizes RNA structural data for use with plotting functions.

    The Sample class stores all of the relevant experimental and computational
    structural data for a single RNA experiment. Between samples, common data
    types should be given a common data keyword so they can be easily compared.

    Parameters
    ----------
    sample : str
        an arbitrary name. This will be used as a label in plot legends and titles
        to differentiate it from other samples
    inherit : Sample or list of these, optional
        Data keywords and associated data from other samples become data keywords
        and associated data of this sample. This does not make additional copies
        of the data: i.e. operations that make changes to inherited data change the
        original sample, and any other samples that inherited that data. This can
        be useful to save time and memory on operations and large data structures
        that are shared between samples.
    keep_inherited_defaults : bool, default = True
        whether to keep inherited default keywords
    **data_keywords:
        There are many built-in data keywords with different expectations and
        behaviors. For a full list with input formats and output behavior, visit:
        https://rnavigate.readthedocs.io/en/latest/loading-data/

    Attributes
    ----------
    sample : str
        the name of the sample
    inputs : dict
        a dictionary of data keywords and their (user-defined) inputs
    data : dict
        a dictionary of data keywords and their associated data
    defaults : dict
        a dictionary of data classes and their default data keywords

    Example
    -------
    >>> sample = rnavigate.Sample(
    ...     sample="My sample",
    ...     shapemap="path/to/shapmapper_profile.txt",
    ...     ss="path/to/structure.ct",
    ...     ringmap="path/to/ringmapper_rings.txt",
    ...     pdb="path/to/pdb.pdb",
    ...     arbitrary_keyword={
    ...         "sites": [10, 20, 30],
    ...         "name": "sites of interest",
    ...         "color": "red",
    ...     },
    ... )
    >>> sample.print_data_keywords()
    My sample data keywords:
      annotations:
        arbitrary_keyword (default)
      profiles:
        shapemap (default)
      structures:
        ss (default)
      interactions:
        ringmap (default)
      pdbs:
        pdb (default)
    """

    def __init__(
        self, sample, inherit=None, keep_inherited_defaults=True, **data_keywords
    ):
        """Creates a Sample."""
        self.sample = sample
        self.inputs = {}
        self.data = {}
        self.defaults = {
            "default_annotation": None,
            "default_profile": None,
            "default_structure": None,
            "default_interactions": None,
            "default_pdb": None,
        }

        # inherit data from other Samples
        self.inherit_data(
            inherit=inherit,
            keep_inherited_defaults=keep_inherited_defaults,
            overwrite=False,
        )

        # get data and store inputs for all data_keywords
        for data_keyword, inputs in data_keywords.items():
            self.set_data(data_keyword=data_keyword, inputs=inputs)
            self.set_as_default(data_keyword=data_keyword, overwrite=False)

    def inherit_data(self, inherit, keep_inherited_defaults, overwrite):
        """retrieves and stores data and data keywords from other samples

        Parameters:
            inherit : Sample or list of Samples
                Other samples from which to inherit data and data keywords
            keep_inherited_defaults : bool
                Use default values from inherited samples
            overwrite : bool
                whether to overwrite any existing keywords with inherited keywords
        """
        if isinstance(inherit, (list, tuple)):
            for inherit_sample in inherit[::-1]:
                self.inherit_data(inherit_sample, keep_inherited_defaults, overwrite)
        elif isinstance(inherit, Sample):
            if overwrite:
                self.data |= inherit.data
                self.inputs |= inherit.inputs
            else:
                self.data = inherit.data | self.data
                self.inputs = inherit.inputs | self.inputs
            if keep_inherited_defaults and overwrite:
                self.defaults |= inherit.defaults
            elif keep_inherited_defaults:
                for k, v in inherit.defaults.items():
                    if self.defaults[k] is None:
                        self.defaults[k] = v
        elif inherit is not None:
            raise ValueError("inherit only accepts rnav.Sample or list of rnav.Sample")

    def set_data(self, data_keyword, inputs, overwrite=False):
        """Add data to Sample using the given data keyword and inputs

        This methods works similarly to the data keywords arguments used
        during Sample initialization:

            my_sample = rnavigate.Sample(
                sample="name",
                data_keyword=inputs)

        is equivalent to:

            my_sample = rnavigate.Sample(
                sample="name")
            my_sample.add_data(
                "data_keyword", inputs)

        Parameters
        ----------
            data_keyword : str
                a data keyword used to store and/or parse the inputs
            inputs : dict or rnavigate.data.Data
                a dictionary used to create the data object or a data object itself
            overwrite : bool, defaults to False
                whether to overwrite a pre-existing data_keyword
        """
        if (data_keyword in self.data) and not overwrite:
            raise ValueError(
                f"'{data_keyword}' is already a data keyword. "
                "Choose a different one."
            )
        try:
            self.data[data_keyword] = create_data(
                sample=self, data_keyword=data_keyword, inputs=inputs
            )
        except BaseException as e:
            raise ValueError(f"issue while loading {data_keyword}") from e
        self.inputs[data_keyword] = inputs

    def get_data(self, data_keyword, data_class=None):
        """Replaces data keyword with data object, even if nested.

        Parameters
        ----------
            data_keyword : rnavigate.data.Data or data keyword or list/dict of these
                If None, returns None.
                If a data keyword, returns associated data from sample
                If Data, returns that data.
                If a list or dictionary, returns list or dictionary with
                    data keyword values replaced with associated Data
            data_class : rnavigate.data.Data class or subclass, optional
                If provided, ensures that returned data is of this type.

        Returns
        -------
            Same type as data_keyword argument, but data keywords are replaced
                with associated data

        Raises
        ------
            ValueError:
                if data is not found in sample
            ValueError:
                if the data retrieved is not of the specified data_class
        """
        # handle special cases
        not_in_sample = ValueError(f"{data_keyword} not in {self.sample}.")
        wrong_class = ValueError(f"{data_keyword} is not {data_class}")
        if data_class is None:
            data_class = data.Sequence
        if isinstance(data_keyword, dict):
            return {k: self.get_data(v, data_class) for k, v in data_keyword.items()}
        elif isinstance(data_keyword, list):
            return [self.get_data(v, data_class) for v in data_keyword]
        elif isinstance(data_keyword, data.Sequence):
            if not isinstance(data_keyword, data_class):
                raise wrong_class
            return data_keyword
        elif data_keyword is None:
            return None
        elif isinstance(data_keyword, str) and data_keyword.startswith("default_"):
            data_keyword = self.defaults[data_keyword]
            if data_keyword is None:
                raise not_in_sample
        # if keyword is in sample.data, retreive Sequence object
        try:
            data_obj = self.data[data_keyword]
        except KeyError as exception:
            raise not_in_sample from exception
        if not isinstance(data_obj, data_class):
            raise wrong_class
        return data_obj

    def set_as_default(self, data_keyword, overwrite=True):
        """Set the given data keyword as the default for its data class

        It's data class is determined automatically. Only one default exists
        per data class and per Sample object.

        Parameters
        ----------
            data_keyword : str
                The data keyword to set as the default
            overwrite : bool, defaults to ``True``
                whether to overwrite a pre-existing default data keyword
        """
        data_object = self.data[data_keyword]
        if isinstance(data_object, data.Annotation):
            default_keyword = "default_annotation"
        elif isinstance(data_object, data.Profile):
            default_keyword = "default_profile"
        elif isinstance(data_object, data.SecondaryStructure):
            default_keyword = "default_structure"
        elif isinstance(data_object, data.Interactions):
            default_keyword = "default_interactions"
        elif isinstance(data_object, data.PDB):
            default_keyword = "default_pdb"
        else:
            return
        if overwrite or self.defaults[default_keyword] is None:
            self.defaults[default_keyword] = data_keyword

    def filter_interactions(
        self,
        interactions,
        metric=None,
        cmap=None,
        normalization=None,
        values=None,
        **kwargs,
    ):
        """sets coloring properties and applies filters to interactions data.

        Parameters
        ----------
            interactions : rnavigate.data.Interactions or data keyword string
                Interactions object to be filtered. If a string, value is
                replaced with self.get_data(interactions)
            metric : str, optional
                column of interactions data to be used as metric for coloring
                interactions.
                "Distance" will compute 3D distance in "pdb", defaulting to
                2'OH atom. "Distance_DMS" or "Distance_[atom id]" will use
                those atoms to compute distance.
            cmap (str | list, optional):
                sets the interactions colormap, used to color interactions
                according to metric values.
            normalization (str, optional):
                `"norm"`: extreme values in colormap are given to the extreme
                    values of interactions metric data
                `"bins"`: data are colored according to which bin they fall in
                    `values` defines bins (list, length = 2 less than cmap)
                `"min_max"`: extreme values in cmap are given to values beyond
                    minimum and maximum, defined by `values`
            values:
                behavior depends on normalization
                `"norm"`: values are not needed
                `"bins"`: list of floats containing the boundaries between bins
                    One fewer than the number of categories
                `"min_max"`: list of floats containing the minimum and maximum
            **kwargs: Other arguments are passed to interactions.filter()
        """
        # check for valid interactions data
        if interactions is None:
            return
        try:
            interactions = self.get_data(interactions, data.Interactions)
        except KeyError as exception:
            raise KeyError(f"{interactions} is not in {self.sample}") from exception

        if metric is None:
            metric = interactions.default_metric
        elif metric.startswith("Distance"):
            if len(metric.split("_")) == 2:
                metric, atom = metric.split("_")
            else:
                atom = "O2'"
            interactions.set_3d_distances(self.get_data("default_pdb"), atom)
        metric = {"metric_column": metric}
        if cmap is not None:
            metric["cmap"] = cmap
        if normalization is not None:
            metric["normalization"] = normalization
        if values is not None:
            metric["values"] = values
        interactions.metric = metric
        if "profile" in kwargs:
            kwargs["profile"] = self.data[kwargs["profile"]]
        else:
            try:
                kwargs["profile"] = self.get_data("default_profile")
            except ValueError:
                kwargs["profile"] = None
        if "structure" in kwargs:
            kwargs["structure"] = self.data[kwargs["structure"]]
        else:
            try:
                kwargs["structure"] = self.get_data("default_structure")
            except ValueError:
                kwargs["structure"] = None
        interactions.filter(**kwargs)

    def print_data_keywords(self, return_dict=False):
        """Print a nicely formatted, organized list of data keywords.

        Returns a dictionary of data keywords, organized by data type, if
        return_dict is True.
        """
        data_keywords = {
            "annotations": [],
            "profiles": [],
            "structures": [],
            "interactions": [],
            "pdbs": [],
        }
        for data_keyword, data_object in self.data.items():
            if data_keyword in self.defaults.values():
                data_keyword += " (default)"
            if isinstance(data_object, data.Annotation):
                data_keywords["annotations"].append(data_keyword)
            if isinstance(data_object, data.Profile):
                data_keywords["profiles"].append(data_keyword)
            if isinstance(data_object, data.SecondaryStructure):
                data_keywords["structures"].append(data_keyword)
            if isinstance(data_object, data.Interactions):
                data_keywords["interactions"].append(data_keyword)
            if isinstance(data_object, data.PDB):
                data_keywords["pdbs"].append(data_keyword)
        if return_dict:
            return data_keywords
        print(f"{self.sample} data keywords:")
        for k, v in data_keywords.items():
            print(f"  {k}:")
            for dkw in v:
                print(f"    {dkw}")
        print()
