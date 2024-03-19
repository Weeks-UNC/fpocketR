from operator import ge, le, gt, lt, eq, ne
import pandas as pd
import matplotlib.colors as mpc
import numpy as np
import seaborn as sns
from rnavigate import data


class Interactions(data.Data):
    """A class for storing and manipulating interactions data.

    Parameters
    ----------
    input_data : string or pandas.DataFrame
        If string, a path to a file containing interactions data.
        If dataframe, the dataframe containing interactions data. The dataframe
        must contain columns "i", "j", and self.metric. Dataframe may also
        include other columns.
    sequence : string or rnavigate.data.Sequence
        The sequence string corresponding to the interactions data.
    metric : string
        The column name to use for visualization.
    metric_defaults : dict
        Keys are metric names and values are dictionaries of metric-specific defaults.
        These defaults include:
            "metric_column" : string
                the column name to use for visualization
            "cmap" : string or matplotlib.colors.Colormap)
                the colormap to use for visualization
            "normalization" : "min_max", "0_1", "none", or "bins"
                The type of normalization to use when mapping values to colors
            "values" : list of float
                The values to used with normalization of the data
            "title" : string
                the title to use for colorbars
            "extend" : "min", "max", "both", or "neither"
                Which ends to extend when drawing the colorbar.
            "tick_labels" : list of string
    read_table_kw : dict
        kwargs passed to pandas.read_table() when reading input_data.
    window : int
        The window size used to generate the interactions data.
    name : str
        The name of the data object.

    Attributes
    ----------
    data : pandas.DataFrame
        The interactions data.
    window : int
        The window size that is being represented by i-j pairs.
    """

    def __init__(
        self,
        input_data,
        sequence,
        metric,
        metric_defaults,
        read_table_kw=None,
        window=1,
        name=None,
    ):
        """Initializes the Interactions object."""
        self.window = window
        super().__init__(
            input_data=input_data,
            sequence=sequence,
            metric=metric,
            metric_defaults=metric_defaults,
            read_table_kw=read_table_kw,
            name=name,
        )
        self.data = self.data.astype({"i": "int", "j": "int"})
        self.reset_mask()

    def mask_on_sequence(self, compliment_only=None, nts=None):
        """Mask interactions based on sequence.

        Parameters
        ----------
        compliment_only : bool, defaults to None
            If True, only keep interactions where i and j are complimentary
            nucleotides.
        nts : str, defaults to None
            If compliment_only is False, only keep interactions where i and j
            are in nts.

        Returns
        -------
        numpy array
            a boolean array of the same length as self.data
        """
        mask = []
        compliment = {"A": "UT", "U": "AG", "T": "AG", "G": "CU", "C": "G"}
        for _, i, j in self.data[["i", "j"]].itertuples():
            keep = []
            for w in range(self.window):
                i_nt = self.sequence[i + w - 1].upper()
                j_nt = self.sequence[j + self.window - w - 2].upper()
                if compliment_only:
                    keep.append(i_nt in compliment[j_nt])
                if nts is not None:
                    keep.append(i_nt in nts and j_nt in nts)
            mask.append(all(keep))
        self.update_mask(mask)
        return mask

    def mask_on_structure(
        self,
        structure,
        min_cd=None,
        max_cd=None,
        ss_only=False,
        ds_only=False,
        paired_only=False,
    ):
        """Masks interactions based on a secondary structure.

        Parameters
        ----------
        structure : rnavigate.data.SecondaryStructure
            The secondary structure to use for masking.
        min_cd : int, defaults to None
            The minimum contact distance to allow.
        max_cd : int, defaults to None
            The maximum contact distance to allow.
        ss_only : bool, defaults to False
            If True, only keep interactions between single-stranded nucleotides.
        ds_only : bool, defaults to False
            If True, only keep interactions between double-stranded nucleotides.
        paired_only : bool, defaults to False
            If True, only keep interactions that are paired in the structure.

        Returns
        -------
        mask : numpy array
            a boolean array of the same length as self.data
        """
        if isinstance(structure, list):
            for each in structure:
                self.mask_on_structure(
                    each, min_cd, max_cd, ss_only, ds_only, paired_only
                )
                return
        alignment = data.SequenceAlignment(self, structure)
        mask = np.full(len(self.data), True)
        i_j_keep = ["i", "j", "mask"]
        for idx, i, j, keep in self.data[i_j_keep].itertuples():
            i = alignment.map_positions(i)
            j = alignment.map_positions(j)
            true_so_far = keep and (-1 not in [i, j])
            if paired_only and true_so_far:
                for w in range(self.window):
                    true_so_far = (
                        structure.pair_nts[i + w - 1] == j + self.window - w - 1
                    )
                    if not true_so_far:
                        break
            if (ss_only or ds_only) and true_so_far:
                percentage = 0
                for w in range(self.window):
                    percentage += int(structure.pair_nts[i - 1 + w] == 0)
                    percentage += int(structure.pair_nts[j - 1 + w] == 0)
                percentage /= self.window * 2
                if ss_only:
                    true_so_far = percentage > 0.501
                if ds_only:
                    true_so_far = percentage < 0.501
            if (min_cd is not None or max_cd is not None) and true_so_far:
                cd = []
                for iw in range(self.window):
                    for jw in range(self.window):
                        cd.append(structure.contact_distance(i + iw, j + jw))
                cd = min(cd)
                if min_cd is not None:
                    true_so_far = cd >= min_cd
                if max_cd is not None:
                    true_so_far = cd <= max_cd
            mask[idx] &= true_so_far
        self.update_mask(mask)
        return mask

    def mask_on_profile(self, profile, min_profile=None, max_profile=None):
        """Masks interactions based on per-nucleotide measurements.

        Parameters
        ----------
        profile : rnavigate.data.Profile
            The profile to use for masking.
        min_profile : float, defaults to None
            The minimum profile value to allow.
        max_profile : float, defaults to None
            The maximum profile value to allow.

        Returns
        -------
        mask : numpy array
            a boolean array of the same length as self.data
        """
        mapping = data.SequenceAlignment(self, profile).mapping
        norm_prof = profile.data["Norm_profile"]
        mask = np.full(len(self.data), True)
        for idx, i, j in self.data[["i", "j"]].itertuples():
            index_i = mapping[i - 1]
            index_j = mapping[j - 1]
            keep_ij = True
            if (index_i != -1) and (index_j != -1):
                prof_i = np.nanmedian(norm_prof[index_i : index_i + self.window])
                prof_j = np.nanmedian(norm_prof[index_j : index_j + self.window])
                if min_profile is not None:
                    keep_ij &= (prof_i >= min_profile) | np.isnan(prof_i)
                    keep_ij &= (prof_j >= min_profile) | np.isnan(prof_j)
                if max_profile is not None:
                    keep_ij &= (prof_i <= max_profile) | np.isnan(prof_i)
                    keep_ij &= (prof_j <= max_profile) | np.isnan(prof_j)
            mask[idx] &= keep_ij
        self.update_mask(mask)
        return mask

    def mask_on_position(self, exclude=None, isolate=None):
        """Mask interactions based on their i and j positions.

        Parameters
        ----------
        exclude : list of int, defaults to None
            A list of positions to exclude.
        isolate : list of int, defaults to None
            A list of positions to isolate.

        Returns
        -------
        mask : numpy array
            a boolean array of the same length as self.data
        """
        mask = np.full(len(self.data), True)
        for index, (i, j) in self.data[["i", "j"]].iterrows():
            if exclude is not None:
                if (i in exclude) or (j in exclude):
                    mask[index] = False
            if isolate is not None:
                if (i not in isolate) and (j not in isolate):
                    mask[index] = False
        self.update_mask(mask)
        return mask

    def mask_on_distance(self, max_dist=None, min_dist=None):
        """Mask interactions based on their distance in sequence space.

        Parameters
        ----------
        max_dist : int, defaults to None
            The maximum distance to allow. If None, no maximum distance is set.
        min_dist : int, defaults to None
            The minimum distance to allow. If None, no minimum distance is set.

        Returns
        -------
        mask : numpy array
            a boolean array of the same length as self.data
        """
        primary_distances = self.data.eval("i - j").abs()
        mask = self.data["mask"]
        if min_dist is not None:
            mask &= primary_distances >= min_dist
        if max_dist is not None:
            mask &= primary_distances <= max_dist
        self.update_mask(mask)
        return mask

    def mask_on_values(self, **kwargs):
        """Mask interactions based on values in self.data.

        Parameters
        ----------
        kwargs : dict
            Each keyword should have the format "column_operator" where column
            is a valid column name of the dataframe and operator is one of:
                "ge": greater than or equal to
                "le": less than or equal to
                "gt": greater than
                "lt": less than
                "eq": equal to
                "ne": not equal to
            The values given to these keywords are then used in the comparison
            and False comparisons are filtered out. e.g.:
                self.mask_on_values(Statistic_ge=23) evaluates to:
                self.update_mask(self.data["Statistic"] >= 23)

        Returns
        -------
        mask : numpy array
            a boolean array of the same length as self.data
        """
        mask = np.full(len(self.data), True)
        for key, value in kwargs.items():
            if key in self.data.keys():
                mask &= self.data[key] > value
            elif "_" in key:
                key2, comparison = key.rsplit("_", 1)
                operators = {"ge": ge, "le": le, "gt": gt, "lt": lt, "eq": eq, "ne": ne}
                if key2 in self.data.columns and comparison in operators:
                    operator = operators[comparison]
                    mask &= operator(self.data[key2], value)
                else:
                    print(f"{key}={value} is not a valid filter.")
            else:
                print(f"{key}={value} is not a valid filter.")
        self.update_mask(mask)
        return mask

    def reset_mask(self):
        """Resets the mask to all True (removes previous filters)"""
        self.data["mask"] = np.ones(len(self.data), dtype=bool)

    def copy(self, apply_filter=False):
        """Returns a copy of the interactions, optionally with masked rows removed.

        Parameters
        ----------
        apply_filter : bool, defaults to False
            If True, masked rows ("mask" == False) are dropped.

        Returns
        -------
        rnavigate.data.Interactions
            A copy of the interactions.
        """
        return self.get_aligned_data(
            alignment=self.null_alignment, apply_filter=apply_filter
        )

    def get_aligned_data(self, alignment, apply_filter=True):
        """Returns a copy mapped to a new sequence with masked rows removed.

        Parameters
        ----------
        alignment : rnavigate.data.SequenceAlignment
            The alignment to use for mapping the interactions.
        apply_filter : bool, defaults to True
            If True, masked rows ("mask" == False) are dropped.

        Returns
        -------
        rnavigate.data.Interactions
            Interactions mapped to a new sequence.
        """
        if apply_filter:
            dataframe = self.data[self.data["mask"]].reset_index(drop=True)
        else:
            dataframe = self.data
        new_data = alignment.map_dataframe(
            dataframe=dataframe, position_columns=["i", "j"]
        )
        return self.__class__(
            input_data=new_data,
            sequence=alignment.target_sequence,
            metric=self._metric,
            metric_defaults=self.metric_defaults,
            window=self.window,
            name=self.name,
        )

    def update_mask(self, mask):
        """Updates the mask by ANDing the current mask with the given mask."""
        self.data["mask"] = self.data["mask"] & mask

    def count_filter(self, **kwargs):
        """Counts the number of interactions that pass the given filters."""
        mask = self.filter(**kwargs)
        return sum(mask)

    def filter(
        self,
        prefiltered=False,
        reset_filter=True,
        # mask on structure
        structure=None,
        min_cd=None,
        max_cd=None,
        paired_only=False,
        ss_only=False,
        ds_only=False,
        # mask on profile
        profile=None,
        min_profile=None,
        max_profile=None,
        # mask on sequence
        compliments_only=False,
        nts=None,
        # mask on position
        max_distance=None,
        min_distance=None,
        exclude_nts=None,
        isolate_nts=None,
        # others
        resolve_conflicts=None,
        **kwargs,
    ):
        """Convenience function that applies the above filters simultaneously.

        Parameters
        ----------
        prefiltered : bool, defaults to False
            If True, the mask is not updated.
        reset_filter : bool, defaults to True
            If True, the mask is reset before applying filters.
        structure : rnavigate.data.SecondaryStructure, defaults to None
            The structure to use for filtering.
        min_cd : int, defaults to None
            The minimum contact distance to allow.
        max_cd : int, defaults to None
            The maximum contact distance to allow.
        paired_only : bool, defaults to False
            If True, only keep interactions that are paired in the structure.
        ss_only : bool, defaults to False
            If True, only keep interactions between single-stranded nucleotides.
        ds_only : bool, defaults to False
            If True, only keep interactions between double-stranded nucleotides.
        profile : rnavigate.data.Profile, defaults to None
            The profile to use for masking.
        min_profile : float, defaults to None
            The minimum profile value to allow.
        max_profile : float, defaults to None
            The maximum profile value to allow.
        compliments_only : bool, defaults to False
            If True, only keep interactions where i and j are complimentary
            nucleotides.
        nts : str, defaults to None
            If compliment_only is False, only keep interactions where i and j
            are in nts.
        max_distance : int, defaults to None
            The maximum distance to allow. If None, no maximum distance is set.
        min_distance : int, defaults to None
            The minimum distance to allow. If None, no minimum distance is set.
        exclude_nts : list of int, defaults to None
            A list of positions to exclude.
        isolate_nts : list of int, defaults to None
            A list of positions to isolate.
        resolve_conflicts : str, defaults to None
            If not None, conflicting windows are resolved using the Maximal
            Weighted Independent Set. The weights are taken from the metric
            value. The graph is first broken into components to speed up the
            identification of the MWIS. Then the mask is updated to only
            include the MWIS.
        **kwargs : dict
            Each keyword should have the format "column_operator" where column
            is a valid column name of the dataframe and operator is one of:
                "ge": greater than or equal to
                "le": less than or equal to
                "gt": greater than
                "lt": less than
                "eq": equal to
                "ne": not equal to
            The values given to these keywords are then used in the comparison
            and False comparisons are filtered out. e.g.:
                self.mask_on_values(Statistic_ge=23) evaluates to:
                self.update_mask(self.data["Statistic"] >= 23)

        Returns
        -------
        mask : numpy array
            a boolean array of the same length as self.data
        """

        def filters_are_on(*filters):
            return any(f not in [None, False] for f in filters)

        if prefiltered:
            return
        if reset_filter:
            self.reset_mask()
        mask = np.full(len(self.data), True)
        if filters_are_on(exclude_nts, isolate_nts):
            mask &= self.mask_on_position(exclude=exclude_nts, isolate=isolate_nts)
        if filters_are_on(max_distance, min_distance):
            mask &= self.mask_on_on_distance(
                max_dist=max_distance, min_dist=min_distance
            )
        if filters_are_on(min_profile, max_profile):
            mask &= self.mask_on_profile(
                profile=profile, min_profile=min_profile, max_profile=max_profile
            )
        if filters_are_on(min_cd, max_cd, ss_only, ds_only, paired_only):
            mask &= self.mask_on_structure(
                structure=structure,
                min_cd=min_cd,
                max_cd=max_cd,
                ss_only=ss_only,
                ds_only=ds_only,
                paired_only=paired_only,
            )
        if filters_are_on(compliments_only, nts):
            mask &= self.mask_on_sequence(compliment_only=compliments_only, nts=nts)
        kwargs, data_specific_mask = self.data_specific_filter(**kwargs)
        mask &= data_specific_mask
        mask &= self.mask_on_values(**kwargs)
        if resolve_conflicts is not None:
            mask &= self.resolve_conflicts(metric=resolve_conflicts)
        return mask

    def data_specific_filter(self, **kwargs):
        """Does nothing for the base Interactions class, can be overwritten in
        subclasses.

        Returns:
            dict: dictionary of keyword argument pairs
        """
        return kwargs, np.full(len(self.data), True)

    def get_ij_colors(self):
        """Gets i, j, and colors lists for plotting interactions.

        i and j are the 5' and 3' ends of each interaction, and colors is the color
        to use for each interaction. Values of self.data[self.metric] are normalized
        to 0 to 1, which correspond to self.min_max values. These are then mapped to
        a color using self.cmap.

        Returns
        -------
        i : list
            5' ends of each interaction
        j : list
            3' ends of each interaction
        colors : list
            colors to use for each interaction
        """
        if len(self.data) == 0:
            return [], [], []

        dataframe = self.get_sorted_data()
        i_list, j_list, colors = [], [], []
        for _, i, j, datum in dataframe[["i", "j", self.metric]].itertuples():
            for w in range(self.window):
                i_list.append(i + w)
                j_list.append(j + self.window - 1 - w)
                colors.append(datum)
        colors = self.cmap.values_to_hexcolors(colors, 0.6)
        return i_list, j_list, colors

    def get_sorted_data(self):
        """Returns a copy of the data sorted by self.metric.

        Returns
        -------
        pandas.DataFrame
            a copy of the data sorted by self.metric
        """
        ascending = self.metric not in ["Distance"]  # high distance is bad
        dataframe = self.data[["i", "j", self.metric]].copy()
        dataframe.sort_values(
            by=self.metric, ascending=ascending, inplace=True, na_position="first"
        )
        return dataframe

    def print_new_file(self, outfile=None):
        """Create a new file with mapped and filtered interactions.

        Parameters
        ----------
        outfile : str, defaults to None
            path to an output file. If None, file string is printed to console.
        """
        data = self.data.copy()
        data = data[data["mask"]]
        if "Sign" in data.columns:
            data.rename({"Sign": "+/-"}, inplace=True)
        csv = data.to_csv(
            columns=self.columns, sep="\t", index=False, line_terminator="\n"
        )
        if outfile is not None:
            with open(outfile, "w") as out:
                out.write(self.header)
                out.write(csv)
        else:
            print(self.header, csv)

    def set_3d_distances(self, pdb, atom):
        """Calculates the distance between atoms in i and j in the PDB structure.

        Parameters
        ----------
        pdb : rnavigate.pdb.PDB
            PDB object to use for calculating distances
        atom : str
            atom id to use for calculating distances
        """
        mapping = data.SequenceAlignment(self, pdb).mapping
        distance_matrix = pdb.get_distance_matrix(atom=atom)
        i = self.data["i"].values
        j = self.data["j"].values
        distances = np.zeros(len(i))
        pairs = self.window**2
        for iw in range(self.window):
            for jw in range(self.window):
                io = mapping[i + iw - 1]
                jo = mapping[j + jw - 1]
                distances += distance_matrix[io, jo] / pairs
        self.data["Distance"] = distances

    def resolve_conflicts(self, metric=None):
        """Uses an experimental method to resolve conflicts.

        Resolves conflicting windows using the Maximal Weighted Independent
        Set. The weights are taken from the metric value. The graph is first
        broken into components to speed up the identification of the MWIS. Then
        the mask is updated to only include the MWIS. This method is computationally
        expensive for large or dense datasets.

        Parameters
        ----------
        metric : str, defaults to None
            The metric to use for weighting the graph. If None, self.metric is used.

        Returns
        -------
        mask : numpy array
            a boolean array of the same length as self.data
        """

        def get_components(graph):
            """Using DFS, returns a list of component subgraphs of graph."""

            def dfs(key, value, component):
                """Depth first seach to identify a component."""
                visited[key] = True
                component[key] = value
                for v in value:
                    if not visited[v]:
                        component = dfs(v, graph[v], component)
                return dict(component)

            components = []
            visited = {key: False for key in graph.keys()}
            for key, value in graph.items():
                if not visited[key]:
                    components.append(dfs(key, value, {}))
            return components

        def graphSets(graph, node_weights):
            """Finds the Maximal Weighted Independent Set of a graph."""
            # Base Case - Given Graph has no nodes
            if len(graph) == 0:
                return []
            # Select a vertex from the graph
            current_vertex = list(graph.keys())[0]
            # Create a copy of the graph and remove the current vertex
            new_graph = dict(graph)
            del new_graph[current_vertex]
            # Case 1 - excluding current Vertex
            case_1 = graphSets(new_graph, node_weights)
            # Case 2 - including current Vertex, and excluding neighbors
            for neighbor in graph[current_vertex]:
                if neighbor in new_graph:
                    del new_graph[neighbor]
            case_2 = [current_vertex] + graphSets(new_graph, node_weights)
            # Our final result is the one which has highest weight, return it
            weight_1 = sum(node_weights[node] for node in case_1)
            weight_2 = sum(node_weights[node] for node in case_2)
            if weight_1 > weight_2:
                return case_1
            return case_2

        # Building the graph ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if metric is None:
            metric = self.metric
        # a list of indices relating back to the original mask
        indices = np.where(self.data["mask"])[0]
        # the correlation data to be considered
        df = self.data.loc[indices].copy()
        window = self.window
        # node_weight[i from indices] = metric value for that correlation
        node_weights = {}
        # graph[i from indices] = list of indices of conflicting correlations
        graph = {}
        # for each considered correlation, record conflicts and node weights
        for v1, i, j, g in df[["i", "j", metric]].itertuples():
            if v1 not in graph:
                graph[v1] = []
            node_weights[v1] = g  # recording the node weight
            # conflicts = indices of all overlapping, but not parallel corrs
            conflicts = np.where(
                (
                    (abs(df["i"] - i) < window)
                    | (abs(df["j"] - i) < window)
                    | (abs(df["i"] - j) < window)
                    | (abs(df["j"] - j) < window)
                )
                & ((df["j"] - j) != (i - df["i"]))
            )[0]
            # adding conflicts to graph dictionary
            for conflict in conflicts:
                v2 = indices[conflict]
                if v2 > v1:  # prevents self-loop and parallel edges
                    if v2 not in graph:
                        graph[v2] = []
                    graph[v1].append(v2)
                    graph[v2].append(v1)
        # Finding the maximum set ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Split the graph into components. Get the max set for each component.
        # The graph's max set is the same as the union of these max sets.
        max_set = []
        for component in get_components(graph):
            max_set += graphSets(component, node_weights)
        # set these indices to true and update the original mask
        new_mask = np.zeros(len(self.data), dtype=bool)
        new_mask[max_set] = True
        self.update_mask(new_mask)
        return new_mask


class SHAPEJuMP(Interactions):
    """A class for storing and manipulating SHAPEJuMP data.

    Parameters
    ----------
    input_data : string or pandas.DataFrame
        If string, a path to a file containing SHAPEJuMP data.
        If dataframe, the dataframe containing SHAPEJuMP data. The dataframe
        must contain columns "i", "j", "Metric" (JuMP rate) and "Percentile"
        (percentile ranking). Dataframe may also include other columns.
    sequence : string or rnavigate.data.Sequence
        The sequence string corresponding to the SHAPEJuMP data.
    metric : string, defaults to "Percentile"
        The column name to use for visualization.
    metric_defaults : dict
        Keys are metric names and values are dictionaries of metric-specific defaults.
        These defaults include:
            "metric_column" : string
                the column name to use for visualization
            "cmap" : string or matplotlib.colors.Colormap)
                the colormap to use for visualization
            "normalization" : "min_max", "0_1", "none", or "bins"
                The type of normalization to use when mapping values to colors
            "values" : list of float
                The values to used with normalization of the data
            "title" : string
                the title to use for colorbars
            "extend" : "min", "max", "both", or "neither"
                Which ends to extend when drawing the colorbar.
            "tick_labels" : list of string
    read_table_kw : dict
        kwargs passed to pandas.read_table() when reading input_data.
    window : int
        The window size used to generate the SHAPEJuMP data.
    name : str
        A name for the interactions object.

    Attributes
    ----------
    data : pandas.DataFrame
        The SHAPEJuMP data.
    """

    def __init__(
        self,
        input_data,
        sequence=None,
        metric="Percentile",
        metric_defaults=None,
        read_table_kw=None,
        window=1,
        name=None,
    ):
        """Constructs an Interactions object from SHAPEJuMP data"""
        if metric_defaults is None:
            metric_defaults = {}
        metric_defaults = {
            "Percentile": {
                "metric_column": "Percentile",
                "cmap": "YlGnBu",
                "normalization": "min_max",
                "values": [0.98, 1.0],
                "title": "SHAPE-JuMP: percentile",
                "extend": "min",
            },
            "Metric": {
                "metric_column": "Metric",
                "cmap": "YlGnBu",
                "normalization": "min_max",
                "values": [0, 0.001],
                "extend": "max",
                "title": "SHAPE-JuMP: rate",
            },
        } | metric_defaults
        super().__init__(
            input_data=input_data,
            sequence=sequence,
            metric=metric,
            metric_defaults=metric_defaults,
            read_table_kw=read_table_kw,
            window=window,
            name=name,
        )

    def read_file(self, input_data, read_table_kw=None):
        """Parses a deletions.txt file and stores it as a dataframe.

        Also calculates a "Percentile" column.

        Parameters
        ----------
        input_data : str
            path to deletions.txt file
        read_table_kw : dict, defaults to {}
            kwargs passed to pandas.read_table().

        Returns
        -------
        pandas.DataFrame
            the SHAPEJuMP data
        """
        column_names = ["Gene", "i", "j", "Metric"]
        data = pd.read_table(input_data, names=column_names, header=0, **read_table_kw)
        data["Percentile"] = data["Metric"].rank(method="max", pct=True)
        return data


class RINGMaP(Interactions):
    """A class for storing and manipulating RINGMaP data.

    Parameters
    ----------
    input_data : string or pandas.DataFrame
        If string, a path to a file containing RINGMaP data.
        If dataframe, the dataframe containing RINGMaP data. The dataframe
        must contain columns "i", "j", "Statistic", and "Zij". Dataframe may
        also include other columns.
    sequence : string or rnavigate.data.Sequence
        The sequence string corresponding to the RINGMaP data.
    metric : string, defaults to "Statistic"
        The column name to use for visualization.
    metric_defaults : dict
        Keys are metric names and values are dictionaries of metric-specific defaults.
        These defaults include:
            "metric_column" : string
                the column name to use for visualization
            "cmap" : string or matplotlib.colors.Colormap)
                the colormap to use for visualization
            "normalization" : "min_max", "0_1", "none", or "bins"
                The type of normalization to use when mapping values to colors
            "values" : list of float
                The values to used with normalization of the data
            "title" : string
                the title to use for colorbars
            "extend" : "min", "max", "both", or "neither"
                Which ends to extend when drawing the colorbar.
            "tick_labels" : list of string
    read_table_kw : dict, optional
        kwargs passed to pandas.read_table() when reading input_data.
    window : int, defaults to 1
        The window size used to generate the RINGMaP data. If an input file is
        provided, this value is overwritten by the value in the header.
    name : str, optional
        A name for the interactions object.

    Attributes
    ----------
    data : pandas.DataFrame
        The RINGMaP data.
    """

    def __init__(
        self,
        input_data,
        sequence=None,
        metric="Statistic",
        metric_defaults=None,
        read_table_kw=None,
        window=1,
        name=None,
    ):
        """Constructs an Interactions object from RINGMaP data"""
        if metric_defaults is None:
            metric_defaults = {}
        metric_defaults = {
            "Statistic": {
                "metric_column": "Statistic",
                "cmap": "bwr",
                "normalization": "min_max",
                "values": [-100, 100],
                "title": "RING-MaP: Gapc",
                "extend": "both",
            },
            "Zij": {
                "metric_column": "Zij",
                "cmap": "bwr",
                "normalization": "min_max",
                "values": [-8, 8],
                "title": "RING-MaP: Zij",
                "extend": "both",
            },
        } | metric_defaults
        super().__init__(
            input_data=input_data,
            sequence=sequence,
            metric=metric,
            metric_defaults=metric_defaults,
            read_table_kw=read_table_kw,
            window=window,
            name=name,
        )

    def read_file(self, filepath, read_table_kw=None):
        """Parses a RINGMaP correlations file and stores data as a dataframe.

        Also sets self.window (usually 1, from header).

        Parameters
        ----------
        filepath : str
            path to correlations file.
        read_table_kw : dict, defaults to {}
            kwargs passed to pandas.read_table().

        Returns
        -------
        pandas.DataFrame
            the RINGMaP data
        """
        with open(filepath, "r") as file:
            self.header = file.readline()
        split_header = self.header.split("\t")
        window = split_header[1].split("=")[1]
        self.window = int(window)
        dataframe = pd.read_table(filepath, header=1, **read_table_kw)
        dataframe.rename(columns={"+/-": "Sign"}, inplace=True)
        return dataframe

    def data_specific_filter(self, positive_only=False, negative_only=False, **kwargs):
        """Adds filters for "Sign" column to parent filter() function

        Parameters
        ----------
        positive_only : bool, defaults to False
            If True, only keep positive correlations.
        negative_only : bool, defaults to False
            If True, only keep negative correlations.

        Returns
        -------
        kwargs : dict
            any additional keyword-argument pairs are returned
        mask : numpy array
            a boolean array of the same length as self.data
        """
        mask = np.full(len(self.data), True)
        if positive_only:
            mask &= self.data["Sign"] == 1
        if negative_only:
            mask &= self.data["Sign"] == -1
        self.update_mask(mask)
        return kwargs, mask

    def get_sorted_data(self):
        """Sorts on the product of self.metric and "Sign" columns.

        Except when self.metric is "Distance".

        Returns
        -------
        pandas.DataFrame
            a copy of the data sorted by (self.metric * "Sign") columns
        """
        if self.metric == "Distance":
            return super().get_sorted_data()
        columns = ["i", "j", self.metric]
        dataframe = self.data[columns + ["Sign"]].copy()
        dataframe.eval(f"{self.metric} = {self.metric} * Sign", inplace=True)
        dataframe.sort_values(by=self.metric, ascending=True, inplace=True)
        return dataframe[columns]


class PAIRMaP(RINGMaP):
    """A class for storing and manipulating PAIRMaP data.

    Parameters
    ----------
    input_data : string or pandas.DataFrame
        If string, a path to a file containing PAIRMaP data.
        If dataframe, the dataframe containing PAIRMaP data. The dataframe
        must contain columns "i", "j", "Statistic", and "Class". Dataframe may
        also include other columns.
    sequence : string or rnavigate.data.Sequence
        The sequence string corresponding to the PAIRMaP data.
    metric : string, defaults to "Class"
        The column name to use for visualization.
    metric_defaults : dict
        Keys are metric names and values are dictionaries of metric-specific defaults.
        These defaults include:
            "metric_column" : string
                the column name to use for visualization
            "cmap" : string or matplotlib.colors.Colormap)
                the colormap to use for visualization
            "normalization" : "min_max", "0_1", "none", or "bins"
                The type of normalization to use when mapping values to colors
            "values" : list of float
                The values to used with normalization of the data
            "title" : string
                the title to use for colorbars
            "extend" : "min", "max", "both", or "neither"
                Which ends to extend when drawing the colorbar.
            "tick_labels" : list of string
    read_table_kw : dict, optional
        kwargs passed to pandas.read_table() when reading input_data.
    window : int, defaults to 1
        The window size used to generate the PAIRMaP data. If an input file is
        provided, this value is overwritten by the value in the header.
    name : str, optional
        A name for the interactions object.

    Attributes
    ----------
    data : pandas.DataFrame
        The PAIRMaP data.
    """

    def __init__(
        self,
        input_data,
        sequence=None,
        metric="Class",
        metric_defaults=None,
        read_table_kw=None,
        window=1,
        name=None,
    ):
        """Constructs an Interactions object from PAIRMaP data"""
        if metric_defaults is None:
            metric_defaults = {}
        metric_defaults = {
            "Class": {
                "metric_column": "Class",
                "cmap": mpc.ListedColormap(
                    [[0.7, 0.7, 0.7], [0.0, 0.0, 0.95], [0.12, 0.76, 1.0]]
                ),
                "normalization": "none",
                "values": None,
                "title": "PAIR-MaP",
                "extend": "both",
                "ticks": [0, 1, 2],
                "tick_labels": ["Complimentary", "Primary", "Secondary"],
            }
        } | metric_defaults
        super().__init__(
            input_data=input_data,
            sequence=sequence,
            metric=metric,
            metric_defaults=metric_defaults,
            read_table_kw=read_table_kw,
            window=window,
            name=name,
        )

    def read_file(self, filepath, read_table_kw=None):
        """Parses a pairmap.txt file and stores data as a dataframe

        Sets self.window (usually 3), from header.

        Parameters
        ----------
        filepath : str
            path to pairmap.txt file
        read_table_kw : dict, defaults to None
            This argument is ignored.
        """
        with open(filepath, "r") as file:
            self.header = file.readline()
        self.window = int(self.header.split("\t")[1].split("=")[1])
        dataframe = pd.read_table(filepath, header=1)
        dataframe.rename(columns={"Sig.": "Statistic"}, inplace=True)
        dataframe["Sign"] = 1
        return dataframe

    def data_specific_filter(self, all_pairs=False, **kwargs):
        """Used by Interactions.filter(). By default, non-primary and
        -secondary pairs are removed. all_pairs=True changes this behavior.

        Parameters
        ----------
        all_pairs : bool, defaults to False
            whether to include all PAIRs.

        Returns
        -------
        kwargs : dict
            any additional keyword-argument pairs are returned
        mask : numpy array
            a boolean array of the same length as self.data
        """
        mask = np.full(len(self.data), True)
        if not all_pairs:
            mask &= self.data["Class"] != 0
        self.update_mask(mask)
        return kwargs, mask

    def get_sorted_data(self):
        """Same as parent function, unless metric is set to "Class", in which
        case ij pairs are returned in a different order.

        Returns
        -------
        pandas.DataFrame
            a copy of the data sorted by self.metric
        """
        if self.metric != "Class":
            return super().get_sorted_data()
        metric = self.metric
        columns = ["i", "j", metric]
        dataframe = self.data[columns].copy()
        # Class data have a weird order
        dataframe = pd.concat(
            [
                dataframe[dataframe[metric] == 0],
                dataframe[dataframe[metric] == 2],
                dataframe[dataframe[metric] == 1],
            ]
        )
        return dataframe


class PairingProbability(Interactions):
    """A class for storing and manipulating pairing probability data.

    Parameters
    ----------
    input_data : string or pandas.DataFrame
        If string, a path to a file containing pairing probability data.
        If dataframe, the dataframe containing pairing probability data. The
        dataframe must contain columns "i", "j", "Probability", and "log10p".
        Dataframe may also include other columns.
    sequence : string or rnavigate.data.Sequence
        The sequence string corresponding to the pairing probability data.
    metric : string, defaults to "Probability"
        The column name to use for visualization.
    metric_defaults : dict
        Keys are metric names and values are dictionaries of metric-specific defaults.
        These defaults include:
            "metric_column" : string
                the column name to use for visualization
            "cmap" : string or matplotlib.colors.Colormap
                the colormap to use for visualization
            "normalization" : "min_max", "0_1", "none", or "bins"
                The type of normalization to use when mapping values to colors
            "values" : list of float
                The values to used with normalization of the data
            "title" : string
                the title to use for colorbars
            "extend" : "min", "max", "both", or "neither"
                Which ends to extend when drawing the colorbar.
            "tick_labels" : list of string
    read_table_kw : dict, optional
        kwargs passed to pandas.read_table() when reading input_data.
    window : int, defaults to 1
        The window size used to generate the pairing probability data.
    name : str, optional
        A name for the PairingProbability object.

    Attributes
    ----------
    data : pandas.DataFrame
        The pairing probability data.
    """

    def __init__(
        self,
        input_data,
        sequence=None,
        metric="Probability",
        metric_defaults=None,
        read_table_kw=None,
        window=1,
        name=None,
    ):
        """Constructs an Interactions object from pairing probability data"""
        if metric_defaults is None:
            metric_defaults = {}
        metric_defaults = {
            "Probability": {
                "metric_column": "Probability",
                "cmap": sns.cubehelix_palette(
                    10, 0.7, 0.9, 1.5, 2.5, 1, 0.4, False, True
                ),
                "normalization": "min_max",
                "values": [0, 1],
                "title": "Pairing probability",
                "extend": "neither",
            },
            "Probability_old": {
                "metric_column": "Probability",
                "cmap": mpc.ListedColormap(
                    [(150, 150, 150), (255, 204, 0), (72, 143, 205), (81, 184, 72)]
                ),
                "normalization": "bins",
                "values": [0.1, 0.3, 0.8],
                "extend": "neither",
                "title": "Pairing probability",
            },
            "Probability_continuous": {
                "metric_column": "Probability",
                "cmap": "plasma_r",
                "normalization": "min_max",
                "values": [0.0, 1.0],
                "extend": "neither",
                "title": "Pairing probability",
            },
        } | metric_defaults
        super().__init__(
            input_data=input_data,
            sequence=sequence,
            metric=metric,
            metric_defaults=metric_defaults,
            read_table_kw=read_table_kw,
            window=window,
            name=name,
        )

    def read_file(self, filepath, read_table_kw=None):
        """Parses a pairing probability file and stores data as a dataframe.

        Parameters
        ----------
        filepath : str
            path to pairing probability file
        read_table_kw : dict, defaults to None
            This argument is ignored.

        Returns
        -------
        pandas.DataFrame
            the pairing probability data
        """
        with open(filepath, "r") as file:
            self.header = file.readline()
        self.window = 1
        dataframe = pd.read_table(filepath, header=1, names=["i", "j", "log10p"])
        dataframe["Probability"] = 10 ** (-dataframe["log10p"])
        return dataframe

    def get_entropy_profile(self, print_out=False, save_file=None):
        """Calculates per-nucleotide Shannon entropy from pairing probabilities.

        Parameters
        ----------
        print_out : bool, defaults to False
            If True, entropy values are printed to console.
        save_file : str, defaults to None
            If not None, entropy values are saved to this file.

        Returns
        -------
        rnavigate.data.Profile
            a Profile object containing the entropy data
        """
        self.data.eval("nlogn = log10p * 10 ** ( - log10p )", inplace=True)
        entropy = np.zeros(self.length)
        args = {"labels": np.arange(self.length) + 1, "fill_value": 0}
        i_sum = self.data[["i", "nlogn"]].groupby("i").sum().reindex(**args)
        j_sum = self.data[["j", "nlogn"]].groupby("j").sum().reindex(**args)
        ij_sum = i_sum + j_sum
        entropy_df = pd.DataFrame(
            {
                "Sequence": list(self.sequence),
                "Nucleotide": np.arange(self.length) + 1,
                "Entropy": ij_sum["nlogn"],
            }
        )
        entropy_df = entropy_df.set_index("Nucleotide")
        # catch rounding errors:
        entropy_df.loc[entropy_df["Entropy"] < 0, "Entropy"] = 0
        if print_out:
            print(*[f"{i+1} {s}" for i, s in enumerate(entropy)], sep="\n")
        if save_file is not None:
            with open(save_file) as outf:
                for i, s in enumerate(entropy):
                    outf.write(f"{i+1}\t{s}\n")
        return data.Profile(
            input_data=entropy_df,
            metric="Entropy",
            metric_defaults={
                "Entropy": {
                    "metric_column": "Entropy",
                    "error_column": None,
                    "color_column": None,
                    "cmap": "rainbow",
                    "normalization": "norm",
                    "values": [0.0, 0.2],
                    "extend": "right",
                    "title": "Shannon entropy",
                    "alpha": 1,
                }
            },
        )

    def data_specific_filter(self, **kwargs):
        """By default, interactions with probabilities less than 0.03 are removed.

        Returns
        -------
        kwargs : dict
            any additional keyword-argument pairs are returned
        mask : numpy array
            a boolean array of the same length as self.data
        """
        mask = self.data["Probability"] >= 0.03
        self.update_mask(mask)
        return kwargs, mask


class AllPossible(Interactions):
    """A class for storing and manipulating all possible interactions.

    Parameters
    ----------
    sequence : string or rnavigate.data.Sequence
        The sequence string corresponding to the pairing probability data.
    metric : string, defaults to "Probability"
        The column name to use for visualization.
    metric_defaults : dict
        Keys are metric names and values are dictionaries of metric-specific defaults.
        These defaults include:
            "metric_column" : string
                the column name to use for visualization
            "cmap" : string or matplotlib.colors.Colormap
                the colormap to use for visualization
            "normalization" : "min_max", "0_1", "none", or "bins"
                The type of normalization to use when mapping values to colors
            "values" : list of float
                The values to used with normalization of the data
            "title" : string
                the title to use for colorbars
            "extend" : "min", "max", "both", or "neither"
                Which ends to extend when drawing the colorbar.
            "tick_labels" : list of string
    read_table_kw : dict, optional
        kwargs passed to pandas.read_table() when reading input_data.
    window : int, defaults to 1
        The window size used to generate the pairing probability data.
    name : str, optional
        A name for the PairingProbability object.

    Attributes
    ----------
    data : pandas.DataFrame
        The pairing probability data.
    """

    def __init__(
        self,
        sequence,
        metric="data",
        input_data=None,
        metric_defaults=None,
        read_table_kw=None,
        window=1,
        name=None,
    ):
        """Constructs an Interactions object from pairing probability data"""
        if isinstance(sequence, data.Sequence):
            sequence = sequence.sequence
        if metric_defaults is None:
            metric_defaults = {}
        metric_defaults = {
            "data": {
                "metric_column": "data",
                "cmap": "magenta",
                "normalization": "none",
                "values": None,
                "title": "Hypthetical pairs",
                "extend": "neither",
            }
        } | metric_defaults
        if input_data is not None:
            dataframe = input_data
        else:
            dataframe = {"i": [], "j": [], "data": []}
            for i in range(1, len(sequence)):
                for j in range(i + 1, len(sequence) + 1):
                    dataframe["i"].append(i)
                    dataframe["j"].append(j)
                    dataframe["data"].append(0)
            dataframe = pd.DataFrame(dataframe)
        super().__init__(
            input_data=dataframe,
            sequence=sequence,
            metric=metric,
            metric_defaults=metric_defaults,
            read_table_kw=read_table_kw,
            window=window,
            name=name,
        )


class StructureAsInteractions(Interactions):
    """A class for storing and manipulating structure data.

    Parameters
    ----------
    input_data : string or pandas.DataFrame
        If string, a path to a file containing structure data.
        If dataframe, the dataframe containing structure data. The dataframe
        must contain columns "i", "j", and "Structure". Dataframe may also
        include other columns.
    sequence : string or rnavigate.data.Sequence
        The sequence string corresponding to the structure data.
    metric : string, defaults to "Structure"
        The column name to use for visualization.
    metric_defaults : dict
        Keys are metric names and values are dictionaries of metric-specific defaults.
        These defaults include:
            "metric_column" : string
                the column name to use for visualization
            "cmap" : string or matplotlib.colors.Colormap
                the colormap to use for visualization
            "normalization" : "min_max", "0_1", "none", or "bins"
                The type of normalization to use when mapping values to colors
            "values" : list of float
                The values to used with normalization of the data
            "title" : string
                the title to use for colorbars
            "extend" : "min", "max", "both", or "neither"
                Which ends to extend when drawing the colorbar.
            "tick_labels" : list of string
    read_table_kw : dict, optional
        kwargs passed to pandas.read_table() when reading input_data.
    window : int, defaults to 1
        The window size used to generate the structure data.
    name : str, optional
        A name for the StructureAsInteractions object.

    Attributes
    ----------
    data : pandas.DataFrame
        The structure data.
    """

    def __init__(
        self,
        input_data,
        sequence,
        metric=None,
        metric_defaults=None,
        window=1,
        name=None,
    ):
        """Constructs an Interactions object from structure data"""
        if metric_defaults is None:
            metric_defaults = {}
        if metric is None:
            metric = "Structure"
        if isinstance(input_data, pd.DataFrame):
            pass
        elif isinstance(input_data, data.SecondaryStructure):
            input_data = input_data.get_interactions_df()
        metric_defaults = {
            "Structure": {
                "metric_column": "Structure",
                "cmap": "grey",
                "normalization": "none",
                "ticks": [],
                "title": "Base-pairs",
                "extend": "neither",
            }
        } | metric_defaults
        super().__init__(
            input_data=input_data,
            sequence=sequence,
            metric=metric,
            metric_defaults=metric_defaults,
            window=window,
            name=name,
        )


class StructureCompareTwo(Interactions):
    """A class for storing and manipulating a comparison of two structures.

    Parameters
    ----------
    input_data : string or pandas.DataFrame
        If string, a path to a file containing structure data.
        If dataframe, the dataframe containing structure data. The dataframe
        must contain columns "i", "j", and "Structure". Dataframe may also
        include other columns.
    sequence : string or rnavigate.data.Sequence
        The sequence string corresponding to the structure data.
    metric : string, defaults to "Structure"
        The column name to use for visualization.
    metric_defaults : dict
        Keys are metric names and values are dictionaries of metric-specific defaults.
        These defaults include:
            "metric_column" : string
                the column name to use for visualization
            "cmap" : string or matplotlib.colors.Colormap
                the colormap to use for visualization
            "normalization" : "min_max", "0_1", "none", or "bins"
                The type of normalization to use when mapping values to colors
            "values" : list of float
                The values to used with normalization of the data
            "title" : string
                the title to use for colorbars
            "extend" : "min", "max", "both", or "neither"
                Which ends to extend when drawing the colorbar.
            "tick_labels" : list of string
    read_table_kw : dict, optional
        kwargs passed to pandas.read_table() when reading input_data.
    window : int, defaults to 1
        The window size used to generate the structure data.
    name : str, optional
        A name for the StructureAsInteractions object.

    Attributes
    ----------
    data : pandas.DataFrame
        The structure data.
    """

    def __init__(
        self,
        input_data,
        sequence,
        metric=None,
        metric_defaults=None,
        window=1,
        name=None,
    ):
        """Constructs an Interactions object from structure data"""
        if metric is None:
            metric = "Which_structure"
        if metric_defaults is None:
            metric_defaults = {}
        if isinstance(input_data, pd.DataFrame):
            pass
        elif isinstance(input_data, list) and len(input_data) == 2:
            ss1, ss2 = [ss.get_interactions_df() for ss in input_data]
            input_data = ss1.merge(
                ss2,
                how="outer",
                on=["i", "j"],
                indicator="Which_structure",
                suffixes=["_left", "_right"],
            )
            categories = {"both": 0, "left_only": 1, "right_only": 2}
            input_data["Which_structure"] = [
                categories[c] for c in input_data["Which_structure"]
            ]
            input_data["Which_structure"].astype(int)
            input_data = input_data.reset_index(drop=True)
        metric_defaults = {
            "Which_structure": {
                "metric_column": "Which_structure",
                "cmap": [
                    # shared
                    (150 / 255.0, 150 / 255.0, 150 / 255.0),
                    # left
                    (38 / 255.0, 202 / 255.0, 145 / 255.0),
                    # right
                    (153 / 255.0, 0.0, 1.0),
                ],
                "normalization": "none",
                "values": [],
                "title": "Base-pairs by structure",
                "extend": "neither",
                "ticks": [0, 1, 2],
                "tick_labels": [
                    "common\npairs",
                    "first\nstructure",
                    "second\nstructure",
                ],
            }
        } | metric_defaults
        super().__init__(
            input_data=input_data,
            sequence=sequence,
            metric=metric,
            metric_defaults=metric_defaults,
            window=window,
            name=name,
        )


class StructureCompareMany(Interactions):
    """A class for storing and manipulating a comparison of many structures.

    Parameters
    ----------
    input_data : string or pandas.DataFrame
        If string, a path to a file containing structure data.
        If dataframe, the dataframe containing structure data. The dataframe
        must contain columns "i", "j", and "Structure". Dataframe may also
        include other columns.
    sequence : string or rnavigate.data.Sequence
        The sequence string corresponding to the structure data.
    metric : string, defaults to "Structure"
        The column name to use for visualization.
    metric_defaults : dict
        Keys are metric names and values are dictionaries of metric-specific defaults.
        These defaults include:
            "metric_column" : string
                the column name to use for visualization
            "cmap" : string or matplotlib.colors.Colormap
                the colormap to use for visualization
            "normalization" : "min_max", "0_1", "none", or "bins"
                The type of normalization to use when mapping values to colors
            "values" : list of float
                The values to used with normalization of the data
            "title" : string
                the title to use for colorbars
            "extend" : "min", "max", "both", or "neither"
                Which ends to extend when drawing the colorbar.
            "tick_labels" : list of string
    read_table_kw : dict, optional
        kwargs passed to pandas.read_table() when reading input_data.
    window : int, defaults to 1
        The window size used to generate the structure data.
    name : str, optional
        A name for the StructureAsInteractions object.

    Attributes
    ----------
    data : pandas.DataFrame
        The structure data.
    """

    def __init__(
        self,
        input_data,
        sequence,
        metric=None,
        metric_defaults=None,
        window=1,
        name=None,
    ):
        """Constructs an Interactions object from structure data"""
        if metric is None:
            metric = "Num_structures"
        if metric_defaults is None:
            metric_defaults = {}
        if isinstance(input_data, pd.DataFrame):
            pass
        elif isinstance(input_data, list) and len(input_data) >= 2:
            input_data = [ss.get_interactions_df() for ss in input_data]
            other_structures = input_data[1:]
            input_data = input_data[0]
            columns = ["Structure_1"]
            input_data = input_data.rename(columns={"Structure": "Structure_1"})
            for i, structure in enumerate(other_structures):
                col = f"Structure_{i+2}"
                columns.append(col)
                structure = structure.rename(columns={"Structure": col})
                input_data = input_data.merge(structure, how="outer", on=["i", "j"])
            input_data[columns] = input_data[columns].fillna(0)
            input_data["Num_structures"] = (
                input_data[columns].sum(axis=1, numeric_only=True) - 1
            )  # to index at 0 for coloring
            input_data = input_data.reset_index(drop=True)
        total_structures = input_data["Num_structures"].max() + 1
        metric_defaults = {
            "Num_structures": {
                "metric_column": "Num_structures",
                "cmap": sns.color_palette("rainbow_r", total_structures),
                "normalization": "none",
                "values": [],
                "title": "Base-pairs by structure",
                "extend": "neither",
                "ticks": [i for i in range(total_structures)],
                "tick_labels": [i + 1 for i in range(total_structures)],
            }
        } | metric_defaults
        super().__init__(
            input_data=input_data,
            sequence=sequence,
            metric=metric,
            metric_defaults=metric_defaults,
            window=window,
            name=name,
        )
