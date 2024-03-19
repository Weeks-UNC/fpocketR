"""annotations.py contains Annotations and subclasses."""
import re
import pandas as pd
from rnavigate import data


class Annotation(data.Sequence):
    """Basic annotation class to store 1D features of an RNA sequence

    Each feature type must be a seperate instance. Feature types include:
        a group of separted nucleotides (e.g. binding pocket)
        regions of interest (e.g. coding sequence, Alu elements)
        sites of interest (e.g. m6A locations)
        primer binding sites.

    Parameters
    ----------
    input_data : list
        List will be treated according to `annotation_type` argument.
        Expected behaviors for each value of `annotation_type`:
        "sites" or "group": 1-indexed location of sites of interest
            example: [1, 10, 20, 30] is four sites, 1, 10, 20, and 30
        "spans": 1-indexed, inclusive locations of spans of interest
            example: [[1, 10], [20, 30]] is two spans, 1 to 10 and 20 to 30
        "primers": Similar to spans, but 5'/3' direction is preserved.
            example: [[1, 10], [30, 20]] forward 1 to 10, reverse 30 to 20
    annotation_type : "group", "sites", "spans", or "primers"
        The type of annotation.
    sequence : str or pandas.DataFrame
        Nucleotide sequence, path to fasta file, or dataframe containing a
        "Sequence" column.
    name : str, defaults to None
        Name of annotation.
    color : matplotlib color-like, defaults to "blue"
        Color to be used for displaying this annotation on plots.

    Attributes
    ----------
    data : pandas.DataFrame
        Stores the list of sites or regions
    name : str
        The label for this annotation for use on plots
    color : valid matplotlib color
        Color to represent annotation on plots
    sequence : str
        The reference sequence string
    """

    def __init__(self, input_data, annotation_type, sequence, name=None, color="blue"):
        """Create an Annotation."""
        super().__init__(sequence, name=name)
        self.color = color

        # make sure input data matches expected format
        valid_types = ["sites", "spans", "group", "primers"]
        self.annotation_type = annotation_type
        if isinstance(input_data, pd.DataFrame):
            self.data = input_data
        elif annotation_type in ["spans", "primers"]:
            self.data = self.from_spans(input_data)
        elif annotation_type in ["sites", "group"]:
            self.data = self.from_sites(input_data)
        else:
            raise ValueError(f"annotation_type not in {valid_types}")

    @classmethod
    def from_boolean_array(
        cls, values, sequence, annotation_type, name, color="blue", window=1
    ):
        """Create an Annotation from an array of boolean values.

        True values are used to create the Annotation.

        Parameters
        ----------
        values : list of True or False
            the boolean array
        sequence : string or rnav.data.Sequence
            the sequence of the Annotation
        annotation_type : "spans", "sites", "primers", or "group"
            the type of the new annotation
            If "spans" or "primers", adjacent True values, or values within
            window are collapse to a region.
        name : string
            a name for labelling the annotation.
        color : string, defaults to "blue"
            a color for plotting the annotation
        window : integer, defaults to 1
            a window around True values to include in the annotation.

        Returns
        -------
        rnavigate.data.Annotation
            the new Annotation
        """
        annotations = []
        current_annotation = None
        pad = window // 2
        for i, value in enumerate(values):
            position = i + 1
            if value and annotation_type in ["sites", "group"]:
                annotations.append(position)
            elif value and annotation_type == "spans":
                start = max(1, position - pad)
                stop = min(len(sequence), position + pad)
                if current_annotation is None:
                    current_annotation = [start, stop]
                elif start <= current_annotation[1] + 1:
                    current_annotation[1] = stop
                elif start > current_annotation[1]:
                    annotations.append(current_annotation)
                    current_annotation = [start, stop]
        if annotation_type == "spans":
            annotations.append(current_annotation)
        return cls(
            input_data=annotations,
            annotation_type=annotation_type,
            color=color,
            sequence=sequence,
            name=name,
        )

    def from_spans(self, spans):
        """Create the self.data dataframe from a list of spans."""
        data_dict = {"start": [], "end": []}
        for span in spans:
            if (len(span) != 2) and any(not isinstance(pos, int) for pos in span):
                raise ValueError(
                    f"{self.annotation_type} must be a list of pairs of "
                    f"integers:\n{spans}"
                )
            start, end = span
            data_dict["start"].append(int(start))
            data_dict["end"].append(int(end))
        return pd.DataFrame(data_dict)

    def from_sites(self, sites):
        """Create the self.data dataframe from a list of sites."""
        if any(not isinstance(site, int) for site in sites):
            raise ValueError(
                f"{self.annotation_type} must be a list of integers:\n{sites}"
            )
        return pd.DataFrame({"site": sites})

    def get_aligned_data(self, alignment):
        """Aligns this Annotation to a new sequence and returns a copy.

        Parameters
        ----------
        alignment : rnavigate.data.Alignment
            Alignment object used to align to a new sequence.

        Returns
        -------
        rnavigate.data.Annotation
            A new Annotation with the same name, color, and annotation
            type, but with the input data aligned to the target sequence.
        """
        new_input_data = alignment.map_dataframe(self.data, self.data.columns)
        return Annotation(
            name=self.name,
            color=self.color,
            sequence=alignment.target_sequence,
            annotation_type=self.annotation_type,
            input_data=new_input_data,
        )

    def get_sites(self):
        """Returns a list of nucleotide positions included in this annotation.

        Returns
        -------
        sites : tuple
            a list of nucleotide positions
        """
        if self.annotation_type in ["spans", "primers"]:
            sites = []
            for _, row in self.data.iterrows():
                sites.extend(list(range(row["start"], row["end"] + 1)))
        elif self.annotation_type in ["sites", "group"]:
            sites = list(self.data["site"].values)
        return sites

    def get_subsequences(self, buffer=0):
        subsequences = []
        if self.annotation_type in ["spans", "primers"]:
            for _, (start, stop) in self.data[["start", "stop"]].iterrows():
                subsequences.append(self.sequence[start - 1 - buffer : stop + buffer])
        elif self.annotation_type in ["sites", "groups"]:
            for _, (site) in self.data[["sites"]].iterrows():
                subsequences.append(self.sequence[site - 1 - buffer : stop + buffer])
        return subsequences

    def __getitem__(self, idx):
        return self.data.loc[idx]

    def __iter__(self):
        for _, row in self.data.iterrows():
            yield row

    def __len__(self):
        return len(self._list)


class Motif(Annotation):
    """Automatically annotates the occurances of a sequence motif as spans.

    Parameters
    ----------
    input_data : str
        sequence motif to search for.
        Uses conventional nucleotide codes.
        e.g. "DRACH" = [AGTU] [AG] A C [ATUC]
    sequence : str or pandas.DataFrame
        Nucleotide sequence, path to fasta file, or dataframe containing a
        "Sequence" column.
    name : str, defaults to None
        Name of annotation.
    color : matplotlib color-like, defaults to "blue"
        Color to be used for displaying this annotation on plots.

    Attributes
    ----------
    data : pandas.DataFrame
        Stores the list of regions that match the motif
    name : str
        The label for this annotation for use on plots
    color : valid matplotlib color
        Color to represent annotation on plots
    sequence : str
        The reference sequence string
    """

    def __init__(self, input_data, sequence, name=None, color="blue"):
        """Creates a Motif annotation"""
        self.motif = input_data
        span_list = self.get_spans_from_motif(sequence, input_data)
        super().__init__(
            name=name,
            sequence=sequence,
            annotation_type="spans",
            input_data=span_list,
            color=color,
        )

    def get_spans_from_motif(self, sequence, motif):
        """Returns a list of spans for each location of motif found within sequence.

        Parameters
        ----------
        sequence : string
            sequence to be searched
        motif : string
            sequence motif to searched for.

        Returns
        -------
        spans : list of lists
            list of [start, end] positions of each motif occurance
        """
        sequence = data.normalize_sequence(sequence=sequence)
        nuc_codes = {
            "A": "A",
            "T": "T",
            "U": "U",
            "G": "G",
            "C": "C",
            "B": "[CGTU]",
            "D": "[ATUG]",
            "H": "[ATUC]",
            "V": "[ACG]",
            "W": "[ATU]",
            "S": "[CG]",  # strong and weak
            "M": "[AC]",
            "K": "[GTU]",  # amino and ketone
            "R": "[AG]",
            "Y": "[CTU]",  # purine and pyrimidine
            "N": "[ATUGC]",
        }  # any nuc
        re_pattern = "".join([nuc_codes[n] for n in motif])
        spans = []
        for match in re.finditer(re_pattern, sequence):
            start, end = match.span()
            spans.append([start + 1, end])
        return spans

    def get_aligned_data(self, alignment):
        """Searches the new sequence for the motif and returns a new Motif annotation.

        Parameters
        ----------
        alignment : rnavigate.data.Alignment
            Alignment object used to align to a new sequence.

        Returns
        -------
        rnavigate.data.Motif
            A new Motif with the same name, color, and motif
            but with the input data aligned to the target sequence.
        """
        return Motif(
            input_data=self.motif,
            name=self.name,
            color=self.color,
            sequence=alignment.target_sequence,
        )


class ORFs(Annotation):
    """Automatically annotations occurances of open-reading frames as spans.

    Parameters
    ----------
    input_data : "longest" or "all"
        which ORFs to annotate. "longest" annotates the longest ORF. "all"
        annotates all potential ORFs.
    sequence : str or pandas.DataFrame
        Nucleotide sequence, path to fasta file, or dataframe containing a
        "Sequence" column.
    name : str, defaults to None
        Name of annotation.
    color : matplotlib color-like, defaults to "blue"
        Color to be used for displaying this annotation on plots.

    Attributes
    ----------
    data : pandas.DataFrame
        Stores the list of regions that match the motif
    name : str
        The label for this annotation for use on plots
    color : valid matplotlib color
        Color to represent annotation on plots
    sequence : str
        The reference sequence string
    """

    def __init__(self, input_data, name=None, sequence=None, color="blue"):
        """Creates an ORF annotation"""
        self.input_data = input_data
        span_list = self.get_spans_from_orf(sequence, which=input_data)
        super().__init__(
            name=name,
            sequence=sequence,
            annotation_type="spans",
            input_data=span_list,
            color=color,
        )

    def get_spans_from_orf(self, sequence, which="all"):
        """Given a sequence string, returns spans for specified ORFs

        Parameters
        ----------
        sequence : string
            RNA nucleotide sequence
        which : "longest" or "all", defaults to "all"
            "all" returns all spans, "longest" returns the longest span

        Returns
        -------
        list of tuples
            (start, end) position of each ORF 1-indexed, inclusive
        """
        if isinstance(sequence, data.Sequence):
            sequence = sequence.sequence
        sequence = sequence.upper().replace("T", "U")
        spans = []
        stop_codons = "UAA|UAG|UGA"
        start_codon = "AUG"
        stop_sites = [m.end() for m in re.finditer(stop_codons, sequence)]
        start_sites = [m.start() + 1 for m in re.finditer(start_codon, sequence)]
        for start in start_sites:
            for stop in stop_sites:
                if (stop - start) % 3 == 2 and start < stop:
                    spans.append([start, stop])
        if which == "all":
            return spans
        if which == "longest":
            lengths = [end - start for start, end in spans]
            index = lengths.index(max(lengths))
            return [spans[index]]

    def get_aligned_data(self, alignment):
        """Searches the new sequence for ORFs and returns a new ORF annotation.

        Parameters
        ----------
        alignment : rnavigate.data.Alignment
            Alignment object used to align to a new sequence.

        Returns
        -------
        rnavigate.data.ORFs
            A new ORFs annotation with the same name, color, and input_data
            but with the input data aligned to the target sequence.
        """
        return ORFs(
            input_data=self.input_data,
            name=self.name,
            color=self.color,
            sequence=alignment.target_sequence,
        )


def domains(input_data, names, colors, sequence):
    """Create a list of Annotations from a list of spans.

    Currently, domains functionality in RNAvigate just uses a list of spans. In the
    future, this should be a dedicated class. Generally, domains should cover an entire
    sequence without overlap, but this is not enforced.
    e.g. [[1, 100], [101, 200]] for a 200 nt sequence.

    Parameters
    ----------
    input_data : list of lists
        list of spans for each domain
    names : list of strings
        list of names for each domain
    colors : list of valid matplotlib colors
        list of colors for each domain
    sequence : string
        sequence to be annotated

    Returns
    -------
    list of rnavigate.data.Annotation
        list of Annotations
    """
    return [
        Annotation(
            input_data=[span],
            annotation_type="spans",
            name=name,
            color=color,
            sequence=sequence,
        )
        for span, name, color in zip(input_data, names, colors)
    ]
