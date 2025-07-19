"""SequenceChecker analysis used to inspect sequence differences.

Given a list of samples, we can inspect which data keywords belong to the
samples, which sequences match up perfectly, and inspect the differences
between sequences.
"""
import pandas as pd
import numpy as np
from rnavigate import data


class SequenceChecker:
    """Check the sequences stored in a list of samples.

    Attributes:
        samples: the list of samples
        sequences: a list of all unique sequence strings stored in the list of
            samples. These are converted to an all uppercase RNA alphabet.
        keywords: a list of unique data keywords stored in the list of samples.
        which_sequences: a dataframe of samples and keywords and which
            sequences each contains.
    """

    def __init__(self, samples):
        """Creates an instance of SequenceChecker given a list of samples

        Arguments:
            samples (list of rnav.Sample)
                samples for which to compare data keywords and sequences.
        """
        self.samples = samples
        self.keywords = self.get_keywords()
        self.sequences = self.get_sequences()
        self.which_sequences = self.get_which_sequences()

    def reset(self):
        """Reset keywords and sequences from sample list in case of changes."""
        self.keywords = self.get_keywords()
        self.sequences = self.get_sequences()
        self.which_sequences = self.get_which_sequences()

    def get_keywords(self):
        """A list of all unique data keywords across samples."""
        return list(set([dkw for s in self.samples for dkw in s.data]))

    def get_sequences(self):
        """A list of all unique sequences (uppercase RNA) across samples."""
        sequences = []
        for sample in self.samples:
            for dkw in self.keywords:
                if dkw not in sample.data:
                    continue
                seq = sample.get_data(dkw).sequence
                seq = seq.upper().replace("T", "U")
                if seq not in sequences:
                    sequences.append(seq)
        return sequences

    def get_which_sequences(self):
        """A DataFrame of sequence IDs (integers) for each data keyword."""
        df = {key: [] for key in ["Sample"] + self.keywords}
        for sample in self.samples:
            df["Sample"].append(sample.sample)
            for dkw in self.keywords:
                if dkw not in sample.data:
                    df[dkw].append(np.nan)
                else:
                    seq = sample.get_data(dkw).sequence
                    seq = seq.upper().replace("T", "U")
                    df[dkw].append(self.sequences.index(seq))
        return pd.DataFrame(df)

    def print_which_sequences(self):
        """Print sequence ID (integer) for each data keyword and sample."""
        print("Sequence IDs")
        for _, row in self.which_sequences.iterrows():
            for column in self.which_sequences.columns:
                if column == "Sample":
                    print(f"    {row[column]}")
                else:
                    if np.isnan(row[column]):
                        continue
                    print(f"        {column:<10} {int(row[column])}")
        print()

    def print_alignments(self, print_format="long", which="all"):
        """Print alignments in the given format for sequence IDs provided.

        Optional arguments:
            print_format (string)
                What format to print the alignments in:
                "cigar" prints the cigar string
                "short" prints the numbers of mismatches and indels
                "long" prints the location and nucleotide identity of all
                    mismatches, insertions and deletions.
                Defaults to "long".
            which (pair of integers)
                two sequence IDs to compare.
                Defaults to every pairwise comparison.
        """
        kwargs = {"print_format": print_format}
        if which == "all":
            num = range(len(self.sequences))
            for s1 in num:
                for s2 in num:
                    if s1 < s2:
                        self.print_alignments(which=(s1, s2), **kwargs)
            return
        i, j = int(which[0]), int(which[1])
        a = self.sequences[i]
        b = self.sequences[j]
        alignment = data.SequenceAlignment(sequence1=a, sequence2=b, full=True)
        print(f"Sequence {i} -> {j} ({len(a)} nts -> {len(b)} nts)")
        alignment.print(print_format=print_format)

    def print_mulitple_sequence_alignment(self, base_sequence):
        """Print the multiple sequence alignment with nice formatting.

        Required arguments:
            base_sequence (string)
                a sequence string that represents the longest common sequence.
                Usually, this is the return value from:
                    rnav.data.set_multiple_sequence_alignment()
        """
        print("Multiple sequence alignment")
        alignments = [
            data.SequenceAlignment(base_sequence, seq).alignment2
            for seq in self.sequences
        ]
        pos = "".join(f"{n:<20}" for n in range(1, len(alignments[0]), 20))
        misses = "".join("X "[len(set(nts)) == 1] for nts in zip(*alignments))
        print("    ID    length   alignment")
        for i, (alignment, seq) in enumerate(zip(alignments, self.sequences)):
            print(f"    {i:<5} {len(seq):<8} {alignment}")
        print(f"    {'mismatches':<14} {misses}")
        print(f"    {'positions':<14} {pos}")
        print()

    def write_fasta(self, filename, which="all"):
        """Write all unique sequences to a fasta file.

        This is very useful for using external multiple sequence aligners such
        as ClustalOmega.
            1) go to https://www.ebi.ac.uk/Tools/msa/clustalo/
            2) upload new fasta file
            3) under STEP 2 output format, select Pearson/FASTA
            4) click 'Submit'
            5) wait for your alignment to finish
            6) download the alignment fasta file
            7) use rnav.data.set_multiple_sequence_alignment()

        Required arguments:
            filename (string)
                path to a new file to which fasta entries are written

        Optional arguments:
            which (list of integers)
                Sequence IDs to write to file.
        """
        if which == "all":
            which = range(len(self.sequences))
        with open(filename, "w") as fasta:
            for i in which:
                fasta.write(f">Sequence_{i}\n")
                fasta.write(f"{self.sequences[i]}\n")
