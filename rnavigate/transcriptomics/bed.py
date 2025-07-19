import pandas as pd
from rnavigate import data


class BedFile:
    def __init__(self, bedfile):
        self.bedfile = bedfile
        self.read_kwargs = {
            "header": None,
            "index_col": False,
            "usecols": [0, 1, 2, 3, 4, 5],
            "names": ["chr", "start", "end", "NA", "score", "strand"],
            "dtype": {
                "chr": str,
                "start": int,
                "end": int,
                "NA": str,
                "score": float,
                "strand": str,
            },
        }
        if self.bedfile.endswith(".gz"):
            self.read_kwargs["compression"] = "gzip"
        self.profile_cols = ["Score"]

    def get_annotation(self, transcript, **kwargs):
        return self.get_annotations([transcript], **kwargs)[transcript]

    def get_annotations(self, transcripts, **kwargs):
        bed_df = pd.read_table(self.bedfile, **self.read_kwargs)
        bed_df["start"] += 1  # bed files are 0-indexed, closed right
        tx_annotations = {}
        for tx in transcripts:
            chrom = tx.chromosome
            strand = tx.strand
            tx_coordinates = tx.coordinate_df["Coordinate"]
            mn = tx_coordinates.min()
            mx = tx_coordinates.max()
            if strand == "-":
                mn, mx = mx, mn
            tx_df = bed_df.query(
                f'(chr == "{chrom}") '
                f'& (strand == "{strand}") '
                f"& ((start > {mn} and start < {mx})"
                f"| (end > {mn} and end < {mx}))"
            )
            spans = []
            for _, row in tx_df.iterrows():
                span = tx.get_tx_range(row["start"], row["end"])
                if span[0] != span[1]:
                    spans.append(span)
            tx_annotations |= {
                tx: data.Annotation(
                    input_data=spans,
                    annotation_type="spans",
                    sequence=tx.sequence,
                    **kwargs,
                )
            }
        return tx_annotations

    def get_profile(self, transcript, **kwargs):
        bed_df = pd.read_table(self.bedfile, **self.read_kwargs)
        chrom = transcript.chromosome
        strand = transcript.strand
        profile = transcript.coordinate_df.copy()
        mn = profile["Coordinate"].min()
        mx = profile["Coordinate"].max()
        bed_df = bed_df.query(
            f'chr == "{chrom}" & strand == "{strand}" & '
            f"((start > {mn} & start < {mx}) | (end > {mn} & end < {mx}))"
        )
        profile[self.profile_cols] = 0
        for _, row in bed_df.iterrows():
            site = profile.eval(
                f'Coordinate > {row["start"]} and Coordinate < {row["end"]}'
            )
            for col in self.profile_cols:
                profile.loc[site, col] = row[col]
        return data.Profile(input_data=profile, metric="Value", **kwargs)

    def get_density_profile(self, transcript, **kwargs):
        bed_df = pd.read_table(self.bedfile, **self.read_kwargs)
        chrom = transcript.chromosome
        strand = transcript.strand
        profile = transcript.coordinate_df.copy()
        mn = profile["Coordinate"].min()
        mx = profile["Coordinate"].max()
        bed_df = bed_df.query(
            f'chr == "{chrom}" & strand == "{strand}" & '
            f"((start > {mn} & start < {mx}) | (end > {mn} & end < {mx}))"
        )
        profile["Density"] = 0
        for _, row in bed_df.iterrows():
            site = profile.eval(
                f'Coordinate > {row["start"]} and Coordinate < {row["end"]}'
            )
            profile.loc[site, "Density"] += 1
            for col in self.profile_cols:
                profile.loc[site, col] += row[col]
        return data.Profile(input_data=profile, metric="Value", **kwargs)


class NarrowPeak(BedFile):
    def __init__(self, bedfile):
        super().__init__(bedfile=bedfile)
        self.read_kwargs["usecols"].extend([6, 7, 8, 9])
        self.read_kwargs["names"].extend(["Value", "P_value", "Q_value", "Peak"])
        self.read_kwargs["dtype"].update(
            {"Value": float, "P_value": float, "Q_value": float, "Peak": int}
        )
        self.profile_cols = ["Score", "Value", "P_value", "Q_value", "Peak"]
