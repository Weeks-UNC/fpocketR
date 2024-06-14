#!/usr/bin/env python3

"""Fragmapper analysis tools.

Description:
FragMapper compares reactivity profile differences between SHAPE-MaP profiles.
The intended application of Fragmapper is to detect fragment or ligand
crosslinking sites in RNA.
"""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from rnavigate import data, Sample

__version__ = "0.1.0"
__author__ = "Seth D. Veenbaas"
__maintainer__ = "Seth D. Veenbaas"
__email__ = "sethv@live.unc.edu"


class FragMaP(data.Profile):
    def __init__(
        self,
        input_data,
        parameters,
        metric="Fragmap_profile",
        metric_defaults=None,
        read_table_kw=None,
        sequence=None,
        name=None,
    ):
        self.parameters = parameters
        if isinstance(input_data, (list, tuple)) and len(input_data) == 2:
            profile1, profile2 = input_data
            input_data = self.get_dataframe(profile1, profile2, **parameters)
        elif not isinstance(input_data, pd.DataFrame):
            raise ValueError(
                "Input data must be a list of 2 shapemap profiles or a "
                "pre-initialized dataframe"
            )
        if metric_defaults is None:
            metric_defaults = {}
        metric_defaults = {
            "Fragmap_profile": {
                "metric_column": "Fragmap_profile",
                "error_column": "Fragmap_err",
                "color_column": "Site",
                "cmap": ["lightgrey", "limegreen"],
                "normalization": "none",
                "values": None,
                "extend": "neither",
                "alpha": 1.0,
                "title": "Frag-MaP sites",
                "ticks": [0, 1],
                "tick_labels": ["other", "sites"],
            }
        } | metric_defaults
        super().__init__(
            input_data=input_data,
            metric=metric,
            metric_defaults=metric_defaults,
            read_table_kw=read_table_kw,
            sequence=sequence,
            name=name,
        )

    @property
    def recreation_kwargs(self):
        return {"parameters": self.parameters}

    def get_dataframe(
        self,
        profile1,
        profile2,
        mutation_rate_threshold,
        depth_threshold,
        p_significant,
        ss_threshold,
        correction_method,
    ):
        columns = [
            "Nucleotide",
            "Sequence",
            "Modified_mutations",
            "Modified_effective_depth",
            "Modified_rate",
            "Std_err",
        ]
        dataframe = pd.merge(
            profile1.data[columns],
            profile2.data[columns],
            how="left",
            on=["Nucleotide", "Sequence"],
            suffixes=("_1", "_2"),
        )

        # Filter data
        dataframe["filter"] = dataframe.eval(
            # Effective read depths must be greater than threshold
            f"Modified_effective_depth_1 > {depth_threshold} "
            f"& Modified_effective_depth_2 > {depth_threshold} "
            # Control mutation rate (profile 2) must be less than threshold
            f"& Modified_rate_2 < {mutation_rate_threshold}"
        ).astype(bool)
        valid = dataframe["filter"]

        # Calculate Z-scores for profile 1 and profile 2
        dataframe[["zscore_1", "zscore_2"]] = np.nan
        # Z-scores per nucleotide
        FragMaP.calc_zscore(
            self,
            valid,
            dataframe,
            incolumn="Modified_rate_1",
            outcolumn="zscore_1",
            base=["A", "U", "C", "G"],
        )
        FragMaP.calc_zscore(
            self,
            valid,
            dataframe,
            incolumn="Modified_rate_2",
            outcolumn="zscore_2",
            base=["A", "U", "C", "G"],
        )

        # Frag-MaP profile is the difference in Z-scores
        dataframe.eval("Fragmap_profile = zscore_1 - zscore_2", inplace=True)
        dataframe.eval(
            "Fragmap_err = sqrt(zscore_1_err ** 2 + zscore_2_err ** 2)", inplace=True
        )

        dataframe["Fragmap_pvalue"] = stats.norm.sf(dataframe["Fragmap_profile"])

        # Calculate T-statistic and P-value
        dataframe["tstat"], dataframe["pvalue"] = stats.ttest_ind_from_stats(
            mean1=dataframe["Modified_rate_1"],
            std1=dataframe["Std_err_1"],
            nobs1=dataframe["Modified_effective_depth_1"],
            mean2=dataframe["Modified_rate_2"],
            std2=dataframe["Std_err_2"],
            nobs2=dataframe["Modified_effective_depth_2"],
            equal_var=False,
        )

        # Bonferoni p-value correction
        p_significant = p_significant / len(dataframe)
        dataframe["Significant"] = dataframe["pvalue"] < p_significant

        # Calculate Delta rate
        dataframe.eval("Delta_rate = Modified_rate_1 - Modified_rate_2", inplace=True)
        dataframe["Delta_rate_zscore"] = stats.zscore(dataframe["Delta_rate"])

        # Perform False Discovery Rate (FDR) correction
        if correction_method is not None:
            p_values = dataframe["Fragmap_pvalue"].values
            reject, pvals_corrected, _, _ = multipletests(
                p_values, alpha=ss_threshold, method=correction_method
            )
            dataframe["Reject"] = reject
            dataframe["Fragmap_pvalue_corrected"] = pvals_corrected
        else:
            dataframe["Reject"] = dataframe["Fragmap_pvalue"] < ss_threshold
            dataframe["Fragmap_pvalue_corrected"] = np.nan

        # Called Frag-MaP sites
        dataframe["Site"] = (dataframe["Reject"]) & (dataframe["Significant"])

        return dataframe

    def calc_zscore(
        self, valid, dataframe, incolumn: str, outcolumn: str, base: list
    ) -> None:
        sele_data = dataframe.loc[dataframe["Sequence"].isin(base)].copy()

        # Calculate outlier nucleotides
        Q1 = sele_data[incolumn].quantile(0.25)
        Q3 = sele_data[incolumn].quantile(0.75)
        IQR = Q3 - Q1
        lower_limit = Q1 - 1.5 * IQR
        upper_limit = Q3 + 1.5 * IQR
        nonoutlier = dataframe[incolumn].between(lower_limit, upper_limit)

        # Calculate zscore
        valid_nt = valid & dataframe["Sequence"].isin(base)

        sem_column = f"Std_err_{incolumn.split('_')[-1]}"
        value = np.log(dataframe.loc[valid_nt, incolumn])
        sem = np.log(dataframe.loc[valid_nt, sem_column])
        mean = np.mean(np.log(dataframe.loc[valid_nt & nonoutlier, incolumn]))
        std = np.std(np.log(dataframe.loc[valid_nt & nonoutlier, incolumn]))

        dataframe.loc[valid_nt, outcolumn] = (value - mean) / std

        # Calculate the partial derivatives
        dz_dy = (mean * np.log(value)) / std
        dz_dx = np.log(value) / std
        dz_ds = -(mean * np.log(value) / (std**2))

        # Propagate error using the error propagation formula
        # Note: The SEMs for "x" and "s" are assumed to be zero.
        fragmap_err = np.sqrt((dz_dy * sem) ** 2 + (dz_dx * 0) ** 2 + (dz_ds * 0) ** 2)

        dataframe.loc[valid_nt, f"{outcolumn}_err"] = fragmap_err

    def get_annotation(self):
        return data.Annotation(
            input_data=self.data.loc[self.data["Site"], "Nucleotide"].to_list(),
            name="Frag-MaP sites",
            sequence=self.sequence,
            annotation_type="sites",
            color="green",
        )


class Fragmapper(Sample):
    def __init__(self, sample1, sample2, parameters=None, profile="shapemap"):
        if sample1.data[profile].sequence != sample2.data[profile].sequence:
            raise ValueError("Profiles must have the same sequence")
        if parameters is None:
            parameters = {}
        self.parameters = {
            "mutation_rate_threshold": 0.025,
            "depth_threshold": 5000,
            "p_significant": 0.001,
            "ss_threshold": 0.05,
            "correction_method": None,
        }
        self.parameters |= parameters
        fragmap = FragMaP(
            input_data=[sample1.get_data(profile), sample2.get_data(profile)],
            parameters=self.parameters,
        )
        super().__init__(
            sample=f"FragMaP: {sample1.sample}/{sample2.sample}",
            inherit=[sample1, sample2],
            fragmap=fragmap,
            fragmap_sites=fragmap.get_annotation(),
        )
        self.sample1 = sample1
        self.sample2 = sample2

    def plot_scatter(self, column="Modified_rate"):
        """Generates scatter plots useful for fragmapper quality control.

        Args:
            column (str, optional):
                Dataframe column containing data to plot (must be avalible for
                the sample and control).
                Defaults to "Modified_rate".

        Returns:
            (matplotlib figure, matplotlib axis)
                Scatter plot with control values on the x-axis, sample values
                on the y-axis, and each point representing a nucleotide not
                filtered out in the fragmapper pipeline.
        """
        fig, ax = plt.subplots(1, 1, figsize=(4, 4))

        columns = [f"{column}_1", f"{column}_2"]
        scatter_data = self.data["fragmap"].data
        scatter_data = scatter_data[scatter_data["filter"]].copy()
        sites = scatter_data[scatter_data["Site"]]
        non_sites = scatter_data[~scatter_data["Site"]]

        for _, (nt, y, x) in sites[["Nucleotide"] + columns].iterrows():
            ax.text(x, y, int(nt))

        ax.scatter(x=non_sites[columns[1]], y=non_sites[columns[0]], s=5, alpha=0.25)
        ax.scatter(x=sites[columns[1]], y=sites[columns[0]], s=25, alpha=1)

        ax.set_xlabel(f"{self.sample2.sample} {column}")
        ax.set_ylabel(f"{self.sample1.sample} {column}")

        plt.legend(["Uncalled Nucleotides", "Frag-MaP Sites"])

        return fig, ax
