"""Performs low SHAPE, low Shannon entropy analysis

Citation:
    Siegfried, N., Busan, S., Rice, G. et al. RNA motif discovery by SHAPE and
        mutational profiling (SHAPE-MaP). Nat Methods 11, 959-965 (2014).
        https://doi.org/10.1038/nmeth.3029

Typical usage example:

import rnavigate as rnav
my_sample = rnav.Sample(
    sample="example sample",
    shapemap="my_shape_profile.txt",
    pairprob="pairing_probabilities.txt",
    ss="MFE_structure.ct"
    )
lowss_sample = rnav.analysis.LowSS(my_sample)
plot = lowss_sample.plot_lowss()
plot.save("lowss_figure.svg")
"""

import numpy as np
from rnavigate import Sample, plots, data


class LowSS(Sample):
    """Creates a new RNAvigate Sample which computes and displays Low SHAPE,
    low Shannon entropy regions (LowSS) given a sample containing SHAPE
    reactivities, pairing probabilities, and MFE structure.

    Methods:
        __init__: performs the analysis
        plot_lowss: displays the result and returns plot object

    Attributes:
        sample (str): the new label for this Sample's data on plots
        parent (rnavigate.Sample): the sample from which data is retrieved
        window (int): size of the windows, must be odd
        median_shape (float): global median SHAPE reactivity
        median_entropy (float): global median Shannon entropy
        data (dictionary): dictionary of data keyword: Data objects, keys are:
            "structure" (rnav.data.SecondaryStructure)
                copy of provided MFE structure
            "shapemap" (rnav.data.SHAPEMaP)
                copy of provided SHAPE-MaP data aligned to "structure"
            "pairprob" (rnav.data.PairingProbability)
                copy of pairing probabilities aligned to "structure"
            "entropies" (rnav.data.Profile)
                Profile of Shannon entropies calculated from "pairprob"
            "lowSS" (rnav.data.Annotations)
                annotations defining low SHAPE, low Shannon entropy regions
    """

    def __init__(
        self,
        sample,
        window=55,
        shapemap="shapemap",
        pairprob="pairprob",
        structure="ss",
    ):
        """Perform Low SHAPE, low Shannon entropy analysis given a sample with
        1) reactivities 2) MFE structure 3) pairing probabilities.

        Required arguments:
            sample (rnav.Sample)
                sample with shapemap, pairing probabilities and MFE structure

        Optional arguments:
            window (integer)
                Window size for calculating median SHAPE and Shannon entropy
                Defaults to 55
            shapemap (data keyword)
                data keyword that points to SHAPE-MaP data
                defaults to "shapemap"
            pairprob (data keyword)
                data keyword that points to pairing probabilities data
                defaults to "pairprob"
            structure (data keyword)
                data keyword that points to MFE structure data
                defaults to "ss"
        """
        # store values
        structure = sample.get_data(structure, data.SecondaryStructure)
        shapemap = sample.get_data(shapemap, data.SHAPEMaP)
        pairprob = sample.get_data(pairprob, data.PairingProbability)
        structure = structure.get_aligned_data(structure.null_alignment)
        shapemap = shapemap.get_aligned_data(
            data.SequenceAlignment(shapemap, structure)
        )
        pairprob = pairprob.get_aligned_data(
            data.SequenceAlignment(pairprob, structure)
        )
        entropy = pairprob.get_entropy_profile()
        super().__init__(
            sample=f"{sample.sample}: lowSS",
            structure=structure,
            shapemap=shapemap,
            pairprob=pairprob,
            entropy=entropy,
        )
        self.parent = sample
        self.window = window
        self.median_shape = shapemap.data["Norm_profile"].median(skipna=True)
        self.median_entropy = entropy.data["Entropy"].median(skipna=True)
        self.reset_window()

    def reset_lowss(self, maximum_shape=None, maximum_entropy=0.08):
        """Generates an annotation of lowSS regions. Stored as self.lowSS

        Args:
            maximum_shape (float, optional): maximum normalized SHAPE
                reactivity to be called lowSS. Defaults to median reactivity.
            maximum_entropy (float, optional): maximum shannon entropy to be
                called lowSS. Defaults to 0.8.
        """
        if maximum_shape is None:
            maximum_shape = self.median_shape
        entropies = self.data["entropy"].data["Windowed_entropy"]
        shapes = self.data["shapemap"].data["Windowed_profile"]
        is_lowss = (entropies < maximum_entropy) & (shapes < maximum_shape)
        self.set_data(
            overwrite=True,
            data_keyword="lowSS",
            inputs=data.Annotation.from_boolean_array(
                values=is_lowss,
                window=self.window,
                annotation_type="spans",
                color="grey",
                sequence=self.data["structure"].sequence,
                name="lowSS",
            ),
        )

    def reset_window(self, window=None):
        """Resets the window size and recalculates windowed SHAPE reactivities
        and shannon entropies and lowSS region annotations.

        Optional arguments:
            window (integer)
                window size, must be an odd integer
                Defaults to value of window attribute
        """
        if window is None:
            window = self.window
        else:
            self.window = window
        # Calculate windowed median SHAPE and entropies
        self.data["shapemap"].calculate_windows(
            column="Norm_profile",
            window=window,
            method="median",
            new_name="Windowed_profile",
            minimum_points=1,
        )
        self.data["entropy"].calculate_windows(
            column="Entropy",
            window=window,
            method="median",
            new_name="Windowed_entropy",
            minimum_points=1,
        )
        self.reset_lowss()

    def plot_lowss(self, region=None, colorbars=True):
        """Visualize LowSS analysis over the given region.

        Optional arguments:
            region (integer or list of 2 integers)
                lowSS region number or start and end positions to plot.
                For region numbers, +/- 150 nts are shown.
                Defaults to entire sequence.
            colorbars (True or False)
                whether to plot colorbars for pairing probability

        Returns:
            rnavigate.plots.AP: LowSS visualization
        """
        shapemap = self.data["shapemap"]
        structure = self.data["structure"]
        pairprob = self.data["pairprob"]
        entropy = self.data["entropy"]
        lowss = self.data["lowSS"]
        # show entire RNA if region is not provided
        if region is None:
            start = 1
            stop = structure.length
            region = [start, stop]
        # if region is an integer, get that lowSS region +/- 150 nts
        elif isinstance(region, int):
            region = lowss[region - 1]
            start, stop = region
            start = max(start - 150, 1)
            stop = min(stop + 150, structure.length)
        # else region should be a list of 2 integers, start and end.
        else:
            start, stop = region
        region_length = stop - start + 1
        # create arc plot instance
        plot = plots.AP(1, region_length, cols=1, rows=1, region=region)
        ax = plot.axes[0, 0]
        x_values = np.arange(start, stop + 1)
        # plot median SHAPE on secondary ax
        prof_ax = ax.twinx()
        prof_ax.set_ylim(-3, 1)
        prof_ax.set_yticks([0.0, 0.4, 0.85])
        prof_ax.fill_between(
            x=x_values,
            y1=[self.median_shape] * region_length,
            y2=shapemap.data["Windowed_profile"][start - 1 : stop],
            fc="0.3",
            lw=0,
        )
        plots.adjust_spines(prof_ax, ["left"])
        plots.clip_spines(prof_ax, ["left"])

        # plot median entropy on second secondary ax
        ent_ax = ax.twinx()
        ent_ax.set_ylim(-1.5, 1.5)
        ent_ax.set_yticks(ticks=[0, 0.08, 0.5])
        ent_ax.fill_between(
            x=x_values,
            y1=[0.08] * region_length,
            y2=entropy.data["Windowed_entropy"][start - 1 : stop],
            fc="C1",
            lw=0,
        )
        plots.adjust_spines(ent_ax, ["left"])
        plots.clip_spines(ent_ax, ["left"])

        # add ss and pairing probabilities track
        plot.plot_data(
            sequence=structure,
            structure=structure,
            interactions=pairprob,
            label=self.sample,
            seqbar=False,
            title=False,
            annotations=[lowss],
            annotation_mode="bar",
        )

        # Place Track Labels
        ax.set_title(
            f"{self.sample}\n{start} - {stop}", loc="left", fontdict={"fontsize": 16}
        )
        ax.text(
            1.002,
            7 / 8,
            "SHAPE\nReactivity",
            fontsize=12,
            transform=ax.transAxes,
            va="center",
        )
        ax.text(
            1.002,
            5 / 8,
            "Shannon\nEntropy",
            fontsize=12,
            transform=ax.transAxes,
            va="center",
        )
        ax.text(
            1.002,
            3 / 8,
            "Secondary\nStructure",
            transform=ax.transAxes,
            fontsize=12,
            va="center",
        )
        ax.text(
            1.002,
            1 / 8,
            "Pairing\nProbability",
            transform=ax.transAxes,
            va="center",
            fontsize=12,
        )

        # Place region labels
        for i, (lssr_start, lssr_stop) in enumerate(lowss):
            middle = (lssr_stop + lssr_start) / 2
            if start < middle < stop:
                ax.text(
                    x=middle, y=550, s=str(i + 1), ha="center", va="center", fontsize=12
                )

        ax.set_ylim([-305, 915])
        ax.set_xticks(ticks=[x for x in range(500, stop + 1, 500) if x > start])
        ax.set_xticks(
            ticks=[x for x in range(100, stop + 1, 100) if x > start], minor=True
        )
        ax.tick_params(axis="x", which="major", labelsize=12)
        plots.adjust_spines(ax, ["bottom"])

        # set figure size so that 200 ax units == 1 inch
        plot.set_figure_size(height_ax_rel=1 / 200, width_ax_rel=1 / 200)
        if colorbars:
            plot.plot_colorbars()
        return plot
