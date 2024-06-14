"""LogCompare compares reactivity profiles for significant differences.

This analysis requires replicates.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from rnavigate.plots import Plot


class LogCompare:
    """Compares 2 experimental samples, given replicates of each sample.
    First, the log10(modified/untreated) rate is calculated. These values a
    scaled to minimize the median of the absolute value of the difference
    between samples. The standard error in these values is computed for each
    replicate. Z-scores between samples are calculated. The results are plotted
    in two panels: (1) the scaled log10(modified/untreated) rate for each
    sample with error bars, and (2) the difference between samples, colored by
    z-score.

    Methods:
        __init__: computes log10(modified/untreated) rates, rescales the data,
            then calls make_plot()
        get_profile_sequence: gets log10(m/u) rate and sequence from sample
        calc_scale_factor: gets scale factor given two profiles
        rescale: rescales a profile to minimize difference to another profile
        load_replicates: calculates average and standard error of replicates
        make_plots: displays the two panels described above.

    Attributes:
        data (str): a key of sample.data to retrieve per-nucleotide data
        groups (dict): a dictionary containing the following key-value pairs:
            1: a dictionary containing these key-value pairs:
                self.data: averaged scaled log10(m/u) across replicates
                "stderr": the standard errors across replicates
                "stacked": 2d array containing each scaled log10(m/u) array
                "seq": the sequence string
            2: same as 1 above, for the second sample
    """

    def __init__(self, samples1, samples2, name1, name2, data="profile", region="all"):
        """Takes replicates of two samples for comparison. Replicates are
        required. Calculates the log division profile
        (log10(modified/untreated)) and minimizes the median of the absolute
        difference between these. Finally, creates a plot.

        Args:
            samples1 (list of Sample objects): Replicates of the first sample
            samples2 (list of Sample objects): Replicates of the second sample
            name1 (string): name of first sample
            name2 (string): name of second sample
            data (str, optional): Datatype to compare. Defaults to "profile".
        """
        if region == "all":
            self.region = [0, samples1[0].data["profile"].length]
        else:
            self.region = region
        self.data = data
        self.groups = {}
        self.load_replicates(*samples1, group=1)
        self.load_replicates(*samples2, group=2)
        sequences_match = self.groups[1]["seq"] == self.groups[2]["seq"]
        assert sequences_match, "Sample sequences do not match."
        self.groups[2][self.data] = self.rescale(
            self.groups[2][self.data], self.groups[1][self.data]
        )
        self.groups[1]["name"] = name1
        self.groups[2]["name"] = name2
        self.make_plots()

    def get_profile_sequence(self, sample):
        """retrieves log10(Modified_rate/Untreated_rate) and sequence from
        sample.data[self.data].

        Args:
            sample (rnavigate.Sample): an rnavigate sample

        Returns:
            np.array, np.array: log10 profile and sequence
        """
        df = sample.data[self.data].data
        plus = df.Modified_rate.values.copy()
        minus = df.Untreated_rate.values.copy()
        profile = plus / minus
        profile = np.log(profile)
        profile[minus > 0.05] = np.nan
        return profile, sample.data[self.data].sequence

    def calc_scale_factor(self, profile, target_profile):
        """Finds the scale factor (x) that minimizes the median value of
        |profile + x - target_profile|.

        Args:
            profile (np.array): log10 profile
            target_profile (np.array): 2nd log10 profile
        """

        def f(offset):
            return np.nanmedian(np.abs(profile + offset - target_profile))

        result = minimize_scalar(f, bounds=[-20, 20])
        offset = result.x
        return offset

    def rescale(self, profile, target_profile):
        """scales profile to minimize difference to target_profile using
        calc_scale_factor

        Args:
            profile (np.array): log10 profile to scale
            target_profile (np.array): 2nd log10 profile

        Returns:
            np.array: scaled profile
        """
        offset = self.calc_scale_factor(profile, target_profile)
        return profile + offset

    def load_replicates(self, *samples, group):
        """calculates average and standard error for a group of replicates,
        stores these values, along with sequence and a np.array containing
        all profiles in self.groups[group].

        Args:
            *samples (list of rnavigate.Sample): replicates to load
            group (int): self.groups key to access replicate data
        """
        profile1, seq1 = self.get_profile_sequence(samples[0])
        primermask = np.array([c.islower() for c in seq1])
        profile1[primermask] = np.nan

        profiles = [profile1]
        for sample in samples[1:]:
            profile, seq = self.get_profile_sequence(sample)
            assert seq == seq1, "Replicate sequences do not match."
            profile[primermask] = np.nan
            profile = self.rescale(profile, profile1)
            profiles.append(profile)
        stacked = np.vstack(profiles)
        avgprofile = np.nanmean(stacked, ax=0)
        stderr = np.std(stacked, ax=0)
        self.groups[group] = {
            self.data: avgprofile,
            "stderr": stderr,
            "stacked": stacked,
            "seq": seq,
        }

    def make_plots(self):
        """Visualize this analysis."""
        # get the replicate data
        prof1, stderr1, stack1, seq, name1 = self.groups[1].values()
        prof2, stderr2, stack2, seq, name2 = self.groups[2].values()

        x = np.array(list(range(1, len(prof1) + 1)))

        # create a two row, one column figure, scaled for RNA length
        _, axes = plt.subplots(2, 1, figsize=(0.1 * len(prof1), 14), sharex=True)
        # first axes contains raw log10 profiles with error bars
        ax = axes[0]
        ax.step(x, prof1, label=name1, color="C0")
        ax.step(x, prof2, label=name2, color="C1")
        ax.fill_between(
            x,
            prof1 - stderr1,
            prof1 + stderr1,
            step="pre",
            color="C0",
            alpha=0.25,
            lw=0,
        )
        ax.fill_between(
            x,
            prof2 - stderr2,
            prof2 + stderr2,
            step="pre",
            color="C1",
            alpha=0.25,
            lw=0,
        )
        ax.legend(loc="upper right")
        ax.axhline(0, color="black", lw=1, zorder=0)
        ax.set_ylabel("ln(Mod/BG)")
        Plot.add_sequence(self, ax=ax, sequence=seq)

        # calculate mean difference and standard error for z-scores
        stack1 = stack1 + 10
        rescale_dms = self.rescale(stack1, stack2)
        stack1 = rescale_dms
        diff = stack2 - stack1
        meandiff = np.nanmean(diff, ax=0)
        std_err = np.std(diff, ax=0)
        # compute z-scores
        z_scores = meandiff / std_err
        # normalize to 0-1 (values > 5 = 1, values < -5 = 0)
        z_scores = (z_scores + 5) / 10
        z_scores[z_scores > 1] = 1
        z_scores[z_scores < 0] = 0
        # get colormapping for differences
        # blue, white, red gradient from -5 to +5 stderrs
        colormap = plt.get_cmap("bwr")

        # create 2nd plot of differences colored by z-score
        ax = axes[1]
        ax.bar(x - 0.5, meandiff, color=colormap(z_scores), lw=0)
        ax.axhline(0, color="black", lw=1)
        ax.set_xlabel("nucleotide")
        ax.set_ylabel(f"{name2}-{name1}")
        ax.errorbar(x - 0.5, meandiff, yerr=std_err, fmt="none", color="black", lw=1)
        Plot.add_sequence(self, ax=ax, sequence=seq)

        # create a color bar scale for z-score differences
        axin1 = ax.inset_axes([0.8, 0.1, 0.15, 0.15])
        cmap = plt.get_cmap("bwr")
        Plot.view_colormap(
            ax=axin1,
            ticks=[0, 5, 10],
            values=[-5, 0, 5],
            title="Z-score",
            cmap=list(cmap(np.arange(cmap.N))),
        )
        plt.tight_layout()
