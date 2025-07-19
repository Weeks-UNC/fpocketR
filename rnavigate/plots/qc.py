import numpy as np
import pandas as pd
import seaborn as sns
from rnavigate import plots


class QC(plots.Plot):
    """Plot QC data from log files.

    Parameters
    ----------
    num_samples : int
        Number of samples to plot.

    Attributes
    ----------
    fig : matplotlib.figure.Figure
        Figure object.
    ax_muts_unt : matplotlib.axes.Axes
        Axes object for the mutations per molecule plot of untreated samples.
    ax_muts_mod : matplotlib.axes.Axes
        Axes object for the mutations per molecule plot of modified samples.
    ax_read_unt : matplotlib.axes.Axes
        Axes object for the read length plot of untreated samples.
    ax_read_mod : matplotlib.axes.Axes
        Axes object for the read length plot of modified samples.
    ax_boxplot : matplotlib.axes.Axes
        Axes object for the boxplot of mutation rates.
    axes : numpy.ndarray of matplotlib.axes.Axes
        Array of axes objects.
    i : int
        Index of the current plot.
    """

    def __init__(self, num_samples):
        """Initialize the plot."""
        super().__init__(num_samples)
        if self.length == 1:
            self.ax_muts_unt = self.axes[0, 0]
            self.ax_muts_mod = self.axes[0, 0]
            self.ax_read_unt = self.axes[0, 1]
            self.ax_read_mod = self.axes[0, 1]
            self.ax_boxplot = self.axes[0, 2]
        elif self.length > 1:
            self.ax_muts_unt = self.axes[0, 0]
            self.ax_muts_mod = self.axes[0, 1]
            self.ax_muts_mod.get_shared_x_axes().join(
                self.ax_muts_mod, self.ax_muts_unt
            )
            self.ax_muts_mod.get_shared_y_axes().join(
                self.ax_muts_mod, self.ax_muts_unt
            )
            self.ax_muts_mod.set_yticklabels([])
            self.ax_read_unt = self.axes[0, 2]
            self.ax_read_mod = self.axes[0, 3]
            self.ax_muts_mod.get_shared_x_axes().join(
                self.ax_read_mod, self.ax_read_unt
            )
            self.ax_muts_mod.get_shared_y_axes().join(
                self.ax_read_mod, self.ax_read_unt
            )
            self.ax_read_mod.set_yticklabels([])
            gs = self.axes[1, 0].get_gridspec()
            for ax in self.axes[1, :]:
                ax.remove()
            self.ax_boxplot = self.fig.add_subplot(gs[1, :])
        self.ax_muts_unt.set(
            xlabel="Mutations per molecule",
            ylabel="Percentage of Reads",
            title="Untreated",
        )
        self.ax_muts_mod.set(xlabel="Mutations per molecule", title="Modified")
        self.ax_read_unt.set(
            xticks=range(12),
            xlabel="Read Length",
            ylabel="Percentage of Reads",
            title="Untreated",
        )
        self.ax_read_mod.set(xticks=range(12), xlabel="Read Length", title="Modified")
        if self.length == 1:
            self.ax_muts_mod.set_title("Mutations per Molecule Distribution")
            self.ax_read_mod.set_title("Read Length Distribution")
            self.ax_boxplot.set_title("Mutation rates")
        xticklabels = [f"{x*50}" for x in range(12)]
        self.ax_read_unt.set_xticklabels(xticklabels, rotation=45)
        self.ax_read_mod.set_xticklabels(xticklabels, rotation=45)
        self.profiles = []

    def set_figure_size(
        self,
        height_ax_rel=None,
        width_ax_rel=None,
        width_ax_in=2,
        height_ax_in=2,
        height_gap_in=1,
        width_gap_in=1,
        top_in=1,
        bottom_in=0.5,
        left_in=1,
        right_in=0.5,
    ):
        """Set the size of the figure.

        Parameters
        ----------
        height_ax_rel : float, optional
            Height of the axes relative to the y-axis limits.
        width_ax_rel : float, optional
            Width of the axes relative to the x-axis limits.
        width_ax_in : float, defaults to 2
            Width of the axes in inches.
        height_ax_in : float, defaults to 2
            Height of the axes in inches.
        height_gap_in : float, defaults to 1
            Height of the gap between axes in inches.
        width_gap_in : float, defaults to 1
            Width of the gap between axes in inches.
        top_in : float, defaults to 1
            Height of the top margin in inches.
        bottom_in : float, defaults to 0.5
            Height of the bottom margin in inches.
        left_in : float, defaults to 0.5
            Width of the left margin in inches.
        right_in : float, defaults to 0.5
            Width of the right margin in inches.
        """
        super().set_figure_size(
            height_ax_rel=height_ax_rel,
            width_ax_rel=width_ax_rel,
            width_ax_in=width_ax_in,
            height_ax_in=height_ax_in,
            height_gap_in=height_gap_in,
            width_gap_in=width_gap_in,
            top_in=top_in,
            bottom_in=bottom_in,
            left_in=left_in,
            right_in=right_in,
        )

    def get_rows_columns(self, rows=None, cols=None):
        """Get the number of rows and columns.

        Parameters
        ----------
        rows : int, optional
            Number of rows. This is ignored.
        cols : int, optional
            Number of columns. This is ignored.

        Returns
        -------
        rows : int
            Number of rows. 1 if there is only one sample, 2 otherwise.
        cols : int
            Number of columns. 3 if there is only one sample, 4 otherwise.
        """
        if self.length == 1:
            return (1, 3)
        else:
            return (2, 4)

    def plot_data(self, profile, label):
        """Plot the data.

        Parameters
        ----------
        profile : rnavigate.profile.Profile
            Profile object.
        label : str
            Label for the sample.
        """
        self.plot_mutations_per_molecule(profile=profile, label=label)
        self.plot_read_lengths(profile=profile, label=label)
        self.profiles.append(profile)
        self.i += 1
        if self.i == self.length:
            handles, labels = self.ax_muts_unt.get_legend_handles_labels()
            if self.length == 1:
                labels = ["Modified", "Untreated"]
            self.ax_muts_unt.legend(handles, labels, title="Samples", loc=1)
            if self.length == 1:
                labels = [label]
            self.make_boxplot(labels)

    def plot_mutations_per_molecule(self, profile, label, upper_limit=12):
        """Plot the mutations per molecule.

        Parameters
        ----------
        profile : rnavigate.profile.Profile
            Profile object.
        label : str
            Label for the sample.
        upper_limit : int, optional
            Upper limit of the x-axis.
        """
        df = profile.mutations_per_molecule
        if df is None:
            raise ValueError("profile is missing QC data from log file")
        x = df.loc[:upper_limit, "Mutation_count"]
        y1 = df.loc[:upper_limit, "Modified_mutations_per_molecule"]
        self.ax_muts_mod.plot(x, y1, label=label)
        y2 = df.loc[:upper_limit, "Untreated_mutations_per_molecule"]
        self.ax_muts_unt.plot(x, y2, label=label)

    def plot_read_lengths(self, profile, label, upper_limit=12):
        """Plot the read lengths.

        Parameters
        ----------
        profile : rnavigate.profile.Profile
            Profile object.
        label : str
            Label for the sample.
        upper_limit : int, optional
            Upper limit of the x-axis.
        """
        df = profile.read_lengths
        if self.length == 1:
            width = 0.4
            x1 = np.arange(upper_limit) - 0.6 + width
            x2 = x1 + width
        else:
            width = 0.8 / self.length
            x1 = np.arange(upper_limit) - 0.4 - (width / 2) + (width * self.i)
            x2 = x1
        y1 = df.loc[: upper_limit - 1, "Modified_read_length"]
        self.ax_read_mod.bar(x1, y1, width, label=label)
        y2 = df.loc[: upper_limit - 1, "Untreated_read_length"]
        self.ax_read_unt.bar(x2, y2, width, label=label)

    def make_boxplot(self, labels):
        """Make the boxplot of mutation rates.

        Parameters
        ----------
        labels : list of str
            Labels for the samples.
        """
        cols = ["Modified_rate", "Untreated_rate"]
        data = []
        for i, profile in enumerate(self.profiles):
            data.append(profile.data[cols].copy().assign(Sample=i + 1))
        data = pd.concat(data)
        data = pd.melt(data, id_vars=["Sample"], var_name=["Rate"])
        ax = sns.violinplot(
            x="Sample",
            y="value",
            hue="Rate",
            data=data,
            ax=self.ax_boxplot,
            split=True,
            scale="width",
        )
        ax.set(
            #    yscale='log',
            #    ylim=(0.00003, 0.3),
            ylabel="Mutation Rate",
            xticklabels=labels,
        )
