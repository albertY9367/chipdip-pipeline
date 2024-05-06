from collections import Counter
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import pysam
import glob
import os
import argparse
from collections import defaultdict

mpl.use("Agg")

"""
Generate maximum representation ecdfs for bead type representation within clusters.
"""


def main():
    args = parse_arguments()
    search = args.directory + "/*" + args.pattern
    files = glob.glob(search)
    ecdf_plot_ax = "None"
    for f in files:
        ecdf_plot_ax = max_representation_ecdf(f, ecdf_plot_ax)
    ecdf_plot_fig = ecdf_plot_ax.get_figure()
    ecdf_plot_fig.savefig(
        args.directory + "/Max_representation_ecdf.pdf", bbox_inches="tight"
    )
    ecdf_counts_ax = "None"
    for f in files:
        ecdf_counts_ax = max_representation_ecdf_counts(f, ecdf_counts_ax, args.xlim)
    ecdf_counts_fig = ecdf_counts_ax.get_figure()
    ecdf_counts_fig.savefig(
        args.directory + "/Max_representation_counts.pdf", bbox_inches="tight"
    )

def max_representation_ecdf(bamfile, ax):
    """
    Plot maximum representation ecdf for bead representation within clusters of a single clusterfile

    Args:
        clusterfile(str): path to clusterfile
        ax(obj): axis to plot on
    """

    results = []
    bamname = bamfile.replace(".bam", "")
    beads_dict = defaultdict(set)
    with pysam.AlignmentFile("bamfile", "rb") as inbam:
        for line in inbam.fetch(until_eof=True):
            barcode = line.get_tag("BC")
            read_type = line.get_tag("RT")
            chromesome = line.reference_name
            if "BEAD" in read_type or "BPM" in read_type:
                beads_dict[barcode].add(chromesome)
    for bc in beads_dict.keys():
        if len(beads_dict[bc]) > 1:
            candidate = Counter(beads_dict[bc]).most_common()[0]
            results.append(candidate[1] / len(beads_dict[bc]))
    if ax == "None":
        ax = plt.figure().subplots()
    ax = sns.ecdfplot(results, linewidth=3, ax=ax, label=bamname)
    ax.set(
        xlabel="Maximum Bead Representation Proportion", ylabel="Proportion of Beads"
    )
    ax.legend()
    return ax

def max_representation_ecdf_counts(bamfile, ax, xlimit):
    """
    Plot maximum representation ecdf for bead representation within clusters of a single clusterfile

    Args:
        clusterfile(str): path to clusterfile
        ax(obj): axis to plot on
    """

    results = []
    bamname = bamfile.replace(".bam", "")
    beads_dict = defaultdict(set)
    with pysam.AlignmentFile("bamfile", "rb") as inbam:
        for line in inbam.fetch(until_eof=True):
            barcode = line.get_tag("BC")
            read_type = line.get_tag("RT")
            chromesome = line.reference_name
            if "BEAD" in read_type or "BPM" in read_type:
                beads_dict[barcode].add(chromesome)
    for bc in beads_dict.keys():
        if len(beads_dict[bc]) > 1:
            candidate = Counter(beads_dict[bc]).most_common()[0]
            results.append(candidate[1])
    if ax == "None":
        ax = plt.figure().subplots()
    ax = sns.ecdfplot(results, linewidth=3, ax=ax, label=bamname)
    ax.set(xlabel="Num Oligos for Most Common Type", ylabel="Proportion of Beads")
    ax.set(xlim=(0, int(xlimit)))
    ax.legend()
    return ax

def max_representation_ecdf_counts(bamfile, ax, xlimit):
    """
    Plot counts of BPM reads with maximum representation

    Args:
        bamfile(str): path to bamfile
        ax(obj): axis to plot on
        xlimit(int): maximum x value to show in plot
    """

    results = []
    bamname = bamfile.replace(".bam", "")
    current_barcode = ""
    current_beads = []
    with pysam.AlignmentFile("bamfile", "rb") as inbam:
        for line in inbam.fetch(until_eof=True):
            barcode = line.get_tag("BC")
            read_type = line.get_tag("RT")
            chromesome = line.reference_name
            if barcode != current_barcode:
                if len(current_beads) > 1:
                    candidate = Counter(current_beads).most_common()[0]
                    results.append(candidate[1])
                current_barcode = barcode
                current_beads = []
            if read_type == "BEAD" or read_type == "BPM":
                current_beads.append(chromesome)
    if ax == "None":
        ax = plt.figure().subplots()
    ax = sns.ecdfplot(results, linewidth=3, ax=ax, label=bamname)
    ax.set(xlabel="Num Oligos for Most Common Type", ylabel="Proportion of Beads")
    ax.set(xlim=(0, int(xlimit)))
    ax.legend()
    return ax

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate the maximum representation ecdf plots to check for bead type uniqueness within clusters."
    )
    parser.add_argument(
        "--directory",
        metavar="FILE",
        action="store",
        required=True,
        help="The directory of clusters file",
    )
    parser.add_argument(
        "--pattern",
        action="store",
        required=True,
        help="The pattern of cluster file names",
    )
    parser.add_argument(
        "--xlim",
        action="store",
        required=True,
        help="The maximum x value on the counts plot",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
