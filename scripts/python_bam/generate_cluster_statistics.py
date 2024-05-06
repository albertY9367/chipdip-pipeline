import tqdm
import glob
import argparse
import pysam

"""
Count 1) number of clusters, 2) number of DPM reads (aligned) and 3) number of  BPM reads within each clusterfile for a directory of clusterfiles.
"""


def main():
    args = parse_arguments()
    search = args.directory + "/*" + args.pattern
    files = glob.glob(search)
    for f in files:
        count_statistics(f)

def count_statistics(bamfile):
    """
    Loop through all clusters within a bamfile, counting DPM and BPM reads

    Args:
        bamfile(str): Path to bamfile
    """
    cluster = 0
    dpm = 0
    bpm = 0
    cluster_barcodes = []
    with pysam.AlignmentFile(bamfile, 'rb') as inbam:
        for read in inbam.fetch(until_eof=True):
            read_type = read.get_tag("RT")
            if "DPM" in read_type:
                dpm += 1
            elif "BEAD" in read_type or "BPM" in read_type:
                bpm += 1
            barcode = read.get_tag("BC")
            if barcode not in cluster_barcodes:
                cluster += 1
                cluster_barcodes.append(barcode)
    print("For bamfile ", bamfile)
    print("Total number of clusters: ", cluster)
    print("Total number of BPM: ", bpm)
    print("Total number of DPM: ", dpm)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate the statistics for all clusterfiles in a directory"
    )
    parser.add_argument(
        "--directory",
        metavar="FILE",
        action="store",
        help="The directory of clusters file",
    )
    parser.add_argument(
        "--pattern", action="store", help="The pattern of cluster file names"
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
