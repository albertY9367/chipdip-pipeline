import pysam
import re
from collections import defaultdict
'''
Python code to add read type and barcode tags to input bam file
'''


def label_bam_file(input_bam, output_bam, num_tags):
    """
    Add antibody label to individual reads of the master DNA bam file

    Args:
        input_bam(str): Path to input master bam file
        output_bam(str): Path to write labeled bam file
        num_tags(int): number of tags in barcode
    """
    count, written, duplicates, skipped = 0, 0, 0, 0
    pattern = re.compile("::" + num_tags * "\[([a-zA-Z0-9_\-]+)\]")
    rt_pattern = re.compile(r"RPM|BPM|DPM|BEAD")
    found = defaultdict(set)
    with pysam.AlignmentFile(input_bam, "rb") as in_bam, \
         pysam.AlignmentFile(output_bam, "wb", template=in_bam) as out_bam:
        for read in in_bam.fetch(until_eof=True):
            count += 1
            if count % 10000000 == 0:
                print(count)
            name = read.query_name
            match = pattern.search(name)
            full_barcode = list(match.groups())
            read_type = rt_pattern.findall(full_barcode[0])[0]
            barcode = full_barcode[1:]
            ref_barcode = ".".join(barcode)
            fc_barcode = ".".join(full_barcode)
            position = read.reference_name + ":" + str(read.reference_start) + '-' + str(read.reference_end)
            if position in found[ref_barcode]:
                duplicates += 1
            else:
                try:
                    found[ref_barcode].add(position)
                    read.set_tag("RT", read_type, replace=True)
                    read.set_tag("BC", ref_barcode, replace=True)
                    read.set_tag("FC", fc_barcode, replace=True)
                    out_bam.write(read)
                    written += 1
                except KeyError:
                    skipped += 1

    print("Total reads:", count)
    print("Reads written:", written)
    print("Duplicate reads:", duplicates)
    print("Reads with an error not written out:", skipped)