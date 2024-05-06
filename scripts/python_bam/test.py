import merge_bam

in_file = "/home/zyang4/chipdip-pipeline/workup/alignments/sample1.DNA.merged.bam"
out_file = "/home/zyang4/chipdip-pipeline/scripts/python_bam/test0503.bam"
merge_bam.label_bam_file(in_file, out_file, 7)