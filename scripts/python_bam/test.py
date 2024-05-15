import os
from tag_bam import *

# input = '/home/zyang4/chipdip-pipeline/workup2/alignments/sample1.merged.BPM.bam'
# output = '/home/zyang4/chipdip-pipeline/workup2/alignments/sample1.tagged.BPM.bam'
# label_bam_file(input, output, 7)

os.system("python '/home/zyang4/chipdip-pipeline/scripts/python_bam/tag_bam.py' --input_bam '/home/zyang4/chipdip-pipeline/workup2/alignments/sample1.merged.BPM.bam' --output_bam '/home/zyang4/chipdip-pipeline/workup2/alignments/sample1.tagged.BPM.bam' --num_tags '7'")