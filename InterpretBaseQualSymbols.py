# This is small python script to interpret plot base quality scores of a single sequence

from Bio import SeqIO
import sys

seqfile=sys.argv[1]
seqname=sys.argv[2]

for record in SeqIO.parse(seqfile, "fastq"):
        if seqname in record.id:
            print record.seq
            print record.letter_annotations["phred_quality"]
            break