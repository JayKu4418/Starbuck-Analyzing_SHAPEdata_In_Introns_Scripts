#!/bin/bash

# This script is going to get featureCounts using STAR alignment

folderToCD=$1 

if [ "$2" = "Transcriptome" ]
then
    annotationFile="SAFformat_GRCh38_latest_genomic_ValidChroms_OnlyGenes.bed"
elif [ "$2" = "Genome" ]
then
    annotationFile="SAFformat_ExtractedGenomicCoordinatesFromGFFfile_OnlyGenes.bed"
elif [ "$2" = "ExonIntron5p" ]
then
    annotationFile="SAFformat_ExtractedExonIntronCoordinatesFromGFFfile_5p_OnlyGenes.bed"
elif [ "$2" = "ExonIntron3p" ]
then
    annotationFile="SAFformat_ExtractedExonIntronCoordinatesFromGFFfile_3p_OnlyGenes.bed"
fi

name=$3

fastafolder="/home/shared/SequencingData/MiSeq_LL_06142018_MadrasinTrial/FASTQ_Generation_2018-06-15_00_25_47Z-103862779"

cd ${folderToCD}
pwd

#STAR --runThreadN 12 --runMode alignReads --outMultimapperOrder Random --outSAMmultNmax 1 --outSAMattributes MD --genomeDir /home/shared/GenomeReference/Long_RefSeqAnnotated_hg38_STARindex --readFilesIn ${fastafolder}/${name}_R1.fastq ${fastafolder}/${name}_R2.fastq --outTmpDir ${name}_templog --outFileNamePrefix ${name}_ --alignEndsType Local --scoreGap -1000000 --scoreDelBase -1 --scoreInsBase -1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 999

#samtools view -b -o ${name}_Aligned.out.bam ${name}_Aligned.out.sam

#samtools sort -@ 4 -n -o ${name}_Aligned.out.sortedByName.bam ${name}_Aligned.out.bam

featureCounts -p --largestOverlap -B -C -M -R SAM -t exon -g gene_id -a ../../data/${annotationFile} -F SAF -o ${name}_${2}_featureCounts.txt ${name}_Aligned.out.sortedByName.bam
cut -f1,6,7 ${name}_${2}_featureCounts.txt > ${name}_${2}_featureCounts_JustCounts.txt
sort -k3,3nr ${name}_${2}_featureCounts_JustCounts.txt > ${name}_${2}_featureCounts_JustCounts.sorted.txt

cd ../../scripts