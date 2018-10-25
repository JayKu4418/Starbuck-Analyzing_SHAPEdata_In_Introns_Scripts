#!/bin/bash

# This script is going to get featureCounts using STAR alignment

folderToCD="/home/jkumar12/Projects/Analyzing_SHAPEdata_In_Introns/tmp/Lela_20uM_MadrasinReadsTest_August2018" 

dataFolder="/home/jkumar12/Projects/Analyzing_SHAPEdata_In_Introns/data"

if [ "$1" = "Transcriptome" ]
then
    annotationFile="SAFformat_GRCh38_latest_genomic_ValidChroms_OnlyGenes.bed"
elif [ "$1" = "Genome" ]
then
    annotationFile="SAFformat_ExtractedGenomicCoordinatesFromGFFfile_OnlyGenes.bed"
elif [ "$1" = "ExonIntron5p" ]
then
    annotationFile="SAFformat_ExtractedExonIntronCoordinatesFromGFFfile_5p_OnlyGenes.bed"
elif [ "$1" = "ExonIntron3p" ]
then
    annotationFile="SAFformat_ExtractedExonIntronCoordinatesFromGFFfile_3p_OnlyGenes.bed"
fi

cd ${folderToCD}
for file in "Madrasin-20uM-4hr-Minus-1_S2" "Madrasin-20uM-4hr-Minus-2_S5" "Madrasin-20uM-4hr-Plus-1_S1" "Madrasin-20uM-4hr-Plus-2_S6" "Madrasin-20uM-8hr-Minus_S4" "Madrasin-20uM-8hr-Plus_S3" "Untreated-Minus_S7" "Untreated-Plus_S8";
do
    echo ${file}
    featureCounts -p --largestOverlap -B -C -M -R SAM -t exon -g gene_id -a ${dataFolder}/${annotationFile} -F SAF -o ${file}_${1}_featureCounts.txt ${file}_MAPQgreaterThan10_SortedByCoords.bam
    cut -f1,6,7 ${file}_${1}_featureCounts.txt > ${file}_${1}_featureCounts_JustCounts.txt
done

cd /home/jkumar12/Projects/Analyzing_SHAPEdata_In_Introns/scripts