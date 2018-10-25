#!/bin/bash

# This script is going to get featureCounts using STAR alignment

folderToCD="/home/jkumar12/Projects/Analyzing_SHAPEdata_In_Introns/tmp/Lela_20uM_MadrasinReadsTest_August2018" 

cd ${folderToCD}
for file in "Madrasin-20uM-4hr-Minus-1_S2" "Madrasin-20uM-4hr-Minus-2_S5" "Madrasin-20uM-4hr-Plus-1_S1" "Madrasin-20uM-4hr-Plus-2_S6" "Madrasin-20uM-8hr-Minus_S4" "Madrasin-20uM-8hr-Plus_S3" "Untreated-Minus_S7" "Untreated-Plus_S8";
do
    echo ${file}
    # Get proper pairs 
    samtools view -b -f 2 -o ${file}_ProperPairs.bam ${file}.bam
    # Sort proper pairs
    samtools sort -o ${file}_ProperPairs_SortedByCoords.bam ${file}.bam
    # Get fragment length distribution using Picard tools
    java -jar $PICARD CollectInsertSizeMetrics I=${file}_ProperPairs_SortedByCoords.bam O=${file}_insert_size_metrics.txt H=${file}_insert_size_histogram.pdf M=0.5
done

cd /home/jkumar12/Projects/Analyzing_SHAPEdata_In_Introns/scripts