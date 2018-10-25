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
    # Grab header for file
    samtools view -H -o ${file}_JustHeader.sam ${file}.bam
    # Grab mapped reads
    samtools view -F 4 -o ${file}_MappedReads.sam ${file}.bam
    # Get multi mapped reads and single mapped reads 
    grep XS:i ${file}_MappedReads.sam > ${file}_MultiMappedReads.sam
    grep -v XS:i ${file}_MappedReads.sam > ${file}_SingleMappedReads.sam
    # Concatenate Header file with multimapped reads and get bam file
    cat ${file}_JustHeader.sam ${file}_MultiMappedReads.sam > ${file}_MultiMappedReadsPlusHeader.sam
    samtools view -b -o ${file}_MultiMappedReads.bam ${file}_MultiMappedReadsPlusHeader.sam
    # Remove sam files for multimapped reads
    rm ${file}_MultiMappedReads.sam
    rm ${file}_MultiMappedReadsPlusHeader.sam
    # Concatenate Header file with single mapped reads and get bam file
    cat ${file}_JustHeader.sam ${file}_SingleMappedReads.sam > ${file}_SingleMappedReadsPlusHeader.sam
    samtools view -b -o ${file}_SingleMappedReads.bam ${file}_SingleMappedReadsPlusHeader.sam
    # Remove sam files for singlemapped reads
    rm ${file}_SingleMappedReads.sam
    rm ${file}_SingleMappedReadsPlusHeader.sam
    # Remove header sam file
    rm ${file}_JustHeader.sam
    
    featureCounts -p --largestOverlap -B -C -M -R SAM -t exon -g gene_id -a ${dataFolder}/${annotationFile} -F SAF -o ${file}_${1}_featureCounts_SingleMappedReads.txt ${file}_SingleMappedReads.bam
    cut -f1,6,7 ${file}_${1}_featureCounts_SingleMappedReads.txt > ${file}_${1}_featureCounts_SingleMappedReads_JustCounts.txt
    
    featureCounts -p --largestOverlap -B -C -M -R SAM -t exon -g gene_id -a ${dataFolder}/${annotationFile} -F SAF -o ${file}_${1}_featureCounts_MultiMappedReads.txt ${file}_MultiMappedReads.bam
    cut -f1,6,7 ${file}_${1}_featureCounts_MultiMappedReads.txt > ${file}_${1}_featureCounts_MultiMappedReads_JustCounts.txt
    
done

cd /home/jkumar12/Projects/Analyzing_SHAPEdata_In_Introns/scripts