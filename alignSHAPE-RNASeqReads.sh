#!/bin/bash

#filename=$1
#echo $filename
#datafolder="../data"
#processeddatafolder="../processed_data"

folderToCD=$1 

datafolder="/home/jkumar12/Projects/Analyzing_SHAPEdata_In_Introns/data"

fastafolder="/home/shared/SequencingData/MiSeq_LL_06142018_MadrasinTrial/FASTQ_Generation_2018-06-15_00_25_47Z-103862779"

folderToWriteTo=$2

cd ${folderToCD}
pwd
for file in *_Aligned.out.sortedByCoord.bam
do
    echo ${file}
    nameFile=$(echo "$file" | cut -d"." -f1)
    realNameFile=${nameFile//Aligned/}
    echo ${realNameFile}
    echo "aligning....."
    #echo ${filename}_R1.fastq
    bowtie2 -p 12 --local -D 20 -R 3 -N 1 -L 15 -i S,1,0.50 --score-min G,20,8 --ma 2 --mp 6,2 --rdg 5,1 --rfg 5,1 --dpad 100 -x ${datafolder}/hg38_bowtie2_index/hg38 -1 ${fastafolder}/${realNameFile}L001_R1_001.fastq -2 ${fastafolder}/${realNameFile}L001_R2_001.fastq -S ${folderToWriteTo}/${realNameFile}.sam
#bowtie2 -p 12 --local -D 20 -R 3 -N 1 -L 15 -i S,1,0.50 --score-min G,20,8 --ma 2 --mp 6,2 --rdg 5,1 --rfg 5,1 --dpad 100 --maxins 1000 -x ${datafolder}/hg38_index/hg38 -1 ${datafolder}/${filename}_R1.fastq -2 $datafolder/${filename}_R2.fastq -S ${processeddatafolder}/${filename}.sam

    echo "converting SAM to BAM........"
# Only get reads that have a MAPQ score >= 1 and are aligned
    samtools view -F 4 -q 1 -b -o ${folderToWriteTo}/${realNameFile}_Aligned_MAPQGreater1.bam ${folderToWriteTo}/${realNameFile}.sam

    echo "sorting and index BAM........"
    samtools sort -@ 4 -o ${processeddatafolder}/${filename}_Aligned_MAPQGreater1_Sorted.bam ${processeddatafolder}/${filename}_Aligned_MAPQGreater1.bam

samtools index ${processeddatafolder}/${filename}_Aligned_MAPQGreater1_Sorted.bam

echo "intersecting BAM with exon-intron boundaries BED........"
bedtools intersect -abam ${processeddatafolder}/${filename}_Aligned_MAPQGreater1_Sorted.bam -b ${processeddatafolder}/Valid_exonIntronBoundaries_5pIntron.bed -wa > ${processeddatafolder}/${filename}_Aligned_MAPQGreater1_Sorted_ExonIntronBoundaries_5p.bam

bedtools intersect -abam ${processeddatafolder}/${filename}_Aligned_MAPQGreater1_Sorted.bam -b ${processeddatafolder}/Valid_exonIntronBoundaries_3pIntron.bed -wa > ${processeddatafolder}/${filename}_Aligned_MAPQGreater1_Sorted_ExonIntronBoundaries_3p.bam

samtools index ${processeddatafolder}/${filename}_Aligned_MAPQGreater1_Sorted_ExonIntronBoundaries_5p.bam

samtools index ${processeddatafolder}/${filename}_Aligned_MAPQGreater1_Sorted_ExonIntronBoundaries_3p.bam

echo "converting BAM to FASTQ.........."
samtools fastq ${processeddatafolder}/${filename}_Aligned_MAPQGreater1_Sorted_ExonIntronBoundaries_5p.bam -1 ${processeddatafolder}/${filename}_Aligned_MAPQGreater1_Sorted_ExonIntronBoundaries_5p_R1.fastq -2 ${processeddatafolder}/${filename}_Aligned_MAPQGreater1_Sorted_ExonIntronBoundaries_5p_R2.fastq

samtools fastq ${processeddatafolder}/${filename}_Aligned_MAPQGreater1_Sorted_ExonIntronBoundaries_3p.bam -1 ${processeddatafolder}/${filename}_Aligned_MAPQGreater1_Sorted_ExonIntronBoundaries_3p_R1.fastq -2 ${processeddatafolder}/${filename}_Aligned_MAPQGreater1_Sorted_ExonIntronBoundaries_3p_R2.fastq
