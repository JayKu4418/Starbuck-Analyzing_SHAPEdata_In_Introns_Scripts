#!/bin/bash

filename=$1
end=$2
echo $filename
datafolder="../data"
processeddatafolder="../processed_data"

echo "aligning....."
echo ${filename}_R1.fastq
bowtie2 -p 4 --local -D 20 -R 3 -N 1 -L 15 -i S,1,0.50 --score-min G,20,8 --ma 2 --mp 6,2 --rdg 5,1 --rfg 5,1 --dpad 100 -x ${processeddatafolder}/Bowtie2Index/Valid_exonIntronBoundaries_${end}pIntron -1 ${datafolder}/${filename}_R1.fastq -2 $datafolder/${filename}_R2.fastq -S ${processeddatafolder}/${filename}_exonIntronBoundaries_${end}p.sam

echo "converting SAM to BAM........"
# Only get reads that have a MAPQ score >= 1 and are aligned
samtools view -F 4 -b -o ${processeddatafolder}/${filename}_exonIntronBoundaries_${end}p_Aligned.bam ${processeddatafolder}/${filename}_exonIntronBoundaries_${end}p.sam

echo "sorting and index BAM........"
samtools sort -@ 4 -o ${processeddatafolder}/${filename}_exonIntronBoundaries_${end}p_Aligned_Sorted.bam ${processeddatafolder}/${filename}_exonIntronBoundaries_${end}p_Aligned.bam

samtools index ${processeddatafolder}/${filename}_exonIntronBoundaries_${end}p_Aligned_Sorted.bam

echo "converting BAM to FASTQ.........."
samtools fastq ${processeddatafolder}/${filename}_exonIntronBoundaries_${end}p_Aligned_Sorted.bam -1 ${processeddatafolder}/${filename}_exonIntronBoundaries_${end}p_Aligned_Sorted_R1.fastq -2 ${processeddatafolder}/${filename}_exonIntronBoundaries_${end}p_Aligned_Sorted_R2.fastq


