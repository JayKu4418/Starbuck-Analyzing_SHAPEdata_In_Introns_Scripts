#!/bin/bash

fastaFolder="/home/shared/SequencingData/HiSeq_YorubanLines" 

datafolder="/home/jkumar12/Projects/Analyzing_SHAPEdata_In_Introns/data"

folderToCD="/home/jkumar12/Projects/Analyzing_SHAPEdata_In_Introns/tmp/HiSeq_19098-Yoruban_cell_lines_SHAPE"

cd ${folderToCD}

for file in "19098DC-Yoruban_cell_lines_SHAPE_TCCTGAG-GCGTAAG" "19098MIN-Yoruban_cell_lines_SHAPE_TAAGGCG-GCGTAAG" "19098PLUS-Yoruban_cell_lines_SHAPE_CGTACTA-GCGTAAG";
do
    echo ${file}
    
    echo "aligning....."
    # The command below is the ShapeMapper 1 alignment parameters
    #bowtie2 -p 12 --local -D 20 -R 3 -N 1 -L 15 -i S,1,0.50 --score-min G,20,8 --ma 2 --mp 6,2 --rdg 5,1 --rfg 5,1 --dpad 100 --maxins 1000 -x ${datafolder}/hg38_bowtie2_index/hg38 -1 ${file}_L001_R1_001.fastq.gz -2 ${file}_L001_R2_001.fastq.gz -S ${folderToWriteTo}/${file}.sam > ${folderToWriteTo}/${file}.log
    # The command below is the ShapeMapper 2 alignment parameters
    bowtie2 -p 12 --local --sensitive-local --mp 3,1 --rdg 5,1 --rfg 5,1 --dpad 30 --maxins 800 --ignore-quals -x ${datafolder}/hg38_bowtie2_index/hg38 -1 ${fastaFolder}/${file}_L001_R1_001.fastq.gz -2 ${fastaFolder}/${file}_L001_R2_001.fastq.gz -S ${file}.sam > ${file}.log
    
    echo "converting SAM to BAM........"
    samtools view -b -o ${file}.bam ${file}.sam
    
    echo "getting proper pairs........"
    # Get proper pairs 
    samtools view -b -f 2 -o ${file}_ProperPairs.bam ${file}.bam
    
    echo "sorting proper pairs........"
    # Sort proper pairs
    samtools sort -o ${file}_ProperPairs_SortedByCoords.bam ${file}.bam
    
    echo "get fragment length distribution........"
    # Get fragment length distribution using Picard tools
    java -jar $PICARD CollectInsertSizeMetrics I=${file}_ProperPairs_SortedByCoords.bam O=${file}_insert_size_metrics.txt H=${file}_insert_size_histogram.pdf M=0.5
done

cd /home/jkumar12/Projects/Analyzing_SHAPEdata_In_Introns/scripts