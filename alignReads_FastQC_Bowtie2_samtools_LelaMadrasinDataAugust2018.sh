#!/bin/bash

folderToCD="/home/shared/SequencingData/Laederach_LL_AM_Madrasin-90563725" 

datafolder="/home/jkumar12/Projects/Analyzing_SHAPEdata_In_Introns/data"

folderToWriteTo="/home/jkumar12/Projects/Analyzing_SHAPEdata_In_Introns/tmp/Lela_20uM_MadrasinReadsTest_August2018"

cd ${folderToCD}
#pwd
for file in "Madrasin-20uM-4hr-Minus-1_S2" "Madrasin-20uM-4hr-Minus-2_S5" "Madrasin-20uM-4hr-Plus-1_S1" "Madrasin-20uM-4hr-Plus-2_S6" "Madrasin-20uM-8hr-Minus_S4" "Madrasin-20uM-8hr-Plus_S3" "Untreated-Minus_S7" "Untreated-Plus_S8";
do
    echo ${file}
    
    echo "unzipping fastq files..."
    gzip -d ${file}_L001_R*_001.fastq.gz
    
    echo "run fastqc..."
    fastqc ${file}_L001_R*_001.fastq
    # Move fastqc files into fastqc folder 
    mv ${file}_L001_R*_001_fastqc fastqc_files
    
    echo "zipping fastq files..."
    gzip ${file}_L001_R*_001.fastq
    
    echo "aligning....."
    # The command below is the ShapeMapper 1 alignment parameters
    #bowtie2 -p 12 --local -D 20 -R 3 -N 1 -L 15 -i S,1,0.50 --score-min G,20,8 --ma 2 --mp 6,2 --rdg 5,1 --rfg 5,1 --dpad 100 --maxins 1000 -x ${datafolder}/hg38_bowtie2_index/hg38 -1 ${file}_L001_R1_001.fastq.gz -2 ${file}_L001_R2_001.fastq.gz -S ${folderToWriteTo}/${file}.sam > ${folderToWriteTo}/${file}.log
    # The command below is the ShapeMapper 2 alignment parameters
    bowtie2 -p 12 --local --sensitive-local --mp 3,1 --rdg 5,1 --rfg 5,1 --dpad 30 --maxins 800 --ignore-quals -x ${datafolder}/hg38_bowtie2_index/hg38 -1 ${file}_L001_R1_001.fastq.gz -2 ${file}_L001_R2_001.fastq.gz -S ${folderToWriteTo}/${file}.sam > ${folderToWriteTo}/${file}.log
    
    echo "converting SAM to BAM........"
    samtools view -b -o ${folderToWriteTo}/${file}.bam ${folderToWriteTo}/${file}.sam
    
    echo "Get unaligned reads....."
    samtools view -b -f 4 -o ${folderToWriteTo}/${file}_UnalignedReads.bam ${folderToWriteTo}/${file}.bam
    bedtools bamtofastq -i ${folderToWriteTo}/${file}_UnalignedReads.bam -fq ${folderToWriteTo}/${file}_UnalignedReads.fastq
    
    echo "Get aligned reads....." 
    samtools view -b -F 4 -o ${folderToWriteTo}/${file}_AlignedReads.bam ${folderToWriteTo}/${file}.bam
    bedtools bamtofastq -i ${folderToWriteTo}/${file}_AlignedReads.bam -fq ${folderToWriteTo}/${file}_AlignedReads.fastq
    
    echo "Obtain reads with MAPQ greater than 10"
    samtools view -b -q 10 -o ${folderToWriteTo}/${file}_MAPQgreaterThan10.bam ${folderToWriteTo}/${file}.bam
    
    echo "sorting and index BAM........"
    samtools sort -@ 4 -o ${folderToWriteTo}/${file}_MAPQgreaterThan10_SortedByCoords.bam ${folderToWriteTo}/${file}_MAPQgreaterThan10.bam

    samtools index ${folderToWriteTo}/${file}_MAPQgreaterThan10_SortedByCoords.bam
done

cd /home/jkumar12/Projects/Analyzing_SHAPEdata_In_Introns/scripts