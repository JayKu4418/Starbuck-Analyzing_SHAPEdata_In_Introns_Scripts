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
#pwd
for file in *_Aligned.out.sortedByCoord.bam
do
    #echo ${file}
    nameFile=$(echo "$file" | cut -d"." -f1)
    realNameFile_withLastUnderscore=${nameFile//Aligned/}
    realNameFile=${realNameFile_withLastUnderscore::-1}
    echo ${realNameFile}
    #echo "aligning....."
    # The command below is the ShapeMapper 1 alignment parameters
    #bowtie2 -p 12 --local -D 20 -R 3 -N 1 -L 15 -i S,1,0.50 --score-min G,20,8 --ma 2 --mp 6,2 --rdg 5,1 --rfg 5,1 --dpad 100 --maxins 1000 -x ${datafolder}/hg38_bowtie2_index/hg38 -1 ${fastafolder}/${realNameFile}_L001_R1_001.fastq -2 ${fastafolder}/${realNameFile}_L001_R2_001.fastq -S ${folderToWriteTo}/${realNameFile}.sam > ${folderToWriteTo}/${realNameFile}.log
    # The command below is the ShapeMapper 2 alignment parameters
    #bowtie2 -p 12 --local --sensitive-local --mp 3,1 --rdg 5,1 --rfg 5,1 --dpad 30 --maxins 800 --ignore-quals -x ${datafolder}/hg38_bowtie2_index/hg38 -1 ${fastafolder}/${realNameFile}_L001_R1_001.fastq -2 ${fastafolder}/${realNameFile}_L001_R2_001.fastq -S ${folderToWriteTo}/${realNameFile}.sam > ${folderToWriteTo}/${realNameFile}.log
    #echo "converting SAM to BAM........"
    #samtools view -b -o ${folderToWriteTo}/${realNameFile}.bam ${folderToWriteTo}/${realNameFile}.sam
    
    #echo "Get unaligned reads....."
    #samtools view -b -f 4 -o ${folderToWriteTo}/${realNameFile}_UnalignedReads.bam ${folderToWriteTo}/${realNameFile}.bam
    #bedtools bamtofastq -i ${folderToWriteTo}/${realNameFile}_UnalignedReads.bam -fq ${folderToWriteTo}/${realNameFile}_UnalignedReads.fastq
    
    #echo "Get aligned reads....." 
    #samtools view -b -F 4 -o ${folderToWriteTo}/${realNameFile}_AlignedReads.bam ${folderToWriteTo}/${realNameFile}.bam
    #bedtools bamtofastq -i ${folderToWriteTo}/${realNameFile}_AlignedReads.bam -fq ${folderToWriteTo}/${realNameFile}_AlignedReads.fastq
    echo "Obtain reads with MAPQ greater than 10"
    samtools view -b -q 10 -o ${folderToWriteTo}/${realNameFile}_MAPQgreaterThan10.bam ${folderToWriteTo}/${realNameFile}.bam
    
    echo "sorting and index BAM........"
    samtools sort -@ 4 -o ${folderToWriteTo}/${realNameFile}_MAPQgreaterThan10_SortedByCoords.bam ${folderToWriteTo}/${realNameFile}_MAPQgreaterThan10.bam

    #samtools index ${folderToWriteTo}/${realNameFile}Aligned_MAPQGreater1_Sorted.bam
done

# samtools merge AllDMS_Madrasin_MAPQgreaterThan10_SortedByCoords.bam Mad-10mm-4hr-DMS_S1_MAPQgreaterThan10_SortedByCoords.bam Mad-10mm-8hr-DMS_S7_MAPQgreaterThan10_SortedByCoords.bam Mad-20mm-4hr-DMS_S2_MAPQgreaterThan10_SortedByCoords.bam Mad-20mm-8hr-DMS_S8_MAPQgreaterThan10_SortedByCoords.bam Mad-30mm-4hr-DMS_S3_MAPQgreaterThan10_SortedByCoords.bam Mad-30mm-8hr-DMS_S9_MAPQgreaterThan10_SortedByCoords.bam
# samtools merge AllMinus_Madrasin_MAPQgreaterThan10_SortedByCoords.bam Mad-10mm-4hr-minus_S4_MAPQgreaterThan10_SortedByCoords.bam Mad-10mm-8hr-minus_S10_MAPQgreaterThan10_SortedByCoords.bam Mad-20mm-4hr-minus_S5_MAPQgreaterThan10_SortedByCoords.bam Mad-20mm-8hr-minus_S11_MAPQgreaterThan10_SortedByCoords.bam Mad-30mm-4hr-minus_S6_MAPQgreaterThan10_SortedByCoords.bam Mad-30mm-8hr-minus_S12_MAPQgreaterThan10_SortedByCoords.bam
cd ../../scripts