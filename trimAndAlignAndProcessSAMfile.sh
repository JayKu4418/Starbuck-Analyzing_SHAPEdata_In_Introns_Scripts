# This script will align the 

folderToRead=$1

name=$2

datafolderToSave=$3

tmpfolderToSave=$4

# First unzip files
#echo "Unzipping files"
#gzip -d -k -c ${folderToRead}/${name}_L001_R1_001.fastq.gz > ${datafolderToSave}/${name}_L001_R1_001.fastq
#gzip -d -k -c ${folderToRead}/${name}_L001_R2_001.fastq.gz > ${datafolderToSave}/${name}_L001_R2_001.fastq


# Align the fastq files
echo "Aligning fastq files using STAR align"
STAR --runThreadN 12 --runMode alignReads --outMultimapperOrder Random --outSAMmultNmax 1 --outSAMattributes MD --genomeDir /home/shared/hg38/Long_RefSeqAnnotated_hg38_STARindex --readFilesIn ${datafolderToSave}/${name}_L001_R1_001.fastq ${datafolderToSave}/${name}_L001_R2_001.fastq --outTmpDir ${tmpfolderToSave}/${name}_templog --outFileNamePrefix ${tmpfolderToSave}/${name}_ --alignEndsType Local --scoreGap -1000000 --scoreDelBase -1 --scoreInsBase -1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 999

# Convert the sam file to bam file
echo "Convert sam file to bam file"
samtools view -b -o ${tmpfolderToSave}/${name}_Aligned.out.bam ${tmpfolderToSave}/${name}_Aligned.out.sam

# Sort the bam file by coordinate 
echo "Sorting bam file"
samtools sort -@ 6 -o ${tmpfolderToSave}/${name}_Aligned.out.sortedByCoord.bam ${tmpfolderToSave}/${name}_Aligned.out.bam

# Index the bam file
echo "Index bam file"
samtools index ${tmpfolderToSave}/${name}_Aligned.out.sortedByCoord.bam