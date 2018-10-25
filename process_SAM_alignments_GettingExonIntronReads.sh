# This script is going to process the output sam file of an alignment file

nameFile=$1

# First convert the sam file to bam file
#echo "Converting sam file to bam file"
#samtools view -b -o ${nameFile}.bam ${nameFile}.sam

# Sort the bam file by read name
echo "Sorting bam file by read name"
samtools sort -@ 6 -n -o ${nameFile}.sorted.bam ${nameFile}.bam

# Only grab the paired reads 
echo "Get paired reads from sorted bam file"
samtools view -bf 2 -o ${nameFile}.sorted.pairedReads.bam ${nameFile}.sorted.bam

# Put not proper pairs into another file
echo "Get unpaired reads from sorted bam file"
samtools view -bF 2 -o ${nameFile}.sorted.unpairedReads.bam ${nameFile}.sorted.bam

# Convert the paired reads into a bed file
echo "Convert paired reads bam file into a bed file with CIGAR string"
bedtools bamtobed -i ${nameFile}.sorted.pairedReads.bam -cigar > ${nameFile}.sorted.pairedReads.bed

# Split bed file into 2 files based on read number 
echo "Splitting the bed file into 2 files based on read number"
awk '{if (substr($4,length($4))==1) print $0}' ${nameFile}.sorted.pairedReads.bed > ${nameFile}.sorted.pairedReads.read1.bed
awk '{if (substr($4,length($4))==2) print $0}' ${nameFile}.sorted.pairedReads.bed > ${nameFile}.sorted.pairedReads.read2.bed

# Paste the 2 files together side by side with a tab delimiter
paste -d '\t' ${nameFile}.sorted.pairedReads.read1.bed ${nameFile}.sorted.pairedReads.read2.bed > ${nameFile}.sorted.pairedReads.bothReads.txt

# Check for reads that both contain N in the CIGAR string 
# N represents the number of bases skipped if read is spliced 
echo "Grab lines where both reads contain N in their CIGAR string"
awk '{if (index($7,"N")!=0 && index($14,"N")!=0) print $0}' ${nameFile}.sorted.pairedReads.bothReads.txt > ${nameFile}.sorted.pairedReads.bothReads.BothN.txt
cut -f 1-7 ${nameFile}.sorted.pairedReads.bothReads.BothN.txt > ${nameFile}.sorted.pairedReads.Read1.BothN.bed
cut -f 8-14 ${nameFile}.sorted.pairedReads.bothReads.BothN.txt > ${nameFile}.sorted.pairedReads.Read2.BothN.bed

# Check for reads where neither contains N in CIGAR string
echo "Grab lines where both reads do NOT contain N in their CIGAR string"
awk '{if (index($7,"N")==0 && index($14,"N")==0) print $0}' ${nameFile}.sorted.pairedReads.bothReads.txt > ${nameFile}.sorted.pairedReads.bothReads.NeitherN.txt
cut -f 1-7 ${nameFile}.sorted.pairedReads.bothReads.NeitherN.txt > ${nameFile}.sorted.pairedReads.Read1.NeitherN.bed
cut -f 8-14 ${nameFile}.sorted.pairedReads.bothReads.NeitherN.txt > ${nameFile}.sorted.pairedReads.Read2.NeitherN.bed

# Check for reads where only 1 read contains N in CIGAR string
echo "Grab lines where one of the reads contains N in their CIGAR string"
awk '{if ((index($7,"N")==0 && index($14,"N")!=0) || (index($7,"N")!=0 && index($14,"N")==0)) print $0}' ${nameFile}.sorted.pairedReads.bothReads.txt > ${nameFile}.sorted.pairedReads.bothReads.SingleN.txt
cut -f 1-7 ${nameFile}.sorted.pairedReads.bothReads.SingleN.txt > ${nameFile}.sorted.pairedReads.Read1.SingleN.bed
cut -f 8-14 ${nameFile}.sorted.pairedReads.bothReads.SingleN.txt > ${nameFile}.sorted.pairedReads.Read2.SingleN.bed

# Sort the bam file by coordinate and indexing
echo "Sorting bam file by coordinates and then index "
samtools sort -@ 6 -o ${nameFile}.sortedByCoord.bam ${nameFile}.bam
samtools index ${nameFile}.sortedByCoord.bam