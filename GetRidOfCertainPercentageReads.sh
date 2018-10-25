#!/bin/bash

# This script is going to get rid of a certain percentage reads in the fastq file
#folderToCD=$1

percentageReads=$1

#datafolder="/home/jkumar12/Projects/Analyzing_SHAPEdata_In_Introns/data"

processeddatafolder="/home/jkumar12/Projects/Analyzing_SHAPEdata_In_Introns/processed_data"

fastafolder="/home/shared/SequencingData/MiSeq_LL_06142018_MadrasinTrial/FASTQ_Generation_2018-06-15_00_25_47Z-103862779"

#numFour=4

#numHundred=100

cd ${fastafolder}
pwd
for file in *.fastq
do
    echo ${file}
    nameFile=$(echo "$file" | cut -d"." -f1)
    echo ${nameFile}
    #realNameFile=${nameFile//Aligned/}
    #echo ${realNameFile}
    numLinesInFile=$(wc -l < "${file}")
    numRecordsInFile=$((numLinesInFile / 4))
    echo ${numRecordsInFile}
    sampleNumber=$((percentageReads * (numRecordsInFile / 100)))
    echo ${sampleNumber}
    #seqtk sample -s100 ${file} ${sampleNumber} > ${processeddatafolder}/${nameFile}_${percentageReads}percent_Sampled.fastq
    #salmon quant -i ${datafolder}/GRCh38_latest_rna_rebuilt_index -l A -1 ${fastafolder}/${realNameFile}L001_R1_001.fastq -2 ${fastafolder}/${realNameFile}L001_R2_001.fastq -p 12 -g ${datafolder}/GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile.tsv -o ${realNameFile}transcripts_quant_genenames
    #sort -k4,4nr ${realNameFile}transcripts_quant_genenames/quant.genes.sf > ${realNameFile}transcripts_quant_genenames/quant.genes.sf.sorted
    #sort -k4,4nr ${realNameFile}transcripts_quant_genenames/quant.sf > ${realNameFile}transcripts_quant_genenames/quant.sf.sorted
done
#cd ${processeddatafolder}