#!/bin/bash

# This script is going to get featureCounts using STAR alignment

folderToCD=$1 

if [ "$2" = "Transcriptome" ]
then
    annotationFile="SAFformat_GRCh38_latest_genomic_ValidChroms_OnlyGenes.bed"
elif [ "$2" = "Genome" ]
then
    annotationFile="SAFformat_ExtractedGenomicCoordinatesFromGFFfile_OnlyGenes.bed"
elif [ "$2" = "ExonIntron5p" ]
then
    annotationFile="SAFformat_ExtractedExonIntronCoordinatesFromGFFfile_5p_OnlyGenes.bed"
elif [ "$2" = "ExonIntron3p" ]
then
    annotationFile="SAFformat_ExtractedExonIntronCoordinatesFromGFFfile_3p_OnlyGenes.bed"
fi

cd ${folderToCD}
pwd
for file in *_Aligned.out.sortedByCoord.bam
do
    echo ${file}
    nameFile=$(echo "$file" | cut -d"." -f1)
    realNameFile=${nameFile//Aligned/}
    echo ${realNameFile}
    #samtools sort -@ 4 -n -o ${realNameFile}Aligned.out.sortedByName.bam ${realNameFile}Aligned.out.bam
    #featureCounts -p --largestOverlap -B -C -M -R SAM -t exon -g gene_id -a ../../data/SAFformat_GRCh38_latest_genomic_ValidChroms_OnlyGenes.bed -F SAF -o ${realNameFile}featureCounts.txt ${realNameFile}Aligned.out.sortedByName.bam
    featureCounts -p --largestOverlap -B -C -M -R SAM -t exon -g gene_id -a ../../data/${annotationFile} -F SAF -o ${realNameFile}${2}_featureCounts.txt ${realNameFile}Aligned.out.sortedByName.bam
    cut -f1,6,7 ${realNameFile}${2}_featureCounts.txt > ${realNameFile}${2}_featureCounts_JustCounts.txt
    sort -k3,3nr ${realNameFile}${2}_featureCounts_JustCounts.txt > ${realNameFile}${2}_featureCounts_JustCounts.sorted.txt
    #salmon quant -i ${datafolder}/GRCh38_latest_rna_rebuilt_index -l A -1 ${fastafolder}/${realNameFile}L001_R1_001.fastq -2 ${fastafolder}/${realNameFile}L001_R2_001.fastq -p 12 -g ${datafolder}/GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile.tsv -o ${realNameFile}transcripts_quant_genenames
    #salmon quant -i ${datafolder}/GRCh38_latest_rna_rebuilt_index -l A -1 ${processeddatafolder}/${realNameFile}L001_R1_001_50percent_Sampled.fastq -2 ${processeddatafolder}/${realNameFile}L001_R2_001_50percent_Sampled.fastq -p 12 -g ${datafolder}/GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile.tsv -o ${realNameFile}transcripts_quant_genenames_Just50percentReads
    #salmon quant -i ${datafolder}/ExtractedGenomicCoordinatesFromGFFfile_sequences_index -l A -1 ${fastafolder}/${realNameFile}L001_R1_001.fastq -2 ${fastafolder}/${realNameFile}L001_R2_001.fastq -p 12 -g ${datafolder}/GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile.tsv -o ${realNameFile}transcripts_quant_genenames_AgainstGenomicCoordinates
    #salmon quant -i ${datafolder}/ExtractedIntronCoordinatesFromGFFfile_sequences_index -l A -1 ${fastafolder}/${realNameFile}L001_R1_001.fastq -2 ${fastafolder}/${realNameFile}L001_R2_001.fastq -p 12 -g ${datafolder}/JustIntrons_GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile.tsv -o ${realNameFile}transcripts_quant_genenames_AgainstIntronicCoordinates
    #salmon quant -i ${datafolder}/Unique_ExtractedExonIntronCoordinates_5p_FromGFFfile_sequences_index -l A -1 ${fastafolder}/${realNameFile}L001_R1_001.fastq -2 ${fastafolder}/${realNameFile}L001_R2_001.fastq -p 12 -g ${datafolder}/ExonIntronCoords_GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile.tsv -o ${realNameFile}transcripts_quant_genenames_AgainstUnique5pExonIntronBoundaryCoordinates
    #salmon quant -i ${datafolder}/Unique_ExtractedExonIntronCoordinates_3p_FromGFFfile_sequences_index -l A -1 ${fastafolder}/${realNameFile}L001_R1_001.fastq -2 ${fastafolder}/${realNameFile}L001_R2_001.fastq -p 12 -g ${datafolder}/ExonIntronCoords_GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile.tsv -o ${realNameFile}transcripts_quant_genenames_AgainstUnique3pExonIntronBoundaryCoordinates
    #sort -k4,4nr ${realNameFile}transcripts_quant_genenames/quant.genes.sf > ${realNameFile}transcripts_quant_genenames/quant.genes.sf.sorted
    #sort -k4,4nr ${realNameFile}transcripts_quant_genenames/quant.sf > ${realNameFile}transcripts_quant_genenames/quant.sf.sorted
    #sort -k4,4nr ${realNameFile}transcripts_quant_genenames_Just50percentReads/quant.genes.sf > ${realNameFile}transcripts_quant_genenames_Just50percentReads/quant.genes.sf.sorted
    #sort -k4,4nr ${realNameFile}transcripts_quant_genenames_Just50percentReads/quant.sf > ${realNameFile}transcripts_quant_genenames_Just50percentReads/quant.sf.sorted
    #sort -k4,4nr ${realNameFile}transcripts_quant_genenames_AgainstGenomicCoordinates/quant.genes.sf > ${realNameFile}transcripts_quant_genenames_AgainstGenomicCoordinates/quant.genes.sf.sorted
    #sort -k4,4nr ${realNameFile}transcripts_quant_genenames_AgainstGenomicCoordinates/quant.sf > ${realNameFile}transcripts_quant_genenames_AgainstGenomicCoordinates/quant.sf.sorted
    #sort -k4,4nr ${realNameFile}transcripts_quant_genenames_AgainstIntronicCoordinates/quant.genes.sf > ${realNameFile}transcripts_quant_genenames_AgainstIntronicCoordinates/quant.genes.sf.sorted
    #sort -k4,4nr ${realNameFile}transcripts_quant_genenames_AgainstIntronicCoordinates/quant.sf > ${realNameFile}transcripts_quant_genenames_AgainstIntronicCoordinates/quant.sf.sorted
done
cd ../../scripts

# featureCounts -p --largestOverlap -B -C -M -R SAM -t exon -g gene_id -a ../../data/SAFformat_ExtractedExonIntronCoordinatesFromGFFfile_3p_OnlyTranscriptIDs.bed -F SAF -o AllMinus_Madrasin_ExonIntron_3p_featureCounts.txt AllMinus_Madrasin_MAPQgreaterThan10_SortedByCoords.bam

# featureCounts -p --largestOverlap -B -C -M -R SAM -t exon -g gene_id -a ../../data/SAFformat_ExtractedExonExonCoordinatesFromGFFfile_3p_OnlyTranscriptIDs.bed -F SAF -o AllMinus_Madrasin_ExonExon_3p_featureCounts.txt AllMinus_Madrasin_MAPQgreaterThan10_SortedByCoords.bam

# featureCounts -p --largestOverlap -B -C -M -R SAM -t exon -g gene_id -a ../../data/SAFformat_ExtractedExonExonCoordinatesFromGFFfile_5p_OnlyTranscriptIDs.bed -F SAF -o AllMinus_Madrasin_ExonExon_5p_featureCounts.txt AllMinus_Madrasin_MAPQgreaterThan10_SortedByCoords.bam

# featureCounts -p --largestOverlap -B -C -M -R SAM -t exon -g gene_id -a ../../data/SAFformat_ExtractedExonIntronCoordinatesFromGFFfile_5p_OnlyTranscriptIDs.bed -F SAF -o AllMinus_Madrasin_ExonIntron_5p_featureCounts.txt AllMinus_Madrasin_MAPQgreaterThan10_SortedByCoords.bam > featureCounts_AllMinus_ExonIntron5p.txt