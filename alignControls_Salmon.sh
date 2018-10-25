#!/bin/bash

# salmon index -t GRCh38_latest_rna.fna -i GRCh38_latest_rna_rebuilt_index --type quasi -k 31

# salmon index -t ExtractedGenomicCoordinatesFromGFFfile_sequences.fa -i ExtractedGenomicCoordinatesFromGFFfile_sequences_index --type quasi -k 31

# salmon index -t ExtractedIntronCoordinatesFromGFFfile_sequences.fa -i ExtractedIntronCoordinatesFromGFFfile_sequences_index --type quasi -k 31


# This script is going to get align to transcriptome and get the transcript count using Salmon

folderToCD=$1 

datafolder="/home/jkumar12/Projects/Analyzing_SHAPEdata_In_Introns/data"

fastafolder="/home/shared/SequencingData/MiSeq_LL_06142018_MadrasinTrial/FASTQ_Generation_2018-06-15_00_25_47Z-103862779"

#processeddatafolder="/home/jkumar12/Projects/Analyzing_SHAPEdata_In_Introns/processed_data"

cd ${folderToCD}
pwd
#for file in *_Aligned.out.sortedByCoord.bam
#do
    #echo ${file}
    #nameFile=$(echo "$file" | cut -d"." -f1)
    #realNameFile=${nameFile//Aligned/}
    #echo ${realNameFile}
    #salmon quant -i ${datafolder}/GRCh38_latest_rna_rebuilt_index -l A -1 ${fastafolder}/${realNameFile}L001_R1_001.fastq -2 ${fastafolder}/${realNameFile}L001_R2_001.fastq -p 12 -g ${datafolder}/GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile.tsv -o ${realNameFile}transcripts_quant_genenames
    #salmon quant -i ${datafolder}/GRCh38_latest_rna_rebuilt_index -l A -1 ${processeddatafolder}/${realNameFile}L001_R1_001_50percent_Sampled.fastq -2 ${processeddatafolder}/${realNameFile}L001_R2_001_50percent_Sampled.fastq -p 12 -g ${datafolder}/GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile.tsv -o ${realNameFile}transcripts_quant_genenames_Just50percentReads
    #salmon quant -i ${datafolder}/ExtractedGenomicCoordinatesFromGFFfile_sequences_index -l A -1 ${fastafolder}/${realNameFile}L001_R1_001.fastq -2 ${fastafolder}/${realNameFile}L001_R2_001.fastq -p 12 -g ${datafolder}/GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile.tsv -o ${realNameFile}transcripts_quant_genenames_AgainstGenomicCoordinates
    #salmon quant -i ${datafolder}/ExtractedIntronCoordinatesFromGFFfile_sequences_index -l A -1 ${fastafolder}/${realNameFile}L001_R1_001.fastq -2 ${fastafolder}/${realNameFile}L001_R2_001.fastq -p 12 -g ${datafolder}/JustIntrons_GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile.tsv -o ${realNameFile}transcripts_quant_genenames_AgainstIntronicCoordinates
    #sort -k4,4nr ${realNameFile}transcripts_quant_genenames/quant.genes.sf > ${realNameFile}transcripts_quant_genenames/quant.genes.sf.sorted
    #sort -k4,4nr ${realNameFile}transcripts_quant_genenames/quant.sf > ${realNameFile}transcripts_quant_genenames/quant.sf.sorted
    #sort -k4,4nr ${realNameFile}transcripts_quant_genenames_Just50percentReads/quant.genes.sf > ${realNameFile}transcripts_quant_genenames_Just50percentReads/quant.genes.sf.sorted
    #sort -k4,4nr ${realNameFile}transcripts_quant_genenames_Just50percentReads/quant.sf > ${realNameFile}transcripts_quant_genenames_Just50percentReads/quant.sf.sorted
    #sort -k4,4nr ${realNameFile}transcripts_quant_genenames_AgainstGenomicCoordinates/quant.genes.sf > ${realNameFile}transcripts_quant_genenames_AgainstGenomicCoordinates/quant.genes.sf.sorted
    #sort -k4,4nr ${realNameFile}transcripts_quant_genenames_AgainstGenomicCoordinates/quant.sf > ${realNameFile}transcripts_quant_genenames_AgainstGenomicCoordinates/quant.sf.sorted
    #sort -k4,4nr ${realNameFile}transcripts_quant_genenames_AgainstIntronicCoordinates/quant.genes.sf > ${realNameFile}transcripts_quant_genenames_AgainstIntronicCoordinates/quant.genes.sf.sorted
    #sort -k4,4nr ${realNameFile}transcripts_quant_genenames_AgainstIntronicCoordinates/quant.sf > ${realNameFile}transcripts_quant_genenames_AgainstIntronicCoordinates/quant.sf.sorted
#done
#cd ../../scripts

#salmon quant -i ${datafolder}/GRCh38_latest_rna_rebuilt_index -l A -1 ${fastafolder}/All_DMSfiles_R1.fastq -2 ${fastafolder}/All_DMSfiles_R2.fastq -p 12 -g ${datafolder}/GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile.tsv -o All_DMSfiles_transcripts_quant_genenames

#salmon quant -i ${datafolder}/GRCh38_latest_rna_rebuilt_index -l A -1 ${fastafolder}/All_Minusfiles_R1.fastq -2 ${fastafolder}/All_Minusfiles_R2.fastq -p 12 -g ${datafolder}/GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile.tsv -o All_Minusfiles_transcripts_quant_genenames

#salmon quant -i ${datafolder}/ExtractedGenomicCoordinatesFromGFFfile_sequences_index -l A -1 ${fastafolder}/All_DMSfiles_R1.fastq -2 ${fastafolder}/All_DMSfiles_R2.fastq -p 12 -g ${datafolder}/GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile.tsv -o All_DMSfiles_transcripts_quant_genenames_AgainstGenomicCoordinates

#salmon quant -i ${datafolder}/ExtractedGenomicCoordinatesFromGFFfile_sequences_index -l A -1 ${fastafolder}/All_Minusfiles_R1.fastq -2 ${fastafolder}/All_Minusfiles_R2.fastq -p 12 -g ${datafolder}/GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile.tsv -o All_Minusfiles_transcripts_quant_genenames_AgainstGenomicCoordinates

salmon quant -i ${datafolder}/Unique_ExtractedExonIntronCoordinates_5p_FromGFFfile_sequences_index_Just5p -l A -1 ${fastafolder}/All_DMSfiles_R1.fastq -2 ${fastafolder}/All_DMSfiles_R2.fastq -p 12 -g ${datafolder}/ExonIntronCoords_GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile_5p.tsv -o All_DMSfiles_transcripts_quant_genenames_AgainstUnique5pExonIntronBoundaryCoordinates

salmon quant -i ${datafolder}/Unique_ExtractedExonIntronCoordinates_5p_FromGFFfile_sequences_index_Just5p -l A -1 ${fastafolder}/All_Minusfiles_R1.fastq -2 ${fastafolder}/All_Minusfiles_R2.fastq -p 12 -g ${datafolder}/ExonIntronCoords_GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile_5p.tsv -o All_Minusfiles_transcripts_quant_genenames_AgainstUnique5pExonIntronBoundaryCoordinates

salmon quant -i ${datafolder}/Unique_ExtractedExonIntronCoordinates_3p_FromGFFfile_sequences_index_Just3p -l A -1 ${fastafolder}/All_DMSfiles_R1.fastq -2 ${fastafolder}/All_DMSfiles_R2.fastq -p 12 -g ${datafolder}/ExonIntronCoords_GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile_3p.tsv -o All_DMSfiles_transcripts_quant_genenames_AgainstUnique3pExonIntronBoundaryCoordinates

salmon quant -i ${datafolder}/Unique_ExtractedExonIntronCoordinates_3p_FromGFFfile_sequences_index_Just3p -l A -1 ${fastafolder}/All_Minusfiles_R1.fastq -2 ${fastafolder}/All_Minusfiles_R2.fastq -p 12 -g ${datafolder}/ExonIntronCoords_GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile_3p.tsv -o All_Minusfiles_transcripts_quant_genenames_AgainstUnique3pExonIntronBoundaryCoordinates