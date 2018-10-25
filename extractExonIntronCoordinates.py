# coding: utf-8

# # Extract introns from GFF file fata
# The purpose of this script is to extract exon-intron boundary coordinates  on the 5' and 3' sitesfrom the GFF file data


import pandas as pd


data_genes = pd.read_csv("../data/GRCh38_latest_genomic_ValidChroms_OnlyGenes.bed",sep="\t",header=None)
print data_genes.shape
print data_genes.head()


# Set column names
data_genes.columns = ["chrm","start","end","strand","gene","transcriptID","RefSeqType","feature"]


# Only grab exons 
data_exons = data_genes[data_genes["feature"]=="exon"]
print data_exons.shape
print data_exons.head()

# Split up exons by plus and minus strand
data_exons_plus = data_exons[data_exons["strand"]=="+"]
print data_exons_plus.head()
data_exons_minus = data_exons[data_exons["strand"]=="-"]
print data_exons_minus.head()


data_exons_plus = data_exons_plus.reset_index(drop=True)
data_exons_minus = data_exons_minus.reset_index(drop=True)

# Get 5' exon-intron boundaries
data_5p = pd.DataFrame({"chrm":pd.concat([data_exons_plus["chrm"],data_exons_minus["chrm"]]),"start":pd.concat([data_exons_plus["end"]-100,data_exons_minus["start"]-100]),
                       "end":pd.concat([data_exons_plus["end"]+100,data_exons_minus["start"]+100]), "strand":pd.concat([data_exons_plus["strand"],data_exons_minus["strand"]]),
                       "transcriptID":pd.concat([data_exons_plus["transcriptID"],data_exons_minus["transcriptID"]]), "gene":pd.concat([data_exons_plus["gene"],data_exons_minus["gene"]])},
                       columns=["chrm","start","end","transcriptID","strand","gene"])
print data_5p.head()
print data_5p.tail()


# Get 3' exon-intron boundaries
data_3p = pd.DataFrame({"chrm":pd.concat([data_exons_plus["chrm"],data_exons_minus["chrm"]]),"start":pd.concat([data_exons_plus["start"]-100,data_exons_minus["end"]-100]),
                       "end":pd.concat([data_exons_plus["start"]+100,data_exons_minus["end"]+100]), "strand":pd.concat([data_exons_plus["strand"],data_exons_minus["strand"]]),
                       "transcriptID":pd.concat([data_exons_plus["transcriptID"],data_exons_minus["transcriptID"]]), "gene":pd.concat([data_exons_plus["gene"],data_exons_minus["gene"]])},
                       columns=["chrm","start","end","transcriptID","strand","gene"])
print data_3p.head()
print data_3p.tail()


# 5p exon-intron boundaries grouped
data_5p_groupedByTranscript = data_5p.groupby(by=["transcriptID"])
# 3p exon-intron boundaries grouped
data_3p_groupedByTranscript = data_3p.groupby(by=["transcriptID"])
# All exons grouped
data_exons_groupedByTranscript = data_exons.groupby(by=["transcriptID"])

all_5p_ExonIntronCoords_ToWrite = pd.DataFrame(columns=["chrm","start","end","transcriptID","strand"])
all_3p_ExonIntronCoords_ToWrite = pd.DataFrame(columns=["chrm","start","end","transcriptID","strand"])
new_geneNameDF = pd.DataFrame(columns=["TranscriptID","Gene"])

for key, item in data_exons_groupedByTranscript:
    
    coords_key = data_exons_groupedByTranscript.get_group(key)
    num_rows = coords_key.shape[0]
    transcriptIDs = [key+"_"+str(i+1) for i in range(num_rows-1)]
    
    exonintron_coords_5p = data_5p_groupedByTranscript.get_group(key)
    exonintron_coords_5p_sorted = exonintron_coords_5p.sort_values(by=["chrm","start"])
    # We are getting rid of the first row for 5p, because that is not actually including the intron, it's the genomic region
    #print list(exonintron_coords_5p_sorted["strand"].values)
    if list(exonintron_coords_5p_sorted["strand"].values)[0]=="+":
        exonintron_coords_5p_sorted = exonintron_coords_5p_sorted.iloc[:-1]
    else:
        exonintron_coords_5p_sorted = exonintron_coords_5p_sorted.iloc[1:]
    new_DF_5p = pd.DataFrame({"chrm":exonintron_coords_5p_sorted["chrm"],"start":exonintron_coords_5p_sorted["start"],"end":exonintron_coords_5p_sorted["end"],"transcriptID":transcriptIDs,"strand":exonintron_coords_5p_sorted["strand"]}
                          ,columns=["chrm","start","end","transcriptID","strand"])
    all_5p_ExonIntronCoords_ToWrite = pd.concat([all_5p_ExonIntronCoords_ToWrite,new_DF_5p],ignore_index=True,axis=0)
    all_5p_ExonIntronCoords_ToWrite.reset_index(drop=True)
    
    exonintron_coords_3p = data_3p_groupedByTranscript.get_group(key)
    exonintron_coords_3p_sorted = exonintron_coords_3p.sort_values(by=["chrm","start"])
    # We are getting rid of the last row for 3p, because that is not actually including the intron, it's the genomic region
    if list(exonintron_coords_3p_sorted["strand"].values)[0]=="+":
        exonintron_coords_3p_sorted = exonintron_coords_3p_sorted.iloc[1:]
    else:
        exonintron_coords_3p_sorted = exonintron_coords_3p_sorted.iloc[:-1]
    new_DF_3p = pd.DataFrame({"chrm":exonintron_coords_3p_sorted["chrm"],"start":exonintron_coords_3p_sorted["start"],"end":exonintron_coords_3p_sorted["end"],"transcriptID":transcriptIDs,"strand":exonintron_coords_3p_sorted["strand"]}
                          ,columns=["chrm","start","end","transcriptID","strand"])
    all_3p_ExonIntronCoords_ToWrite = pd.concat([all_3p_ExonIntronCoords_ToWrite,new_DF_3p],ignore_index=True,axis=0)
    all_3p_ExonIntronCoords_ToWrite.reset_index(drop=True)
    
    new_geneNameDF = pd.concat([new_geneNameDF,pd.DataFrame({"Gene":list(coords_key["gene"].values)[1:],"TranscriptID":transcriptIDs},columns=["TranscriptID","Gene"])],ignore_index=True,axis=0)
    new_geneNameDF.reset_index(drop=True)
    
new_geneNameDF.to_csv("../data/ExonIntronCoords_GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile.tsv",sep="\t",header=False,index=False)

all_5p_ExonIntronCoords_ToWrite["chrm"] = all_5p_ExonIntronCoords_ToWrite["chrm"].replace("chr23","chrX")
all_5p_ExonIntronCoords_ToWrite["chrm"] = all_5p_ExonIntronCoords_ToWrite["chrm"].replace("chr24","chrY")


all_5p_ExonIntronCoords_ToWrite_Sorted = all_5p_ExonIntronCoords_ToWrite.sort_values(by=["chrm","start"])
all_5p_ExonIntronCoords_ToWrite_Sorted.to_csv("../data/ExtractedExonIntronCoordinates_5p_FromGFFfile_ForSequence.bed",sep="\t",header=False,index=False)

all_3p_ExonIntronCoords_ToWrite["chrm"] = all_3p_ExonIntronCoords_ToWrite["chrm"].replace("chr23","chrX")
all_3p_ExonIntronCoords_ToWrite["chrm"] = all_3p_ExonIntronCoords_ToWrite["chrm"].replace("chr24","chrY")


all_3p_ExonIntronCoords_ToWrite_Sorted = all_3p_ExonIntronCoords_ToWrite.sort_values(by=["chrm","start"])
all_3p_ExonIntronCoords_ToWrite_Sorted.to_csv("../data/ExtractedExonIntronCoordinates_3p_FromGFFfile_ForSequence.bed",sep="\t",header=False,index=False)