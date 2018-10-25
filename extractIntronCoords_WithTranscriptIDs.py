import pandas as pd

allintrons_DF = pd.read_csv("../data/ExtractedIntronCoordinatesFromGFFfile.bed",sep="\t",header=None)
#print allintrons_DF.shape
#allintrons_DF.head()


allintrons_DF_groupedByTranscript = allintrons_DF.groupby(by=[5])

allIntronicCoords_ToWrite = pd.DataFrame(columns=["chrm","start","end","transcriptID","strand"])
new_geneNameDF = pd.DataFrame(columns=["TranscriptID","Gene"])
for key, item in allintrons_DF_groupedByTranscript:
    introns = allintrons_DF_groupedByTranscript.get_group(key)
    introns_sorted = introns.sort_values(by=[0,1])
    num_rows = introns.shape[0]
    transcriptIDs = [key+"_"+str(i+1) for i in range(num_rows)]
    new_DF = pd.DataFrame({"chrm":introns_sorted[0],"start":introns_sorted[1],"end":introns_sorted[2],"transcriptID":transcriptIDs,"strand":introns_sorted[3]},columns=["chrm","start","end","transcriptID","strand"])
    allIntronicCoords_ToWrite = pd.concat([allIntronicCoords_ToWrite,new_DF],ignore_index=True,axis=0)
    allIntronicCoords_ToWrite.reset_index(drop=True)
    new_geneNameDF = pd.concat([new_geneNameDF,pd.DataFrame({"Gene":introns_sorted[4],"TranscriptID":transcriptIDs},columns=["TranscriptID","Gene"])],ignore_index=True,axis=0)
    new_geneNameDF.reset_index(drop=True)
    
new_geneNameDF.to_csv("../data/JustIntrons_GeneNameTranscriptID_NCBI_RefSeq_hg38_FromGFFfile.tsv",sep="\t",header=False,index=False)

allIntronicCoords_ToWrite["chrm"] = allIntronicCoords_ToWrite["chrm"].replace("chr23","chrX")
allIntronicCoords_ToWrite["chrm"] = allIntronicCoords_ToWrite["chrm"].replace("chr24","chrY")


allIntronicCoords_ToWrite_Sorted = allIntronicCoords_ToWrite.sort_values(by=["chrm","start"])

allIntronicCoords_ToWrite_Sorted.to_csv("../data/ExtractedIntronCoordinatesFromGFFfile_ForSequence.bed",sep="\t",header=False,index=False)


