# coding: utf-8

# # Extract introns from GFF file fata
# The purpose of this script is to extract intron coordinates from the GFF file data
# We will look at features labeled exons and grab the end coordinates and start coordinates

# In[1]:


import pandas as pd


# In[5]:


data_genes = pd.read_csv("../data/GRCh38_latest_genomic_ValidChroms_OnlyGenes.bed",sep="\t",header=None)
print data_genes.shape
print data_genes.head()


# In[7]:


# Set column names
data_genes.columns = ["chrm","start","end","strand","gene","transcriptID","RefSeqType","feature"]


# In[8]:


# Only grab exons 
data_exons = data_genes[data_genes["feature"]=="exon"]
print data_exons.shape
print data_exons.head()


# In[10]:


# Let's group by the transcript ID
data_exons_GrpedByFeature = data_exons.groupby(by=["transcriptID"])


# In[40]:

# Let's extract the entire genomic coordinate
allGenomicCoords_DF = pd.DataFrame(columns=["chrm","start","end","strand","gene","transcriptID"])
for key, item in data_exons_GrpedByFeature:
    exons = data_exons_GrpedByFeature.get_group(key)
    chrm = [list(exons["chrm"].values)[0]]
    gene = [list(exons["gene"].values)[0]]
    strand = [list(exons["strand"].values)[0]]
    transcriptID = [key]
    clump_start_end = sorted(list(exons["start"].values) + list(exons["end"].values))
    start = [clump_start_end[0]]
    end = [clump_start_end[len(clump_start_end)-1]]
    # Only do this if there are more than 1 exon
    newGenomicCoord_DF = pd.DataFrame({"chrm":chrm,"start":start,"end":end,"strand":strand,"transcriptID":transcriptID,"gene":gene},columns=["chrm","start","end","strand","gene","transcriptID"])
    allGenomicCoords_DF = pd.concat([allGenomicCoords_DF,newGenomicCoord_DF],ignore_index=True,axis=0)
    allGenomicCoords_DF = allGenomicCoords_DF.reset_index(drop=True)


allGenomicCoords_DF.to_csv("../data/ExtractedGenomicCoordinatesFromGFFfile.bed",sep="\t",header=False,index=False)