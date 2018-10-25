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


allintrons_DF = pd.DataFrame(columns=["chrm","start","end","strand","gene","transcriptID"])
for key, item in data_exons_GrpedByFeature:
    exons = data_exons_GrpedByFeature.get_group(key)
    # Only do this if there are more than 1 exon
    if exons.shape[0] > 1:
        clump_start_end = sorted(list(exons["start"].values) + list(exons["end"].values))
        # Get start coordinates of introns 
        start_coords_intron = [clump_start_end[i] for i in range(1,len(clump_start_end)-2,2)]
        # Get end coordinates of introns
        end_coords_intron = [clump_start_end[i]-1 for i in range(2,len(clump_start_end)-1,2)]
        if len(start_coords_intron)==len(end_coords_intron):
            chrm = [list(exons["chrm"].values)[0]]*len(start_coords_intron)
            gene = [list(exons["gene"].values)[0]]*len(start_coords_intron)
            strand = [list(exons["strand"].values)[0]]*len(start_coords_intron)
            transcriptID = [key]*len(start_coords_intron)
            newintron_DF = pd.DataFrame({"chrm":chrm,"start":start_coords_intron,"end":end_coords_intron,"strand":strand,"transcriptID":transcriptID,"gene":gene},columns=["chrm","start","end","strand","gene","transcriptID"])
            allintrons_DF = pd.concat([allintrons_DF,newintron_DF],ignore_index=True,axis=0)
            allintrons_DF = allintrons_DF.reset_index(drop=True)
        else:
            print key


# In[41]:


#allintrons_DF.head()


allintrons_DF.to_csv("../data/ExtractedIntronCoordinatesFromGFFfile.bed",sep="\t",header=False,index=False)