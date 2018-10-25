#!/bin/bash

folderToCD=$1 

cd ${folderToCD}
pwd
for bamfile in *_Aligned.out.sortedByCoord.bam
do
    echo ${bamfile}
    nameFile=$(echo "$bamfile" | cut -d"." -f1)
    bamCoverage --bam ${bamfile} --binSize 10 -p max --effectiveGenomeSize 2913022398 --normalizeUsing CPM -o ${nameFile}.out.sortedByCoord.bedgraph -of bedgraph
    bedtools coverage -a ${nameFile}.out.sortedByCoord.bedgraph -b ../../data/NCBI_RefSeq_Curated_hg38_HighlyExpressedGenes_Introns.bed -counts > ${nameFile}.out.Counts.HighlyExpressedGenesIntrons.txt
done
cd ../../scripts