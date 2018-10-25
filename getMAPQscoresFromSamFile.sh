#!/bin/bash

samFileToExtract=$1
mapqFileToWrite=$2

#while read line;
#do
#    if [ ${line:0:1} != '@' ]
#    then
        #echo ${line}
#        A="$(cut -f5 <<<"$line")"
#        echo ${A} >> ${mapqFileToWrite}
        #cut -f5 $line >> ${mapqFileToWrite}
#    fi
#done < ${samFileToExtract}


samtools view ${samFileToExtract} | cut -f5 > ${mapqFileToWrite}