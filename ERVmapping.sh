#!/bin/bash

logMsg() {
    if [[ $# -lt 2 ]];then
        printf "\nUsage: logMsg <error_type> msg\nwhere <error_type> is one of 'INFO', 'WARN', 'DEBUG', 'CMD', or 'ERROR'.\n"
        exit 1
    fi
    if [[ "$1" == CMD ]];then 
        printf "+ [%s]:\t[CMD]\t%s\n" "$(date)" "$2"
    else
        printf "[%s]:\t[%s]\t%s\n" "$(date)" "$1" "$2"
    fi
    if [[ "$1" == ERROR ]];then
        exit 1
    fi
} 

if [ -z ${CPUS+x} ];then export CPUS=1;fi
if [ -z ${LIMIT_RAM+x} ];then export LIMIT_RAM=35129075129;fi
if [ -z ${READ_1+x} ];then 
    logMsg "ERROR" "Read 1 not set"
fi
if [ -z ${READ_2+x} ];then 
    logMsg "ERROR" "Read 2 not set"
fi
if [[ ! -e "$READ_1" ]];then
    logMsg "ERROR" "Read_1 file does not exist: ($READ_1)"
fi
if [[ ! -e "$READ_2" ]];then
    logMsg "ERROR" "Read_2 file does not exist: ($READ_2)"
fi

logMsg "INFO" "---- Alignment"
STAR --genomeDir /resources --runThreadN $CPUS --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM $LIMIT_RAM --outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.02 --alignIntronMin 20 --alignIntronMax 1000000 --alignMateGapMax 1000000 --readFilesIn $READ_1 $READ_2 --readFilesCommand zcat
logMsg "INFO" "---- Alignment Complete"

logMsg "INFO" "---- Indexing"
samtools index Aligned.sortedByCoord.out.bam
logMsg "INFO" "---- Indexing Complete"

logMsg "INFO" "---- Finding ERVs"
coverageBed -nonamecheck -a /scripts/ERVmap.bed -b Aligned.sortedByCoord.out.bam -counts > ERVresults.txt
logMsg "INFO" "---- Finding ERVs complete"
#~/project/genome/hg38cut_L1_ERV.bed -counts > $i.tx
