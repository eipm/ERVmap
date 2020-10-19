#!/bin/bash

_DEBUG="off"

function DEBUG() {
 [ "$_DEBUG" == "on" ] &&  $@
}

function usage() {
    echo "Usage: ERVmapping.sh <-r1|--read1> SAMPLE_1.fastq.gz <-r2|--read2> SAMPLE_1.fastq.gz [-o|--output] results/SAMPLE <-m|--mode> {STAR|BED|ALL} [-c|--cpus] Ncpus [-l|--limit-ram] 35129075129 [-d|--debug]"
}


MODE="!{mode}"
OUT_PREFIX="!{outPrefix}"
CPUS=!{cpus}
LIMIT_RAM=!{limitMemory}
_DEBUG="!{debug}"

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
if [ -z ${MODE+x} ];then 
    usage
    logMsg "ERROR" "MODE not set."
fi

if [ -z ${OUT_PREFIX+x} ];then
    OUT_PREFIX="$RANDOM""_"
    logMsg "WARN" "This prefix will be used as output: $OUT_PREFIX"
fi 

logMsg "DEBUG" "CPUs:($CPUS)"
logMsg "DEBUG" "Limit RAM:($LIMIT_RAM)"
logMsg "DEBUG" "OUT_PREFIX:($OUT_PREFIX)"
logMsg "DEBUG" "MODE:($MODE)"

logMsg "INFO" "-------- START (mode: $MODE) ---------"

if [[ "$MODE" == "STAR" || "$MODE" == "ALL" ]]; then
    READS="!{reads}"
    logMsg "DEBUG" "Reads: ($READS)"
    BAM="./results/$OUT_PREFIX""Aligned.sortedByCoord.out.bam"
    if [[ ! -e "$BAM" ]];then
    logMsg "INFO" "---- Alignment ----"
    STAR --genomeDir /genome --runThreadN $CPUS --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM $LIMIT_RAM --outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.02 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesIn $READS --readFilesCommand zcat --outFileNamePrefix ./results/$OUT_PREFIX
    logMsg "INFO" "---- Alignment Complete ----"
    if [[ ! -e  "$BAM" ]];then
        logMsg "DEBUG" "PWD: $(pwd);$(ls -la ./)"
        logMsg "DEBUG" "RESULTS: $(ls -l results)"
        logMsg "ERROR" "BAM files not available:($BAM)"
    fi
    logMsg "INFO" "---- Indexing"
    samtools index -@ $CPUS "$BAM"
    logMsg "INFO" "---- Indexing Complete"
    else
        # BAM file already exists, skipping
        logMsg "WARN" "BAM file already exists; skipping this step"
    fi
fi

if [[ "$MODE" == "BED" || "$MODE" == "ALL" ]]; then
    if [[ -z $BAM ]] BAM="!{bam}"
    if [[ ! -e "$BAM" ]];then
        logMsg "ERROR" "Cannot find BAM file: ( $BAM )"
    fi
    logMsg "INFO" "---- Finding ERVs ----"
    coverageBed -nonamecheck -a /resources/ERVmap.bed -b "data/$OUT_PREFIX""Aligned.sortedByCoord.out.bam" -counts -sorted > ./results/"$OUT_PREFIX""ERVresults.txt"
    logMsg "INFO" "---- Finding ERVs complete ----"
fi

logMsg "INFO" "-------- END (mode: $MODE) ---------"
