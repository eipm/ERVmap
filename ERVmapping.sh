#!/bin/bash
_DEBUG="off"
echo $1

function DEBUG() {
 [ "$_DEBUG" == "on" ] &&  $@
}

function usage() {
    echo "Usage: ERVmapping.sh <-r1|--read1> SAMPLE_1.fastq.gz <-r2|--read2> SAMPLE_1.fastq.gz [-o|--output] results/SAMPLE <-m|--mode> {STAR|BED|ALL} [-c|--cpus] Ncpus [-l|--limit-ram] 35129075129 [-d|--debug]"
}

while [ "$1" != "" ]; do
    case $1 in
        -r1 | --read1 )         shift
                                READ1=$1
                                ;;
        -r2 | --read2 )         shift
                                READ2=$1
                                ;;
        -o | --output )         shift
                                OUT_PREFIX=$1
                                ;;
        -m | --mode )           shift
                                MODE=$1
                                ;;
        -c | --cpus )           shift
                                CPUS=$1
                                ;;
        -l | --limit-ram )      shift
                                LIMIT_RAM=$1
                                ;;
        -d | --debug  )         _DEBUG="on"
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

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
if [ -z ${READ1+x} ];then 
    logMsg "ERROR" "Read 1 not set"
fi
if [ -z ${READ2+x} ];then 
    logMsg "ERROR" "Read 2 not set"
fi
if [[ ! -e "$READ1" ]];then
    logMsg "ERROR" "Read_1 file does not exist: ($READ1)"
fi
if [[ ! -e "$READ2" ]];then
    logMsg "ERROR" "Read_2 file does not exist: ($READ2)"
fi
if [ -z ${OUT_PREFIX+x} ];then
    OUT_PREFIX="./RESULTS/Alignment_$RANDOM"
    logMsg "WARN" "This prefix will be used as output: $OUT_PREFIX"
fi 
BAM="/results/$OUT_PREFIX""Aligned.sortedByCoord.out.bam"

logMsg "DEBUG" "CPUs:($CPUS)"
logMsg "DEBUG" "Limit RAM:($LIMIT_RAM)"
logMsg "DEBUG" "Read1:($READ1)"
logMsg "DEBUG" "Read2:($READ2)"
logMsg "DEBUG" "OUT_PREFIX:($OUT_PREFIX)"
logMsg "DEBUG" "MODE:($MODE)"

logMsg "INFO" "-------- START (mode: $MODE) ---------"

if [[ "$MODE" == "STAR" || "$MODE" == "ALL" ]]; then
    logMsg "INFO" "---- Alignment ----"
    logMsg "DEBUG" "$(ls -la /results)"
    STAR --genomeDir /resources --runThreadN $CPUS --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM $LIMIT_RAM --outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.02 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesIn $READ1 $READ_2 --readFilesCommand zcat --outFileNamePrefix /results/$OUT_PREFIX
    logMsg "INFO" "---- Alignment Complete ----"
    if [[ ! -e  "$BAM" ]];then
        logMsg "ERROR" "BAM files not available:($BAM)"
    fi
    logMsg "INFO" "---- Indexing"
    samtools index "$BAM"
    logMsg "INFO" "---- Indexing Complete"
fi

if [[ "$MODE" == "BED" || "$MODE" == "ALL" ]]; then
    logMsg "INFO" "---- Finding ERVs ----"
    coverageBed -nonamecheck -a /resources/ERVmap.bed -b "$BAM" -counts -sorted > $OUT_PREFIX.ERVresults.txt
    logMsg "INFO" "---- Finding ERVs complete ----"
fi

logMsg "INFO" "-------- END (mode: $MODE) ---------"
