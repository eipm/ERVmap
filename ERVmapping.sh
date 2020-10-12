#!/bin/bash
_DEBUG="off"

function DEBUG()
{
 [ "$_DEBUG" == "on" ] &&  $@
}

function usage()
{
    echo "Usage: ERVmapping.sh <-r1|--read1> SAMPLE_1.fastq.gz <-r2|--read2> SAMPLE_1.fastq.gz <-m|--mode> {STAR|BED} [-c|--cpus] Ncpus [-l|--limit-ram] 35129075129 [-d|--debug]"
}

while [ "$1" != "" ]; do
    case $1 in
        -r1 | --read1 )         shift
                                READ_1=$1
                                ;;
        -r2 | --read2 )         shift
                                READ_2=$1
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
if [ -z ${OUT_PREFIX+x} ];then
    OUT_PREFIX="RESULTS/Alignment_$RANDOM"
    logMsg "WARN" "This prefix will be used as output: $OUT_PREFIX"
fi

logMsg "INFO" "-------- START (mode: $MODE) ---------"

if [[ $MODE == "STAR" ]]; then
    logMsg "INFO" "---- Alignment"
    STAR --genomeDir /resources --runThreadN $CPUS --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM $LIMIT_RAM --outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.02 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesIn $READ_1 $READ_2 --readFilesCommand zcat --outFileNamePrefix $OUT_PREFIX
    logMsg "INFO" "---- Alignment Complete"

    logMsg "INFO" "---- Indexing"
    samtools index Aligned.sortedByCoord.out.bam
    logMsg "INFO" "---- Indexing Complete"
fi

if [[ $MODE == "BED" ]]; then
    logMsg "INFO" "---- Finding ERVs"
    coverageBed -nonamecheck -a /scripts/ERVmap.bed -b Aligned.sortedByCoord.out.bam -counts -sorted > ERVresults.txt
    logMsg "INFO" "---- Finding ERVs complete"
fi

logMsg "INFO" "-------- END (mode: $MODE) ---------"
