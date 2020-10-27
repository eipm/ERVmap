#!/bin/bash
# set -x

_DEBUG=!{debug}

logMsg() {
    if [[ $# -lt 2 ]];then
        printf "\nUsage: logMsg <error_type> msg\nwhere <error_type> is one of 'INFO', 'WARN', 'DEBUG', 'CMD', or 'ERROR'.\n"
        exit 1
    fi
    case "$1" in
        'CMD' )
            printf "+ [%s]:\t[CMD]\t%s\n" "$(date)" "$2"
            ;;
        'DEBUG' | 'true' )
            if [ $_DEBUG == 'on' ];then printf "[%s]:\t[%s]\t%s\n" "$(date)" "$1" "$2";fi
            ;;
        'WARN' | 'INFO' | 'ERROR' ) 
            printf "[%s]:\t[%s]\t%s\n" "$(date)" "$1" "$2"
            ;;
    esac
    if [[ "$1" == ERROR ]];then
        exit 1
    fi
} 

# initializing parameters for STAR
READS="!{reads}"
CPUS=!{task.cpus}
LIMIT_RAM=!{limitMemory}
OUT_PREFIX="!{sample}."
LOCAL_OUTDIR="!{localOutputDir}"

# checking the prefix of the output BAM
if [ -z ${OUT_PREFIX+x} ];then
    OUT_PREFIX="$RANDOM""_"
    logMsg "WARN" "This prefix will be used as output: $OUT_PREFIX"
fi 
# setting default parameters if not defined
if [ -z ${CPUS+x} ];then export CPUS=1;fi
if [ -z ${LIMIT_RAM+x} ];then export LIMIT_RAM=35129075129;fi
[ -e "/genome/genomeParameters.txt" ] || logMsg "ERROR" "The indexed genome cannot be found. Check that it is present and you have read permissions."

logMsg "DEBUG" "OUT_PREFIX:($OUT_PREFIX)"
logMsg "DEBUG" "Local OutDir: $(pwd)/$LOCAL_OUTDIR"
logMsg "DEBUG" "Reads: ($READS)"
logMsg "DEBUG" "CPUs:($CPUS)"
logMsg "DEBUG" "Limit RAM:($LIMIT_RAM)"

logMsg "INFO" "-------- START ERValign ---------"

BAM="$LOCAL_OUTDIR/$OUT_PREFIX""Aligned.sortedByCoord.out.bam"

logMsg "INFO" "---- Alignment ----"
STAR --genomeDir /genome --runThreadN $CPUS --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM $LIMIT_RAM --outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.02 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesIn $READS --readFilesCommand zcat --outFileNamePrefix $LOCAL_OUTDIR/$OUT_PREFIX --outTmpDir /STAR_tmp/$OUT_PREFIX""tmp
[ $? == 0 ] || logMsg  "ERROR" "The alignment didn't complete succesfully. Check the logs."

logMsg "INFO" "---- Alignment Complete ----"
[ -e  "$BAM" ] || logMsg "ERROR" "BAM files not available:($BAM)"

logMsg "INFO" "---- Indexing"
samtools index -@ $CPUS "$BAM"
[ $? == 0 ] || logMsg  "ERROR" "The indexing didn't complete succesfully. Check the logs."
logMsg "INFO" "---- Indexing Complete"

logMsg "INFO" "-------- END ERValign ---------"
