#!/bin/bash

function usage() {
    echo "Usage: ERVcount.sh <-b|--bam> Aligned.bam [-o|--output] results/SAMPLE [-d|--debug {off|on}]"
}

BAM="!{bam}"
_DEBUG="!{debug}"

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
        'INFO' | 'WARN' | 'ERROR' ) 
            printf "[%s]:\t[%s]\t%s\n" "$(date)" "$1" "$2"
            ;;
    esac
    if [[ "$1" == ERROR ]];then
        exit 1
    fi
}  

OUT_PREFIX=!{outPrefix}
if [ -z ${OUT_PREFIX+x} ];then
    OUT_PREFIX="$RANDOM""_"
    logMsg "WARN" "This prefix will be used as output: $OUT_PREFIX"
fi 

logMsg "DEBUG" "OUT_PREFIX:($OUT_PREFIX)"

logMsg "INFO" "-------- START ERVcount ---------"

if [[ ! -e "$BAM" ]];then
    logMsg "ERROR" "Cannot find BAM file: ( $BAM )"
fi

logMsg "INFO" "---- Finding ERVs ----"
coverageBed -nonamecheck -a /resources/ERVmap.bed -b "$BAM" -counts -sorted > "$OUT_PREFIX""ERVresults.txt"
logMsg "INFO" "---- Finding ERVs complete ----"
logMsg "INFO" "-------- END ERVcount ---------"
