#!/usr/bin/env nextflow


// check parameters
// if (!params.read1) {
//     exit 1, "read1 parameter is missing."
// }
// if (!params.read2) {
//     exit 1, "read2 parameter is missing."
// }
// if (!new File(params.read1).exists()) {
//     exit 1, "The read1 file does not exists".params.read1."\n"
// }
// if (!new File(params.read2).exists()) {
//     exit 1, "The read2 file does not exists.".params.read2."\n"
// }
// if (!params.genome) {
//     exit 1, "genome parameter is missing."
// }
// if (!new File(params.genome).exists()) {
//     exit 1, "The indexed genome does not exists.".params.genome."\n"
// // }
// if (!params.ervbed) {
//     exit 1, "ERVbed parameter is missing."
// }
// if (!new File(params.ervbed).exists()) {
//     exit 1, "The ERV bed file does not exists.".params.ervbed."\n"
// }
if (!params.inputDir) {
    exit 1, "inputDir parameter is missing."
}
if (!params.outputDir) {
    exit 1, "outputDir parameter is missing."
}
if ( !new File(params.inputDir).exists()) {
    exit 1, "The input folder does not exists."+params.inputDir+"\n"
}
if (!params.outPrefix) {
    exit 1, "Output prefix parameter is missing."
}
if (!params.debug) {
    exit 1, "Debug prefix parameter is missing."
}
pairFiles_ch = Channel.fromFilePairs( params.inputDir+"/*{1,2}.fastq.gz", size: 2, checkIfExists: true )
bam_align_ch = Channel.fromPath( params.outputDir )
bam_count_ch = Channel.fromPath( params.outputDir )

process ERValign {
    // tag ${sample}
    
    // executor configuration
    time '8h'
    // memory '35.GB'
    scratch true
    
    // other configuration
    echo true
    errorStrategy 'terminate'
    
    mode='STAR'

    input:
    val(outPrefix) from params.outPrefix    
    val(mode) from mode
    val(cpus) from params.cpus
    val(limitMemory) from params.limitMemory
    val(debug) from params.debug
    tuple val(sample), file(reads) from pairFiles_ch
    path results from bam_align_ch

    output:
    path ( "results/${outPrefix}Aligned.sortedByCoord.out.bam" ) into star_bam_ch

    // """
    // echo 'STAR:' "$mode - $sample - $reads"
    // """
    shell:
    template 'ERValign.sh'
}

process ERVcount {
//     tag ${sample}

    // executor configuration
    time '3h'
//     memory '8.GB'
    scratch true

    // other configuration
    echo true
    errorStrategy 'terminate'
    mode = 'BED'
    stageInMode 'symlink'
    publishDir 'results/nextflow'

    input:
    path (bam) from star_bam_ch
    val(outPrefix) from params.outPrefix
    val(debug) from params.debug
    path bam from bam_count_ch
    
    output: 
    path ("results") into final_results_ch

    shell:
    template "ERVcount.sh"
}

final_results_ch.subscribe { println "File: ${it.name}" }

// ~~~~~~~~~~~~~~~ PIPELINE COMPLETION EVENTS ~~~~~~~~~~~~~~~~~~~ //

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${workflow.success ? 'OK' : 'FAILED'}"

    // getting the default data from the workflow
    // PipelineMessage.completed(workflow).forTopic(pipeline_topic_completed).data(runIDKey, "${runID}").send()
}

workflow.onError {
    log.info "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
    log.info "Pipeline execution stopped with the following report: ${workflow.errorReport}"

    // getting the default data from the workflow
    // PipelineMessage.progress().forTopic(pipeline_topic_failure).data(runIDKey, "${runID}").data('onError', 'inside').send()
    // PipelineMessage.failed(workflow).forTopic(pipeline_topic_failure).data(runIDKey, "${runID}").send()
}


// ~~~~~~~~~~~~~~~ HELPER FUNCTIONS ~~~~~~~~~~~~~~~~~~~ //

/**
 * Create a dir.
 *
 * @param dir
 * @return
 */
def mkDir(def dir) {
    "mkdir -p ${dir}".execute()
    "chmod -R a+r+w ${dir}".execute()
}