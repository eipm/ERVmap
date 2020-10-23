#!/usr/bin/env nextflow


// check parameters
if (!params.inputDir) {
    exit 1, 'inputDir parameter is missing.'
}
if ( !new File(params.inputDir).exists()) {
    exit 1, "The input folder does not exists."+params.inputDir+"\n"
}
if (!params.outputDir) {
    exit 1, "outputDir parameter is missing."
}
if (!params.starTmpDir) {
    exit 1, "starTmpDir parameter is missing."
}

if ( !new File(params.starTmpDir).exists()) {
    exit 1, 'The STAR temporary folder does not exists. ('+params.starTmpDir+')\n'
}
if (!params.outPrefix) {
    exit 1, "Output prefix parameter is missing."
}
if (!params.localOutDir) {
    params.localOutDir='bam'
}
if (!params.debug) {
    exit 1, "Debug prefix parameter is missing."
}


pairFiles_ch = Channel.fromFilePairs( params.inputDir+"/*{1,2}.fastq.gz", size: 2, checkIfExists: true )

process ERValign {
    tag "${sample}"
    
    // executor configuration
    time '8h'
    // memory '35.GB'
    scratch true
    cpus params.cpus
    storeDir params.outputDir

    // other configuration
    echo true
    errorStrategy 'terminate'

    input:
    val(outPrefix) from params.outPrefix
    val(localOutDir) from params.localOutDir
    val(limitMemory) from params.limitMemory
    val(debug) from params.debug
    tuple val(sample), file(reads) from pairFiles_ch
    
    output:
    path ( "${localOutDir}/${outPrefix}Aligned.sortedByCoord.out.bam" ) into bam_ch
    path ( "${localOutDir}/${outPrefix}Aligned.sortedByCoord.out.bam.bai" ) into bai_ch

    shell:
    template 'ERValign.sh'
}

process ERVcount {
    tag "${sample}"

    // executor configuration
    time '3h'
    memory '8.GB'
    scratch true
    storeDir params.outputDir

    // other configuration
    echo true
    errorStrategy 'terminate'
    mode = 'BED'
    stageInMode 'symlink'
    
    input:
    val(outPrefix) from params.outPrefix
    val(debug) from params.debug
    path (bam) from bam_ch
    path (bai) from bai_ch
    
    output: 
    path (params.outPrefix+'ERVresults.txt') into final_results_ch

    shell:
    template "ERVcount.sh"
}

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
