#!/usr/bin/env nextflow


// check parameters
if (!params.inputDir) {
    exit 1, 'inputDir parameter is missing.'
}
if (!params.skipAlignment) {
    params.skipAlignment=false
}
if (!new File(params.inputDir).exists()) {
    exit 1, "The input folder does not exists."+params.inputDir+"\n"
}
if (!params.outputDir) {
    exit 1, "outputDir parameter is missing."
}
if (!params.starTmpDir) {
    exit 1, "starTmpDir parameter is missing."
}

if (!new File(params.starTmpDir).exists()) {
    exit 1, 'The STAR temporary folder does not exists. ('+params.starTmpDir+')\n'
}
if (!params.localOutputDir) {
    params.localOutputDir='bam'
}
if (!params.debug) {
    exit 1, "Debug prefix parameter is missing."
}

// Adding the option to skip the alignment if skipAlignment is set to true.
(fastqPairs_ch, bamPairs_ch) = ( params.skipAlignment
                 ? [Channel.empty(), Channel.fromFilePairs( params.inputDir+'/'+params.inputPattern, checkIfExists: true )  ]
                 : [Channel.fromFilePairs( params.inputDir+'/'+params.inputPattern, size: 2, checkIfExists: true ), Channel.empty() ] )

bamPairs_ch.flatten().collate(3).set { newBamPairs_ch }

process ERValign {
    tag "${sample}"
    
    // executor configuration
    time '8h'
    memory '35 GB'
    scratch true
    cpus params.cpus
    publishDir params.outputDir, mode: 'copy'
    stageOutMode 'rsync'

    // other configuration
    echo true
    errorStrategy 'terminate'

    input:
    tuple val(sample), path (reads) from fastqPairs_ch
    val (localOutputDir) from params.localOutputDir
    val (limitMemory) from params.limitMemory
    val (debug) from params.debug
    
    output:
    tuple sample, path ("${localOutputDir}/${sample}.Aligned.sortedByCoord.out.bam*") into star_bam_ch  
    // path ("${localOutputDir}/${sample}.Aligned.sortedByCoord.out.bam.bai") into star_bai_ch
    // val (sample) into prefix_ch

    when:
    !params.skipAlignment

    shell:
    template 'ERValign.sh'
}
star_bam_ch.subscribe{ println "File: ${it.name} => ${it.text}" }

// process ERVcount {
//     tag "${sample}"

//     // executor configuration
//     time '3h'
//     memory '8.GB'
//     scratch true
//     storeDir params.outputDir

//     // other configuration
//     echo true
//     errorStrategy 'terminate'
//     stageInMode 'symlink'
    
//     input:
//     tuple val(sample), path (bam), path (bai) from newBamPairs_ch
//     val (debug) from params.debug
//     // path (bam) from bam_ch.mix(star_bam_ch) // mixing with the channel from ERValign if started from FASTQs and skipAlignment=false
//     // path (bai) from bai_ch.mix(star_bai_ch)  // mixing with the channel from ERValign if started from FASTQs and skipAlignment=false
    
//     output: 
//     path ("${sample}"+'.ERVresults.txt') into final_results_ch

//     shell:
//     template 'ERVcount.sh'
// }

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
