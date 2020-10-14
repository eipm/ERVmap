#!/usr/bin/env nextflow

process STAR-Alignment {
    tag ${sample}
    echo true
    scratch true
    errorStrategy 'terminate'

    input:
    

    output:

    shell:
    template 'ERVmapping.sh'
}

process ERV-counting {
        tag ${sample}
    echo true
    scratch true
    errorStrategy 'terminate'

    input:
    

    output:

    shell:
    template 'ERVmapping.sh'
}

process cleanUp {
    
}