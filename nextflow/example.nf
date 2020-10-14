#!/usr/bin/env nextflow

// help message
def helpMessage() {
    log.info"""
    Usage:

        nextflow run example.nf --ref ref.fasta[.gz] --left path/to/reads_1.fastq --right path/to/reads_2.fastq --outdir path/to/resultdir

    --ref       Path to reference transcript (fasta[.gz])
    --left      Path to left strand (fastq)
    --right     Path to right strand (fastq)
    --outdir    Output directory (where output results will be saved)  
    --help      Displays usage message  

    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


// Create reference transcriptome index using Salmom
process salmonIndex {

    output:
    file index into txome_index

    """
    salmon index -t '${params.ref}' -i index
    """
}


// Transcriptome alignment and quantification using Salmon
process salmonAlignQuant {

    input:
    file index from txome_index

    """
    salmon quant -i $index -l A -1 '${params.left}' -2 '${params.right}' --validateMappings -o '${params.outdir}'
    """
}

// FastQC
process fastQC {

    output:
    stdout result

    """
    fastqc --quiet '${params.left}' '${params.right}' --outdir '${params.outdir}'
    """
}
