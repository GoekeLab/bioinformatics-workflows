#!/usr/bin/env nextflow
params.ref = 'path/to/ref.fasta'
params.left = 'path/to/reads_1.fastq'
params.right = 'path/to/reads_2.fastq'
params.outdir = 'results'


// Create reference transcriptome index using Salmom
process salmonIndex {
    input: 
    path ref from params.ref
    output:
    path index into txome_index

    """
    salmon index -t '${ref}' -i index
    """
}


// Transcriptome alignment and quantification using Salmon
process salmonAlignQuant {
    publishDir params.outdir

    input:
      path index from txome_index
      path left from params.left
      path right from params.right
    output:
      path 'quant'

    """
    salmon quant -i $index -l A -1 '${left}' -2 '${right}' --validateMappings -o quant
    """
}

// FastQC
process fastQC {
    publishDir params.outdir

    input:
      path index from txome_index
      path left from params.left
      path right from params.right
    output:
      path 'qc'

    """
    mkdir qc && fastqc --quiet '${params.left}' '${params.right}' --outdir qc
    """
}
