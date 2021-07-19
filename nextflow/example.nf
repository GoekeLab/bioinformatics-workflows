#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define pipeline inpit parameretes 
// note: input files requires the use of absolute paths
params.ref = '/path/to/ref.fasta'
params.left = '/path/to/reads_1.fastq'
params.right = '/path/to/reads_2.fastq'
params.outdir = 'results'


// Create reference transcriptome index using Salmom
process SALMON_INDEX {
  input: 
    path ref
  output:
    path index

  """
    salmon index -t '${ref}' -i index
  """
}


// Transcriptome alignment and quantification using Salmon
process SALMON_ALIGN_QUANT {
  publishDir params.outdir

  input:
    path index
    path left 
    path right
  output:
    path 'quant'

  """
    salmon quant -i $index -l A -1 '${left}' -2 '${right}' --validateMappings -o quant
  """
}

// FastQC
process FASTQC {
  publishDir params.outdir

  input:
    path index
    path left
    path right
  output:
    path 'qc'

  """
    mkdir qc && fastqc --quiet '${params.left}' '${params.right}' --outdir qc
  """
}

workflow {
  SALMON_INDEX(params.ref)
  SALMON_ALIGN_QUANT( SALMON_INDEX.out, params.left, params.right )
  FASTQC( SALMON_INDEX.out, params.left, params.right )
}
