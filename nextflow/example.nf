#!/usr/bin/env nextflow


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
