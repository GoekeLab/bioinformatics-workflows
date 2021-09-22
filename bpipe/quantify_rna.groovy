
title 'Bpipe Example - RNASeq Transcript Quantification'

REF='../test_data/transcriptome.fa'

index_reference = {

    doc 'Indexes the reference transcript fasta file'

    output.dir = 'index'

    produce('refseq.bin') {
        exec """
            salmon index -t $REF -i $output.dir
        """
    }
}

fastqc = {

    doc 'Calculates aggregate statistics from FASTQ reads for quality control'

    output.dir='qc'

    exec """
        fastqc --quiet $inputs.fq.gz  --outdir $output.dir
    """
}

quantify = {

    doc 'Counts transcripts present in reference fasta'

    output.dir = 'quant'

    exec """
        salmon quant -i index -l A -1 $input1.fq.gz -2 $input2.fq.gz --validateMappings -o $output.dir
    """
}

run {
   index_reference + [ fastqc, quantify ]
}




