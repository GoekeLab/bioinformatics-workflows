# FastQC
rule fastqc:
    input:
        "data/reads_1.fastq",
        "data/reads_2.fastq"
    output:
        directory("results/fastqc")
    params:
        "--quiet --outdir"
    shell:
        "mkdir results/fastqc; fastqc {input} {params} {output}"

# Create reference transcriptome index using Salmon
rule salmonIndex:
    input:
        "data/transcripts.fasta"
    output:
        directory("results/salmon/index")
    shell:
        "salmon index -t {input} -i {output}"

# Transcriptome alignment and quantification using Salmon
rule salmonAlignQuant:
    input:
        index = "results/salmon/index",
        left = "data/reads_1.fastq",
        right = "data/reads_2.fastq"
    output:
        directory("results/salmon/quant")
    shell:
        "salmon quant -i {input.index} -l A -1 {input.left} -2 {input.right} --validateMappings -o {output}"
