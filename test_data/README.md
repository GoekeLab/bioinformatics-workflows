

## Test data for the RNA-Seq workflow

This folder contains a test data set for the RNA-Seq workflow.

- [reads_1.fq.gz](reads_1.fq.gz) and [reads_2.fq.gz](reads_2.fq.gz): Fastq files containing the simulated paired-end RNA-Seq data
- [transcriptome.fa](transcriptome.fa): Transcriptome reference fasta file (representing 582 transcripts)
- [truth.tsv](truth.tsv): True read counts for the simulated data

## Test data generation

The test data set was generated to represent by the following procedure (thanks to @rob-p):

1. The sample [ERR188297](https://www.ebi.ac.uk/ena/browser/view/ERR188297) was downloaded from ENA (this is an experimental sample from [GEUVADIS](https://www.ebi.ac.uk/ena/browser/view/PRJEB3366)).

2. The sample was quantified against the [Gencode v38 human transcriptome](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz).

3. The results were loaded in R with [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html) and aggregating to the gene level.

4. Expressed genes were randomly pulled out until the sum of their estimated read counts exceeded 100,000 (resulting in 66 genes).

5. All transcripts from these genes were selected to generate the [test data transcriptome reference file](transcriptome.fa) (582 transcripts).

6.The estimated transcript level counts were then used to simulate the test data with [polyester](https://bioconductor.org/packages/release/bioc/html/polyester.html) using `simulate_experiment_countmat`.

7. The reads were shuffled (while maintaining the pairing) using [bbmap](https://sourceforge.net/projects/bbmap/).

8. Fake quality scores were added to the reads, using [bbmap](https://sourceforge.net/projects/bbmap/).
