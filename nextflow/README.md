## Running the proof of concept Nextflow pipeline

1. Install `nextflow` by following the instructions [here](https://www.nextflow.io/)
2. Optionally move the `nextflow` directory to a directory in your path
   i.e.  `echo $PATH` and then `mv nextflow usr/bin`
3. Navigate to the directory where you saved `example.nf`
4. Run the pipeline with the following command:
   - `nextflow run example.nf --ref ref.fasta[.gz] --left path/to/reads_1.fastq --right path/to/reads_2.fastq --outdir path/to/resultdir`
5. Or run `nextflow run example.nf --help` for more information

