## Running the proof of concept Snakemake pipeline

1. Install `snakemake` by following the instructions [here](https://snakemake.readthedocs.io/en/stable/)
2. In the directory where your snakefile is located, create a directory called `data`.
3. Copy your data into the `data` folder with the following file names (these are the original file names of the salmon sample_data):
   - `reads_1.fastq`
   - `reads_2.fastq`
   - `transcripts.fasta`
4. Run the a dry-run of the pipeline with the following command:
   - `snakemake -n results/salmon/quant/ results/fastqc/`
   - to run snakemake we tell the program what output files we want rather than passing input parameters
5. Run with the following command, specifying the number of cores
    - `snakemake --cores [numcores] results/salmon/quant/ results/fastqc/`
    - For example, to run with 2 cores the command would be:
    - `snakemake --cores 2 results/salmon/quant/ results/fastqc/`
