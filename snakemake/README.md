## About snakemake
[snakemake](https://snakemake.readthedocs.io/en/stable/#) is a DSL workflow manager based on Python.

## Training material and documentation
- [snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html)
- [snakemake documentation](https://snakemake.readthedocs.io/en/stable/index.html)

## Community-developed workflows in snakemake
Best practices workflows for snakemake can be found here: https://github.com/snakemake-workflows/docs

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

## Notes and Contribution
This pipeline is a minimal example of using snakemake. We welcome contributions to the documentation and workflow, please create an issue or submit a pull request!
