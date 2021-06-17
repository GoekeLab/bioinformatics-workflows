## About snakemake
[snakemake](https://snakemake.github.io) is a workflow management system relying on a Python based DSL.

## Training material and documentation
- [snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html)
- [snakemake documentation](https://snakemake.readthedocs.io)

## Community-developed workflows in snakemake
Best practices workflows for snakemake can be found here: https://github.com/snakemake-workflows/docs

## Running the proof of concept Snakemake pipeline

1. Install `snakemake` by following the instructions [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
2. In the directory where your snakefile is located, create a directory called `data`.
3. Copy your data into the `data` folder with the following file names (these are the original file names of the salmon sample_data; using different paths or names would just require a modification of the configuration file `config/config.yaml`):
   - `reads_1.fastq`
   - `reads_2.fastq`
4. Run the a dry-run of the pipeline with the following command:
   - `snakemake -n --use-conda`
   - Snakemake will make use of the configuration file `config/config.yaml` provided in this example to determine which samples to run.
5. Run with the following command, specifying the number of cores
    - `snakemake --cores [numcores] --use-conda`
    - For example, to run with 2 cores the command would be:
    - `snakemake --cores 2 --use-conda`
    - All required software will be automatically deployed by Snakemake via the [mamba package manager](https://github.com/mamba-org/mamba)

## Notes and Contribution
We welcome contributions to the documentation and workflow, please create an issue or submit a pull request!

## How to cite Snakemake
Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33, [DOI:10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2).
