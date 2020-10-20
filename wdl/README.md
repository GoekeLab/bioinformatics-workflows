## Running the proof of concept WDL pipeline

1. Install `cromwell` by following the instructions [here](https://snakemake.readthedocs.io/en/stable/)
2. Install the files from this repo locally
   - `example.wdl` is the workflow specification
   - `wdl_input.json` is the input data - change this file to point to your desired paths
   - `options.json` allows you to specify extra options for `cromwell`, such as an output directory.
3. Run the pipeline with the following command:
   - `snakemake -n results/salmon/quant/ results/fastqc/`
   - to run snakemake we tell the program what output files we want rather than passing input parameters
4. Run with the following command, specifying the number of cores
    - `java -jar path/to/cromwell-53.1.jar run example.wdl --inputs wdl_input.json -o options.json`