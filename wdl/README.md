## Running the proof of concept WDL pipeline

1. Install `cromwell` by following the instructions [here](https://support.terra.bio/hc/en-us/articles/360037487871)
2. Install the files from this repo locally
   - `example.wdl` is the workflow specification
   - `wdl_input.json` is the input data - change this file to point to your desired paths
   - `options.json` allows you to specify extra options for `cromwell`, such as an output directory. More information [here](https://cromwell.readthedocs.io/en/stable/wf_options/Overview/)
3. Run the pipeline with the following command:
    - `java -jar path/to/cromwell-53.1.jar run example.wdl --inputs wdl_input.json -o options.json`
