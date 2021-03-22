## About WDL
The Workflow Description Language ([WDL](https://openwdl.org/)) is a workflow specification language. Workflow written in WDL can be executed with [cromwell](https://cromwell.readthedocs.io/en/stable/) or [Toil(experimental)](https://toil.readthedocs.io/en/latest/gettingStarted/quickStart.html#running-a-basic-wdl-workflow).

## Training material and documentation
- [WDL tutorials](https://support.terra.bio/hc/en-us/sections/360007347652-WDL-Tutorials)
- [WDL documentation](https://support.terra.bio/hc/en-us/sections/360007274612-WDL-Documentation)

## Community-developed workflows in snakemake
[BioWDL](https://github.com/biowdl) contains a collection of community developed pipelines written in WDL.

## Running the proof of concept WDL pipeline

1. Install `cromwell` by following the instructions [here](https://support.terra.bio/hc/en-us/articles/360037487871)
2. Install the files from this repo locally
   - `example.wdl` is the workflow specification
   - `wdl_input.json` is the input data - change this file to point to your desired paths
   - `options.json` allows you to specify extra options for `cromwell`, such as an output directory. More information [here](https://cromwell.readthedocs.io/en/stable/wf_options/Overview/)
3. Run the pipeline with the following command:
    - `java -jar path/to/cromwell-53.1.jar run example.wdl --inputs wdl_input.json -o options.json`

## Notes and Contribution
This pipeline is a minimal example of using WDL. We welcome contributions to the documentation and workflow, please create an issue or submit a pull request!
