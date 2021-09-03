## About Toil
Toil is a scalable, efficient, cross-platform pipleline management system written entirely in Python and designed around the principles of functional programming.
https://toil.ucsc-cgl.org/

## Training material and documentation
https://toil.readthedocs.io/en/latest/

Tutorial for using Toil Python API for workflows 
https://toil.readthedocs.io/en/latest/developingWorkflows/developing.html

Tutorial for using Toil with WDL
https://toil.readthedocs.io/en/latest/running/wdl.html

Tutorial for using Toil with CWL
https://toil.readthedocs.io/en/latest/running/cwl.html

## Community-developed workflows in Toil
Add description and links to community-developed pipelines.

## Running the proof of concept Toil pipeline

# General

1. Install Toil by following the instructions listed [here](https://toil.readthedocs.io/en/latest/gettingStarted/install.html)

2. (Optional) Install Docker by following the instructions listed [here](https://docs.docker.com/get-docker/)
Docker isn't a necessity however we offer the option of using docker to increase the portability of the workflows

# WDL

1. Install the files under `toil/wdl` locally
    - `example.wdl` is the workflow specification 
    - `wdl_input.json` is the input data, change this file to point to desired paths
2. Run the pipeline with the following command:
    `toil-wdl-runner example.wdl wdl_input.json`

# CWL

1. Install the files under `toil/cwl` locally
    - `example.cwl` is the workflow specification
    - `fastqc.cwl` is a command line tool used by the workflow
    - `salmonIndex.cwl` is a command line tool used by the workflow
    - `salmonQuant.cwl` is a command line tool used by the workflow
    - `cwl_input.yaml` is the input data, change this file to point to desired paths
2. Run the pipeline with the following command:
    `toil-cwl-runner example.cwl cwl_input.yaml`

# Python

1. Install the files under `toil/python` locally
    - `example.py` is both workflow specification and input
2. Run the pipeline with the following command:
    `python example.py file:my-job-store`


## Notes and Contribution
This pipeline is a minimal example of using Toil. We welcome contributions to the documentation and workflow, please create an issue or submit a pull request!

## How to cite X

