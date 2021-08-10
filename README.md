# Reproducible, scalable, and shareable analysis pipelines with bioinformatics workflow managers

Workflow managers provide an easy and intuitive way to simplify pipeline development. Here we provide basic proof-of-concept implementations for selected workflow managers. The analysis workflow is based on a small portion of an RNA-seq pipeline, using fastqc for quality controls and salmon for transcript quantification.
These implementations are designed for basic illustrations. Workflow managers provide many more powerful features than what we use here, please visit the official documentations to explore those in detail.

![](docs/imgs/workflow_online_description_small.png)

## The RNA-Seq workflow

The RNA-Seq analysis workflow performs quality controls with fastqc and quantifies transcripts expression using Salmon. Here we will use local installation (see documentation for [salmon](https://github.com/COMBINE-lab/salmon/) and [fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)). For the local installations you can add a symbolic link to the executables to your $PATH:

`sudo ln -s /absolute/path/salmon/bin/salmon /usr/local/bin/salmon`

`sudo ln -s /absolute/path/FastQC/fastqc /usr/local/bin/fastqc`

You can test the installation using the help function of the two tools (i.e. `salmon -h` and `fastqc -h`). 
   
## Test Data
This repository contains a [simulated test data set](test_data) which can be used to run the example implementations. The test data contains RNA-Seq reads ([reads_1.fq.gz](test_data/reads_1.fq.gz) and [reads_2.fq.gz](test_data/reads_2.fq.gz)), a transcriptome reference file ([transcriptome.fa](test_data/transcriptome.fa)) and the true counts from the simulation experiments ([truth.tsv](test_data/truth.tsv))

## Basic proof-of-concept implementations
Each workflow manager folder in this repository has a README detailing how to run the proof-of-concept pipeline:

- [Galaxy](galaxy)
- [Nextflow](nextflow)
- [Snakemake](snakemake)
- [SciPipe](scipipe)
- [WDL](wdl)
- [GenPipes](genpipes)

## Online Documentation for Workflow managers
Workflow managers have many more features which are not used in these implementations, and there are many additional workflow managers. You can read more about each workflow manager in their official documentation:

- [Galaxy](https://docs.galaxyproject.org/en/master/)
- [KNIME](https://docs.knime.com/)
- [Nextflow](https://www.nextflow.io/docs/latest/index.html)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [GenPipes](https://genpipes.readthedocs.io/en/genpipes-v-3.3.0/)
- [Bpipe](http://docs.bpipe.org/)
- [Pachyderm](https://docs.pachyderm.com/latest/)
- [SciPipe](https://scipipe.org/)
- [Luigi](https://luigi.readthedocs.io/en/stable/)
- [WDL](https://openwdl.org/)
- [CWL](https://www.commonwl.org/user_guide/index.html)
- [Toil](https://toil.readthedocs.io/en/latest/)

## Contact and Call for Contribution
This repository was created by [Laura Wratten](https://github.com/lwratten). We very much encourage contributions by users of these workflows. If you would like to add an implementation for any of these workflow managers you can follow the [template](template/README.md). If you would like to suggest changes to any of the existing implementations, please raise an issue and submit a pull request.

## Acknowledgements

We would like to thank the following people for their contribution to this repository:
- [Rob Patro](https://github.com/rob-p) for creating the [RNA-Seq test data set](test_data)
- [Paolo Di Tommaso](https://github.com/pditommaso) for the [Nextflow workflow](nextflow)
- [Johannes Köster](https://github.com/johanneskoester) for the [Snakemake workflow](snakemake)
- [Samuel Lampa](https://github.com/samuell) for the [SciPipe workflow](scipipe)
