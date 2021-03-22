# bioinformatics-workflows

## Introduction
These pipelines are intended as a proof-of-concept demonstrate the distinctions between bioinformatics workflow managers.
The analysis workflow used for the proof-of-concept is a small portion of an RNA-seq workflow, using the tool `salmon` for transcriptome alignment and quantification and `fastqc` for quality control.

![](docs/imgs/workflow_online_description_small.png)

## Tools
To run locally, install [salmon](https://github.com/COMBINE-lab/salmon/) and [fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
   -  You will want to add a symbolic link to the executables to your `$PATH`
   -  `sudo ln -s /absolute/path/salmon/bin/salmon /usr/local/bin/salmon`
   -  `sudo ln -s /absolute/path/FastQC/fastqc /usr/local/bin/fastqc`
   -  You can test this works by using the help function of the two tools (i.e. `salmon -h` and `fastqc -h`)
   
## Test Data
Installing `salmon` as also installs a `sample_data` directory inside the `salmon` directory.
The files `reads_1.fastq`, `reads_2.fastq` and `transcripts.fasta` were as the sample data set to test these pipelines.

## Documentation
Each workflow manager folder in the repository has a README detailing how to run the proof-of-concept pipeline.
