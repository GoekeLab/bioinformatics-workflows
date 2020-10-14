## Running the proof of concept Nextflow pipeline

1. Install `nextflow` by following the instructions [here](https://www.nextflow.io/)
2. Optionally move the `nextflow` directory to a directory in your path
   i.e.  `echo $PATH` and then `mv nextflow usr/bin`
3. To run locally, install [salmon](https://github.com/COMBINE-lab/salmon/) and [fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
   -  You will want to add a symbolic link to the executables to your `$PATH`
   -  `sudo ln -s /absolute/path/salmon/bin/salmon /usr/local/bin/salmon`
   -  `sudo ln -s /absolute/path/FastQC/fastqc /usr/local/bin/fastqc`
   -  You can test this works by using the help function of the two tools (i.e. `salmon -h` and `fastqc -h`)
4. Navigate to the directory where you saved `example.nf`
5. Run the pipeline with the folling command:
   - `nextflow run example.nf --ref ref.fasta[.gz] --left path/to/reads_1.fastq --right path/to/reads_2.fastq --outdir path/to/resultdir`
6. Or run `nextflow run example.nf --help` for more information

