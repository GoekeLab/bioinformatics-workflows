1. Install `nextflow` by following the instructions [here](https://www.nextflow.io/)
2. Optionally move the `nextflow` directory to a directory in your path
   i.e.  `echo $PATH` and then `mv nextflow usr/bin`
3. To run locally, install `salmon` and `fastqc`
   1. You will want to add a symbolic link to the executables to your file path
   2. `sudo ln -s /absolute/path/salmon/bin/salmon /usr/local/bin/salmon`
   3. `sudo ln -s /absolute/path/FastQC/fastqc /usr/local/bin/fastqc`
   4. You can test this works by using the help function of the two tools
      1. `salmon -h` `fastqc -h`

