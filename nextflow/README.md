## About nextflow
Nextflow is a DSL workflow manager (https://www.nextflow.io/).

## Training material and documentation
- [nextflow tutorials](https://nf-co.re/usage/nextflow)
- [nextflow documentation](https://www.nextflow.io/docs/latest/index.html)

## Community-developed workflows in nextflow
The [nf-core projects](https://nf-co.re/) hosts community curated nextflow pipelines.

## Running the proof of concept Nextflow pipeline

1. Install `nextflow` by following the instructions [here](https://www.nextflow.io/)
2. Optionally move the `nextflow` directory to a directory in your path
   i.e.  `echo $PATH` and then `mv nextflow usr/bin`
3. Navigate to the directory where you saved `example.nf`
4. Run the pipeline with the following command:
   - `nextflow run example.nf --ref ref.fasta[.gz] --left path/to/reads_1.fastq --right path/to/reads_2.fastq --outdir path/to/resultdir`
5. Or run `nextflow run example.nf --help` for more information

## Notes and Contribution
This pipeline is a minimal example of using nextflow. We welcome contributions to the documentation and workflow, please create an issue or submit a pull request!
