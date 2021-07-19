## About Nextflow
Nextflow is a DSL workflow manager (https://www.nextflow.io/).

## Training material and documentation
- [Nextflow tutorials](https://nf-co.re/usage/nextflow)
- [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html)

## Community-developed workflows in Nextflow
The [nf-core projects](https://nf-co.re/) hosts community curated Nextflow pipelines.

## Running the proof of concept Nextflow pipeline

1. Install `nextflow` by following the instructions [here](https://www.nextflow.io/)
2. Optionally move the `nextflow` directory to a directory in your path
   i.e.  `echo $PATH` and then `mv nextflow usr/bin`
3. Navigate to the directory where you saved `example.nf`
4. Run the pipeline with the following command:

      ```
      nextflow run nextflow/example.nf \
            --left $PWD/test_data/reads_1.fq.gz \
            --right $PWD/test_data/reads_2.fq.gz \
            --ref $PWD/test_data/transcriptome.fa
      ```
5. Enable the use Docker containers adding the following option to the above command line: `-with-docker quay.io/nextflow/rnaseq-nf:v1.0`

Note: input files requires the use of absolute paths

## Notes and Contribution
This pipeline is a minimal example of using Nextflow. We welcome contributions to the documentation and workflow, please create an issue or submit a pull request!

## How to cite Nextflow
Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316–319. doi:10.1038/nbt.3820
