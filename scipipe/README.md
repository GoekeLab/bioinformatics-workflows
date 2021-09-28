## About SciPipe

SciPipe is a workflow library implemented as a programming library in Go.
See [scipipe.org](https://scipipe.org) for more information!

## Training material and documentation

- [Link to tutorial](https://scipipe.org/writing_workflows/)
- [Link to online documentation](https://scipipe.org)
- [A basic tutorial as a screencast](https://www.youtube.com/watch?v=hi0Uqwddrtg)
- [Case study workflows for the SciPipe paper](https://github.com/pharmbio/scipipe-demo)

## Community-developed workflows in SciPipe

- "Reproducible Probabilistic Target Profiles"
  - Repository: https://github.com/pharmbio/ptp-project
  - Description: A workflow developed by the workflow tool authors, for building
    target binding profiles, using QSAR and SVM (see the [paper](https://doi.org/10.3389/fphar.2018.01256) for more details).
- Kleuren SciPipe workflow
  - Repository: https://github.com/Colelyman/kleuren-scipipe-workflow
  - Description: A scipipe workflow for building the necessary data structures
    for running kleuren.

## Running the proof of concept SciPipe pipeline

1. Install the Go toolchain by following [this page](https://go.dev/learn/) for
   the specific instructions for your operating system.
   - Note that on Windows, you need to have a POSIX-compliant bash environment
     to run Scipipe, such as Windows Subsystem for Linux (WSL). Find more
     information about [installing WSL here](https://docs.microsoft.com/en-us/windows/wsl/install-win10)
2. Go into the scipipe directory:
   ```bash
   cd scipipe
   ```
3. Run the workflow using the `run.sh` bash script:
   ```bash
   ./run.sh
   ```
   Alternatively you can run the full command directly:
   ```bash
   go run example.go -ref ../test_data/transcriptome.fa -left ../test_data/reads_1.fq.gz -right ../test_data/reads_2.fq.gz -outdir results
   ```

Done!

## Notes and Contribution

This pipeline is a minimal example of using SciPipe. We welcome contributions to the
documentation and workflow, please create an issue or submit a pull request!

## How to cite SciPipe

Please cite SciPipe as:

> Lampa, S., Dahlö, M., Alvarsson, J., & Spjuth, O. (2019). SciPipe: A workflow
> library for agile development of complex and dynamic bioinformatics pipelines.
> GigaScience, 8(5), giz044. DOI: [10.1093/gigascience/giz044](https://doi.org/10.1093/gigascience/giz044)
