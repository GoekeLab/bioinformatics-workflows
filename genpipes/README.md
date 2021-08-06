## About GenPipes

[GenPipes](https://www.computationalgenomics.ca/genpipes/) is a flexible Python-based framework that facilitates the 
development and deployment of multi-step workflows optimized for High-Performance Computing clusters and the cloud. 
GenPipes comes with 13 validated and scalable pipelines for various genomics applications.

## Training material and documentation

- [Link to tutorial](https://genpipes.readthedocs.io/en/master/tutorials/list_tutorials.html)
- [Link to online documentation](https://genpipes.readthedocs.io/en/master/index.html)
- [Link to developer guide](https://genpipes.readthedocs.io/en/latest/development/dev_guide.html)

## Community-developed workflows in GenPipes

- [Amplicon Sequencing Pipeline](https://genpipes.readthedocs.io/en/latest/user_guide/pipelines/gp_ampliconseq.html)
- [ChIP Sequencing Pipeline](https://genpipes.readthedocs.io/en/latest/user_guide/pipelines/gp_chipseq.html)
- [CoV Sequencing Pipeline](https://genpipes.readthedocs.io/en/latest/user_guide/pipelines/gp_covseq.html)
- [DNA Sequencing Pipeline](https://genpipes.readthedocs.io/en/latest/user_guide/pipelines/gp_dnaseq.html)
- [DNA Sequencing High Coverage Pipeline](https://genpipes.readthedocs.io/en/latest/user_guide/pipelines/gp_dnaseq_highcov.html)
- [HiC Sequencing Pipeline](https://genpipes.readthedocs.io/en/latest/user_guide/pipelines/gp_hicseq.html)
- [Illumina Run Processing Pipeline](https://genpipes.readthedocs.io/en/latest/user_guide/pipelines/gp_illumina.html)
- [Methylation Sequencing Pipeline](https://genpipes.readthedocs.io/en/latest/user_guide/pipelines/gp_wgs_methylseq.html)
- [Nanopore Pipeline](https://genpipes.readthedocs.io/en/latest/user_guide/pipelines/gp_nanopore.html)
- [RNA Sequencing Pipeline](https://genpipes.readthedocs.io/en/latest/user_guide/pipelines/gp_rnaseq.html)
- [De-Novo RNA Sequencing Pipeline](https://genpipes.readthedocs.io/en/latest/user_guide/pipelines/gp_rnaseq_denovo.html)
- [RNA Sequencing Light Pipeline](https://genpipes.readthedocs.io/en/latest/user_guide/pipelines/gp_rnaseq_light.html)
- [Tumor Pair Sequencing Pipeline](https://genpipes.readthedocs.io/en/latest/user_guide/pipelines/gp_tumourpair.html)

## Running the proof of concept GenPipes pipeline

This proof of concept pipeline is implemented in the following [branch](https://bitbucket.org/mugqic/genpipes/src/goekelab_manuscript). 
To run this proof of concept, you need to clone the repository and change to the appropriate branch. The location of the cloned directory
will be referred to as `$GENPIPES` in the following instructions. 

1. From the cloned location, launch the pipeline using the attached [readset.txt](./readset.txt) file: 

```
python $GENPIPES/pipelines/goeke_rnaseq.py -c $GENPIPES/pipelines/goeke_rnaseq.slurm.ini --steps 1-3 --output pipeline_output --readsets readset.txt --job-scheduler slurm > run_commands.sh 
	
```
 
2. To submit jobs to the scheduler just run the following command: 
```
bash run_commands.sh
```

*Notes:*
 
- Please note that this proof of concept pipeline requires de-compressed fastqs to work properly
- When running from Compute Canada or other supported GenPipes clusters, the pipeline will load software using linux modules, if running from 
another location, use the [containerized version of GenPipes](https://genpipes.readthedocs.io/en/master/tutorials/genpipes_in_the_container.html).
- GenPipes currently runs using python v2.7, a python v3 release is coming soon. 

## Notes and Contribution

This pipeline is a minimal example of a Salmon-based pipeline using GenPipes. We welcome contributions to the documentation 
and workflow, please create an issue or submit a pull request!

## How to cite GenPipes

Bourgey M, Dali R, Eveleigh R, Chen KC, Letourneau L, Fillon J, Michaud M, Caron M, Sandoval J, Lefebvre F, Leveque G, Mercier E, 
Bujold D, Marquis P, Van PT, Anderson de Lima Morais D, Tremblay J, Shao X, Henrion E, Gonzalez E, Quirion PO, Caron B, Bourque G. 
**GenPipes: an open-source framework for distributed and scalable genomic analyses.** *Gigascience.* 2019 Jun 1;8(6):giz037. 
doi: 10.1093/gigascience/giz037. PMID: 31185495; PMCID: PMC6559338.
