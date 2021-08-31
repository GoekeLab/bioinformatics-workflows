import os
import subprocess

from toil.common import Toil
from toil.job import Job
from toil.lib.docker import apiDockerCall
'''
align = Job.wrapJobFn(apiDockerCall,
                      image='ubuntu',
                      working_dir=os.getcwd(),
                      parameters=['ls', '-lha'])
'''
def fastQC():
    subprocess.run(["fastqc", "data/reads_1.fq", "data/reads_2.fq"])
    #os.system("fastqc data/reads_1.fq")
#"/home/hexotical/bioinformatics-workflows/test_data/reads_1.fq.gz"
fastQCone = Job.wrapJobFn(apiDockerCall,
                        image='pegi3s/fastqc',
                        working_dir='/home/hexotical/bioninformatics-workflows/toil/python',
                        parameters=['fastqc', 'data/reads_1.fq', 'data/reads_2.fq'])
'''
salmonIndex = Job.wrapJobFn(apiDockerCall,
                        image='combinelab/salmon',
                        working_dir='/home/hexotical/bioninformatics-workflows/toil/python',
                        parameters=['salmon', 'index', '-t', '/home/hexotical/bioinformatics-workflows/toil/python/data/transcriptome.fa'])
'''
def salmonIndex():
    subprocess.run(["salmon", "index", "-t", "data/transcriptome.fa"])
    #subprocess.run(['tar', '-cvzf', 'index.tar.gz', 'index'])
#almon index -t "${ref_txome}" -i ${index}; tar -cvzf index.tar.gz ${index}
#/home/hexotical/bioinformatics-workflows/toil/python/data/transcriptome.fa
#salmonAlignQuaint = Job.wrapJobFn
def salmonAlignQuaint():
    subprocess.run(["salmon", "quant", "-i", "index", "-l", "A", "-1", "data/reads_1.fq", "-2", "data/reads_2.fq", "--validateMappings", "-o", "quant"])
#salmon quant -i "${indexdir}" -l A -1 "${reads1}" -2 "${reads2}" --validateMappings -o quant;
#AddChildren


if __name__=="__main__":
    options = Job.Runner.getDefaultOptions("./toilWorkflowRun")
    options.logLevel = "DEBUG"
    options.clean = "always"

    with Toil(options) as toil:
       toil.start(Job.wrapFn(fastQC))
       toil.start(Job.wrapFn(salmonIndex))
       toil.start(Job.wrapFn(salmonAlignQuaint))
