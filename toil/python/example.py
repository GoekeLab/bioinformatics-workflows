import os
import subprocess

from toil.common import Toil
from toil.job import Job
from toil.lib.docker import apiDockerCall

#Docker flag to change for using docker images
DOCKER_FLAG = True

#Necessary file paths
fastQ_one = 'data/reads_1.fq'
fastQ_two = 'data/reads_2.fq'
transcript_fastA = 'data/transciptome.fa'


fastQCDock = Job.wrapJobFn(apiDockerCall,
                        image='pegi3s/fastqc',
                        working_dir= os.getcwd(),
                        parameters=['fastqc', fastQ_one, fastQ_two])

salmonIndexDock = Job.wrapJobFn(apiDockerCall,
                        image='combinelab/salmon',
                        working_dir='/home/hexotical/bioninformatics-workflows/toil/python',
                        parameters=['salmon', 'index', '-t', transcript_fastA, '-i', 'index'])

salmonQuantDock = Job.wrapJobFn(apiDockerCall,
                        image='combinelab/salmon',
                        working_dir='/home/hexotical/bioninformatics-workflows/toil/python',
                        parameters=['salmon', 'quant', '-i', 'index', '-l', 'A', '-1', fastQ_one, '-2', fastQ_two, '--validateMappings', '-o', 'quant'])

def fastQC(self):
    subprocess.run(["fastqc", fastQ_one, fastQ_two])


def salmonIndex(self):
    subprocess.run(["salmon", "index", "-t", transcript_fastA, '-i', 'index'])

def salmonAlignQuant(self):
    subprocess.run(["salmon", "quant", "-i", "index", "-l", "A", "-1", fastQ_one, "-2", fastQ_two, "--validateMappings", "-o", "quant"])

fastQ = Job.wrapJobFn(fastQC)
salmonI = Job.wrapJobFn(salmonIndex)
salmonQ = Job.wrapJobFn(salmonAlignQuant)



if __name__=="__main__":
    options = Job.Runner.getDefaultOptions("./toilWorkflowRun")
    options.logLevel = "DEBUG"
    options.clean = "always"


    # We add a child to ensure that the necessary files exist by the time Toil attempts to run the job
    salmonIndexDock.addChild(salmonQuantDock)
    salmonI.addChild(salmonQ)

    with Toil(options) as toil:
        if(DOCKER_FLAG):
            toil.start(fastQCDock)
            toil.start(salmonIndexDock)
            toil.start(salmonQuantDock)
        else:
            toil.start(fastQ)
            toil.start(salmonI)
            toil.start(salmonQ)
