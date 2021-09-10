import os
import errno
import logging
from toil import customDockerInitCmd

from toil.job import Job
from toil.common import Toil
from toil.lib.docker import apiDockerCall
from toil.wdl.wdl_functions import generate_docker_bashscript_file

current_working_dir = os.getcwd()

logger = logging.getLogger(__name__)


class FastQConeCls(Job):
    def __init__(self, reads=None, *args, **kwargs):
        super(FastQConeCls, self).__init__(*args, **kwargs)
        Job.__init__(self)

        self.reads = reads

    def run(self, fileStore):
        fileStore.logToMaster("FastQCone")
        tempDir = fileStore.getLocalTempDir()

        try:
            os.makedirs(os.path.join(tempDir, 'execution'))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        fpath = fileStore.readGlobalFile(self.reads, userPath=os.path.join(tempDir, os.path.basename(self.reads)))
        cmd = f'zcat "{fpath}" | fastqc stdin:readsone'

        generate_docker_bashscript_file(temp_dir=tempDir, docker_dir=tempDir, globs=[], cmd=cmd, job_name='FastQCone')

        # apiDockerCall() with demux=True returns a tuple of bytes objects (stdout, stderr).
        stdout, stderr = \
            apiDockerCall(self,
                          image='pegi3s/fastqc',
                          working_dir=tempDir,
                          parameters=[os.path.join(tempDir, "FastQCone_script.sh")],
                          entrypoint="/bin/bash",
                          user='root',
                          stderr=True,
                          demux=True,
                          volumes={tempDir: {"bind": tempDir}})

        with open(os.path.join(current_working_dir, 'FastQCone.log'), 'wb') as f:
            if stdout:
                f.write(stdout)
            if stderr:
                f.write(stderr)

        output_file_id = fileStore.writeGlobalFile(os.path.join(tempDir, 'execution', 'readsone_fastqc.html'))
        fastqc_output_path = os.path.join(os.path.abspath(current_working_dir), 'readsone_fastqc.html')
        fileStore.exportFile(output_file_id, f'file://{fastqc_output_path}')

        return {"fastqc_output_path": output_file_id}


class FastQCtwoCls(Job):
    def __init__(self, reads=None, *args, **kwargs):
        super(FastQCtwoCls, self).__init__(*args, **kwargs)
        Job.__init__(self)

        self.reads = reads

    def run(self, fileStore):
        fileStore.logToMaster("FastQCtwo")
        tempDir = fileStore.getLocalTempDir()

        try:
            os.makedirs(os.path.join(tempDir, 'execution'))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        fpath = fileStore.readGlobalFile(self.reads, userPath=os.path.join(tempDir, os.path.basename(self.reads)))
        cmd = f'zcat "{fpath}" | fastqc stdin:readstwo'

        generate_docker_bashscript_file(temp_dir=tempDir, docker_dir=tempDir, globs=[], cmd=cmd, job_name='FastQCtwo')

        # apiDockerCall() with demux=True returns a tuple of bytes objects (stdout, stderr).
        stdout, stderr = \
            apiDockerCall(self,
                          image='pegi3s/fastqc',
                          working_dir=tempDir,
                          parameters=[os.path.join(tempDir, "FastQCtwo_script.sh")],
                          entrypoint="/bin/bash",
                          user='root',
                          stderr=True,
                          demux=True,
                          volumes={tempDir: {"bind": tempDir}})

        with open(os.path.join(current_working_dir, 'FastQCtwo.log'), 'wb') as f:
            if stdout:
                f.write(stdout)
            if stderr:
                f.write(stderr)

        output_file_id = fileStore.writeGlobalFile(os.path.join(tempDir, 'execution', 'readstwo_fastqc.html'))
        fastqc_output_path = os.path.join(os.path.abspath(current_working_dir), 'readstwo_fastqc.html')
        fileStore.exportFile(output_file_id, f'file://{fastqc_output_path}')

        return {"fastqc_output_path": output_file_id}


class SalmonIndexCls(Job):
    def __init__(self, ref_txome=None, index="", *args, **kwargs):
        super(SalmonIndexCls, self).__init__(*args, **kwargs)
        Job.__init__(self)

        self.ref_txome = ref_txome
        self.index = index

    def run(self, fileStore):
        fileStore.logToMaster("SalmonIndex")
        tempDir = fileStore.getLocalTempDir()

        try:
            os.makedirs(os.path.join(tempDir, 'execution'))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        fpath = fileStore.readGlobalFile(self.ref_txome, userPath=os.path.join(tempDir, os.path.basename(self.ref_txome)))

        index_fpath = os.makedirs(os.path.join(tempDir, 'execution', 'index'))


        cmd = f'salmon index -t {fpath} --index "{index_fpath}";  tar -cvzf index.tar.gz {index_fpath}'#; tar -cvzf index.tar.gz "{index_fpath}"'#-i index; tar -cvzf index.tar.gz index'

        generate_docker_bashscript_file(temp_dir=tempDir, docker_dir=tempDir, globs=[], cmd=cmd, job_name='SalmonIndex')

        # apiDockerCall() with demux=True returns a tuple of bytes objects (stdout, stderr).
        stdout, stderr = \
            apiDockerCall(self,
                          image='combinelab/salmon',
                          working_dir=tempDir,
                          parameters=[os.path.join(tempDir, "SalmonIndex_script.sh")],
                          entrypoint="/bin/bash",
                          user='root',
                          stderr=True,
                          demux=True,
                          volumes={tempDir: {"bind": tempDir}})

        with open(os.path.join(current_working_dir, 'SalmonIndex.log'), 'wb') as f:
            if stdout:
                f.write(stdout)
            if stderr:
                f.write(stderr)


        output_file_id = fileStore.writeGlobalFile(os.path.join(tempDir, 'execution', 'index.tar.gz'))
        index_output_path = os.path.join(os.path.abspath(current_working_dir), 'index.tar.gz')
        fileStore.exportFile(output_file_id, f'file://{index_output_path}')

        return {"index": output_file_id}


class SalmonAlignQuantCls(Job):
    def __init__(self, reads1=None, reads2=None, index=None, indexdir="", quantex="", *args, **kwargs):
        super(SalmonAlignQuantCls, self).__init__(*args, **kwargs)
        Job.__init__(self)

        self.reads1 = reads1
        self.reads2 = reads2
        self.index = index
        self.indexdir = indexdir
        self.quantex = quantex

    def run(self, fileStore):
        fileStore.logToMaster("SalmonAlignQuant")
        tempDir = fileStore.getLocalTempDir()

        try:
            os.makedirs(os.path.join(tempDir, 'execution'))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        
        fpath_reads1 = fileStore.readGlobalFile(self.reads1, userPath=os.path.join(tempDir, os.path.basename(self.reads1)))

        fpath_reads2 = fileStore.readGlobalFile(self.reads2, userPath=os.path.join(tempDir, os.path.basename(self.reads2)))

        fpath_index = fileStore.readGlobalFile(self.index, userPath=os.path.join(tempDir, os.path.basename(self.index)))


        cmd = f'tar -xvf "{fpath_index}"; salmon quant -i index -l A -1 "{fpath_reads1}" -2 "{fpath_reads2}" --validateMappings -o quant; tar -cvzf quant.tar.gz quant'

        generate_docker_bashscript_file(temp_dir=tempDir, docker_dir=tempDir, globs=[], cmd=cmd,
                                        job_name='SalmonAlignQuant')

        # apiDockerCall() with demux=True returns a tuple of bytes objects (stdout, stderr).
        stdout, stderr = \
            apiDockerCall(self,
                          image='combinelab/salmon',
                          working_dir=tempDir,
                          parameters=[os.path.join(tempDir, "SalmonAlignQuant_script.sh")],
                          entrypoint="/bin/bash",
                          user='root',
                          stderr=True,
                          demux=True,
                          volumes={tempDir: {"bind": tempDir}})
        with open(os.path.join(current_working_dir, 'SalmonAlignQuant.log'), 'wb') as f:
            if stdout:
                f.write(stdout)
            if stderr:
                f.write(stderr)

        output_file_id = fileStore.writeGlobalFile(os.path.join(tempDir, 'execution', 'quant.tar.gz'))
        quant_output_path = os.path.join(os.path.abspath(current_working_dir), 'quant.tar.gz')
        fileStore.exportFile(output_file_id, f'file://{quant_output_path}')

        return {"quant": output_file_id}


if __name__ == "__main__":
    options = Job.Runner.getDefaultOptions("./toilWorkflowRun")
    options.clean = 'always'

    pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

    with Toil(options) as fileStore:
        # import all files into the jobstore and retrieve their file ID references
        # that way jobs run remotely can fetch from the centralized jobstore without sharing a filesystem
        reads1_file_id = fileStore.importFile(f'file://{os.path.join(pkg_root, "test_data/reads_1.fq.gz")}')
        reads2_file_id = fileStore.importFile(f'file://{os.path.join(pkg_root, "test_data/reads_2.fq.gz")}')
        ref_transcriptome_file_id = fileStore.importFile(f'file://{os.path.join(pkg_root, "test_data/transcriptome.fa")}')
        
        fastqc_job_1 = FastQConeCls(reads=reads1_file_id)  # this is our root job, which runs first
        fastqc_output_report_path = fastqc_job_1.rv("fastqc_res")  # "rv" == return value

        fastqc_job_2 = FastQCtwoCls(reads=reads2_file_id)
        FastQCtwo_fastqc_res = fastqc_job_2.rv("fastqc_res")
        fastqc_job_1.addChild(fastqc_job_2)  # fastqc_job_2 will run after our root job
        
        salmon_index_job = SalmonIndexCls(ref_txome=ref_transcriptome_file_id)
        index_file_id = salmon_index_job.rv("index")
        fastqc_job_1.addChild(salmon_index_job)  # salmon_index_job will run after our root job

        salmon_align_quant_job = SalmonAlignQuantCls(reads1=reads1_file_id, reads2=reads2_file_id, index=index_file_id)
        fastqc_job_1.addFollowOn(salmon_align_quant_job)
        # we don't do anything our results, but we could
        salmon_align_quant_file_id = salmon_align_quant_job.rv("quant")

        fileStore.start(fastqc_job_1)

#/home/hexotical/bioinformatics-workflows/toil/python/index.tar.gz