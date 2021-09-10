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

        return {"fastqc_output_path": fastqc_output_path}


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

        return {"fastqc_output_path": fastqc_output_path}


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
        cmd = f'salmon index -t "{fpath}" -i index; tar -cvzf index.tar.gz index'

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

        
        return {"index_output_path" : index_output_path}


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
        
        fpath_reads1 = fileStore.readGlobalFile(self.reads, userPath=os.path.join(tempDir, os.path.basename(self.reads)))

        fpath_reads2 = fileStore.readGlobalFile(self.reads, userPath=os.path.join(tempDir, os.path.basename(self.reads)))

        fpath_index = fileStore.readGlobalFile(self.reads, userPath=os.path.join(tempDir, os.path.basename(self.reads)))


        cmd = f'tar -xvf "{fpath_index}"; salmon quant -i index -l A -1 "{fpath_reads1}" -2 "{fpath_reads2}" --validateMappins -o quant; tar -cvzf quant.tar.gz quant'

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

        return {"quant": quant_output_path}


if __name__ == "__main__":
    options = Job.Runner.getDefaultOptions("./toilWorkflowRun")
    options.clean = 'always'

    pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

    with Toil(options) as fileStore:
        reads1 = fileStore.importFile('file://' + os.path.join(pkg_root, "test_data/reads_1.fq.gz"))
        reads2 = fileStore.importFile('file://' + os.path.join(pkg_root, "test_data/reads_2.fq.gz"))
        ref_txome = fileStore.importFile('file://' + os.path.join(pkg_root, "test_data/transcriptome.fa"))
        
        FastQCone = FastQConeCls(reads=reads1)
        fastqc_output_report_path = FastQCone.rv("fastqc_res")

        FastQCtwo = FastQCtwoCls(reads=reads2)
        FastQCtwo_fastqc_res = FastQCtwo.rv("fastqc_res")
        FastQCone.addChild(FastQCtwo)
        
        SalmonIndex = FastQCone.addChild(SalmonIndexCls(ref_txome=ref_txome))
        SalmonIndex_index = SalmonIndex.rv("index")
        
        SalmonAlignQuant = SalmonIndex.addChild(SalmonAlignQuantCls(reads1=reads1, reads2=reads2, index=(SalmonIndex_index)))
        SalmonAlignQuant_quant = SalmonAlignQuant.rv("quant")

        fileStore.start(FastQCone)
