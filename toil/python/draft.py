import os
import errno
import logging

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


# class FastQCtwoCls(Job):
#     def __init__(self, reads=None, *args, **kwargs):
#         super(FastQCtwoCls, self).__init__(*args, **kwargs)
#         Job.__init__(self)
#
#         self.id_reads = reads
#
#     def run(self, fileStore):
#         fileStore.logToMaster("FastQCtwo")
#         tempDir = fileStore.getLocalTempDir()
#
#         _toil_wdl_internal__stdout_file = os.path.join(tempDir, 'stdout')
#         _toil_wdl_internal__stderr_file = os.path.join(tempDir, 'stderr')
#
#         try:
#             os.makedirs(os.path.join(tempDir, 'execution'))
#         except OSError as e:
#             if e.errno != errno.EEXIST:
#                 raise
#
#         reads = process_and_read_file(abspath_file(self.id_reads, current_working_dir), tempDir,
#                                       fileStore, docker=True)
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command3 = r'''
#              zcat "'''
#         except:
#             command3 = ''
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command4 = str(
#                 reads if not isinstance(reads, WDLFile) else process_and_read_file(reads, tempDir, fileStore)).strip(
#                 "\n")
#         except:
#             command4 = ''
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command5 = r'''" | fastqc stdin:readstwo
#           '''
#         except:
#             command5 = ''
#
#         cmd = command3 + command4 + command5
#         cmd = textwrap.dedent(cmd.strip("\n"))
#         generate_docker_bashscript_file(temp_dir=tempDir, docker_dir=tempDir, globs=[], cmd=cmd, job_name='FastQCtwo')
#
#         # apiDockerCall() with demux=True returns a tuple of bytes objects (stdout, stderr).
#         _toil_wdl_internal__stdout, _toil_wdl_internal__stderr = \
#             apiDockerCall(self,
#                           image='pegi3s/fastqc',
#                           working_dir=tempDir,
#                           parameters=[os.path.join(tempDir, "FastQCtwo_script.sh")],
#                           entrypoint="/bin/bash",
#                           user='root',
#                           stderr=True,
#                           demux=True,
#                           volumes={tempDir: {"bind": tempDir}})
#         with open(os.path.join(current_working_dir, 'FastQCtwo.log'), 'wb') as f:
#             if _toil_wdl_internal__stdout:
#                 f.write(_toil_wdl_internal__stdout)
#             if _toil_wdl_internal__stderr:
#                 f.write(_toil_wdl_internal__stderr)
#
#         _toil_wdl_internal__stdout_file = generate_stdout_file(_toil_wdl_internal__stdout,
#                                                                tempDir,
#                                                                fileStore=fileStore)
#         _toil_wdl_internal__stderr_file = generate_stdout_file(_toil_wdl_internal__stderr,
#                                                                tempDir,
#                                                                fileStore=fileStore,
#                                                                stderr=True)
#
#         fastqc_res = WDLFileType().create(
#             'readstwo_fastqc.html', output=True)
#         fastqc_res = process_outfile(fastqc_res, fileStore, tempDir, '/home/quokka/git/bioinformatics-workflows')
#
#         rvDict = {"fastqc_res": fastqc_res}
#         return rvDict
#
#
# class SalmonIndexCls(Job):
#     def __init__(self, ref_txome=None, index="", *args, **kwargs):
#         super(SalmonIndexCls, self).__init__(*args, **kwargs)
#         Job.__init__(self)
#
#         self.id_ref_txome = ref_txome
#         self.id_index = WDLStringType().create(
#             'index')
#
#     def run(self, fileStore):
#         fileStore.logToMaster("SalmonIndex")
#         tempDir = fileStore.getLocalTempDir()
#
#         _toil_wdl_internal__stdout_file = os.path.join(tempDir, 'stdout')
#         _toil_wdl_internal__stderr_file = os.path.join(tempDir, 'stderr')
#
#         try:
#             os.makedirs(os.path.join(tempDir, 'execution'))
#         except OSError as e:
#             if e.errno != errno.EEXIST:
#                 raise
#
#         ref_txome = process_and_read_file(abspath_file(self.id_ref_txome, current_working_dir),
#                                           tempDir, fileStore, docker=True)
#         index = self.id_index
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command6 = r'''
#              salmon index -t "'''
#         except:
#             command6 = ''
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command7 = str(
#                 ref_txome if not isinstance(ref_txome, WDLFile) else process_and_read_file(ref_txome, tempDir,
#                                                                                            fileStore)).strip("\n")
#         except:
#             command7 = ''
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command8 = r'''" -i '''
#         except:
#             command8 = ''
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command9 = str(
#                 index if not isinstance(index, WDLFile) else process_and_read_file(index, tempDir, fileStore)).strip(
#                 "\n")
#         except:
#             command9 = ''
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command10 = r'''; tar -cvzf index.tar.gz '''
#         except:
#             command10 = ''
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command11 = str(
#                 index if not isinstance(index, WDLFile) else process_and_read_file(index, tempDir, fileStore)).strip(
#                 "\n")
#         except:
#             command11 = ''
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command12 = r'''
#           '''
#         except:
#             command12 = ''
#
#         cmd = command6 + command7 + command8 + command9 + command10 + command11 + command12
#         cmd = textwrap.dedent(cmd.strip("\n"))
#         generate_docker_bashscript_file(temp_dir=tempDir, docker_dir=tempDir, globs=[], cmd=cmd, job_name='SalmonIndex')
#
#         # apiDockerCall() with demux=True returns a tuple of bytes objects (stdout, stderr).
#         _toil_wdl_internal__stdout, _toil_wdl_internal__stderr = \
#             apiDockerCall(self,
#                           image='combinelab/salmon',
#                           working_dir=tempDir,
#                           parameters=[os.path.join(tempDir, "SalmonIndex_script.sh")],
#                           entrypoint="/bin/bash",
#                           user='root',
#                           stderr=True,
#                           demux=True,
#                           volumes={tempDir: {"bind": tempDir}})
#         with open(os.path.join(current_working_dir, 'SalmonIndex.log'), 'wb') as f:
#             if _toil_wdl_internal__stdout:
#                 f.write(_toil_wdl_internal__stdout)
#             if _toil_wdl_internal__stderr:
#                 f.write(_toil_wdl_internal__stderr)
#
#         _toil_wdl_internal__stdout_file = generate_stdout_file(_toil_wdl_internal__stdout,
#                                                                tempDir,
#                                                                fileStore=fileStore)
#         _toil_wdl_internal__stderr_file = generate_stdout_file(_toil_wdl_internal__stderr,
#                                                                tempDir,
#                                                                fileStore=fileStore,
#                                                                stderr=True)
#
#         index = WDLFileType().create(
#             'index.tar.gz', output=True)
#         index = process_outfile(index, fileStore, tempDir, '/home/quokka/git/bioinformatics-workflows')
#
#         rvDict = {"index": index}
#         return rvDict
#
#
# class SalmonAlignQuantCls(Job):
#     def __init__(self, reads1=None, reads2=None, index=None, indexdir="", quantex="", *args, **kwargs):
#         super(SalmonAlignQuantCls, self).__init__(*args, **kwargs)
#         Job.__init__(self)
#
#         self.id_reads1 = reads1
#         self.id_reads2 = reads2
#         self.id_index = index
#         self.id_indexdir = WDLStringType().create(
#             'index')
#         self.id_quantex = WDLStringType().create(
#             'quant')
#
#     def run(self, fileStore):
#         fileStore.logToMaster("SalmonAlignQuant")
#         tempDir = fileStore.getLocalTempDir()
#
#         _toil_wdl_internal__stdout_file = os.path.join(tempDir, 'stdout')
#         _toil_wdl_internal__stderr_file = os.path.join(tempDir, 'stderr')
#
#         try:
#             os.makedirs(os.path.join(tempDir, 'execution'))
#         except OSError as e:
#             if e.errno != errno.EEXIST:
#                 raise
#
#         reads1 = process_and_read_file(abspath_file(self.id_reads1, current_working_dir), tempDir,
#                                        fileStore, docker=True)
#         reads2 = process_and_read_file(abspath_file(self.id_reads2, current_working_dir), tempDir,
#                                        fileStore, docker=True)
#         index = process_and_read_file(abspath_file(self.id_index, current_working_dir), tempDir,
#                                       fileStore, docker=True)
#         indexdir = self.id_indexdir
#         quantex = self.id_quantex
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command13 = r'''
#              tar -xzf '''
#         except:
#             command13 = ''
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command14 = str(
#                 index if not isinstance(index, WDLFile) else process_and_read_file(index, tempDir, fileStore)).strip(
#                 "\n")
#         except:
#             command14 = ''
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command15 = r''';
#              salmon quant -i "'''
#         except:
#             command15 = ''
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command16 = str(indexdir if not isinstance(indexdir, WDLFile) else process_and_read_file(indexdir, tempDir,
#                                                                                                      fileStore)).strip(
#                 "\n")
#         except:
#             command16 = ''
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command17 = r'''" -l A -1 "'''
#         except:
#             command17 = ''
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command18 = str(
#                 reads1 if not isinstance(reads1, WDLFile) else process_and_read_file(reads1, tempDir, fileStore)).strip(
#                 "\n")
#         except:
#             command18 = ''
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command19 = r'''" -2 "'''
#         except:
#             command19 = ''
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command20 = str(
#                 reads2 if not isinstance(reads2, WDLFile) else process_and_read_file(reads2, tempDir, fileStore)).strip(
#                 "\n")
#         except:
#             command20 = ''
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command21 = r'''" --validateMappings -o quant;
#              tar -cvzf quant.tar.gz '''
#         except:
#             command21 = ''
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command22 = str(quantex if not isinstance(quantex, WDLFile) else process_and_read_file(quantex, tempDir,
#                                                                                                    fileStore)).strip(
#                 "\n")
#         except:
#             command22 = ''
#
#         try:
#             # Intended to deal with "optional" inputs that may not exist
#             # TODO: handle this better
#             command23 = r'''
#           '''
#         except:
#             command23 = ''
#
#         cmd = command13 + command14 + command15 + command16 + command17 + command18 + command19 + command20 + command21 + command22 + command23
#         cmd = textwrap.dedent(cmd.strip("\n"))
#         generate_docker_bashscript_file(temp_dir=tempDir, docker_dir=tempDir, globs=[], cmd=cmd,
#                                         job_name='SalmonAlignQuant')
#
#         # apiDockerCall() with demux=True returns a tuple of bytes objects (stdout, stderr).
#         _toil_wdl_internal__stdout, _toil_wdl_internal__stderr = \
#             apiDockerCall(self,
#                           image='combinelab/salmon',
#                           working_dir=tempDir,
#                           parameters=[os.path.join(tempDir, "SalmonAlignQuant_script.sh")],
#                           entrypoint="/bin/bash",
#                           user='root',
#                           stderr=True,
#                           demux=True,
#                           volumes={tempDir: {"bind": tempDir}})
#         with open(os.path.join(current_working_dir, 'SalmonAlignQuant.log'), 'wb') as f:
#             if _toil_wdl_internal__stdout:
#                 f.write(_toil_wdl_internal__stdout)
#             if _toil_wdl_internal__stderr:
#                 f.write(_toil_wdl_internal__stderr)
#
#         _toil_wdl_internal__stdout_file = generate_stdout_file(_toil_wdl_internal__stdout,
#                                                                tempDir,
#                                                                fileStore=fileStore)
#         _toil_wdl_internal__stderr_file = generate_stdout_file(_toil_wdl_internal__stderr,
#                                                                tempDir,
#                                                                fileStore=fileStore,
#                                                                stderr=True)
#
#         quant = WDLFileType().create(
#             'quant.tar.gz', output=True)
#         quant = process_outfile(quant, fileStore, tempDir, '/home/quokka/git/bioinformatics-workflows')
#
#         rvDict = {"quant": quant}
#         return rvDict


if __name__ == "__main__":
    options = Job.Runner.getDefaultOptions("./toilWorkflowRun")
    options.clean = 'always'

    pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

    with Toil(options) as fileStore:
        reads1 = fileStore.importFile(os.path.join(pkg_root, "test_data/reads_1.fq.gz"))
        reads2 = fileStore.importFile(os.path.join(pkg_root, "test_data/reads_2.fq.gz"))
        ref_txome = fileStore.importFile(os.path.join(pkg_root, "test_data/transcriptome.fa"))

        FastQCone = FastQConeCls(reads=reads1)
        fastqc_output_report_path = FastQCone.rv("fastqc_res")

        # FastQCtwo = job0.addChild(FastQCtwoCls(reads=reads2))
        # FastQCtwo_fastqc_res = FastQCtwo.rv("fastqc_res")
        #
        # SalmonIndex = job0.addChild(SalmonIndexCls(ref_txome=ref_txome))
        # SalmonIndex_index = SalmonIndex.rv("index")
        #
        # SalmonAlignQuant = job0.addChild(SalmonAlignQuantCls(reads1=reads1, reads2=reads2, index=(SalmonIndex_index)))
        # SalmonAlignQuant_quant = SalmonAlignQuant.rv("quant")

        fileStore.start(FastQCone)
