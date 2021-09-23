cwlVersion: v1.2
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: pegi3s/fastqc

baseCommand: "fastqc"
arguments: 
  - valueFrom: $(runtime.outdir)
    prefix: "-o"
  - valueFrom: "--noextract"


inputs:
  reads_file:
    type:
      - File
    inputBinding:
      position: 50
    doc: |

outputs:

  zipped_file:
    type:
      - File
    outputBinding:
      glob: '*.zip'
  html_file:
    type:
      - File
    outputBinding:
      glob: '*.html'
