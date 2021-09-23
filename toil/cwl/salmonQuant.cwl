cwlVersion: v1.2
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: combinelab/salmon

baseCommand: ["salmon", "quant"]

inputs:
  salmonindex:
    type: Directory
    inputBinding:
      position: 1
      prefix: --index
  libtype:
    type: string
    default: A
    inputBinding:
      position: 2
      prefix: --libType
  fq1:
    type: File
    inputBinding:
      position: 3
      prefix: --mates1
  fq2:
    type: File
    inputBinding:
      position: 4
      prefix: --mates2
  quant:
    type: string
    inputBinding:
      position: 5
      prefix: --output
  
outputs:
  quant:
    type: Directory
    outputBinding:
      glob: "$(inputs.quant)"