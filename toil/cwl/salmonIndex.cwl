cwlVersion: v1.2
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: combinelab/salmon
    
baseCommand: ["salmon", "index"]


inputs:
  transcripts:
    type: File
    inputBinding:
      prefix: -t
  index:
    type: string
    inputBinding:
      prefix: -i

outputs:
  salmon_index:
    type: Directory
    outputBinding:
      glob: "$(inputs.index)"