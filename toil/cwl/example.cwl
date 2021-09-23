cwlVersion: v1.2
class: Workflow
label: Example CWL workflow using Toil

inputs:
  fqOne: File
  fqTwo: File
  transcriptome: File
  index: string
  quant: string

steps:
    fastqc:
        run: fastqc.cwl
        in:
            reads_file: fqOne
        out: [html_file]
    fastqcTwo:
        run: fastqc.cwl
        in:
            reads_file: fqTwo
        out: [html_file]
    salmonIndex:
        run: salmonIndex.cwl
        in:
            index: index
            transcripts: transcriptome
        out: [salmon_index]
    salmonQuant:
        run: salmonQuant.cwl
        in:
            salmonindex: salmonIndex/salmon_index
            fq1: fqOne
            fq2: fqTwo
            quant: quant
        out: [quant]

outputs:
    qc_html:
        type: File
        outputSource: fastqc/html_file
    qc_html_two:
        type: File
        outputSource: fastqcTwo/html_file
    index:
        type: Directory
        outputSource: salmonIndex/salmon_index
    quant:
        type: Directory
        outputSource: salmonQuant/quant