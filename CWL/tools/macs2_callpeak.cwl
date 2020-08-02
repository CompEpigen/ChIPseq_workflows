cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
  DockerRequirement:
    dockerPull: genomicpariscentre/macs2:2.1.0.20140616
 
### BASE COMMAND AND ARGUMENTS:
##################################################
baseCommand: ["macs2", "callpeak"]
arguments:  
  - valueFrom: |
      ${
        if ( inputs.is_paired_end ){
           return "BAMPE";
        }
        else {
          return "BAM";
        }
      }
    prefix: --format
    position: 1
  - valueFrom: "--nomodel"
    position: 2
  - valueFrom: "all"
    prefix: "--keep-dup"
    position: 2
  - valueFrom: $(inputs.bam.nameroot + ".macs2")
    prefix: "--name"
    position: 100

inputs:
  bam:
    type: File
    inputBinding:
        position: 101
        prefix: "--treatment"
  genome_size:
    type: long
    inputBinding:
        position: 3
        prefix: "--gsize"
  broad:
    type: boolean
    inputBinding:
        position: 3
        prefix: "--broad"
  qvalue:
    type: float
    inputBinding:
        position: 3
        prefix: "--qvalue"
  is_paired_end:
    type: boolean
    default: true
 
outputs:
  peaks_bed:    
    type: 
      type: array
      items: File
    outputBinding:
      glob: "*Peak"
  peaks_xls:
    type: File
    outputBinding:
      glob: "*_peaks.xls"