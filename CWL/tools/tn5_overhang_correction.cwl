cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing: 
      - entryname: tn5_overhang_correction.sh
        entry: |
          BAM="$1"
          OUTPUT_SUFFIX="$2"
          IS_PAIRED_END="$3"
          if [[ "$IS_PAIRED_END" == TRUE ]]
          then
            samtools view -h -f 3 -F 16 "$BAM" | awk -f shifting_paired_end_plus_reads.awk > correcting.sam
            samtools view -f 19 "$BAM" | awk -f shifting_paired_end_minus_reads.awk >> correcting.sam
          else
            samtools view -h -F 16 "$BAM" | awk -f shifting_single_end_plus_reads.awk > correcting.sam
            samtools view -f 16 "$BAM" | awk -f shifting_single_end_minus_reads.awk >> correcting.sam
          fi
          OUTFILE="\${BAM##*/}"
          OUTFILE="\${OUTFILE%.bam}_\${OUTPUT_SUFFIX}.bam"
          samtools sort -@ $(runtime.cores) -O bam -T sorting.bam -o "$OUTFILE" correcting.sam
      - entryname: shifting_paired_end_plus_reads.awk
        entry: |
          BEGIN {OFS="\t"}
          {
            if ( $1 ~ /^@/) { print }
            else if ($9>=38) {
              $4=$4+4; $8=$8-5; $9=$9-9; $11="*";
              if ($8>0){ print } else {$8=1; print}
            }
            else if ($9<=-38) {
              $4=$4+4; $8=$8-5; $9=$9+9; $11="*";
              if ($8>0){ print } else {$8=1; print}
            }
          }
      - entryname: shifting_paired_end_minus_reads.awk
        entry: |
          BEGIN {OFS="\t"}
          {
            if ($9>=38) {
              $4=$4-5; $8=$8+4; $9=$9-9; $11="*";  
              if ($4>0){ print } else {$4=1; print}
            }
            else if ($9<=-38) {
              $4=$4-5; $8=$8+4; $9=$9+9; $11="*"; 
              if ($4>0){ print } else {$4=1; print}
            }
          }
      - entryname: shifting_single_end_plus_reads.awk
        entry: |
          BEGIN {OFS="\t"}
          {
            if ( $1 ~ /^@/) { print }
            else {
              $4=$4+4; $7="*"; $8=0; $9=0; $11="*"; print
            }
          }
      - entryname: shifting_single_end_minus_reads.awk
        entry: |
          BEGIN {OFS="\t"}
          {
            $4=$4-5; $7="*"; $8=0; $9=0; $11="*"; print
          }

hints:
  ResourceRequirement:
    coresMin: 4
    ramMin: 20000
    #tmpdirMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/samtools:1.7
  
baseCommand: ["bash", "tn5_overhang_correction.sh"]

inputs:
  is_paired_end:
    type: boolean
    inputBinding:
      prefix: "TRUE"
      position: 3
  bam:
    type: File
    inputBinding:
      position: 1
  out_suffix:
    type: string
    default: "tn5correct"
    inputBinding:
      position: 2
 
outputs:
  bam_tn5_corrected:
    type: File
    outputBinding:
      glob: $("*_" + inputs.out_suffix + ".bam")