cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 4
    ramMin: 20000
    #tmpdirMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/samtools:1.7
  
baseCommand: ["bash", "-c"]
arguments:
  - valueFrom: |
      ${
        var cmd_line = "";
        
        if ( inputs.is_paired_end ){ // for paired end data
                                      // unpaired will be removed
        
          ////// shift + strand reads
          cmd_line += "samtools view -h -f 3 -F 16 " + inputs.bam.path; // only properly paired reads with the
                                                // first read on the + strand are output;
                                                // the header is included
          
          cmd_line += " | awk  \'BEGIN {OFS=\"\\t\"} " +
                      "{ " + 
                      "if ( $1 ~ /^@/) { print }" + //header lines are printed unmodified
                      "else if ($9>=38) { " + 
                      "$4=$4+4; $8=$8-5; $9=$9-9; $11=\"*\"; if ($8>0){ print } else {$8=1; print}}" + // read start positions are shifted
                      "else if ($9<=-38) { " + 
                      "$4=$4+4; $8=$8-5; $9=$9+9; $11=\"*\"; if ($8>0){ print } else {$8=1; print}}" + // read start positions are shifted
                      "}\' > correcting.sam";
          ///// shift - strand reads      
          cmd_line += " ; samtools view -f 19 " + inputs.bam.path; // only properly paired reads with the
                                                // first read on the - strand are output;
                                                // the header is excluded
          
          cmd_line += " | awk  \'BEGIN {OFS=\"\\t\"} " +
                      "{ " + 
                      "if ($9>=38) { " + 
                      "$4=$4-5; $8=$8+4; $9=$9-9; $11=\"*\";  if ($4>0){ print } else {$4=1; print}}" + // read start positions are shifted
                      "else if ($9<=-38) { " + 
                      "$4=$4-5; $8=$8+4; $9=$9+9; $11=\"*\"; if ($4>0){ print } else {$4=1; print}}" + // read start positions are shifted
                      "}\' >> correcting.sam";
                      
          
        }
        else { // for single end data
        
          ////// shift + strand reads
          cmd_line += "samtools view -h -F 16 " + inputs.bam.path; // paired end as well as
                                                // - strand reads are excluded
                                                // the header is included
          
          cmd_line += " | awk  \'BEGIN {OFS=\"\\t\"} " +
                      "{ " + 
                      "if ( $1 ~ /^@/ ) { print }" + //header lines are printed unmodified
                      "else { $4=$4+4; $7=\"*\"; $8=0; $9=0; $11=\"*\"; print}" + // read start positions are shifted
                      "}\' > correcting.sam";
          ///// shift - strand reads      
          cmd_line += " ; samtools view -f 16 " + inputs.bam.path; // paired end as well as
                                                // + strand reads are excluded
                                                // the header is included
          
          cmd_line += " | awk  \'BEGIN {OFS=\"\\t\"} " +
                      "{ " + 
                      "$4=$4-5; $7=\"*\"; $8=0; $9=0; $11=\"*\"; print " + // read start positions are shifted
                      "}\' >> correcting.sam";
        }
        
        cmd_line += " ; samtools sort -@ " + runtime.cores + " -O bam -T sorting.bam -o " + inputs.bam.nameroot + "_" + inputs.out_suffix + ".bam correcting.sam";
        
        return cmd_line;
      }
        
inputs:
  is_paired_end:
    type: boolean
  bam:
    type: File
  out_suffix:
    type: string
    default: "tn5correct"
 
outputs:
  bam_tn5_corrected:
    type: File
    outputBinding:
      glob: $(inputs.bam.nameroot + "_" + inputs.out_suffix + ".bam")