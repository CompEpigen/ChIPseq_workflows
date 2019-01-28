cwlVersion: v1.0
class: Workflow

requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}

### INPUT PART:
##################################################
inputs:
  sample_id:
    type: string
  fastq1:
    type: 
      type: array
      items: [File, "null"]
  fastq2: 
    type:
      type: array
      items: [File, "null"]
  reference:
    type: File
    secondaryFiles:
      - .fai
      - ^.1.bt2
      - ^.2.bt2
      - ^.3.bt2
      - ^.4.bt2
      - ^.rev.1.bt2
      - ^.rev.2.bt2
  adapter1: 
    type: 
      type: array
      items: [string, "null"]
  adapter2:
    type: 
      type: array
      items: [string, "null"]
  fragment_size:
    type: long?
  is_paired_end:
    type: boolean
  effective_genome_size:
    type: long
  bin_size:
    type: int?
    default: 10
  ignoreForNormalization:
    type:
      type: array
      items: string
    default: ["chrX", "chrY", "chrM"]
        
### WORKFLOW STEPS:
##################################################
steps:
  trim_and_map:
    run: "../workflow_modules/trim_and_map.cwl"
    scatter: [fastq1, fastq2, adapter1, adapter2]
    scatterMethod: 'dotproduct'
    in:
      fastq1:
        source: fastq1
      fastq2: 
        source: fastq2
      reference:
        source: reference
      adapter1: 
        source: adapter1
      adapter2:
        source: adapter2
      is_paired_end:
        source: is_paired_end
    out:
      - pre_trim_fastqc_zip
      - pre_trim_fastqc_html
      - fastq1_trimmed
      - fastq2_trimmed
      - trim_galore_log
      - post_trim_fastqc_html
      - post_trim_fastqc_zip
      - bam
      - bowtie2_log

  merge_duprem_filter:
    run: "../workflow_modules/merge_duprem_filter.cwl"
    in:
      sample_id:
        source: sample_id
      bams: 
        source: trim_and_map/bam
      is_paired_end:
        source: is_paired_end
    out:
      - post_filter_fastqc_zip
      - post_filter_fastqc_html
      - picard_markdup_stdout
      - bam

  tn5_overhang_correction:
    run: "../tools/tn5_overhang_correction.cwl"
    in:
      bam:
        source: merge_duprem_filter/bam
      is_paired_end:
        source: is_paired_end
    out:
      - bam_tn5_corrected

  indexing_shifted_bam:
    run: "../tools/samtools_index_hack.cwl"
    in:
      bam_sorted:
        source: tn5_overhang_correction/bam_tn5_corrected
    out:
       - bam_sorted_indexed
  
  generate_coverage_tracks:
    run: "../tools/deeptools_bamCoverage.cwl"
    in:
      bam:
        source: indexing_shifted_bam/bam_sorted_indexed
      is_paired_end:
        source: is_paired_end
      fragment_size:
        source: fragment_size
      effective_genome_size:
        source: effective_genome_size
      bin_size:
        source: bin_size
      ignoreForNormalization:
        source: ignoreForNormalization
    out:
      - bigwig
  
  chip_qc:
    run: "../workflow_modules/chip_qc.cwl"
    in:
      sample_id:
        source: sample_id
      bam:
        source: merge_duprem_filter/bam
      is_paired_end:
        source: is_paired_end
      fragment_size:
        source: fragment_size
    out:
      - qc_plot_coverage_plot
      - qc_plot_coverage_tsv
      - qc_plot_fingerprint_plot
      - qc_plot_fingerprint_tsv
      - qc_plot_fingerprint_stderr
      - qc_crosscorr_summary
      - qc_crosscorr_plot
      - qc_phantompeakqualtools_stderr
      - qc_phantompeakqualtools_stdout

  create_summary_qc_report:
    doc: |
      multiqc summarizes the qc results from fastqc 
      and other tools
    run: "../tools/multiqc_hack.cwl"
    in:
      qc_files_array_of_array:
        source:
          - trim_and_map/pre_trim_fastqc_zip
          - trim_and_map/pre_trim_fastqc_html
          - trim_and_map/post_trim_fastqc_html
          - trim_and_map/post_trim_fastqc_zip
          - trim_and_map/trim_galore_log
        linkMerge: merge_flattened
      qc_files_array:
        source:
          - trim_and_map/bowtie2_log
          - merge_duprem_filter/post_filter_fastqc_zip
          - merge_duprem_filter/post_filter_fastqc_html
          - chip_qc/qc_plot_coverage_tsv
          - chip_qc/qc_plot_coverage_plot
          - chip_qc/qc_plot_fingerprint_tsv
          - chip_qc/qc_plot_fingerprint_plot
          - chip_qc/qc_crosscorr_summary
          - chip_qc/qc_crosscorr_plot
          - chip_qc/qc_phantompeakqualtools_stdout
          - chip_qc/qc_crosscorr_summary
          - merge_duprem_filter/picard_markdup_stdout
        linkMerge: merge_flattened
      report_name:
        source: sample_id
    out:
      - multiqc_zip
      - multiqc_html


### OUTPUTS:
##################################################
outputs:
  pre_trim_fastqc_zip:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_and_map/pre_trim_fastqc_zip
  pre_trim_fastqc_html:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_and_map/pre_trim_fastqc_html
  trim_galore_log:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_and_map/pre_trim_fastqc_zip
  post_trim_fastqc_html:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_and_map/post_trim_fastqc_html
  post_trim_fastqc_zip:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_and_map/post_trim_fastqc_zip
  bowtie2_log:
    type:
      type: array
      items: File
    outputSource: trim_and_map/bowtie2_log

  post_filter_fastqc_zip:
    type:
      type: array
      items: File
    outputSource: merge_duprem_filter/post_filter_fastqc_zip
  post_filter_fastqc_html:
    type:
      type: array
      items: File
    outputSource: merge_duprem_filter/post_filter_fastqc_html
  bam:
    type: File
    secondaryFiles: .bai
    outputSource: merge_duprem_filter/bam
  picard_markdup_stdout:
    type: File
    outputSource: merge_duprem_filter/picard_markdup_stdout

  bam_tn5_corrected:
    type: File
    outputSource: indexing_shifted_bam/bam_sorted_indexed

  bigwig:
    type: File
    outputSource: generate_coverage_tracks/bigwig

  qc_plot_coverage_plot:
    type: File
    outputSource: chip_qc/qc_plot_coverage_plot
  qc_plot_coverage_tsv:
    type: File
    outputSource: chip_qc/qc_plot_coverage_tsv
  qc_plot_fingerprint_plot:
    type: File?
    outputSource: chip_qc/qc_plot_fingerprint_plot
  qc_plot_fingerprint_tsv:
    type: File?
    outputSource: chip_qc/qc_plot_fingerprint_tsv
  qc_plot_fingerprint_stderr:
    type: File
    outputSource: chip_qc/qc_plot_fingerprint_stderr
  qc_crosscorr_summary:
    type: File?
    outputSource: chip_qc/qc_crosscorr_summary
  qc_crosscorr_plot:
    type: File?
    outputSource: chip_qc/qc_crosscorr_plot
  qc_phantompeakqualtools_stderr:
    type: File?
    outputSource: chip_qc/qc_phantompeakqualtools_stderr

  multiqc_zip:
    type: File
    outputSource: create_summary_qc_report/multiqc_zip
  multiqc_html:
    type: File
    outputSource: create_summary_qc_report/multiqc_html