cwlVersion: v1.0
class: Workflow

requirements:
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  sample_id:
    type: string
  bam:
    type: File
    secondaryFiles: .bai
  fragment_size:
    type: int?
  is_paired_end:
    type: boolean

steps:
  qc_plot_coverage:
    doc: |
      deeptools plotCoverage - plots how many times a certain fraction of the 
      genome was covered (consideres the complete fragment between a reads pair).
    run: "../tools/deeptools_plotCoverage.cwl"
    in:
      bam:
        source: bam
      sample_id:
        source: sample_id
      is_paired_end:
        source: is_paired_end
      fragment_size:
        source: fragment_size
    out:
      - qc_plot_coverage_plot  
      - qc_plot_coverage_tsv

  qc_plot_fingerprint:
    doc: |
      Applies deeptools plotFingerprint to generate meaningful plots for comparing IP and
      control samples in terms of quality control. The main question which can be answered is:
      Did the antibody lead to enough enrichment that the IP can be distinguished from the control? 
    run: "../tools/deeptools_plotFingerprint.cwl"
    in:
      bam:
        source: bam
      sample_id:
        source: sample_id
      is_paired_end:
        source: is_paired_end
      fragment_size:
        source: fragment_size
    out:
      - qc_plot_fingerprint_plot  
      - qc_plot_fingerprint_tsv
      - qc_plot_fingerprint_stderr

  qc_phantompeakqualtools:
    run: "../tools/phantompeakqualtools.cwl"
    in:
      bam:
        source: bam
    out:
      - qc_crosscorr_summary  
      - qc_crosscorr_plot
      - qc_phantompeakqualtools_stderr
      - qc_phantompeakqualtools_stdout
      
outputs:
  qc_plot_coverage_plot:
    type: File
    outputSource: qc_plot_coverage/qc_plot_coverage_plot
  qc_plot_coverage_tsv:
    type: File
    outputSource: qc_plot_coverage/qc_plot_coverage_tsv
  qc_plot_fingerprint_plot:
    type: File?
    outputSource: qc_plot_fingerprint/qc_plot_fingerprint_plot
  qc_plot_fingerprint_tsv:
    type: File?
    outputSource: qc_plot_fingerprint/qc_plot_fingerprint_tsv
  qc_plot_fingerprint_stderr:
    type: File
    outputSource: qc_plot_fingerprint/qc_plot_fingerprint_stderr
  qc_crosscorr_summary:
    type: File?
    outputSource: qc_phantompeakqualtools/qc_crosscorr_summary
  qc_crosscorr_plot:
    type: File?
    outputSource: qc_phantompeakqualtools/qc_crosscorr_plot
  qc_phantompeakqualtools_stderr:
    type: File?
    outputSource: qc_phantompeakqualtools/qc_phantompeakqualtools_stderr
  qc_phantompeakqualtools_stdout:
    type: File?
    outputSource: qc_phantompeakqualtools/qc_phantompeakqualtools_stdout
