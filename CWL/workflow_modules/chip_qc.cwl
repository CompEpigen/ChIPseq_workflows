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
  qc_phantompeakqualtools:
    run: "../tools/phantompeakqualtools.cwl"
    in:
      bam:
        source: bam
    out:
      - qc_crosscorr_summary  
      - qc_crosscorr_plot
      - qc_crosscorr_fragment_size
      - qc_phantompeakqualtools_stderr
      - qc_phantompeakqualtools_stdout

  fragment_size_decision_maker:
    doc: |
      If no user-defined fragment size was set,
      the fragment size infered from cross-correlation analysis 
      will be used.
    run: "../tools/frag_size_decision_maker.cwl"
    in:
      user_def_fragment_size:
        source: fragment_size
      cc_fragment_size:
        source: qc_phantompeakqualtools/qc_crosscorr_fragment_size
    out:
      - fragment_size

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
        source: fragment_size_decision_maker/fragment_size
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
        source: fragment_size_decision_maker/fragment_size
    out:
      - qc_plot_fingerprint_plot  
      - qc_plot_fingerprint_tsv
      - qc_plot_fingerprint_stderr
      
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
  fragment_size:
    type: int?
    outputSource: fragment_size_decision_maker/fragment_size
