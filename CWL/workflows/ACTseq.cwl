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
    doc: |
      Sample ID used for naming the output files.
    type: string
  fastq1:
    doc: |
      List of fastq files containing the first mate of raw reads.
      Muliple files are provided if multiplexing of the same library has been done
      on multiple lanes. The reads comming from different fastq files are pooled
      after alignment. Also see parameter "fastq2".
    type: 
      type: array
      items: [File]
  fastq2: 
    doc: |
      List of fastq files containing the second mate of raw reads in case of paired end
      (also see parameter "fastq1").
      Important: this list has to be of same length as parameter "fastq1" no matter if paired or single end is used.
      In case of single end data specify "null" for every entry of fastq1.
    type:
      type: array
      items: [File, "null"]
  is_paired_end:
    doc: |
      If paired end data is used set to true, else set to false.
    type: boolean
  adapter1: 
    doc: |
      Adapter sequence for first reads.
      If not specified (set to "null"), trim_galore will try to autodetect whether ...
      - Illumina universal adapter (AGATCGGAAGAGC)
      - Nextera adapter (CTGTCTCTTATA)
      - Illumina Small RNA 3-prime Adapter (TGGAATTCTCGG)
      ... was used.
      You can directly choose one of the above configurations
      by setting the string to "illumina", "nextera", or "small_rna".
      Or you specify the adaptor string manually (e.g. "AGATCGGAAGAGC").
    type: string?
  adapter2: 
    doc: |
      Adapter sequence for second reads (only relevant for paired end data).
      If it is not specified (set to "null"), trim_galore will try to autodetect whether ...
      - Illumina universal adapter (AGATCGGAAGAGC)
      - Nextera adapter (CTGTCTCTTATA)
      - Illumina Small RNA 3-prime Adapter (TGGAATTCTCGG)
      ... was used.
      You can directly choose one of the above configurations
      by setting the string to "illumina", "nextera", or "small_rna".
      Or you specify the adaptor string manually (e.g. "AGATCGGAAGAGC").
    type: string?
  genome:
    doc: |
      Path to reference genome in fasta format.
      Bowtie2 index files (".1.bt2", ".2.bt2", ...) as well as a samtools index (".fai")
      has to be located in the same directory.
      All of these files can be downloaded for the most common genome builds at 
      https://support.illumina.com/sequencing/sequencing_software/igenome.html.
      Alternatively, you can use "bowtie2-build" or "samtools index" to create them yourself.
    type: File
    secondaryFiles:
      - .fai
      - ^.1.bt2
      - ^.2.bt2
      - ^.3.bt2
      - ^.4.bt2
      - ^.rev.1.bt2
      - ^.rev.2.bt2
  fragment_size:
    doc: |
      Mean library fragment size, used to reconstruct entire 
      fragments from single end reads. Not relevant in case of paired end data.
    type: int?
  effective_genome_size:
    doc: |
      The effectively mappable genome size, please see: 
      https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
    type: long
  bin_size:
    doc: |
      Bin size used for generation of coverage tracks.
      The larger the bin size the smaller are the coverage tracks, however,
      the less precise is the signal. For single bp resolution set to 1.
    type: int
    default: 10
  ignoreForNormalization:
    doc: |
      List of space-delimited chromosome names that shall be ignored
      when calculating the scaling factor. 
    type: string?
    default: "chrX chrY chrM"
        
### WORKFLOW STEPS:
##################################################
steps:
  trim_and_map:
    run: "../workflow_modules/trim_and_map.cwl"
    scatter: [fastq1, fastq2]
    scatterMethod: 'dotproduct'
    in:
      fastq1:
        source: fastq1
      fastq2: 
        source: fastq2
      genome:
        source: genome
      adapter1: 
        source: adapter1
      adapter2:
        source: adapter2
      is_paired_end:
        source: is_paired_end
    out:
      - raw_fastqc_zip
      - raw_fastqc_html
      - fastq1_trimmed
      - fastq2_trimmed
      - trim_galore_log
      - trimmed_fastqc_html
      - trimmed_fastqc_zip
      - bam
      - bowtie2_log

  merge_filter:
    run: "../workflow_modules/merge_filter.cwl"
    in:
      sample_id:
        source: sample_id
      bams: 
        source: trim_and_map/bam
      is_paired_end:
        source: is_paired_end
    out:
      - fastqc_zip
      - fastqc_html
      - bam
  

  tn5_overhang_correction:
    run: "../tools/tn5_overhang_correction.cwl"
    in:
      bam:
        source: merge_filter/bam
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
  
  chip_qc:
    run: "../workflow_modules/chip_qc.cwl"
    in:
      sample_id:
        source: sample_id
      bam:
        source: merge_filter/bam
      is_paired_end:
        source: is_paired_end
      user_def_fragment_size:
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
      - fragment_size

  generate_coverage_tracks:
    run: "../tools/deeptools_bamCoverage.cwl"
    in:
      bam:
        source: indexing_shifted_bam/bam_sorted_indexed
      is_paired_end:
        source: is_paired_end
      fragment_size:
        source: chip_qc/fragment_size
      effective_genome_size:
        source: effective_genome_size
      bin_size:
        source: bin_size
      ignoreForNormalization:
        source: ignoreForNormalization
    out:
      - bigwig

  create_summary_qc_report:
    doc: |
      multiqc summarizes the qc results from fastqc 
      and other tools
    run: "../tools/multiqc_hack.cwl"
    in:
      qc_files_array_of_array:
        source:
          - trim_and_map/raw_fastqc_zip
          - trim_and_map/raw_fastqc_html
          - trim_and_map/trimmed_fastqc_html
          - trim_and_map/trimmed_fastqc_zip
          - trim_and_map/trim_galore_log
        linkMerge: merge_flattened
      qc_files_array:
        source:
          - trim_and_map/bowtie2_log
          - merge_filter/fastqc_zip
          - merge_filter/fastqc_html
          - chip_qc/qc_plot_coverage_tsv
          - chip_qc/qc_plot_coverage_plot
          - chip_qc/qc_plot_fingerprint_tsv
          - chip_qc/qc_plot_fingerprint_plot
          - chip_qc/qc_crosscorr_summary
          - chip_qc/qc_crosscorr_plot
          - chip_qc/qc_phantompeakqualtools_stdout
          - chip_qc/qc_crosscorr_summary
        linkMerge: merge_flattened
      report_name:
        source: sample_id
    out:
      - multiqc_zip
      - multiqc_html


### OUTPUTS:
##################################################
outputs:
  raw_fastqc_zip:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_and_map/raw_fastqc_zip
  raw_fastqc_html:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_and_map/raw_fastqc_html
  trim_galore_log:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_and_map/raw_fastqc_zip
  trimmed_fastqc_html:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_and_map/trimmed_fastqc_html
  trimmed_fastqc_zip:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_and_map/trimmed_fastqc_zip
  bowtie2_log:
    type:
      type: array
      items: File
    outputSource: trim_and_map/bowtie2_log

  fastqc_zip:
    type:
      type: array
      items: File
    outputSource: merge_filter/fastqc_zip
  fastqc_html:
    type:
      type: array
      items: File
    outputSource: merge_filter/fastqc_html
  bam:
    type: File
    secondaryFiles: .bai
    outputSource: merge_filter/bam

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