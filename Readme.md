# CWL-based ChIP-seq and Cut&Run Workflows  
  
This repository includes workflows for upstream processing of ChIP-seq and Cut&Run data. The workflows are written in CWL and they are ready to run with a containerization solution like docker or sigularity.

## CWL Workflows:  
Following workflows are available in the directory `./CWL/workflows`:
  * `ChIPseq_pipeline.cwl` - works with single or paired end standard ChIPseq data
  * `ChIPseq_pipeline_spike_in.cwl` - same as above but with support for spike-in normalization
  * `ChIPmentation_pipeline.cwl` - works with single or paired end ChIPmentation data  
  
They depend on subworkflows located in `./CWL/workflow_modules` and CommandLineTool wrapper located in `./CWL/tools`.

## How to use:
### 0. How to set up:
Clone this repository, or download and extract a specific release.

You need to install a cwlrunner. We recommend cwltool for execution on a single machine or toil if you would like to use an HPC cluster or cloud infrastructure. A snapshot of a dedicated conda environment can be found at `conda_toil_env.yml`. However, every other runner that is fully compatible with CWL v1.0 should work, too.

To use the docker containers specified in the workflows, you need to set up docker or singularity.

Moreover, you need to install node.js for evaluation of the built-in java script expressions.

### 1. How to run:
For workflow execution with cwltool and toil, you can find example shell scripts in the top directory (`run_*.sh`). They can be use in the following fashion (however, you likely need to slightly adjust the cwlrunner command line specified in the scripts):  
`bash ./run_cwltool.sh <cwl_workflow> <input_file_or_dir> <output_dir>`  
Or: `bash ./run_cwltool.sh <cwl_workflow> <input_file_or_dir> <output_dir>`   
  
  - `<cwl_workflow>` is a path to a cwl workflow.  
  - `<input_file_or_dir>` is the path to a YAML file containing input parameters, or the path to a directory with multiple YAML jobs that should be executed
  - `<output_dir>` is the directory where final output files will be collected
  
Monitoring:  
A log file (the name is ending with `.yaml.log`) will be created for each sample/run. You can open them to check the progress.  
  
Example/test run:  
Please try to run the following tiny test example to check whether your system is correctly configured.
`bash ./run_cwltool.sh ./CWL/workflows/ChIPseq_pipeline.cwl ./tiny_test/test_main.chip.yml ~/cwl_test`  
  
  
Terminating a job:  
If you like to terminate a specific run and all its subprocessed, you can find the corresponding job-ids in the file `cwl_background_job_ids.log`. Run: `kill <job_id>`.  

### 2. Check the output
Once the processing finished, you can find one QC report per workflow iteration (name ends with `*_report.html`). These reports were created using multiqc.

They can give you a deep understanding of your data's quality. However, you should also keep in mind that unexpected results can come from bugs in the custom scripts or the tools which are applied by the workflow. Therefore, check the report carefully. Especially, examine the "read flow" through the pipeline using the fastqc output: pay attention whether you encounter an unexpected loss of reads.
  
## Testing and debugging
There are tiny synthetic tests located in `./tiny_test/`. They are perfect for initial testing and debugging.  