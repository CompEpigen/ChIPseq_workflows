#!/bin/bash

## input parameters:
CWLSCRIPT="$1" # Path to CWL workflow
INPUT="$2" # Path to folder containing multiple yaml input files
		   # Or path to a single yaml input file
OUTDIR="$3" # Path to output dir (will be created if doesn't exist)

## working and temp dirs - please adapt accordingly:
BASEDIR=${HOME}/cwl_working_dir/base
WORKDIR=${HOME}/cwl_working_dir/work
TMPDIR=${HOME}/cwl_working_dir/tmp
TMPOUTDIR=${HOME}/cwl_working_dir/tmp_out

## create dirs if not exist
if [ ! -d "$OUTDIR" ]
then
	mkdir -p "$OUTDIR"
fi

if [ ! -d "$BASEDIR" ]
then
	mkdir -p "$BASEDIR"
fi

if [ ! -d "$WORKDIR" ]
then
	mkdir -p "$WORKDIR"
fi

if [ ! -d "$TMPDIR" ]
then
	mkdir -p "$TMPDIR"
fi

if [ ! -d "$TMPOUTDIR" ]
then
	mkdir -p "$TMPOUTDIR"
fi

## get input files:
if [ -d "$INPUT" ]
then
	## input is a directory containing multiple yml files
	YAML_INPUTS=($( ls "${INPUT}"/*.y*ml ))
else
	## a single yml file was given as input
	YAML_INPUTS=( "$INPUT" )
fi


## Start running:
echo ">>> Working with $CWLSCRIPT"
echo "> Starting jobs:"

for((i=0; i<"${#YAML_INPUTS[@]}"; i++))
do
	echo "- ${YAML_INPUTS[i]}"
	cwltool --parallel --debug --singularity \
		--tmp-outdir-prefix "$TMPOUTDIR" --tmpdir-prefix "$TMPDIR" \
		--basedir "$BASEDIR" --outdir "$OUTDIR" \
		"$CWLSCRIPT" "${YAML_INPUTS[i]}" > "${OUTDIR}/${YAML_INPUTS[i]##*/}".log 2>&1 &
	echo -e "${YAML_INPUTS[i]}\t$!" >> "${OUTDIR}/cwl_background_job_ids.log"
	echo "  job id: $!"
	disown %-
done

## packed copy of the workflow to output dir
cwltool --pack "$CWLSCRIPT" > "${OUTDIR}/workflow.cwl" 2> /dev/null

echo "> done"
echo "> jobs are running in the background and are detached from this terminal session"
echo "> for checking the progress please have a look at the \".yaml.log\" files in the output directory"
echo "> OUTPUT DIR: $OUTDIR"
echo "> You can close the terminal now."
