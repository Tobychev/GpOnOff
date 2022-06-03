#!/usr/bin/env bash

rootdir=$1
condadir=$2
condaenv=$3
srcname=$4

startdate=$(date --iso)

cd $rootdir
source $condadir
conda activate $condaenv
conffile=$(printf "%s/gp_%s.yaml" "$srcname" "$srcname")
outdir=$(printf "%s/GP_Result" "$srcname")
reportname=$(printf "%s-%s" "$srcname" "$startdate")
printf "Running analysis %s" "$conffile"
python gplib/OnOff-Analysis.py "$conffile" --debug
python gplib/SpectrumFit.py "$conffile"
python gplib/MakeGPReport.py "$outdir" "$reportname"
