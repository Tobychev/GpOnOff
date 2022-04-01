#!/usr/bin/env bash

rootdir=$1
condadir=$2
condaenv=$3
srcname=$4

cd $rootdir
source $condadir
conda activate $condaenv
conffile=$(printf "%s/gp_%s.yaml" "$srcname" "$srcname")
printf "Running analysis %s" "$conffile"
python gplib/OnOff-Analysis.py $conffile
