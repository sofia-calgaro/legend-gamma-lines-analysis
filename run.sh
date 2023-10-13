#!/usr/bin/env bash

# Check if the folder is not empty; then, remove content
if [ -n "$(ls -A tmp/)" ]; then
  rm tmp/*
fi

gamma_src_code="$(grep -oP '(?<="gamma-src-code": ")[^"]*' $1)"
output_folder="$(grep -oP '(?<="output": ")[^"]*' $1)"

gamma_lines="$(grep -oP '(?<="gamma-lines": )(true|false)' $1)"
peak_search="$(grep -oP '(?<="run": )(true|false)' $1)"

# check if you selected at least one analysis
if [ "$gamma_lines" == "$peak_search" ]
then
  echo "You have to select at least one analysis (gamma lines or peak search). Exit here!"
  exit 1
fi

# get cut and detector to inspect
cut="$(grep -oP '(?<="cut": ")[^"]*' $1)"
cut=$(echo "${cut//\ }")
detector_type="$(grep -oP '(?<="detectors": ")[^"]*' $1)"

# perform gamma lines analysis, if ture
if [ "$detector_type" == "single" ]; then
  if [ "$peak_search" == "true" ]; then
    echo "Peak search is still not available for single detectors. Exit here!"
      exit 1
  fi
  file="src/settings/list_detectors.json"
  detector_list=$(python-3.9.6 -c \
"import json;
with open('$file', 'r') as f:
  print(json.load(f)['detectors'])")
  for d in $detector_list
   do
    d_name=$(echo "${d//\'}")
    d_name=$(echo "${d_name//\[}")
    d_name=$(echo "${d_name//\]}")
    d_name=$(echo "${d_name//\,}")
    analysis="gamma-lines"
    JobName="line_""$cut"_"$d_name""_bat"
    qsub -N $JobName run-bat.qsub $1 $gamma_src_code $output_folder $detector_type $d_name $cut $analysis
    sleep 10
    done
else
    if [ "$gamma_lines" == "true" ]; then
      analysis="gamma-lines"
      JobName="line_""$cut"_"$detector_type""_bat"
      qsub -N $JobName run-bat.qsub $1 $gamma_src_code $output_folder $detector_type $cut $analysis
    fi
    if [ "$peak_search" == "true" ]; then
      analysis="peak-search"
      JobName="peak_""$cut"_"$detector_type""_bat"
      qsub -N $JobName run-bat.qsub $1 $gamma_src_code $output_folder $detector_type $cut $analysis
    fi
fi


