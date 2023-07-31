#!/usr/bin/env bash

detector_type="$(grep -oP '(?<="detectors": ")[^"]*' $1)"
if [ "$detector_type" == "single" ]; then
   file="./list_detectors_p3p4.json"
   detector_list=$(python-3.9.6 -c \
"import json;
with open('$file', 'r') as f:
  print(json.load(f)['detector_p3p4'])")
  for d in $detector_list
   do
    d_name=$(echo "${d//\'}")
    d_name=$(echo "${d_name//\[}")
    d_name=$(echo "${d_name//\]}")
    d_name=$(echo "${d_name//\,}")
    JobName="$d_name""_bat"
    echo $d_name
    echo $d
    #./run_test.sh $1 $detector_type $d_name
    qsub -N $JobName run-bat.qsub $1 $detector_type $d_name
    done
   else
    JobName="$detector_type""_bat"
    qsub -N $JobName run-bat.qsub $1 $detector_type;
fi

