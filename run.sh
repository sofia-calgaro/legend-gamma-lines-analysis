#!/usr/bin/env bash

detector_type="$(grep -oP '(?<="detectors": ")[^"]*' $1)"
if [ "$detector_type" == "single" ]; then
   detector_list=$(cat ./list_detectors_p3p4.json | jq -r '.[]')
   for d in $detector_list
   do
    if [ "$d" == ']' ] || [ "$d" == '[' ]; then
       continue
    fi
    d_name=$(echo $d| sed 's/"//g')
    d_name=$(echo $d_name| sed 's/,//g')
    JobName="$d""_bat"
    #qsub -N $JobName run-bat.qsub $1 $d_name $d
    done
else
    JobName="$detector""_bat"
    #qsub -N $JobName run-bat.qsub $1 $detector;
fi


