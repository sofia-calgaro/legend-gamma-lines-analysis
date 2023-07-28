#!/usr/bin/env bash

detector="$(grep -oP '(?<="detectors": ")[^"]*' $1)"
JobName="$detector""_bat"
qsub -N $JobName run-bat.qsub $1 $detector;
