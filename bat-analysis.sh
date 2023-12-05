#!/bin/bash
#
#$ -j y     ### have logoutput and erroroutput in the same file
#$ -o $PWD
#$ -e $PWD

LOGFILE=$4".out"
OUT_FOLDER=$3

# Check if the folders already exist
if [ ! -d "$OUT_FOLDER" ]; then
    mkdir -p "$OUT_FOLDER"
fi
if [ ! -d "$OUT_FOLDER/log" ]; then
    mkdir -p "$OUT_FOLDER/log"
fi
if [ ! -d "$OUT_FOLDER/$4" ]; then
    mkdir -p "$OUT_FOLDER/$4"
fi

cd $2

if [ $4 == 'All' ] || [ $4 == 'BEGe' ] || [ $4 == 'ICPC' ] || [ $4 == 'COAX' ] || [ $4 == 'PPC' ]; then
    singularity exec /lfs/l1/legend/software/singularity/legendexp_legend-software_latest.sif bash -c "source /usr/local/bin/thisroot.sh && cd $2 && python src/main.py --config $1 --det $4" &> $OUT_FOLDER/log/$LOGFILE
    source ~/.bashrc
    #cd src/gamma-fitter/runGammaAnalysis
    #make runGammaAnalysis &&
    if [ "$6" == "gamma-lines" ]; then
        src/gamma-fitter/runGammaAnalysis $5 $4 false $6 &>> $OUT_FOLDER/log/$LOGFILE
    fi
    if [ "$6" == "peak-search" ]; then
        # get additional parameters
        E_min="$(grep -oP '"energy-min-keV": \K\d+' $1)"
        E_max="$(grep -oP '"energy-max-keV": \K\d+' $1)"
        E_step="$(grep -oP '"energy-step-keV": \K\d+' $1)"
        E_window="$(grep -oP '"fit-window-width-keV": \K\d+' $1)"
        src/gamma-fitter/runGammaAnalysis $5 $4 false $6 $E_min $E_max $E_step $E_window $OUT_FOLDER/$4 &>> $OUT_FOLDER/log/$LOGFILE
    fi
else
    singularity exec /lfs/l1/legend/software/singularity/legendexp_legend-software_latest.sif bash -c "source /usr/local/bin/thisroot.sh && cd $2 && python src/main.py --config $1 --det $5" &> $OUT_FOLDER/log/$LOGFILE
    source ~/.bashrc
    #cd src/gamma-fitter/runGammaAnalysis
    #make runGammaAnalysis &&
    src/gamma-fitter/runGammaAnalysis $6 ${5: -7} true $7 &> $OUT_FOLDER/log/$LOGFILE

fi
