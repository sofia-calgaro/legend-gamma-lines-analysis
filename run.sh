#!/bin/bash

singularity exec /lfs/l1/legend/software/singularity/legendexp_legend-software_latest.sif bash -c "source /usr/local/bin/thisroot.sh && python main.py --config config.json"
#cd gamma-fitter-BATv100 && make 
./gamma-fitter-BATv100/runGammaAnalysis 