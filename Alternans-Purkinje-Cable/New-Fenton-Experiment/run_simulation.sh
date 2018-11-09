#!/bin/bash

run_simulation () {
    echo "****************************************************************************************************"
    echo "[!] Running simulation with $1ms pacing ..."
    echo "****************************************************************************************************"
    # Steady-State 
    ./bin/FentonExperiment -c example_configs/feedback-380cm-s/sst_noble_5cm_$1ms.ini
    cp output.sst steady_state/cable-5cm-$1ms-Gap.sst
    ./clear_results.sh
    # Experiment
    ./bin/FentonExperiment -c example_configs/feedback-380cm-s/simple_noble_5cm_$1ms.ini
    mkdir scripts/$1ms
    cp output/*.dat scripts/$1ms
    echo "****************************************************************************************************"
}

#PACINGS=( 200 210 220 230 240 250 260 270 280 290 300 )
PACINGS=($(seq 280 -5 250))

for PACING in "${PACINGS[@]}"; do
    run_simulation $PACING
done