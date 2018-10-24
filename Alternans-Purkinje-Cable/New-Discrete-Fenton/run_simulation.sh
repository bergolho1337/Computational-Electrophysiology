#!/bin/bash

run_simulation () {
    echo "*********************************************************************************************"
    echo "[!] Running simulation with $1ms pacing ..."
    echo "*********************************************************************************************"
    # Steady-State 
    ./purkinje example_configs/sst_simple_purkinje_5cm_$1ms.in
    cp output.sst steady_state/cable-5cm-$1ms-Gap.sst
    make clcResults
    # Experiment
    ./purkinje example_configs/simple_purkinje_5cm_$1ms.in
    mkdir scripts/$1ms
    cp output/*.dat scripts/$1ms
    echo "*********************************************************************************************"
}

if [ ! -f purkinje ]; then
    make
fi

#PACINGS=( 200 210 220 230 240 250 260 270 280 290 300 )
PACINGS=($(seq 190 -10 100))

for PACING in "${PACINGS[@]}"; do
    run_simulation $PACING
done
