#!/bin/bash

run_simulation () {
    echo "****************************************************************************************************"
    echo "[!] Running simulation with $1ms pacing and a fiber with $2um of diameter and a size of $3cm ..."
    echo "****************************************************************************************************"
    # Steady-State 
    ./purkinje example_configs/$2um/$3cm/sst_simple_purkinje_$3cm_$1ms.in
    cp output.sst steady_state/$2um/$3cm/cable-$3cm-$1ms-Gap.sst
    make clcResults
    # Experiment
    ./purkinje example_configs/$2um/$3cm/simple_purkinje_$3cm_$1ms.in
    mkdir scripts/$2um/$3cm/$1ms
    cp output/*.dat scripts/$2um/$3cm/$1ms
    echo "****************************************************************************************************"
}

if [ ! -f purkinje ]; then
    make
fi

#PACINGS=( 200 210 220 230 240 250 260 270 280 290 300 )
PACINGS=($(seq 190 -10 100))
DIAMETER=3000
CABLE_SIZE=10

for PACING in "${PACINGS[@]}"; do
    run_simulation $PACING $DIAMETER $CABLE_SIZE
done
