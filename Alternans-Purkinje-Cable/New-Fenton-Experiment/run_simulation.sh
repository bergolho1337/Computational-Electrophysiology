#!/bin/bash

run_simulation_Noble () {
    echo "****************************************************************************************************"
    echo "[!] Running Noble simulation with $1ms pacing ..."
    echo "****************************************************************************************************"
    # Steady-State 
    ./bin/FentonExperiment -c example_configs/feedback-NoGap-Noble/sst_noble_5cm_NoGap_$1ms.ini
    cp output.sst steady_state/cable-5cm-$1ms-NoGap.sst 
    ./clear_results.sh
    # Experiment
    ./bin/FentonExperiment -c example_configs/feedback-NoGap/simple_noble_5cm_NoGap_$1ms.ini
    mkdir scripts/$1ms
    cp output/*.dat scripts/$1ms
    echo "****************************************************************************************************"
}

run_simulation_BR () {
    echo "****************************************************************************************************"
    echo "[!] Running Beeler-Reuter simulation with $1ms pacing ..."
    echo "****************************************************************************************************"
    # Steady-State 
    ./bin/FentonExperiment -c example_configs/feedback-NoGap-BR/sst_beeler_reuter_8cm_NoGap_$1ms.ini
    cp output.sst steady_state/cable-BR-8cm-$1ms-NoGap.sst 
    ./clear_results.sh
    # Experiment
    ./bin/FentonExperiment -c example_configs/feedback-NoGap-BR/simple_beeler_reuter_8cm_NoGap_$1ms.ini
    mkdir scripts/$1ms
    cp output/*.dat scripts/$1ms
    echo "****************************************************************************************************"
}

#PACINGS=( 200 210 220 230 240 250 260 270 280 290 300 )
PACINGS=($(seq 295 -5 200))

for PACING in "${PACINGS[@]}"; do
    #run_simulation_Noble $PACING
    run_simulation_BR $PACING
done
