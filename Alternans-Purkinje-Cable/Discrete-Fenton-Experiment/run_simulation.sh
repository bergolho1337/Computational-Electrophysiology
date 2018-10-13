#!/bin/bash

echo "******************** SIMULATION **********************************"
./purkinje -s 0.1 900 Meshes/cable5cm-dog.msh steady_state/cable-5cm-300ms-Gap.sst Plot/cable5cm-dog.plt
echo "******************** SIMULATION **********************************"