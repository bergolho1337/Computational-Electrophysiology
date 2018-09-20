#!/bin/bash

if test "$#" -ne 1; then
    echo "[!] ERROR! Illegal number of parameters"
    exit 1
fi

FILE_PATH="$1/V_t" 

# Call the plot script
python plot_ap.py $FILE_PATH