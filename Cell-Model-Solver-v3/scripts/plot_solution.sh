#!/bin/bash

if test "$#" -ne 1; then
    echo "[!] ERROR! Illegal number of parameters"
    exit 1
fi

FOLDER_PATH=$1

# Get all the files and sorted in numerical order
FILES_ARRAY=$(ls $1V_t_* | sort -n -t _ -k 5)

# Get all the solution files from each timestep
LIST_FILES=""
for ITEM in ${FILES_ARRAY[*]}
do
    LIST_FILES="$LIST_FILES $ITEM"
done 

# Merge them together
cat $LIST_FILES > solution.dat

# Call the plot script
python plot_ap.py solution.dat