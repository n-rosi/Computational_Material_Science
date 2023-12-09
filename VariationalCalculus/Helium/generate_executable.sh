#!/bin/bash

# Set your GSL include and library paths
GSL_INCLUDE_DIR="/home/nrosi/C-lib/gsl/include"
GSL_LIB_DIR="/home/nrosi/C-lib/gsl/lib"

# Input and output files
MAIN_OBJECTIVE_FILE="main.o"
SOURCE_OBJECTIVE_FILE="source.o"
LINALG_OBJECTIVE_FILE="linalg.o"

# Compilation command
gcc "$LINALG_OBJECTIVE_FILE" "$SOURCE_OBJECTIVE_FILE" "$MAIN_OBJECTIVE_FILE" -I"$GSL_INCLUDE_DIR" -L"$GSL_LIB_DIR" -lgsl -lgslcblas -lm

# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo "Executable generated succesfully."
else
    echo "Compilation failed."
fi

