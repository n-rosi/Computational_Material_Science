#!/bin/bash

# Set your GSL include and library paths
GSL_INCLUDE_DIR="/home/nrosi/C-lib/gsl/include"
GSL_LIB_DIR="/home/nrosi/C-lib/gsl/lib"

# Input and output files
SOURCE_FILE="main.c"
OBJECT_FILE="main.o"

# Compilation command
gcc -c "$SOURCE_FILE" -o "$OBJECT_FILE" -I"$GSL_INCLUDE_DIR" -L"$GSL_LIB_DIR" -lgsl -lgslcblas -lm

# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful."
else
    echo "Compilation failed."
fi

