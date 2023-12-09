#!/bin/bash

# Input and output files
SOURCE_FILE="source.c"
OBJECT_FILE="source.o"

# Compilation command
gcc -c "$SOURCE_FILE" -o "$OBJECT_FILE" -lm

# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful."
else
    echo "Compilation failed."
fi

