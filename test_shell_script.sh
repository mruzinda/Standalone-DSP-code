#!/bin/bash

#ls -Art | head -n 1 # Oldest file
#ls -Art | tail -n 1 # Newest file

# Path to file in first argument of execution
path_to_raw_files=$1

# Finds the oldest file in the directory
oldest_file=$(ls $path_to_raw_files -Art | head -n 1)

# Removes everything after the first period in the RAW file name
basefile=${oldest_file%%.*}

echo $oldest_file
echo $basefile

