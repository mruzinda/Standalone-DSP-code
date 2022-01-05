#!/bin/bash

#ls -Art | head -n 1 # Oldest file
#ls -Art | tail -n 1 # Newest file
#oldest_file=$(ls $path_to_raw_files -Art | head -n 1) # Finds the oldest file in the directory based on timestamp
#basefile=$path_to_raw_files${oldest_file%%.*}

# Path to file in first argument of execution
path_to_raw_files=$1

# Finds the first RAW file in the directory
first_file=$(find $path_to_raw_files -type f -iname "*0000.raw")

# Removes everything after the first period in the RAW file name
basefile=${first_file%%.*}

echo $first_file
echo $basefile

