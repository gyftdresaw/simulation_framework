#!/bin/bash

# submit all jobs in current directory to queue

for filename in l2d*.py; do
    echo "exporting $filename"
    export filename;
    sbatch submit.sh
done
