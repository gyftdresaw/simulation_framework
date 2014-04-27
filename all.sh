#!/bin/bash

# submit all jobs in current directory to queue

for f in l2d*.py; do
    echo $f
    export f;
done
