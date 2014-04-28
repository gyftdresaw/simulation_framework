#!/bin/bash

# move completed simulations into done folder
pyext=".py"

for fname in *.p; do
    filebase=$(echo $fname | awk '{split($0,a,"."); print a[1]}')
    mv $fname done
    mv $filebase$pyext done
done
