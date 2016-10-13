#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Wrong number of parameters"
    exit
fi

nbRuns=$1


for i in `seq 1 $1`
do
    ./evaluateRandomMultisynapse 0 > logsRuns/$i.log
    ./evaluateRandomMultisynapse 1 > logsRuns/$i.multi.log
done
