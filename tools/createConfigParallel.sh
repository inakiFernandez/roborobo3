#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Wrong number of parameters"
    exit
fi

RESULT=""
DIRPATH="${1%/*}"
COMMAND="./roborobo -l "

BASENAME="${1##*/}"
EXTENSION="${BASENAME##*.}"
BASENAME="${BASENAME%.*}"
mkdir "logs/$BASENAME" 2> /dev/null
name=$1
#printf "$name \n $BASENAME $EXTENSION \n $DIRPATH \n $RESULT \n"

for i in `seq 1 $2`
do
    printf "$i\n"
    #if different parameters needed
    name="$BASENAME-run$i.properties" 
    cp $1 $name
    

    RESULT="$RESULT$COMMAND $DIRPATH/$BASENAME.$EXTENSION > logs/$BASENAME/$BASENAME-run$i.log\n"

done

printf "$RESULT"

#launch from roborobo root folder
#prepare first the parameter file.Check it is in batch mode. 
#launch with ./createConfigParallel.sh config/tColl2/test.properties 30 > parallel_file

#launch parallel with: parallel -j 8               -a parallel_file
#                                 Number of cores

#join output files with paste logs/$BASENAME/$BASENAME-run*.log > logs/$BASENAME-allRuns.log

#copy file to remote PC: scp createConfigParallel.sh adminmaia@maia1-is:~/fernandi_local/roborobo3/roborobo3/tools
#python3 tools/plotItemsByColor.py logs/items.log logs/itemsIter.log logs/colorChanges.log logs/graphs
# scp  adminmaia@maia1-is:~/fernandi_local/roborobo3/roborobo3/logs/$BASENAME-allRuns.log ./results