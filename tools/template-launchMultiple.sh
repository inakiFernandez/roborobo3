#!/bin/bash

if [ $# -ne 4 ]; then
    echo "Wrong number of parameters. [progname] template outconfigbasename outlogbasename nbRuns"
    exit
fi
tmplate=$1
outbasename=$2
outlogbasename=$3
nbRuns=$4

nbRob='50 100 200 300 400'
taskSeq='1,-1 2,-1'
ctrlSetup='3'

listProp=`parallel --header : echo R{1}.T{2}.B{3} gInitialNumberOfRobots={f1} gTaskSeq={f2} gBrait={f3} ::: f1 $nbRob ::: f2 $taskSeq ::: f3 $ctrlSetup`

#echo "$listProp"

program="./roborobo -l "

mkdir $outlogbasename
commandFile=$tmplate.parallel

rm $commandFile
touch $commandFile

while read -r line
do
suffix=`echo $line | cut -f1 -d" " | sed -e 's/,-1//' |sed -e 's/,/T/g'`
#echo $suffix 

cp $tmplate $outbasename-$suffix.properties
mkdir $outlogbasename/$suffix
#echo $line
echo $line | perl -p -e  's/^.*? //' | sed -e 's/ /\n/g' >> $outbasename-$suffix.properties

    for (( j=1; j<=$nbRuns; j++))
    do
	echo "$program $outbasename-$suffix.properties > $outlogbasename/$suffix/run-$j.log" >> $commandFile    
    done

done <<< "$listProp"


parallel -j8 -a $commandFile



#./tools/template-launchMultiple.sh config/tColl2/template-medea config/tColl2/medea-colors logs/expColors 30


