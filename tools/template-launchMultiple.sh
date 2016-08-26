#!/bin/bash

if [ $# -ne 4 ]; then
    echo "Wrong number of parameters. [progname] template outconfigbasename outlogbasename nbRuns"
    exit
fi
tmplate=$1
outbasename=$2
outlogbasename=$3
nbRuns=$4

#properties=('gInitialNumberOfRobots' 'gTaskSeq') #(p1,p2,...)
#values1=('50' '100' '200' '300' '400')
#values2=('1,-1' '2,-1')
#values=(('50','100','200','300','400'),('1,-1','2,-1'))
#values=($values1 $values2)
#parallel echo ::: 'gInitialNumberOfRobots=' ::: '50' '100' '200' '300' '400'
#parallel echo ::: 'gInitialNumberOfRobots=' ::: '50' '100' '200' '300' '400' | parallel echo ::: 'gTaskSeq=' ::: '1,-1' '2,-1' :::: -
#parallel --header :  echo gInitialNumber={gInitialNumberOfRobots} gTaskSeq={gTaskSeq} ::: gInitialNumberOfRobots 50 100 200 300 400 ::: gTaskSeq 1,-1 2,-1 | sed -e 's/ /\n/g'
#rep=`seq 1 1 30`
#listProp=`parallel --header : echo  gInitialNumberOfRobots={f1} gTaskSeq={f2} it={f3} ::: f1 $nbRob ::: f2 $taskSeq ::: f3 $rep`

nbRob='50 100 200 300 400'
taskSeq='1,-1 2,-1'
ctrlSetup='1 2'

listProp=`parallel --header : echo R{1}.T{2}.B{3} gInitialNumberOfRobots={f1} gTaskSeq={f2} gBrait=f{3} ::: f1 $nbRob ::: f2 $taskSeq ::: f3 $ctrlSetup`

#echo "$listProp"

program="./roborobo -l "

mkdir $outlogbasename

#i=1
while read -r line
do
suffix=`echo $line | cut -f1 -d" " | sed -e 's/,-1//' |sed -e 's/,/T/g'`
#echo $suffix 

cp $tmplate $outbasename-$suffix.properties
mkdir $outlogbasename/$suffix

commandFile=$tmplate.parallel

rm $commandFile
touch $commandFile

echo $line | perl -p -e  's/^.*? //' | sed -e 's/ /\n/g' >> $outbasename-$suffix.properties


    for (( j=1; j<=$nbRuns; j++))
    do
	echo "$program $outbasename-$suffix.properties > $outlogbasename/$suffix/run-$j.log" >> $commandFile    
    done

#i=$(($i+1))
done <<< "$listProp"


parallel --dry-run -j8 -a $commandFile



#./tools/template-launchMultiple.sh config/tColl2/template-medea config/tColl2/medea-colors logs/expColors 30




#indexProp=0
#for i in ${properties[@]}; do
#    for j in ${val[@]}; do	
#	echo "" # "$i=$j"
#    done
#    indexProp=$(($indexProp+1))
#done

