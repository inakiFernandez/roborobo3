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
nbRob='200'

taskSeq='1,-1 2,-1'
taskSeq='2,-1'
taskSeq='1,2,1,2,-1'
taskSeq='1,-1'

#taskTimeChange='0,150000,400000,450000,-1'
taskTimeChange='0,-1'


ctrlSetup='1 2 3' #'3'
ctrlSetup='5 6' 
ctrlSetup='0' 

evotop='3 0' #3=evotopo, 0,1,2= fixed topo (MLP, Elman, Perceptron)
evotop='3'


sigma='0.5'
sigma='0.5 0.05 0.1'

selpres='0.0 0.25 0.5 0.75 1.0'
selpres='1.0'

multisynapses='false true'
#multisynapses='false'


listProp=`parallel --header : echo R{1}.T{2}.Top{6}.B{3}.S{4}.SP{7}.Time{5}.M{8} gInitialNumberOfRobots={f1} gTaskSeq={f2} gBrait={f3} gSigmaRef={f4} gTimeChange={f5} gControllerType={f6} gSelPressure={f7} allowMultisynapses={f8} ::: f1 $nbRob ::: f2 $taskSeq ::: f3 $ctrlSetup ::: f4 $sigma ::: f5 $taskTimeChange ::: f6 $evotop ::: f7 $selpres ::: f8 $multisynapses`

#listProp=`parallel --header : echo R{1}.T{2}.B{3}.S{4} gInitialNumberOfRobots={f1} gTaskSeq={f2} gBrait={f3} gSigmaRef={f4} ::: f1 $nbRob ::: f2 $taskSeq ::: f3 $ctrlSetup ::: f4 $sigma`

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

#FOR SEQUENCE ONLY
#echo "gTimeChange=$taskTimeChange" >> $outbasename-$suffix.properties


    for (( j=1; j<=$nbRuns; j++))
    do
	echo "$program $outbasename-$suffix.properties > $outlogbasename/$suffix/run-$j.log" >> $commandFile    
    done

done <<< "$listProp"


parallel -j15 -a $commandFile



#./tools/template-launchMultiple.sh config/tColl2/template-medea config/tColl2/medea-colors logs/expColors 30


#for i in * ; do paste $i/*.log > $i.all.log; done

#for i in *.log ; do tail -n +21 $i > tmp ; rm $i ; mv tmp $i ; done ;

#for i in *.log ; do head -n -6 $i > tmp ; rm $i ; mv tmp $i ; done ;