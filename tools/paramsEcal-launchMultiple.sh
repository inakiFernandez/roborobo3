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


#taskSeq='1,2,1,2,-1'

#taskTimeChange='0,150000,400000,450000,-1'

#ctrlSetup='1 2 3' #'3'
#ctrlSetup='5 6' 
#ctrlSetup='0' 

#evotop='3 0' #3=evotopo, 0,1,2= fixed topo (MLP, Elman, Perceptron)
evotop='0'


#sigma='0.5 0.05 0.1'
sigma='0.1'

#selpres='0.0 0.25 0.5 0.75 1.0'
selpres='1.0 0.0'

#multisynapses='false true'
#multisynapses='false'


listProp=`parallel --header : echo R{1}.Top{3}.S{2}.SP{4} gInitialNumberOfRobots={f1} gSigmaRef={f2} gControllerType={f3} gSelPressure={f4} ::: f1 $nbRob ::: f2 $sigma ::: f3 $evotop ::: f4 $selpres`

#echo "$listProp"

#exit

program="./roborobo -l "

mkdir $outlogbasename
commandFile=$tmplate.parallel

rm $commandFile
touch $commandFile

while read -r line
do
    suffix=`echo $line | cut -f1 -d" " | sed -e 's/,-1//' |sed -e 's/,/T/g'`
    #echo $suffix 
    
    #same file for all runs
    #cp $tmplate $outbasename-$suffix.properties

    mkdir $outlogbasename/$suffix
    #echo $line
    echo $line | perl -p -e  's/^.*? //' | sed -e 's/ /\n/g' >> $outbasename-$suffix.properties


    for (( j=1; j<=$nbRuns; j++))
    do
	#Different file for different runs
	cp $tmplate $outbasename-$suffix-run$j.properties
	echo "logItemName = $outlogbasename/$suffix/items-run-$j.log" >>  $outbasename-$suffix-run$j.properties
	echo "logItGatheredName = $outlogbasename/$suffix/itemsIter-run-$j.log" >>  $outbasename-$suffix-run$j.properties
	echo "logColorChangesName = $outlogbasename/$suffix/colorChanges-run-$j.log" >>  $outbasename-$suffix-run$j.properties
#	echo "$program $outbasename-$suffix.properties > $outlogbasename/$suffix/run-$j.log" >> $commandFile    
	echo "$program $outbasename-$suffix-run$j.properties > $outlogbasename/$suffix/run-$j.log" >> $commandFile    
    done

done <<< "$listProp"



#parallel -j15 -a $commandFile
parallel -j4 -a $commandFile

for (( j=1; j<=$nbRuns; j++))
do
    echo "python3 tools/plotItemsByColor.py $outlogbasename/$suffix/items-run-$j.log $outlogbasename/$suffix/itemsIter-run-$j.log $outlogbasename/$suffix/colorChanges-run-$j.log $outlogbasename/$suffix/graph-run$j --png"
done

#./tools/paramsEcal-launchMultiple.sh config/tColl2/C2-8colors-sigma0.1-800steps.properties config/tColl2/C2-8colors logs/exp8Colors 30


#for i in * ; do paste $i/run*.log > $i.all.log; done

#for i in * ; do paste $i/itemsIter-run*.log > $i.itemsIter.all.log; done

#for i in * ; do paste $i/items-run*.log > $i.items.all.log; done

#for i in * ; do paste $i/colorChanges-run*.log > $i.colorChanges.all.log; done



#for i in *.log ; do tail -n +21 $i > tmp ; rm $i ; mv tmp $i ; done ;

#for i in *.log ; do head -n -6 $i > tmp ; rm $i ; mv tmp $i ; done ;