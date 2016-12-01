#!/bin/bash

if [ $# -ne 5 ]; then
    echo "Wrong number of parameters. [progname] templateNotMulti templateMulti outconfigbasename nbIter nbRuns"
    exit
fi
template=$1
templateMulti=$2
outbasename=$3
iter=$4
nbRuns=$5


##TODO ATTENTION USE template and templateMulti variables
program="./neatLearnFunction "
rm ../config/template-params-multi
cp ../config/template-params ../config/template-params-multi
#Add allowMulti 
echo "allowMultisynapses=true" >> ../config/template-params-multi

rm -rf $outbasename
mkdir $outbasename

mkdir $outbasename/noMulti
echo "No multi"
for (( j=1; j<=$nbRuns; j++))
do
echo $j
$program $template $iter > $outbasename/noMulti/run-$j.log
done

mkdir $outbasename/multi
echo "Multi"
for (( j=1; j<=$nbRuns; j++))
do
echo $j
$program $templateMulti $iter > $outbasename/multi/run-$j.log
done

#./tools/template-launchMultiple.sh config/tColl2/template-medea config/tColl2/medea-colors logs/expColors 30
cp $template $outbasename
cp $templateMulti $outbasename

for i in $outbasename/multi/* ; 
do 
awk '{ print $3 }' $i > $i.ave.log
awk '{ print $4 }' $i > $i.best.log
awk '{ print $5 }' $i > $i.species.log
awk '{ print $6 }' $i > $i.speciesSize.log
awk '{ print $7 }' $i > $i.depthBest.log
awk '{ print $8 }' $i > $i.linksBest.log
awk '{ print $9 }' $i > $i.nodesBest.log
awk '{ print $10 }' $i > $i.avgDepth.log
done
for i in $outbasename/noMulti/* ; 
do 
awk '{ print $3 }' $i > $i.ave.log
awk '{ print $4 }' $i > $i.best.log
awk '{ print $5 }' $i > $i.species.log
awk '{ print $6 }' $i > $i.speciesSize.log
awk '{ print $7 }' $i > $i.depthBest.log
awk '{ print $8 }' $i > $i.linksBest.log
awk '{ print $9 }' $i > $i.nodesBest.log
awk '{ print $10 }' $i > $i.avgDepth.log
done
paste $outbasename/multi/*ave.log >> $outbasename/multiAverageAll.log; 
paste $outbasename/multi/*best.log >> $outbasename/multiBestAll.log; 
paste $outbasename/noMulti/*ave.log >> $outbasename/noMultiAverageAll.log; 
paste $outbasename/noMulti/*best.log >> $outbasename/noMultiBestAll.log; 

paste $outbasename/multi/*species.log >> $outbasename/multiSpecies.log; 
paste $outbasename/multi/*speciesSize.log >> $outbasename/multiSpeciesSize.log; 
paste $outbasename/noMulti/*species.log >> $outbasename/noMultiSpecies.log; 
paste $outbasename/noMulti/*speciesSize.log >> $outbasename/noMultiSpeciesSize.log; 

paste $outbasename/multi/*depthBest.log >> $outbasename/multiDepthBest.log; 
paste $outbasename/multi/*linksBest.log >> $outbasename/multiLinksBest.log; 
paste $outbasename/multi/*nodesBest.log >> $outbasename/multiNodesBest.log; 
paste $outbasename/noMulti/*depthBest.log >> $outbasename/noMultiDepthBest.log; 
paste $outbasename/noMulti/*linksBest.log >> $outbasename/noMultiLinksBest.log; 
paste $outbasename/noMulti/*nodesBest.log >> $outbasename/noMultiNodesBest.log; 

paste $outbasename/multi/*avgDepth.log >> $outbasename/multiAvgDepth.log; 
paste $outbasename/noMulti/*avgDepth.log >> $outbasename/noMultiAvgDepth.log; 

filenames=`realpath $outbasename/*.log` ;
echo "fnames += ["
for i in $filenames ;
do
echo \"$i\",
done
echo "]"
