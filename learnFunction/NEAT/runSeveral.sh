#!/bin/bash
oneRun()
{
    local j=$1
    echo $j
    $program $template $iter > $outbasename/noMulti/run-$j.log
}
if [ $# -ne 6 ]; then
    echo "Wrong number of parameters. [progname] templateNotMulti templateMulti outconfigbasename nbIter nbRuns playVideo"
    exit
fi
template=$1
templateMulti=$2
outbasename=$3
iter=$4
nbRuns=$5
playVideo=$6
program="./neatLearnFunction "
nbExp=1
doMulti="false"
if [ $doMulti = "true" ]; then
    echo "Do multi"
    nbExp=2
    ##TODO ATTENTION USE template and templateMulti variables    
    rm ../config/template-params-multi # $2
    cp ../config/template-params ../config/template-params-multi # $1 $2
    #Add allowMulti 
    echo "allowMultisynapses=true" >> ../config/template-params-multi #$2
fi

rm -rf $outbasename
mkdir $outbasename

mkdir $outbasename/noMulti
echo "No multi"
for (( j=1; j<=$nbRuns; j++))
do
echo $j
$program $template $iter > $outbasename/noMulti/run-$j.log
#oneRun $j &
done

folderFiles="./sandbox"
folderHeatmapFiles="./heatmapData"
if [ $doMulti = "true" ]; then
rm $folderFiles/*.dat   
cd $folderFiles
find -maxdepth 1 -type f -name 'g*.dat-ind*' -exec rm {} \;
cd ..
rm $folderHeatmapFiles/*.csv

mkdir $outbasename/multi
echo "Multi"
for (( j=1; j<=$nbRuns; j++))
do
echo $j
$program $templateMulti $iter > $outbasename/multi/run-$j.log 
done
cp $templateMulti $outbasename
fi

cp $template $outbasename


sleep 2

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
awk '{ print $11 }' $i > $i.otherBest.log
awk '{ print $12 }' $i > $i.otherAve.log
awk '{$1=$2=$3=$4=$5=$6=$7=$8=$9=$10=$11=$12=""; print $0 }' $i > $i.allFitnessAndOther.log
done

if [ $doMulti = "true" ]; then
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
	awk '{ print $11 }' $i > $i.otherBest.log
	awk '{ print $12 }' $i > $i.otherAve.log
	awk '{$1=$2=$3=$4=$5=$6=$7=$8=$9=$10=$11=$12=""; print $0 }' $i > $i.allFitnessAndOther.log
    done

paste $outbasename/multi/*ave.log >> $outbasename/multiAverageAll.log; 
paste $outbasename/multi/*best.log >> $outbasename/multiBestAll.log; 
paste $outbasename/multi/*species.log >> $outbasename/multiSpecies.log; 
paste $outbasename/multi/*speciesSize.log >> $outbasename/multiSpeciesSize.log; 
paste $outbasename/multi/*depthBest.log >> $outbasename/multiDepthBest.log; 
paste $outbasename/multi/*linksBest.log >> $outbasename/multiLinksBest.log; 
paste $outbasename/multi/*nodesBest.log >> $outbasename/multiNodesBest.log; 
paste $outbasename/multi/*avgDepth.log >> $outbasename/multiAvgDepth.log; 
paste $outbasename/multi/*otherBest.log >> $outbasename/multiOtherBestAll.log; 
paste $outbasename/multi/*otherAve.log >> $outbasename/multiOtherAverageAll.log;

#TOFIX Paste is not good with several columns per run
paste $outbasename/multi/*allFitnessAndOther.log >> $outbasename/multiallFitnessAndOtherAll.log; 

fi
 
paste $outbasename/noMulti/*species.log >> $outbasename/noMultiSpecies.log; 
paste $outbasename/noMulti/*speciesSize.log >> $outbasename/noMultiSpeciesSize.log; 
paste $outbasename/noMulti/*ave.log >> $outbasename/noMultiAverageAll.log; 
paste $outbasename/noMulti/*best.log >> $outbasename/noMultiBestAll.log; 
paste $outbasename/noMulti/*depthBest.log >> $outbasename/noMultiDepthBest.log; 
paste $outbasename/noMulti/*linksBest.log >> $outbasename/noMultiLinksBest.log; 
paste $outbasename/noMulti/*nodesBest.log >> $outbasename/noMultiNodesBest.log; 
paste $outbasename/noMulti/*avgDepth.log >> $outbasename/noMultiAvgDepth.log; 
paste $outbasename/noMulti/*otherBest.log >> $outbasename/noMultiOtherBestAll.log; 
paste $outbasename/noMulti/*otherAve.log >> $outbasename/noMultiOtherAverageAll.log; 

#TOFIX Paste is not good with several columns per run
paste $outbasename/noMulti/*allFitnessAndOther.log >> $outbasename/noMultiallFitnessAndOtherAll.log; 

filenames=`realpath $outbasename/*.log` ;

buildVideo=true

if [ "$buildVideo" = true ] ; then
    rm ../tmpDataFilesExp.py
    touch ../tmpDataFilesExp.py

    printf "fnames = [" >> ../tmpDataFilesExp.py
    for i in $filenames ;
    do
	printf \"$i\", >> ../tmpDataFilesExp.py
    done
    printf "]" >> ../tmpDataFilesExp.py
    sleep 2
    #python fitness
    fitnessDataScript="../plotFitnessInMultiExp.py"
    
    #Extract duration per function
    durPerFunction=`grep "durationPerFunction = .*$" ../config/template-params | sed -e "s/durationPerFunction = \(.*\)$/\1/"`
    python3 $fitnessDataScript $folderFiles/$outbasename.png $nbExp --png $durPerFunction
    
    echo "Fitness Img done"
    
    sleep 2

    videoApproxScript="../gnuplotFolder.sh"
    videoName="$($videoApproxScript $folderFiles true true)"
    
    echo "Video approx done"
    
    sleep 2
    echo $videoName

    heatmapDataScript="../heatmap.py"

    #python heatmap data
    python3 $heatmapDataScript $folderHeatmapFiles $iter --norm --rank
    
    echo "Heatmap images done"
    
    rm $folderHeatmapFiles/*.csv
    sleep 2 

    avconv -loglevel quiet -y -r 4 -start_number 1 -i $folderHeatmapFiles/phenotypicDistance%d.png -b:v 1000k $folderHeatmapFiles/beh.mp4
avconv -loglevel quiet -y -r 4 -start_number 1 -i $folderHeatmapFiles/genotypicDistance%d.png -b:v 1000k $folderHeatmapFiles/gen.mp4
    rm $folderHeatmapFiles/*istance*.png
    
    echo "Heatmap video done"
    allVideoScript="../heatmapsRun.sh"
   
    $allVideoScript $folderHeatmapFiles $folderFiles $outbasename.png $videoName $iter $playVideo
    
    mkdir $outbasename/stats
    mv $folderHeatmapFiles/scatterRelGenoPheno.flv $outbasename/stats
    mv $folderHeatmapFiles/binsScatterRelGenoPheno.flv $outbasename/stats
    mv $folderHeatmapFiles/*.data $outbasename/stats
    mv $folderHeatmapFiles/*.png $outbasename/stats
    mv $folderFiles/$outbasename*.flv $outbasename/stats
    

    rm $folderFiles/*.png
    rm $folderHeatmapFiles/*.mp4

    rm $folderFiles/*.nn

    echo "Full video done"
    rm $folderFiles/*.dat
else
    echo "fnames += ["
    for i in $filenames ;
    do
	echo \"$i\",
    done
    echo "]"
fi