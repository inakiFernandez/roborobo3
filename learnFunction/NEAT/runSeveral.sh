#!/bin/bash
start=`date +%s`

if [ $# -ne 8 ]; then
    echo "Wrong number of parameters. [progname] templateNotMulti templateMulti outconfigbasename nbIter nbRuns playVideo nbCores doBuildVideo"
    exit
fi
template=$1
templateMulti=$2
outbasename=$3
iter=$4
nbRuns=$5
playVideo=$6
nbCores=$7 

program="./neatLearnFunction "
nbExp=1

scriptGenerate="../cmdDatasets.sh"
doGenerate="true"

#Generate dataset every time (before experiments)
if [ $doGenerate == "true" ]; then
    echo "Generate dataset $outbasename"
    $scriptGenerate $template $outbasename
fi



doMulti="false"
if [ $doMulti = "true" ]; then
    echo "Do multi"
    nbExp=2
   
    rm $templateMulti 
    cp $template $templateMulti 

    #Add allowMulti 
    echo "allowMultisynapses=true" >> $templateMulti
fi

rm -rf $outbasename
rm -rf sandbox/$outbasename
rm -rf headmapData/$outbasename
mkdir -p $outbasename

mkdir $outbasename/noMulti
mkdir -p sandbox/$outbasename
mkdir -p heatmapData/$outbasename
echo "No multi"
touch $outbasename/noMulti.parallel
for (( j=1; j<=$nbRuns; j++))
do
    #Non paralelized version
    #$program $template $iter > $outbasename/noMulti/run-$j.log #oneRun $j &   
    echo "$program $template $iter $j-runFolder $outbasename > $outbasename/noMulti/run-$j.log; echo $j" >> $outbasename/noMulti.parallel
    
    mkdir -p sandbox/$outbasename/$j-runFolder
    mkdir -p heatmapData/$outbasename/$j-runFolder
done

echo "Running $nbRuns runs of $outbasename experiments in $nbCores processors for $iter generations"
#--results $outbasename/outExp
parallel -j $nbCores -a $outbasename/noMulti.parallel

folderFiles="./sandbox/$outbasename"
folderHeatmapFiles="./heatmapData/$outbasename"
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
        #TODO parallelized version
	$program $templateMulti $iter $j-runFolder $outbasename > $outbasename/multi/run-$j.log 
    done
    cp $templateMulti $outbasename
fi

cp $template $outbasename


sleep 2
echo "End run exp in parallel"
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

#TODO TOFIX Paste is not good with several columns per run
paste $outbasename/noMulti/*allFitnessAndOther.log >> $outbasename/noMultiallFitnessAndOtherAll.log; 

echo "End paste log files"

filenames=`realpath $outbasename/*.log` ;

#In case of variable-length experiments (e.g. stop condition on performance bound)
wcOut=($(wc -l $outbasename/noMultiBestAll.log))
iter=${wcOut[0]}
echo "Number of iterations: $iter"

    echo "Before last scripts"
    #python fitness
    fitnessDataScript="../plotFitnessInMultiExp.py"
    videoApproxScript="../gnuplotFolder.sh"
    heatmapDataScript="../heatmap.py"
    allVideoScript="../heatmapsRun.sh"
    #Extract duration per function
    durPerFunction=`grep "durationPerFunction= .*$" $template | sed -e "s/durationPerFunction= \(.*\)$/\1/"`
   
    echo "Before loop create data treatment parallel script"
    #rm  $outbasename/dataplots.parallel
    touch $outbasename/dataplots.parallel
    sleep 2;
    mkdir -p $folderHeatmapFiles/progressbar/
    convert -size 815x100 xc:white $folderHeatmapFiles/progressbar/empty.png
    xmax=800
    xinterval=$(($xmax / $iter))
    for ((j=1; j<=$iter; j++)); 
    do  
	xend=$(($xinterval * ($j + 1) + 15 )); 
#echo $xend; 
	it=$(($j))
	convert $folderHeatmapFiles/progressbar/empty.png -strokewidth 0 -fill "rgba( 0, 255, 50 )" -draw "rectangle 0,0 $xend,100 " -fill "rgba( 0, 0, 0 )" -gravity Center -weight 2200 -pointsize 50 -annotate 0 "It. $it"  -append $folderHeatmapFiles/progressbar/progressbar$j.jpg; 
    done;
    avconv -loglevel quiet  -y -r 4 -start_number 1 -i $folderHeatmapFiles/progressbar/progressbar%d.jpg -b:v 1000k  -vcodec mpeg4  $folderHeatmapFiles/progressbar/progressbar.mp4 
    
    rm $folderHeatmapFiles/progressbar/progressbar*.jpg
    isWithSplit="false"
    for (( j=1; j<=$nbRuns; j++))
    do
	rm ../tmpDataFilesExp$j.py
	touch ../tmpDataFilesExp$j.py
	printf "fnames = [" >> ../tmpDataFilesExp$j.py
	for i in $filenames ;
	do
	    printf \"$i\", >> ../tmpDataFilesExp$j.py
	done
	printf "]" >> ../tmpDataFilesExp$j.py
	
	parallelCommand="python3 $fitnessDataScript $folderFiles/$j-runFolder/$outbasename.png $nbExp --png $durPerFunction $j; echo 'Fitness Img done: run '$j; sleep 2; $videoApproxScript $folderFiles/$j-runFolder $isWithSplit true $j $outbasename; echo 'Video approx done: run '$j' ';  sleep 2;  python3 $heatmapDataScript $folderHeatmapFiles/$j-runFolder $iter --norm --rank $j; echo 'Heatmap images done: run '$j; rm $folderHeatmapFiles/$j-runFolder/*.csv; sleep 2; avconv -loglevel quiet -y -r 4 -start_number 1 -i $folderHeatmapFiles/$j-runFolder/phenotypicDistance%d.png -b:v 1000k $folderHeatmapFiles/$j-runFolder/beh.mp4; avconv -loglevel quiet -y -r 4 -start_number 1 -i $folderHeatmapFiles/$j-runFolder/genotypicDistance%d.png -b:v 1000k $folderHeatmapFiles/$j-runFolder/gen.mp4; rm $folderHeatmapFiles/$j-runFolder/*istance*.png; echo 'Heatmap videos done: run '$j; $allVideoScript $folderHeatmapFiles $folderFiles/$j-runFolder $outbasename.png $folderFiles/$j-runFolder/vidApprox.mp4 $iter $playVideo $j; sleep 3; mv $folderHeatmapFiles/$j-runFolder/scatterRelGenoPheno$j.flv $outbasename/stats; mv $folderHeatmapFiles/$j-runFolder/binsScatterRelGenoPheno$j.flv $outbasename/stats; mv $folderHeatmapFiles/$j-runFolder/*.data $outbasename/stats; mv $folderHeatmapFiles/$j-runFolder/*.png $outbasename/stats; mv $folderFiles/$j-runFolder/$outbasename.pngAll$j.flv $outbasename/stats; rm $folderFiles/$j-runFolder/*.png; rm $folderHeatmapFiles/$j-runFolder/*.mp4; rm $folderFiles/$j-runFolder/*.nn; echo 'Full video done: run '$j; rm $folderFiles/$j-runFolder/*.dat; rm -rf $folderFiles/$j-runFolder; rm -rf $folderHeatmapFiles/$j-runFolder; " 
	echo $parallelCommand >> $outbasename/dataplots.parallel  
    done
    sleep 2   

    mkdir $outbasename/stats
    echo "Before data treatment parallel"
    buildVideo=$8

if [ "$buildVideo" = true ] ; then
    #--results $outbasename/outData
    parallel -j $nbCores -a $outbasename/dataplots.parallel
 
else
    echo "TODO"
    echo "parallel -j $nbCores -a $outbasename/dataplots.parallel"
fi

end=`date +%s`

runtime=$((end-start))
echo "Done: $runtime seconds" 
