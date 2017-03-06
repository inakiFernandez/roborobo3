#!/bin/bash

if [ $# -ne 5 ]; then
    echo "Wrong number of parameters. [progname] folder withSplit withAllPop idRun nameExp"
    exit
fi
export LC_ALL=C
folder=$1
folderScript=`dirname $0`
withSplit=$2 # attention to gnu script
withAllPop=$3
numRun=$4
outbasename=$5

echo "Start $0"
gnuScript="$folderScript/plotdatapoints.gnu"
if [ "$withAllPop" = true ] ; then
    gnuScript="$folderScript/plotdataAllPop.gnu"
fi

f1="$folderScript/datasets/in-o1-$outbasename.sorted.dat"  
f2="$folderScript/datasets/in-o2-$outbasename.sorted.dat"
#rm $folder/*.png
#rm $folder/*.sorted
#rm $folder/*.dataa
#rm $folder/*.datab
#rm $folder/*ind*aa
#rm $folder/*ind*ab

filesBest=`ls $folder | grep '\.dat$' | grep '^g' | grep -v '\-ind'`
#(>&2 echo $filesBest)

# in $folder/g*.dat ; 
for i in $filesBest ; 
do
   othername="$folder/other$(basename "$i")"
   nbLines=`wc -l < $folder/$i`
   halfLines=$(($nbLines / 2))
   if [ "$withSplit" = true ] ; then
       split -l $halfLines $folder/$i $folder/$i   
       rm $folder/$i
       sort -g -t$'\t' -k 1,2 $folder/${i}aa > $folder/${i}aa.sorted
       sort -g -t$'\t' -k 1,2 $folder/${i}ab > $folder/${i}ab.sorted
       rm $folder/${i}aa
       rm $folder/${i}ab
       split -l $halfLines $othername $othername   
       sort -g -t$'\t' -k 1,2 ${othername}aa > ${othername}aa.sorted
       sort -g -t$'\t' -k 1,2 ${othername}ab > ${othername}ab.sorted
       rm ${othername}aa
       rm ${othername}ab
       if [ "$withAllPop" = true ] ; then
	   for indivFilename in $folder/${i}-ind* ;
	   do
	       split -l $halfLines $indivFilename $indivFilename
	       sort -g -t$'\t' -k 1,2 ${indivFilename}aa > ${indivFilename}aa.sorted
	       sort -g -t$'\t' -k 1,2 ${indivFilename}ab > ${indivFilename}ab.sorted
	       rm ${indivFilename}aa
	       rm ${indivFilename}ab
	   done
	   numberInd=`ls -F $folder/${i}-ind*aa.sorted | wc -l`
	   #echo $i

	   gnuplot -e "datafile='$folder/${i}aa.sorted'; datafile2='$folder/${i}ab.sorted'; dataother='${othername}aa.sorted' ; dataother2='${othername}ab.sorted'; f1='$f1' ; f2='$f2' ; sizePop='$numberInd' ; database='$folder/${i}" $gnuScript > $folder/$i.png  #2> /dev/null
       else
	   gnuplot -e "datafile='$folder/${i}aa.sorted'; datafile2='$folder/${i}ab.sorted'; dataother='${othername}aa.sorted' ; dataother2='${othername}ab.sorted'; f1='$f1' ; f2='$f2'" $gnuScript > $folder/$i.png #2> /dev/null
       fi
   else
       sort -g -t$'\t' -k 1,2 $folder/${i} > $folder/${i}.sorted
       rm $folder/${i}
       sort -g -t$'\t' -k 1,2 ${othername} > ${othername}.sorted
       rm ${othername}
       if [ "$withAllPop" = true ] ; then
	   for indivFilename in $folder/${i}-ind* ;
	   do
	        sort -g -t$'\t' -k 1,2 ${indivFilename} > ${indivFilename}.sorted
		rm ${indivFilename}
	   done
	   numberInd=`ls -F $folder/${i}-ind*.sorted | wc -l`
	   #echo $numberInd

	   gnuplot -e "datafile='$folder/${i}.sorted'; dataother='${othername}.sorted'; f1='$f1' ; f2='$f2' ; sizePop='$numberInd' ; database='$folder/${i}" $gnuScript > $folder/$i.png  #2> /dev/null
       else
	   gnuplot -e "datafile='$folder/${i}.sorted'; dataother='${othername}.sorted'; f1='$f1' ; f2='$f2'" $gnuScript > $folder/$i.png  #2> /dev/null  
       fi
   fi 
   #rm $folder/*.sorted
done

sleep 2;
listFiles=`ls  $folder/g*-t-0.dat.png`
for j in $listFiles ; #$folder/g*-t-0.dat.png ; 
do 
#echo $j
index=`echo $j | sed -e 's/g\(.*\)-t-0.dat.png/\1/'` ; 
mv "$j" "$index.png"; 
done

listFiles=`ls  $folder/g*-t-1.dat.png`
for j in $listFiles ; #$folder/g*-t-1.dat.png ; 
do 
#echo $j
index=`echo $j | sed -e 's/g\(.*\)-t-1.dat.png/\1/'` ; 
mv "$j" "$index.png"; 
done

cd $folder
find -maxdepth 1 -type f -name 'g*-t-*.dat-ind*' -delete # -exec rm {} \;
find -maxdepth 1 -type f -name '*.sorted' -delete #-exec rm {} \;
cd - > /dev/null
#rm $folder/*.sorted
#rm $folder/g*-t-*.dat-ind*


#DATE=$(date +"%Y%m%d%H%M")

avconv -loglevel quiet  -y -r 4 -start_number 1 -i $folder/%d.png -b:v 1000k $folder/vidApprox.mp4
echo $folder/vidApprox.mp4 

#avconv -i beh.mp4 -i gen.mp4 -i ../sandbox/201701181431.mp4 -filter_complex "[0:v]pad=iw*2:ih*2[a]; [a][1:v]overlay=w; [2:v]overlay=w:ih[bg]" allInOne2.mp4

#-b 800000 
 #NOTE: -b is the bitrate, if not -b Number, no compression

#mogrify -crop 1736x1048+190+100 +repage ~/sandbox/test50foo3.png

# avconv -i beh.mp4 -i gen.mp4 -i ../sandbox/201701181431.mp4 -i ../../../../../sandbox/test50foo3.png -i progressbar/progressbar.mp4 -filter_complex "[0:v]pad=iw*2:ih*2[a]; [a][1:v]overlay=w[x]; [2:v]scale=iw*1.23:ih*1.23[y]; [x][y]overlay=0:h+10[z]; [3:v] scale=iw*.4685:ih*.4685[w]; [z][w] overlay=w-26:h+110[v]; [v][4:v] overlay=main_w-w:main_h-h-10 [out]" -map "[out]" allInOne7.mp4
