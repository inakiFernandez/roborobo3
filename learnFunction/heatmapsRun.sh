#!/bin/bash

if [ $# -ne 6 ]; then
    echo "Wrong number of parameters. $0 folderH folderS errorCurveImage approxSolutionsVideo nbIterations"
    exit
fi
folderH=$1
folderS=$2
folderScript=`dirname $0`
errorCurveFile=$3
approxSolutionsVideo=$4
numberGen=$5
playVideo=$6

#avconv -loglevel quiet -y -r 4 -start_number 1 -i $folderH/behavioralDistance%d.png -b:v 1000k $folderH/beh.mp4
#avconv -loglevel quiet -y -r 4 -start_number 1 -i $folderH/genotypicDistance%d.png -b:v 1000k $folderH/gen.mp4
#avconv -i beh.mp4 -i gen.mp4 -filter_complex "[0:v:0]pad=iw*2:ih[bg]; [bg][1:v:0]overlay=w" allInOne.mp4
#avconv -i beh.mp4 -i gen.mp4 -i ../sandbox/201701181431.mp4 -filter_complex "[0:v]pad=iw*2:ih*2[a]; [a][1:v]overlay=w; [2:v]overlay=w:ih[bg]" allInOne2.mp4
#NOTE: -b is the bitrate, if not -b Number, no compression

mogrify -crop 1736x1048+190+100 +repage $folderS/$errorCurveFile
mkdir -p $folderH/progressbar/

convert -size 815x100 xc:white $folderH/progressbar/empty.png

xmax=800

xinterval=$(($xmax / $numberGen))
for ((j=1; j<=$numberGen; j++)); 
do  
xend=$(($xinterval * ($j + 1) + 15 )); 
#echo $xend; 
it=$(($j))
convert $folderH/progressbar/empty.png -strokewidth 0 -fill "rgba( 0, 255, 50 )" -draw "rectangle 0,0 $xend,100 " -fill "rgba( 0, 0, 0 )" -gravity Center -weight 2200 -pointsize 50 -annotate 0 "It. $it"  -append $folderH/progressbar/progressbar$j.jpg; 
done;
avconv -loglevel quiet  -y -r 4 -start_number 1 -i $folderH/progressbar/progressbar%d.jpg -b:v 1000k  -vcodec mpeg4 -t 100  $folderH/progressbar/progressbar.mp4 #TOFIX maybe remove -t 100
rm $folderH/progressbar/progressbar*.jpg


avconv -loglevel quiet -y -i $folderH/beh.mp4 -i $folderH/gen.mp4 -i $approxSolutionsVideo -i $folderS/$errorCurveFile -i $folderH/progressbar/progressbar.mp4 -filter_complex "[0:v]pad=iw*2:ih*2[a]; [a][1:v]overlay=w[x]; [2:v]scale=iw*1.23:ih*1.23[y]; [x][y]overlay=0:h+10[z]; [3:v] scale=iw*.4685:ih*.4685[w]; [z][w] overlay=w-26:h+110[v]; [v][4:v] overlay=main_w-w:main_h-h-10 [out]" -map "[out]" $folderS/${errorCurveFile}All.mp4 #-b 800000
avconv -loglevel quiet -y -i $folderS/${errorCurveFile}All.mp4 $folderS/${errorCurveFile}All.flv

avconv -loglevel quiet -y -r 4 -start_number 1 -i $folderH/scatterRelGenoPheno%d.png -b:v 1000k -vcodec mpeg4 -t 100 $folderH/scatterRelGenoPheno.mp4
avconv -loglevel quiet -y -i $folderH/scatterRelGenoPheno.mp4 $folderH/scatterRelGenoPheno.flv 

rm $folderH/scatterRelGenoPheno*.png

avconv -loglevel quiet -y -r 4 -start_number 1 -i $folderH/binsScatterRelGenoPheno%d.png -b:v 1000k -vcodec mpeg4 -t 100 $folderH/binsScatterRelGenoPheno.mp4
avconv -loglevel quiet -y -i $folderH/binsScatterRelGenoPheno.mp4 $folderH/binsScatterRelGenoPheno.flv 

rm $folderH/binsScatterRelGenoPheno*.png

if [ $playVideo = "true" ]; then
    vlc --fullscreen $folderS/${errorCurveFile}All.flv
fi
