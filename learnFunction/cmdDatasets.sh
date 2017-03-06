#!/bin/bash
#TOTEST add ../ to run from learnFunction/NEAT,  get name-files from argument
if [ $# -ne 2 ]; then
    echo "Wrong number of parameters. [progname] templateNotMulti outconfigbasename"
   exit
fi
template=$1
outbasename=$2

echo "Generating $outbasename dataset from $template config file"
g++ -I../../include/contrib/ -I../../include/contrib/odneatgc/ -I.. -I../../include/core/Utilities/ -I../../include/contrib/zsu ../generateDataset.cpp -o ../generateDataset -std=c++0x ../../src/contrib/odneatgc/*.o ../../src/core/ExtendedProperties.o ../../src/contrib/zsu/*.o ../../src/core/Misc.o -L ../build/ -l auxlearnfunction -g -DODNEAT_FUNCTIONS

#generatedataset with outbasename

../generateDataset $template $outbasename

 python ../../tools/puregenome_draw.py -f ../datasets/o1-$outbasename.nn -d ../datasets/o1-$outbasename.dot -p ../datasets/o1-$outbasename.png
python ../../tools/puregenome_draw.py -f ../datasets/o2-$outbasename.nn -d ../datasets/o2-$outbasename.dot -p ../datasets/o2-$outbasename.png
#eog ../datasets/o1-$outbasename.png &
#eog ../datasets/o2-$outbasename.png &

rm ../datasets/in-o1-$outbasename.dat ../datasets/in-o2-$outbasename.dat 

paste ../datasets/in-$outbasename.dat ../datasets/o1-$outbasename.dat > ../datasets/in-o1-$outbasename.dat
paste ../datasets/in2-$outbasename.dat ../datasets/o2-$outbasename.dat > ../datasets/in-o2-$outbasename.dat

sort -g -t$'\t' -k 1,2 ../datasets/in-o1-$outbasename.dat > ../datasets/in-o1-$outbasename.sorted.dat
sort -g -t$'\t' -k 1,2 ../datasets/in-o2-$outbasename.dat > ../datasets/in-o2-$outbasename.sorted.dat

rm ../datasets/in-o1-$outbasename.png 
rm ../datasets/in-o2-$outbasename.png 
rm ../datasets/in-both-$outbasename.png

gnuplot -e "filename='../datasets/in-o1-$outbasename.sorted.dat'" ../plot_function.plg > ../datasets/in-o1-$outbasename.png
gnuplot -e "filename='../datasets/in-o2-$outbasename.sorted.dat'" ../plot_function.plg > ../datasets/in-o2-$outbasename.png
gnuplot -e "filename1='../datasets/in-o1-$outbasename.sorted.dat'; filename2='../datasets/in-o2-$outbasename.sorted.dat'" ../plot_2_functions.plg > ../datasets/in-both-$outbasename.png

###################################

#eog ../datasets/in-o1-$outbasename.png &
#eog ../datasets/in-o2-$outbasename.png &
#eog ../datasets/in-both-$outbasename.png &

###################################

rm ../datasets/oAll-$outbasename.dat
rm ../datasets/inAll-$outbasename.data
cp ../datasets/in-o1-$outbasename.dat ../datasets/in-o-$outbasename.dat
cat ../datasets/in-o2-$outbasename.dat >> ../datasets/in-o-$outbasename.dat
cp ../datasets/in-$outbasename.dat ../datasets/inAll-$outbasename.data
cat ../datasets/in2-$outbasename.dat >> ../datasets/inAll-$outbasename.data
cp ../datasets/o1-$outbasename.dat ../datasets/oAll-$outbasename.dat
cat ../datasets/o2-$outbasename.dat >> ../datasets/oAll-$outbasename.dat



#Plot all networks in folder# To adapt to current folder
#for i in *.nn; do python ../../../tools/puregenome_draw.py -f $i -d $i.dot -p $i.png; done

#Invert dataset 1D (symmetric over x axis)
#while read in; 
#do 
#echo "-1*$in" | bc |  awk '{printf "%f\n", $0}' ; 
#done < datasets/o1-template-params.inverted.dat > datasets/o1-template-params.inverted.dat2
#sandbox plot of approximated output
#fname='NEAT/sandbox/g80-1.dat'
#fSname'NEAT sandbox/g80-1.sorted.dat'
#rm splotCmds.txt
#echo "splot 'datasets/in-o1.sorted.dat', 'datasets/in-o2.sorted.dat', "
#for i in NEAT/sandbox/*.dat ;
#do
   #sort -g -t$'\t' -k 1,2 $i > $i.sorted.dats
#echo "splot 'datasets/in-o1.sorted.dat', '$i.sorted.dats' "
#echo "'$i.sorted.dats', "
#done > splotCmds.txt
