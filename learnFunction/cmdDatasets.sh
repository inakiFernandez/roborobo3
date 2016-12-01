g++ -I../include/contrib/ -I. -I../include/core/Utilities/ -I../include/contrib/zsu generateDataset.cpp -o generateDataset -std=c++0x ../src/contrib/odneatgc/*.o ../src/core/ExtendedProperties.o ../src/contrib/zsu/*.o ../src/core/Misc.o -L build/ -l auxlearnfunction -g
./generateDataset config/template-params
 python ../tools/puregenome_draw.py -f datasets/o1-template-params.nn -d datasets/o1-template-params.dot -p datasets/o1-template-params.png
python ../tools/puregenome_draw.py -f datasets/o2-template-params.nn -d datasets/o2-template-params.dot -p datasets/o2-template-params.png
eog datasets/o1-template-params.png &
eog datasets/o2-template-params.png &
rm datasets/in-o*
paste datasets/in-template-params.dat datasets/o1-template-params.dat > datasets/in-o1.dat
paste datasets/in2-template-params.dat datasets/o2-template-params.dat > datasets/in-o2.dat

 sort -g -t$'\t' -k 1,2 datasets/in-o1.dat >datasets/in-o1.sorted.dat
 sort -g -t$'\t' -k 1,2 datasets/in-o2.dat >datasets/in-o2.sorted.dat

rm datasets/in-o1.png 
rm datasets/in-o2.png 
rm datasets/in-both.png
gnuplot -e "filename='datasets/in-o1.sorted.dat'" plot_function.plg > datasets/in-o1.png

gnuplot -e "filename='datasets/in-o2.sorted.dat'" plot_function.plg > datasets/in-o2.png


gnuplot -e "filename1='datasets/in-o1.sorted.dat'; filename2='datasets/in-o2.sorted.dat'" plot_2_functions.plg > datasets/in-both.png

eog datasets/in-o1.png &
eog datasets/in-o2.png &
eog datasets/in-both.png &

###################################
rm datasets/oAll.dat
rm datasets/inAll.data
cp datasets/in-o1.dat datasets/in-o.dat
cat datasets/in-o2.dat >> datasets/in-o.dat
cp datasets/in-template-params.dat datasets/inAll.data
cat datasets/in2-template-params.dat >> datasets/inAll.data
cp datasets/o1-template-params.dat datasets/oAll.dat
cat datasets/o2-template-params.dat >> datasets/oAll.dat



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


