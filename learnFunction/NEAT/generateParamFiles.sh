#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Wrong number of parameters. [progname] template-file  outbasename "
    exit
fi
#Generate parameter files for combinations of parameters from a template neat file

template=$1
outbasename=$2

# run with: ./generateParamFiles.sh ../config/template-params params-basename 


#"disjoint_coeff" "excess_coeff" "mutdiff_coeff" 
#"recur_prob" "recur_only_prob" "mutate_only_prob"
paramNames=("weigh_mut_power" "compat_thresh" "mutate_link_weights_prob" "mutate_toggle_enable_prob" "mutate_add_node_prob" "mutate_add_link_prob" "pop_size")


listVal=("0.1" "0.75" "0.6" "0.0" "0.0" "0.0" "100") #"40 100") 

paramEqNames=("durationPerFunction=" "loadGenome=" "randomTopoInit=" "injectGenome=" "initTopoN" "initTopoL" "probbiggaussian")
listEqVal=("100000" "false" "true" "false" "50" "150" "0.0")

commandParallel="parallel echo "
for((i=0;i<${#paramNames[@]};i++))
do
    pName=${paramNames[$i]}
    #echo "Name: "$pName
    pValues=${listVal[$i]}
    commandParallel="$commandParallel ::: $pName ::: $pValues"
done 


for((i=0;i<${#paramEqNames[@]};i++))
do
    pName=${paramEqNames[$i]}
    #echo "Name: "$pName
    pValues=${listEqVal[$i]}
    commandParallel="$commandParallel ::: $pName ::: $pValues"
done 

#echo $commandParallel

eval $commandParallel > $outbasename.paramLists

tmpFile="tmpFile.param"
tmpFile2="tmpFile2.param"
folderParams="paramFiles"
rm -rf $folderParams
mkdir -p $folderParams
listFilenames=""

while read line; do
    touch $tmpFile
    touch $tmpFile2
    cp $template $tmpFile
    #echo "Line: " $line
    nbWords=`echo $line | wc -w`
    IFS=', ' read -r -a array <<< "$line"
    filename="$folderParams/$outbasename"
    for((i=0;i<$nbWords;i++))
    do
	#echo "Param: ("${array[i]} ", "${array[$(( $i + 1 ))]}")"
	sedRegEx="s/${array[i]} .*$/${array[$i]} ${array[$(( i + 1 ))]}/"
	#echo "Sed regex:" $sedRegEx

	cat $tmpFile | sed -e "$sedRegEx" >  $tmpFile2
	cp $tmpFile2 $tmpFile
	#echo "File:" $tmpFile2
	
	initial="$(echo ${array[i]} | head -c 1)"
	filename="$filename-$initial${array[$(( $i + 1 ))]}"
	i=$(( i + 1 ))	
    done    

    filename="$filename.params"
    cp $tmpFile2 $filename
    listFilenames="${listFilenames}\n${filename}"
    rm $tmpFile
    rm $tmpFile2
done < $outbasename.paramLists

echo -e $listFilenames > $outbasename.paramFiles
