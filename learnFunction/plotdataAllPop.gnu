set terminal png 
#enhanced background rgb '#DCDCDC'
set xrange [-1:1]
set yrange [-1:1]
set style line 13 lc rgb 'white' lt 1 lw 1
set label datafile at -0.8,0.8
set grid ls 13
unset key
set object 1 rectangle from graph 0, graph 0 to graph 1, graph 1 behind fc rgbcolor '#DCDCDC' fs noborder

numberApprox = 10
 #(sizePop-1)
#print numberApprox

if (!exists("datafile2")){
plot for [i=0:numberApprox] database.'-ind'.i.'.sorted' using 1:3 with lines, datafile using 1:3 lt rgb "#AA3300" lw 3 with lines, dataother using 1:3 lt rgb "#DD8800" with lines, f1 using 1:3 lt rgb "#FF0000" lw 4 with lines, f2 using 1:3 lt rgb "#0000FF" lw 4 with lines, f1 using 1:3 lt rgb "#FF0000" lw 4 with lines
} else {
plot for [i=0:numberApprox] database.'-ind'.i.'aa.sorted' using 1:3 with lines, for [i=0:(sizePop-1)] database.'-ind'.i.'ab.sorted' using 1:3 with lines, datafile using 1:3 lt rgb "#AA3300" lw 3 with lines, datafile2 using 1:3 lt rgb "#0033AA" lw 3 with lines, dataother using 1:3 lt rgb "#DD8800" with lines, dataother2 using 1:3 lt rgb "#0088DD" with lines, f2 using 1:3 lt rgb "#0000FF" lw 4 with lines, f1 using 1:3 lt rgb "#FF0000" lw 4 with lines
}

set print
#print datafile
#print database
#print dataother
