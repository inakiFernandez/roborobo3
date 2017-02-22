set terminal png
set xrange [-1:1]
set yrange [-1:1]
set label datafile at -0.8,0.8


if (!exists("datafile2"))
{
plot datafile using 1:3 with lines, dataother using 1:3 with lines, f1 using 1:3 lt rgb "#FF0000" lw 3 with lines, f2 using 1:3 lt rgb "#0000FF" lw 3 with lines;
print "Single"
}
else
{
plot datafile using 1:3 lt rgb "#AA3300" with lines, datafile2 using 1:3 lt rgb "#0033AA" with lines, dataother using 1:3 lt rgb "#DD8800" with lines, dataother2 using 1:3 lt rgb "#0088DD" with lines, f1 using 1:3 lt rgb "#FF0000" lw 2 with lines, f2 using 1:3 lt rgb "#0000FF" lw 2 with lines 
print "Split"
}
set print
print datafile
