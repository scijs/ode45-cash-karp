set terminal postscript
set output "plot.ps"

set xtics font "Verdana,8"
set ytics font "Verdana,8"

set multiplot

set origin 0,0
set yrange [-1.1:1.1]
set size 1, 0.5
set notitle
set grid
set nokey
set xtics 1
set ytics 1
set style line 1 lc rgb '#0060ad' lt 1 lw 1 pt 6 pi -1 ps 0.5
set pointintervalbox 0.4
plot "rkck.dat" using 1:2 with linespoints ls 1

set origin 0,0.5

set size 1, 0.5
set yrange [-1.1:1.1]
set notitle
set grid
set nokey
set xtics 1
set ytics 1
set style line 1 lc rgb '#0060ad' lt 1 lw 1 pt 6 pi -1 ps 0.5
set pointintervalbox 0.4
plot "rk4.dat" using 1:2 with linespoints ls 1


set nomultiplot
