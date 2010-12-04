set key left top 
set border 15 lw 1.5
set style line 1 lt 1 lw 8 pt 1 ps 2.5 
set size ratio 1 
set xlabel "Magentic Field (T)"
set ylabel "Moment (mb/Nd)" 
plot "moment.clc"  using 1:2 title "moment.rtplot"  with linesp ls 1 
pause 2
set terminal postscript color eps
set output "moment.eps" 
replot