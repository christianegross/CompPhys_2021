#Christiane Gross, Nico Dichter
#Exercise 5 Computational Physics WS20/21
#script for plotting theoretical and experimental results of magnetization and energy

set style line 1 lc 7 lt 7 pt 7 #ps 0.2
set style line 2 lc 1 lt 7 pt 8
set style line 3 lc 2 lt 7 pt 11
set style line 4 lc 3 lt 7 pt 3
set style line 5 lc 4 lt 7 pt 5#4 gut
set style line 6 lc 5 lt 7 pt 4#4 gut

set style line 8 lc 2 lt 7 pt 11

vfile="../data/vcycle.dat"
wfile="../data/wcycle.dat"

set terminal pdfcairo size 4in, 4in
set out 'test.pdf'

set xlabel "MC time"
set ylabel "Gamma"
plot vfile using 1  with linespoints linestyle 1 title "v-cycle", wfile using 1  with linespoints linestyle 2 title "w-cycle"
set xrange [0:50]
plot vfile using 1  with linespoints linestyle 1 title "v-cycle", wfile using 1  with linespoints linestyle 2 title "w-cycle"
set xrange [0:200]
plot vfile using 1  with linespoints linestyle 1 title "v-cycle", wfile using 1  with linespoints linestyle 2 title "w-cycle"
set xrange [500:600]
plot vfile using 1  with linespoints linestyle 1 title "v-cycle", wfile using 1  with linespoints linestyle 2 title "w-cycle"



set ter epslatex size 15 cm, 10 cm color colortext
unset title
set output 'autocorrelationzoomedin.tex'
set xlabel 'MC time $\tau$'
set ylabel '$\Gamma(\tau)$'
set xrange [0:60]
plot vfile using 1  with linespoints linestyle 1 title 'v-cycle', wfile using 1  with linespoints linestyle 2 title 'w-cycle'
set output
