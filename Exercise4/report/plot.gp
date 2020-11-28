#Christiane Gross, Nico Dichter
#Exercise 4 Computational Physics WS20/21
#script for plotting theoretical and experimental results of magnetization and energy

set style line 1 lc 7 lt 7 pt 7 #ps 0.2
set style line 2 lc 1 lt 7 pt 8
set style line 3 lc 2 lt 7 pt 11
set style line 4 lc 3 lt 7 pt 3
set style line 5 lc 4 lt 7 pt 5#4 gut
set style line 6 lc 5 lt 7 pt 4#4 gut

set style line 8 lc 2 lt 7 pt 11

chainfile4='../data/magnetizationnmd4.dat'
correlationfile4='../data/magnetizationcorrelationnmd4.dat'
datafile4='../data/raw4.dat'
errorstabilityfile4='../data/magnetizationerrorstabilitynmd100.dat'

chainfile100='../data/magnetizationnmd100.dat'
correlationfile100='../data/magnetizationcorrelationnmd100.dat'
datafile100='../data/raw100.dat'
errorstabilityfile100='../data/magnetizationerrorstabilitynmd100.dat'

set terminal pdfcairo size 4in, 4in
set out 'test.pdf'

set xlabel 'step'
set ylabel 'magnetization'
set xrange [0:500]
#plot chainfile4 u 0:1 ls 1 ps 0.1 title 'N_{md}=4' ,chainfile100 u 0:1 ls 2 ps 0.1 title 'N_{md}=100'
unset xrange
set ylabel 'Correlation'
#plot correlationfile4 u 0:1 ls 1 ps 0.1 title 'N_{md}=4', correlationfile100 u 0:1 ls 2 ps 0.1 title 'N_{md}=100'

plot 'test.txt' u 1 with lines dt 3 title 'everything', 'test.txt' every :::0::0 u 1 with lines dt 2 title 'zeroeth block', 'test.txt' every :::1::1 u 1 with lines title 'first block'

do for [n=0:6]{
set title sprintf('binlength=%d', 2**n)
#plot correlationfile100 every :::n::n u 1 ls 1 title ''
}


set xrange [0:100]
set title 'N_md=100'
plot for [n=0:6] correlationfile100 every :::n::n u 1 ls n w lp title sprintf('binlength=%d', 2**n)

set title 'N_md=4'
plot for [n=0:6] correlationfile4 every :::n::n u 1 ls n w lp title sprintf('binlength=%d', 2**n)

set ter epslatex size 15 cm, 10 cm color colortext
unset title
set out 'markovchaincomparison.tex'
set xlabel 'time'
set ylabel 'magnetization'
set xrange [0:500]
plot chainfile4 u 0:1 ls 1 title '$N_{md}=4$' ,chainfile100 u 0:1 ls 2 title '$N_{md}=100$'

set out 'simplecorrelation.tex'
set ylabel '$\Gamma(\tau)$'
plot correlationfile4 every :::0::0 u 1 ls 1 w lp title '$N_{md}=4$', correlationfile100 every :::0::0 u 1 ls 2 w lp title '$N_{md}=100$'

set xrange [0:100]
set out 'correlationbinnmd4.tex'
plot for [n=1:6] correlationfile4 every :::n::n u 1 ls n w lp title sprintf('binlength=%d', 2**n)

set xrange [0:50]
set out 'correlationbinnmd100.tex'
plot for [n=1:6] correlationfile100 every :::n::n u 1 ls n w lp title sprintf('binlength=%d', 2**n)

set out 'naiveerrorbinned.tex'
set key top left
set xtics 2,2,64
unset mxtics
set xlabel 'binlength'
set ylabel '$\langle m\rangle$'
set logscale x
set xrange [2:64]
plot datafile4 using 7:6 ls 1 w lp title '$N_{md}=4$', datafile100 using 7:6 ls 2 w lp title '$N_{md}=100$'

set out 'bootstraperrorbinned.tex'
plot datafile4 using 7:4 ls 1 w lp title '$N_{md}=4$', datafile100 using 7:4 ls 2 w lp title '$N_{md}=100$'


set out 'errorstabilitynmd100.tex'
set xlabel '$N_{bs}$'
set xrange [2:65536]
set xtics 2, 4, 65536
plot for [n=0:5] errorstabilityfile100 every :::n::n u 1:3 ls (n+1) w lp title sprintf('binlength=%d', 2**(n+1))
set output
