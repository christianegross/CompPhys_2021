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

datafile4='../data/magnetizationnmd4.dat'
correlationfile4='../data/magnetizationcorrelationnmd4.dat'

datafile100='../data/magnetizationnmd100.dat'
correlationfile100='../data/magnetizationcorrelationnmd100.dat'

set terminal pdfcairo size 4in, 4in
set out 'test.pdf'

set xlabel 'step'
set ylabel 'magnetization'
set xrange [0:500]
#plot datafile4 u 0:1 ls 1 ps 0.1 title 'N_{md}=4' ,datafile100 u 0:1 ls 2 ps 0.1 title 'N_{md}=100'
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
plot for [n=0:6] correlationfile100 every :::n::n u 1 ls n w lines title sprintf('binlength=%d', 2**n)

set title 'N_md=4'
plot for [n=0:6] correlationfile4 every :::n::n u 1 ls n w lines title sprintf('binlength=%d', 2**n)