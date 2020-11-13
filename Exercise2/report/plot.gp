#Christiane Gross, Nico Dichter
#Exercise 2 Computational Physics WS20/21
#script for plotting theoretical and experimental results of magnetization and energy

set style line 1 lc 7 lt 7 pt 7 #ps 0.2
set style line 2 lc 1 lt 7 pt 8
set style line 3 lc 2 lt 7 pt 11
set style line 4 lc 3 lt 7 pt 3
set style line 5 lc 4 lt 7 pt 5#4 gut
set style line 6 lc 5 lt 7 pt 4#4 gut

set style line 8 lc 2 lt 7 pt 11

magnetization(J)=(1-1/(sinh(2*J))**4)**(1.0/8.0)
epsilonfile="epsilonexpected.txt"

set ter pdfcairo size 4in, 4in
#pdf for easy viewing

set out "energy.pdf"
set title "<E>"
set xlabel "J"
set ylabel "<e>"
set xrange [0.25:2]

plot epsilonfile u 1:2 w lines ls 1 title "Literatur"

set out "magnetization.pdf"
set title "J fix"
set xlabel "h"
set ylabel "<m>"
set xrange [-1:1]


set title "h=0"
set xlabel "J"
set ylabel "<|m|>"
set y2label "<m>"
set xrange [0.25:1]
set key bottom right

plot magnetization(x) ls 1 title "Literatur"

set ter epslatex size 15 cm, 10.6cm color colortext
#.tex for easy inputting in report

set output
