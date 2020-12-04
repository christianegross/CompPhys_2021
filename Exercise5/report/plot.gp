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


set terminal pdfcairo size 4in, 4in
set out 'test.pdf'

set ter epslatex size 15 cm, 10 cm color colortext
unset title
set output