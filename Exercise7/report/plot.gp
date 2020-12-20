#Christiane Gross, Nico Dichter
#Exercise 6 Computational Physics WS20/21

set style line 1 lc 7 lt 7 pt 7 #ps 0.2
set style line 2 lc 1 lt 7 pt 8
set style line 3 lc 2 lt 7 pt 11
set style line 4 lc 3 lt 7 pt 3
set style line 5 lc 4 lt 7 pt 5#4 gut
set style line 6 lc 5 lt 7 pt 4#4 gut
set style line 7 lc 6 lt 7 pt 7#4 gut

set style line 8 lc 2 lt 7 pt 9
set style line 9 lc 3 lt 7 pt 10
set style line 10 lc 5 lt 7 pt 12

result3="../data/result3.dat"
result4="../data/result4.dat"

set terminal pdfcairo size 4in, 4in
set out 'test.pdf'

do for [angsize=4:60:8]{
set title sprintf("with pmax angsize=%d", angsize)
plot for [size=4:56:4] result3 u (($2==angsize&&$1==size)?$3:1/0):4 title sprintf("%d", size)
}
do for [angsize=4:60:8]{
set title sprintf("wo pmax angsize=%d", angsize)
plot for [size=4:56:4] result3 u (($2==angsize&&$1==size)?$3:1/0):7 title sprintf("%d", size)
}


set ter epslatex size 15 cm, 10 cm color colortext
unset title

set output
