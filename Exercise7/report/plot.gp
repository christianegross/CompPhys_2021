#Christiane Gross, Nico Dichter
#Exercise 7 Computational Physics WS20/21

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

PI=3.14159

result3="../data/result3.dat"
result4="../data/result4.dat"
result5="../data/result5.dat"

set terminal pdfcairo size 4in, 4in
set out 'test.pdf'

#do for [angsize=4:60:8]{
#set title sprintf("with pmax angsize=%d", angsize)
#plot for [size=4:56:4] result3 u (($2==angsize&&$1==size)?$3:1/0):6 title sprintf("%d", size)
#}
#do for [angsize=4:60:8]{
#set title sprintf("wo pmax angsize=%d", angsize)
#plot for [size=4:36:4] result3 u (($2==angsize&&$1==size)?$3:1/0):9 title sprintf("%d", size)
#}

#plot result4 u 2:4, result4 u 2:6

set ter epslatex size 15 cm, 10 cm color colortext
unset title

set xlabel '$p_{max}/\si{\per\femto\meter}$
set ylabel '$|t_0(q,q)|/\si{\femto\meter^2}$'
set key top rmargin title 'Grid points'

set out 'tnnsmallangsize.tex'
plot for [size=4:56:8] result3 u (($2==4&&$1==size)?$3:1/0):6 ls ((size+4)/8) title sprintf("%d", size)

set out 'tnnbigangsize.tex'
plot for [size=4:56:8] result3 u (($2==60&&$1==size)?$3:1/0):6 ls ((size+4)/8) title sprintf("%d", size)

set out 'tnnsmallangsizewopm.tex'
plot for [size=4:56:8] result3 u (($2==4&&$1==size)?$3:1/0):9 ls ((size+4)/8) title sprintf("%d", size)

set out 'tnnbigangsizewopm.tex'
plot for [size=4:56:8] result3 u (($2==60&&$1==size)?$3:1/0):9 ls ((size+4)/8) title sprintf("%d", size)

set key inside top right title ''

set output 'delta.tex'
set xlabel '$q/\si{\per\femto\meter}$
set ylabel '$\delta(q)$'
set ytics ('$\pi/8$' PI/8,'$\pi/4$' PI/4,'$3\pi/8$' 3*PI/8, '$\pi/2$' PI/2)
set yrange [0:3*PI/8]
plot result4 u 2:($4/2) ls 1 title 'with variation of $p_{max}$'

unset yrange
unset ytics
set output 'crossect.tex'
set ylabel '$\derivative{\sigma}{\hat{q}_f}/\si{\femto\meter^2\mega\electronvolt}$
set xlabel '$x$'
set key inside top left title 'l='
plot for [l=0:6:1] result5 u (($2==l)?$1:1/0):3 ls ((l+1)) title sprintf("%d", l)