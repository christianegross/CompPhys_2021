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

magnetization(J)=0*(J<0.4406)+(1-1/(sinh(2*J))**4)**(1.0/8.0)*(J>0.4406)
epsilonfile='epsilonexpected.txt'
datafile='../data/N_J_h.dat'
#datafile N J h m m_err e e_err
absfile='../data/abs_m.dat'
#absfile N J m m_err 

set ter pdfcairo size 4in, 4in
#pdf for easy viewing

set out 'energy.pdf'
set title '<E>'
set xlabel 'J'
set ylabel '<e>'
set xrange [0.25:2]

plot epsilonfile u 1:2 w lines ls 1 title 'Literatur', for [N=4:20:4] datafile u (($1==N&&$3==0)?$2:1/0):6:7 w yerrorbars ls (N/4+1) title sprintf('N=%d', N)

set out 'magnetization.pdf'
set key top left
set xlabel 'h'
set ylabel '<m>'
set xrange [-1:1]
set title 'J =0.3'

plot for [N=4:20:4] datafile u (($1==N&&$2==0.3)?$3:1/0):4:5 w yerrorbars ls (N/4+1) title sprintf('N=%d', N)

set title 'J =0.8'

plot for [N=4:20:4] datafile u (($1==N&&$2==0.8)?$3:1/0):4:5 w yerrorbars ls (N/4+1) title sprintf('N=%d', N)

set title 'h=0'
set xlabel 'J'
set ylabel '<|m|>'
set y2label '<m>'
set link y2
set xrange [0.25:1]
set key bottom right

plot magnetization(x) ls 1 title 'Literatur', for [N=4:20:4] datafile u (($1==N&&$3==0)?$2:1/0):4:5 w yerrorbars ls (N/4+1) title sprintf('N=%d', N), for [N=4:20:4] absfile u (($1==N)?$2:1/0):3:4 w yerrorbars lc (N/4+7) title sprintf('N=%d, m_abs', N)
set title 'magnetization'
plot magnetization(x) ls 1 title 'Literatur', for [N=4:20:4] datafile u (($1==N&&$3==0)?$2:1/0):4:5 w yerrorbars ls (N/4+1) title sprintf('N=%d', N)

set title 'abs(magnetization)'
plot magnetization(x) ls 1 title 'Literatur', for [N=4:20:4] absfile u (($1==N)?$2:1/0):3:4 w yerrorbars lc (N/4+7) title sprintf('N=%d, m_abs', N)

unset y2label

clear
set ter epslatex size 15 cm, 10.6cm color colortext
#.tex for easy inputting in report

unset title
set out 'avenergy.tex'
set title ''
set xlabel '$J$'
set ylabel '$\langle\epsilon\rangle$'
set xrange [0.25:2]
set key top right

plot epsilonfile u 1:2 w lines ls 1 title 'Literatur', for [N=4:20:4] datafile u (($1==N&&$3==0)?$2:1/0):6:7 w yerrorbars ls (N/4+1) ps 2 title sprintf('$N=$%d', N)

set out 'magnetizationfixJ.tex'

set key top left
set xlabel '$h$'
set ylabel '$\langle m\rangle$'

#set multiplot layout 1,2
set xrange [-1:1]
set title '$J=0.3$'

plot for [N=4:20:4] datafile u (($1==N&&$2==0.3)?$3:1/0):4:5 w yerrorbars ls (N/4+1) ps 2 title sprintf('$N=$%d', N)

#set title '$J=0.8$'

#plot for [N=4:20:4] datafile u (($1==N&&$2==0.8)?$3:1/0):4:5 w yerrorbars ls (N/4+1) title sprintf('$N=$%d', N)

unset title
#unset multiplot

set xrange [0.25:1]
set key bottom right
set out 'magnetizationh0.tex'
set xlabel '$J$'
set xrange [0.25:1]
set key bottom right

plot magnetization(x) ls 1 title 'Literatur', for [N=4:20:4] datafile u (($1==N&&$3==0)?$2:1/0):4:5 w yerrorbars ls (N/4+1) ps 2 title sprintf('$N=$%d', N)

set out 'absmagnetization.tex'
set ylabel '$\langle |m|\rangle$'

plot magnetization(x) ls 1 title 'Literatur', for [N=4:20:4] absfile u (($1==N)?$2:1/0):3:4 w yerrorbars ls (N/4+1) ps 2 title sprintf('$N=$%d', N)


set output
