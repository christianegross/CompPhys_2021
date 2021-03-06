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

plot epsilonfile u 1:2 w lines ls 1 title 'theory', for [N=4:20:4] datafile u (($1==N&&$3==0)?$2:1/0):6:7 w yerrorbars ls (N/4+1) title sprintf('N=%d', N)

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

plot magnetization(x) ls 1 title 'theory', for [N=4:20:4] datafile u (($1==N&&$3==0)?$2:1/0):4:5 w yerrorbars ls (N/4+1) title sprintf('N=%d', N), for [N=4:20:4] absfile u (($1==N)?$2:1/0):3:4 w yerrorbars lc (N/4+7) title sprintf('N=%d, m_abs', N)
set title 'magnetization'
plot magnetization(x) ls 1 title 'theory', for [N=4:20:4] datafile u (($1==N&&$3==0)?$2:1/0):4:5 w yerrorbars ls (N/4+1) title sprintf('N=%d', N)

set title 'abs(magnetization)'
plot magnetization(x) ls 1 title 'theory', for [N=4:20:4] absfile u (($1==N)?$2:1/0):3:4 w yerrorbars lc (N/4+7) title sprintf('N=%d, m_abs', N)

unset y2label

set out 'heatcap.pdf'
set ylabel 'C/J^2'
set xlabel 'J^-1'
set key top right
set xrange [1:4]
plot epsilonfile u (1/$1):3 w lines ls 1 title 'theory', for [N=4:20:4] datafile u (($1==N&&$3==0)?1/$2:1/0):($7*$7*N*N/$2/$2) ls (N/4+1) title sprintf('N=%d', N)

clear
set ter epslatex size 15 cm, 10cm color colortext
#.tex for easy inputting in report

unset title
set out 'avenergy.tex'
set title ''
set xlabel '$J$'
set ylabel '$\langle\epsilon\rangle$'
set xrange [0.25:2]
set key top right

plot epsilonfile u 1:2 w lines ls 1 title 'theory', for [N=4:20:4] datafile u (($1==N&&$3==0)?$2:1/0):6:7 w yerrorbars ls (N/4+1) title sprintf('$N=$%d', N)

set xrange [0.25:1]
set key top right
set out 'magnetizationh0.tex'
set xlabel '$J^{-1}$'
set xrange [1:4]
set key top right

plot magnetization(1/x) ls 1 title 'theory', for [N=4:20:4] datafile u (($1==N&&$3==0)?(1/$2):1/0):4:5 w yerrorbars ls (N/4+1) title sprintf('$N=$%d', N)

set out 'absmagnetization.tex'
set ylabel '$\langle |m|\rangle$'

plot magnetization(1/x) ls 1 title 'theory', for [N=4:20:4] absfile u (($1==N)?(1/$2):1/0):3:4 w yerrorlines ls (N/4+1) title sprintf('$N=$%d', N)

set ter epslatex size 15 cm, 21cm color colortext

set out 'magnetizationfixJ.tex'

set key top left
set xlabel '$h$'
set ylabel '$\langle m\rangle$'

set multiplot layout 2,1
set xrange [-1:1]
set title '$J=0.3$'

plot for [N=4:20:4] datafile u (($1==N&&$2==0.3)?$3:1/0):4:5 w yerrorbars ls (N/4+1) title sprintf('$N=$%d', N)

set title '$J=0.8$'

plot for [N=4:20:4] datafile u (($1==N&&$2==0.8)?$3:1/0):4:5 w yerrorbars ls (N/4+1) title sprintf('$N=$%d', N)

unset title
unset multiplot

set ter epslatex size 15 cm, 9.5cm color colortext

set tmargin at screen 0.99
set out 'heatcapacity.tex'
set ylabel '$C/J^2$'
set xlabel '$J^{-1}$'
set key top right
set xrange [1:4]
plot epsilonfile u (1/$1):3 w lines ls 1 title 'theory', for [N=4:20:4] datafile u (($1==N&&$3==0)?1/$2:1/0):($7*$7*N*N/$2/$2) ls (N/4+1) title sprintf('$N=$%d', N)

set output
