#Christiane Gross, Nico Dichter
#Project Computational Physics WS20/21

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


beta=4.2
dim=2

wilsonloopfile=sprintf("../data/wilsonsu%dbeta%.3f.dat", dim, beta)
potentialfile=sprintf("../data/potentialsu%dbeta%.3f.dat", dim, beta)
potentialtotalfile="../data/potentialtotal.dat"
plaquettetotalfile="../data/plaquettetotal.dat"

set terminal pdfcairo size 4in, 4in
set out "test.pdf"

set mytics 4
set mxtics 2

set xrange [0:50]

set out "comparisoncreutz.pdf"

set title "plaquette, measured during sweep"
plot for [betap=12:32:4] sprintf("../data/beta%d.dat", betap) using 1:($3) linestyle ((betap-8)/4) title sprintf("%.1f", betap/10.0)

set title "1-plaquette, measured during sweep"
plot for [betap=12:32:4] sprintf("../data/beta%d.dat", betap) using 1:(1-$3) linestyle ((betap-8)/4) title sprintf("%.1f", betap/10.0) 
	 
set title "plaquette, measured after sweep"
plot for [betap=12:32:4] sprintf("../data/beta%d.dat", betap) using 1:($4) linestyle ((betap-8)/4) title sprintf("%.1f", betap/10.0) 

set title "1-plaquette, measured after sweep"
plot for [betap=12:32:4] sprintf("../data/beta%d.dat", betap) using 1:(1-$4) linestyle ((betap-8)/4) title sprintf("%.1f", betap/10.0)	 

set out "su3plaquette.pdf"
set title "hot start, epsilon=0.3"
set xlabel "steps"
set ylabel "<P>"
unset xrange
set yrange[0.4:1]
plot "../data/su3beta5p5fixed.dat" u 1:3 ls 1 w lines title "beta=5.5", "../data/su3beta5p5fixed.dat" u 1:4 w lines ls 2 title "beta=5.5"
set xrange [0:50]
plot "../data/su3beta5p5fixed.dat" u 1:3 ls 1 w lines title "beta=5.5", "../data/su3beta5p5fixed.dat" u 1:4 w lines ls 2 title "beta=5.5"


v1(x)=a1*exp(-b1*x)
v2(x)=a2*exp(-b2*x)
v3(x)=a3*exp(-b3*x)
v4(x)=a4*exp(-b4*x)

#a3=0.9
#a3_err=0.1
#b3=5
#b3_err=0.1

#a4=0.9
#b4=6
#a4_err=0.1
#b4_err=0.1

fit v1(x) wilsonloopfile u (($1==1)&&($5==32)?$2:1/0):3:4 yerrors via a1, b1
fit v2(x) wilsonloopfile u (($1==2)&&($5==32)?$2:1/0):3:4 yerrors via a2, b2
fit v3(x) wilsonloopfile u (($1==3)&&($5==32)?$2:1/0):3:4 yerrors via a3, b3
fit v4(x) wilsonloopfile u (($1==4)&&($5==32)?$2:1/0):3:4 yerrors via a4, b4


set print potentialfile
print(sprintf("R\ta\ta_err\tb\tb_err"))
print(sprintf("1\t%f\t%f\t%f\t%f",a1, a1_err, b1, b1_err))
print(sprintf("2\t%f\t%f\t%f\t%f",a2, a2_err, b2, b2_err))
print(sprintf("3\t%f\t%f\t%f\t%f",a3, a3_err, b3, b3_err))
print(sprintf("4\t%f\t%f\t%f\t%f",a4, a4_err, b4, b4_err))

potential(x)=asquaresigma*x-b/x+c
#b=0.1
fit potential(x) potentialfile u 1:4:5 yerrors via asquaresigma, b, c

print(sprintf("#result potentialfit\tasquaresigma=\t%e +/-\t%e\tb= \t%e+/-\t%e\tc=\t%e +/-\t%e\tbeta=\t%f\tdim=\t%d\t",asquaresigma, asquaresigma_err, b, b_err, c, c_err, beta, dim))
print(sprintf("#a%d%d=%f\n#b%d%d=%f\n#c%d%d=%f\n", beta*10, dim, asquaresigma, beta*10, dim, b, beta*10, dim, c))
unset print


set out sprintf("wilsonloopsu%dbeta%.3f.pdf", dim, beta)
set title sprintf("beta=%.1f, SU(%d)", beta, dim)
unset xrange 
unset yrange 
set xlabel "T/a"
set ylabel "W(R/a, T/a)"
plot wilsonloopfile u (($1==1)&&($5==32)?$2:1/0):3:4 w yerrorbars ls 1 title "R=1a", v1(x) ls 1 title "",\
	 wilsonloopfile u (($1==2)&&($5==32)?$2:1/0):3:4 w yerrorbars  ls 2 title "R=2a", v2(x) ls 2 title "",\
	 wilsonloopfile u (($1==3)&&($5==32)?$2:1/0):3:4 w yerrorbars  ls 3 title "R=3a", v3(x) ls 3 title "",\
	 wilsonloopfile u (($1==4)&&($5==32)?$2:1/0):3:4 w yerrorbars  ls 4 title "R=4a", v4(x) ls 4 title ""
	 
	 
set xlabel "R/a"
set ylabel "aV(R/a)"
#set yrange [-15:25]
set xrange [0.5:5]
set key top left
plot potentialfile using 1:4:5 w yerrorbars ls 1 title "V(R)", potential(x) title "" #:5 w yerrorbars

	 
set ter epslatex size 9 cm, 8 cm color colortext
unset title
unset xrange
unset yrange
set output 'comparisoncreutzreport.tex'
set key horizontal title '$\beta=$'
set yrange [0.1:1.1]
set xlabel 'steps'
set ylabel '1-$\langle P\rangle$'
plot for [betap=12:32:4] sprintf('../data/beta%d.dat', betap) using 1:(1-$3) linestyle ((betap-8)/4) w linespoints title sprintf('%.1f', betap/10.0)
#use betap instead of beta so wilsonloops still work
unset yrange
set key vertical top right title ''

	 
set out 'su3plaquettereport.tex'
set xrange [0:500]
set yrange[0.4:0.9]
set ylabel '$\langle P\rangle$'
plot '../data/su3beta5p5fixed.dat' u 1:3 ls 1 w linespoints title '$\beta=5.5$'
unset xrange 
unset yrange


set out sprintf('wilsonloopsu%dbeta%2dreport.tex', dim, beta*10)
print(sprintf('wilsonloopsu%dbeta%2dreport.tex', dim, beta*10))
set xlabel '$T/a$'
set ylabel '$W(R/a, T/a)$'
plot wilsonloopfile u (($1==1)&&($5==32)?$2:1/0):3:4 w yerrorbars ls 1 title '$R=1a$', v1(x) ls 1 title '',\
	 wilsonloopfile u (($1==2)&&($5==32)?$2:1/0):3:4 w yerrorbars  ls 2 title '$R=2a$', v2(x) ls 2 title '',\
	 wilsonloopfile u (($1==3)&&($5==32)?$2:1/0):3:4 w yerrorbars  ls 3 title '$R=3a$', v3(x) ls 3 title '',\
	 wilsonloopfile u (($1==4)&&($5==32)?$2:1/0):3:4 w yerrorbars  ls 4 title '$R=4a$', v4(x) ls 4 title ''


set out sprintf('potentialsu%dbeta%2dreport.tex', dim, beta*10)	 
set xlabel '$R/a$'
set ylabel '$aV(R/a)$'
set xrange [0.5:5]
set key top left
plot potentialfile using 1:4:5 w yerrorbars ls 1 title 'V(R)', potential(x) title '' #:5 w yerrorbars

potentialtot(x, a, b, c)=a*x-b/x+c

load "potentialparameters.cfg"

set out 'potentialsu2report.tex'
set yrange [0:4]
plot '../data/potentialsu2beta1.600.dat' using 1:4:5 w yerrorbars ls 1 title '$\beta=$1.6', potentialtot(x, a162, b162, c162) ls 1 title '',\
	 '../data/potentialsu2beta2.000.dat' using 1:4:5 w yerrorbars ls 2 title '$\beta=$2.0', potentialtot(x, a202, b202, c202) ls 2 title '',\
	 '../data/potentialsu2beta3.200.dat' using 1:4:5 w yerrorbars ls 3 title '$\beta=$3.2', potentialtot(x, a322, b322, c322) ls 3 title '',\
	 '../data/potentialsu2beta4.200.dat' using 1:4:5 w yerrorbars ls 4 title '$\beta=$4.2', potentialtot(x, a422, b422, c422) ls 4 title '',\
	 '../data/potentialsu2beta5.500.dat' using 1:4:5 w yerrorbars ls 5 title '$\beta=$5.5', potentialtot(x, a552, b552, c552) ls 5 title ''
unset yrange

set out 'potentialsu3report.tex'
plot '../data/potentialsu3beta2.000.dat' using 1:4:($5/100.0) w yerrorbars ls 1 title '$\beta=$2.0', potentialtot(x, a203, b203, c203) ls 1 title '',\
	 '../data/potentialsu3beta2.300.dat' using 1:4:($5/100.0) w yerrorbars ls 2 title '$\beta=$2.3', potentialtot(x, a233, b233, c233) ls 2 title '',\
	 '../data/potentialsu3beta5.500.dat' using 1:4:($5/100.0) w yerrorbars ls 4 title '$\beta=$5.5', potentialtot(x, a553, b553, c553) ls 4 title ''
unset xrange

set out 'stringtensionreport.tex'
set xlabel '$\beta$'
set ylabel '$a^2\sigma$'
set yrange [0:2]
set xrange [1.5:6]
plot potentialtotalfile u (($13==2)?$11:1/0):2:3 w yerrorbars ls 1 title 'SU(2)',\
	 potentialtotalfile u (($13==3)?$11:1/0):2:3 w yerrorbars ls 2 title 'SU(3)'
unset yrange
unset xrange

set out 'plaquettereport.tex'
set ylabel '$\langle P\rangle$'
set xrange [1.5:6]
set yrange[0:1]
set key bottom right 

strong(x)=x/4
weak(x)=1-3/(4*x)

plot plaquettetotalfile using (($1==2&&$3==32)?$2:1/0):4:5 w yerrorbars ls 1 title 'SU(2)',\
	 plaquettetotalfile using (($1==3&&$3==32)?$2:1/0):4:5 w yerrorbars ls 2 title 'SU(3)',\
	 strong(x) ls 3 dt 1 title 'strong coupling SU(2)',\
	 weak(x) ls 4 dt 2 title	 'weak coupling SU(2)'	 
set output
