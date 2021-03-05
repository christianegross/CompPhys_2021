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

set terminal pdfcairo size 4in, 4in
set out 'test.pdf'

set mytics 4
set mxtics 2

set xrange [0:50]

set out 'comparisoncreutz.pdf'

set title "plaquette, measured during sweep"
plot for [beta=12:32:4] sprintf("../data/beta%d.dat", beta) using 1:($3) linestyle ((beta-8)/4) title sprintf('%.1f', beta/10.0)

set title "1-plaquette, measured during sweep"
plot for [beta=12:32:4] sprintf("../data/beta%d.dat", beta) using 1:(1-$3) linestyle ((beta-8)/4) title sprintf('%.1f', beta/10.0) 
	 
set title "plaquette, measured after sweep"
plot for [beta=12:32:4] sprintf("../data/beta%d.dat", beta) using 1:($4) linestyle ((beta-8)/4) title sprintf('%.1f', beta/10.0) 

set title "1-plaquette, measured after sweep"
plot for [beta=12:32:4] sprintf("../data/beta%d.dat", beta) using 1:(1-$4) linestyle ((beta-8)/4) title sprintf('%.1f', beta/10.0)	 

set out "su3plaquette.pdf"
set title "hot start, epsilon=0.3"
set xlabel "steps"
set ylabel "<P>"
unset xrange
set yrange[0.4:1]
plot "../data/su3beta5p5fixed.dat" u 1:3 ls 1 w lines title "beta=5.5", "../data/su3beta5p5fixed.dat" u 1:4 w lines ls 2 title "beta=5.5"
set xrange [0:50]
plot "../data/su3beta5p5fixed.dat" u 1:3 ls 1 w lines title "beta=5.5", "../data/su3beta5p5fixed.dat" u 1:4 w lines ls 2 title "beta=5.5"

beta=2.0
dim=2

wilsonloopfile=sprintf("../data/wilsonsu%dbeta%.3f.dat", dim, beta)
potentialfile=sprintf("../data/potentialsu%dbeta%3f.dat", dim, beta)

v1(x)=a1*exp(-b1*x)
v2(x)=a2*exp(-b2*x)
v3(x)=a3*exp(-b3*x)
v4(x)=a4*exp(-b4*x)

#a3=0.9
#a3_err=0.1
#b3=1.8
#b3_err=0.1

#a4=0.9
#b4=2.5
#a4_err=0.1
#b4_err=0.1

fit v1(x) wilsonloopfile u (($1==1)&&($5==16)?$2:1/0):3:4 yerrors via a1, b1
fit v2(x) wilsonloopfile u (($1==2)&&($5==16)?$2:1/0):3:4 yerrors via a2, b2
fit v3(x) wilsonloopfile u (($1==3)&&($5==16)?$2:1/0):3:4 yerrors via a3, b3
fit v4(x) wilsonloopfile u (($1==4)&&($5==16)?$2:1/0):3:4 yerrors via a4, b4


set print potentialfile
print(sprintf("R\ta\ta_err\tb\tb_err"))
print(sprintf("1\t%f\t%f\t%f\t%f",a1, a1_err, b1, b1_err))
print(sprintf("2\t%f\t%f\t%f\t%f",a2, a2_err, b2, b2_err))
print(sprintf("3\t%f\t%f\t%f\t%f",a3, a3_err, b3, b3_err))
print(sprintf("4\t%f\t%f\t%f\t%f",a4, a4_err, b4, b4_err))

potential(x)=asquaresigma*x-b/x+c

fit potential(x) potentialfile u 1:4:5 yerrors via asquaresigma, b, c

print(sprintf("#result potentialfit\t%e\t%e\t%e\t%e\t%e\t%e\n",asquaresigma, asquaresigma_err, b, b_err, c, c_err))

set out sprintf("wilsonloopsu%dbeta%.3f.pdf", dim, beta)
set title sprintf("beta=%.1f, SU(%d)", dim, beta)
unset xrange 
unset yrange 
set xlabel "T/a"
set ylabel "W(R/a, T/a)"
plot wilsonloopfile u (($1==1)&&($5==32)?$2:1/0):3:4 w yerrorbars ls 1 title "R=1", v1(x) ls 1 title "",\
	 wilsonloopfile u (($1==2)&&($5==32)?$2:1/0):3:4 w yerrorbars  ls 2 title "R=2", v2(x) ls 2 title "",\
	 wilsonloopfile u (($1==3)&&($5==32)?$2:1/0):3:4 w yerrorbars  ls 3 title "R=3", v3(x) ls 3 title "",\
	 wilsonloopfile u (($1==4)&&($5==32)?$2:1/0):3:4 w yerrorbars  ls 4 title "R=4", v4(x) ls 4 title ""
	 
	 
set xlabel "R/a"
set ylabel "aV(R/a)"
#set yrange [-15:25]
set xrange [0:5]
set key top left
plot potentialfile using 1:4 ls 1 title "V(R)", potential(x) title "" #:5 w yerrorbars

	 
set ter epslatex size 9 cm, 10 cm color colortext
unset title
unset xrange
unset yrange
set output 'comparisoncreutzreport.tex'
set key title '$\beta=$'
set yrange [0:1.1]
plot for [beta=12:32:4] sprintf("../data/beta%d.dat", beta) using 1:(1-$3) linestyle ((beta-8)/4) title sprintf('%.1f', beta/10.0)

set output
