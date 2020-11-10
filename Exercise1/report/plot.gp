set ter epslatex size 15 cm, 10.6cm color colortext

set style line 1 lc 7 lt 7 pt 7 #ps 0.2
set style line 2 lc 1 lt 7 pt 8
set style line 3 lc 2 lt 7 pt 11
set style line 4 lc 3 lt 7 pt 3
set style line 5 lc 4 lt 7 pt 5#4 gut
set style line 6 lc 5 lt 7 pt 4#4 gut

set style line 8 lc 2 lt 7 pt 11

#do not consider J because J=1, set T=1 in plot commands
lambdaplus(T, h) =exp(1.0/T)*(cosh(h/T)+sqrt(sinh(h/T)*sinh(h/T)+exp(-4.0/T)))
lambdaminus(T, h)=exp(1.0/T)*(cosh(h/T)-sqrt(sinh(h/T)*sinh(h/T)+exp(-4.0/T)))
Z(T, h, N)=(lambdaplus(T,h))**N+(lambdaminus(T,h))**N

fracsinh(T, h)=sinh(h/T)*cosh(h/T)/sqrt(sinh(h/T)*sinh(h/T)+exp(-4.0/T))

magnetization(T,h,N)=(((lambdaplus(T,h))**(N-1))*exp(1/T)*(sinh(h/T)+fracsinh(T, h))+((lambdaminus(T,h))**(N-1))*exp(1/T)*(sinh(h/T)-fracsinh(T, h)))/Z(T,h,N)

set ter pdfcairo size 4in, 4in
#pdf for easy viewing

datafile="../data/N_h.dat"
Nfix=16
hfix=0.4
set yrange [-1.5:1.5]
set ylabel "magnetization"


set out "varyingh.pdf"
set xlabel "h"
set xrange [-1:1]
set key top left
do for [N in "2 4 8 12 16 20"]{
set title sprintf("N=%s", N)
plot magnetization(1,x,N) title "expectation", datafile u (($1==N)?$2:1/0):3:4 w yerrorbars title "measurement"#sprintf("N=%s", N)
}


set out "varyingN.pdf"
set xlabel "N"
set xrange [1:21]
set key top right
do for [a=-10:10:2]{
h=a/10.0
set title sprintf("h=%.1f", h)
plot magnetization(1,h,x) title "expectation", datafile u (($2==h)?$1:1/0):3:4 w yerrorbars title "measurement"#sprintf("h=%.1f", h)
}


set ter epslatex size 15 cm, 10.6cm color colortext
#.tex for easy inputting in report

unset title

set out "magnetizationvaryingh.tex"
set xrange [-1:1]

set key top left
plot magnetization(1,x,Nfix) ls 1 lw 2 title "expectation", datafile u (($1==Nfix)?$2:1/0):3:4 w yerrorbars ls 2 ps 2 title "measurement"#sprintf("N=%s", Nfix)

set out "magnetizationvaryingN.tex"
set xrange [1:21]
set yrange [-0.5:1.5]
set key bottom right
plot magnetization(1,hfix,x) ls 1 lw 2 title "expectation", datafile u (($2==hfix)?$1:1/0):3:4 w yerrorbars ls 2 ps 2 title "measurement"#sprintf("h=%.1f", hfix)

set output
