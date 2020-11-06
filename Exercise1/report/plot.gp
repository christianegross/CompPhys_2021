set ter epslatex size 15 cm, 10.6cm color colortext

set style line 1 lc 7 lt 7 pt 7 #ps 0.2
set style line 2 lc 1 lt 7 pt 8
set style line 3 lc 2 lt 7 pt 11
set style line 4 lc 3 lt 7 pt 3
set style line 5 lc 4 lt 7 pt 5#4 gut
set style line 6 lc 5 lt 7 pt 4#4 gut

#do not consider J because J=1
lambdaplus(T, h) =exp(1.0/T)*(cosh(h/T)+sqrt(sinh(h/T)*sinh(h/T)+exp(-4.0/T)))
lambdaminus(T, h)=exp(1.0/T)*(cosh(h/T)-sqrt(sinh(h/T)*sinh(h/T)+exp(-4.0/T)))
Z(T, h, N)=(lambdaplus(T,h))**N+(lambdaminus(T,h))**N

fracsinh(T, h)=sinh(h/T)*cosh(h/T)/sqrt(sinh(h/T)*sinh(h/T)+exp(-4.0/T))

magnetization(T,h,N)=-(((lambdaplus(T,h))**(N-1))*exp(1/T)*(sinh(h/T)+fracsinh(T, h))+((lambdaminus(T,h))**(N-1))*exp(1/T)*(sinh(h/T)-fracsinh(T, h)))/Z(T,h,N)

set ter pdfcairo size 4in, 4in
#pdf for easy viewing
set out "firsttry.pdf"

set xrange [0:5]
#plot magnetization(x,0.2,5) lt 2 title "h=0.2, N=5"

set xlabel "temperature"
set ylabel "magnetization"
do for [N=1:8:2]{
set title sprintf("N=%d", N)
#plot for [h=-2:2] magnetization (x, h/2.0, N) title sprintf("h=%.1f", h/2.0)
}

set title "h=1"
plot for [N=1:6:1] magnetization (x, 1.0, N) title sprintf("N=%d", N)


set ter epslatex size 15 cm, 10.6cm color colortext
#.tex for easy inputting in report

unset title
set out "magnetizationvaryingh.tex"
plot for [h=-2:2] magnetization (x, h/2.0, N) title sprintf("$h=%.1f$", h/2.0)

set out "magnetizationvaryingN.tex"
plot for [N=1:6:1] magnetization (x, 1.0, N) title sprintf("$N=%d$", N)
