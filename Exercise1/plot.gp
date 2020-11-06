set ter epslatex size 15 cm, 10.6cm color colortext

set style line 1 lc 7 lt 7 pt 7 #ps 0.2
set style line 2 lc 1 lt 7 pt 8
set style line 3 lc 2 lt 7 pt 11
set style line 4 lc 3 lt 7 pt 3
set style line 5 lc 4 lt 7 pt 5#4 gut
set style line 6 lc 5 lt 7 pt 4#4 gut

#do not consider J because J=1
lambdaplus(T, h)=exp(1/T)*(cosh(h/T)+sqrt(sinh(h/T)*sinh(h/T)+exp(-4/T)))
lambdaminus(T, h)=exp(1/T)*(cosh(h/T)+sqrt(sinh(h/T)*sinh(h/T)+exp(-4/T)))
Z(T, h, N)=(lambdaplus(T,h))**N+(lambdaminus(T,h))**N

fracsinh(T, h)=sinh(h/T)*cosh(h/T)/sqrt(sinh(h/T)*sinh(h/T)+exp(-4/T))

magnetization(T,h,N)=(((lambdaplus(T,h))**(N-1))*exp(1/T)*(sinh(h/T)+fracsinh(T, h))+((lambdaminus(T,h))**(N-1))*exp(1/T)*(sinh(h/T)-fracsinh(T, h)))/Z(T,h,N)

set ter pdfcairo size 4in, 4in
set out "firsttry.pdf"
set xrange [0:10]
#plot magnetization(x,0.2,5) lt 2 title "h=0.2, N=5"

set xlabel "Temperatur"
set ylabel "magnetization"
do for [N=4:20:4]{
set title sprintf("N=%d", N)
plot for [h=-2:2] magnetization (x, h/2.0, N) title sprintf("h=%.1f", h/2.0)
}

plot for [N=4:20:4] magnetization (x, 1, N) title sprintf("N=%d", N)
