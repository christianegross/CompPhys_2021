#Christiane Gross, Nico Dichter
#Exercise 3 Computational Physics WS20/21
#script for plotting theoretical and experimental results of magnetization and energy

set style line 1 lc 7 lt 7 pt 7 #ps 0.2
set style line 2 lc 1 lt 7 pt 8
set style line 3 lc 2 lt 7 pt 11
set style line 4 lc 3 lt 7 pt 3
set style line 5 lc 4 lt 7 pt 5#4 gut
set style line 6 lc 5 lt 7 pt 4#4 gut

set style line 8 lc 2 lt 7 pt 11

#used definitons for expected values from sheet
binomial(n,k)			=n!/k!/(n-k)!
smallf(J,h,x)			=exp(0.5*J*x**2+h*x)
Z(N, J, h)				=(sum [i=0:N] binomial(N, i)*smallf(J,h,(N-2*i)))
betaepsilon(N,J,h)		=-1.0/N/Z(N,J,h)*(sum [i=0:N] binomial(N, i)*(0.5*J*(N-2*i)**2+h*(N-2*i))*smallf(J,h,N-2*i))
magnetization(N,J,h)	=1.0/N/Z(N,J,h)*(sum [i=0:N] binomial(N, i)*(N-2*i)*smallf(J,h,N-2*i))

convergefile='../data/converge.dat'
datafile='../data/raw.dat'

set ter pdfcairo size 4in, 4in
#pdf for easy viewing

set out 'converge.pdf'
set title 'Convergence check'
set xlabel 'N_{md}'
set ylabel '|delta|'
set logscale y
set logscale x
set grid

plot convergefile u 1:2 ls 2 title 'theory'
unset logscale x
unset logscale y

set out 'test.pdf'
binwidth=0.01
bin(x,width)=width*floor(x/width)+width/2.0
set boxwidth binwidth

f(x,J,h,N)=1/sqrt(2*3.1415*J/N)*exp(-x**2/(2*J/N)+N*log(2*cosh(h+x)))/5000000
#plot '../data/raw.dat' using (bin($1,binwidth)):(0.0002) smooth freq with boxes,f(x,1.,0.5,10.)
plot f(x,1.,0.5,10.)

set out 'bootstrap.pdf'
set logscale x
set xlabel 'length of bin'
set ylabel 'error
set key top left
do for [size=15:20:5]{
set title sprintf("magnetization, size=%d", size)
plot for [n=4:40:1] '../data/bootstrapbinlength.dat' u ($2==n/20.0&&$1==size?$7:1/0):4 w linespoints lc n title sprintf("J=%f",n/20.0)
}
do for [size=15:20:5]{
set title sprintf("energy, size=%d", size)
plot for [n=4:40:1] '../data/bootstrapbinlength.dat' u ($2==n/20.0&&$1==size?$7:1/0):6 w linespoints lc n title sprintf("J=%f",n/20.0)
}

unset logscale x

set xrange [0.2:2]
set yrange[0:1]
set out 'magnetization.pdf'
set title 'magnetization'
set xlabel 'J'
set ylabel 'magnetization'
plot for [N=5:20:5] magnetization(N,x/N,0.5) ls (N/5) title sprintf("N=%d",N), for [N=5:20:5] datafile using ($1==N&&$7==64?$2:1/0):3:4 w yerrorbars ls (N/5) title ''
unset yrange

set out 'energy.pdf'
set title 'energy'
set xlabel 'J'
set ylabel 'energy'
plot for [N=5:20:5] betaepsilon(N,x/N,0.5) ls (N/5) title sprintf("N=%d",N), for [N=5:20:5] datafile using ($1==N&&$7==64?$2:1/0):($5+1):6 w yerrorbars ls (N/5) title ''

unset xrange
unset yrange

set ter epslatex size 15cm, 10cm color colortext

set out 'converge_t.tex'
set title 'Convergence check'
set xlabel '$N_{md}$'
set ylabel '$|\Delta|$'
set logscale y
set grid

plot convergefile u 1:2 ls 2 title 'theory'
unset logscale y

set xrange [0.2:2]
set yrange[0:1]
set out 'magnetization_t.tex'
set title ''
set xlabel '$J$'
set ylabel '$\langle m \rangle$'
set key bottom right
plot for [N=5:20:5] magnetization(N,x/N,0.5) ls (N/5) title '', for [N=5:20:5] datafile using ($1==N&&$7==64?$2:1/0):3:4 w yerrorbars ls (N/5) title sprintf("$N=$%d",N)
unset yrange


set out 'energy_t.tex'
set title ''
set xlabel '$J$'
set ylabel '$\langle \beta\epsilon \rangle$'
set key top right
plot for [N=5:20:5] betaepsilon(N,x/N,0.5) ls (N/5) title '', for [N=5:20:5] datafile using ($1==N&&$7==64?$2:1/0):($5+1):6 w yerrorbars ls (N/5) title sprintf("$N=$%d",N)
