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

file="../data/results.dat"
radiusfile="../data/radius.dat"

set fit errorvariables

F300(x)=a300+b300*x+c300*x**2+d300*x**3
fit [0:1] F300(x) file using ((($1==300)&&($2==20))?$4:1/0):5 via a300,b300, c300, d300


F400(x)=a400+b400*x+c400*x**2+d400*x**3
fit [0:1] F400(x) file using ((($1==400)&&($2==20))?$4:1/0):5 via a400,b400, c400, d400


F500(x)=a500+b500*x+c500*x**2+d500*x**3
fit [0:1] F500(x) file using ((($1==500)&&($2==20))?$4:1/0):5 via a500,b500, c500, d500


F600(x)=a600+b600*x+c600*x**2+d600*x**3
fit [0:1] F600(x) file using ((($1==600)&&($2==20))?$4:1/0):5 via a600,b600, c600, d600


F700(x)=a700+b700*x+c700*x**2+d700*x**3
fit [0:1] F700(x) file using ((($1==700)&&($2==20))?$4:1/0):5 via a700,b700, c700, d700


F800(x)=a800+b800*x+c800*x**2+d800*x**3
fit [0:1] F800(x) file using ((($1==800)&&($2==20))?$4:1/0):5 via a800,b800, c800, d800


F900(x)=a900+b900*x+c900*x**2+d900*x**3
fit [0:1] F900(x) file using ((($1==900)&&($2==20))?$4:1/0):5 via a900,b900, c900, d900


F1000(x)=a1000+b1000*x+c1000*x**2+d1000*x**3
fit [0:1] F1000(x) file using ((($1==1000)&&($2==20))?$4:1/0):5 via a1000,b1000, c1000, d1000


F1100(x)=a1100+b1100*x+c1100*x**2+d1100*x**3
fit [0:1] F1100(x) file using ((($1==1100)&&($2==20))?$4:1/0):5 via a1100,b1100, c1100, d1100


F1200(x)=a1200+b1200*x+c1200*x**2+d1200*x**3
fit [0:1] F1200(x) file using ((($1==1200)&&($2==20))?$4:1/0):5 via a1200,b1200, c1200, d1200


set print radiusfile
print(sprintf("300\t%f\t%f\t%f\t%f\n",b300, b300_err, sqrt(-6*b300), 3*b300_err/sqrt(-6*b300)))
print(sprintf("400\t%f\t%f\t%f\t%f\n",b400, b400_err, sqrt(-6*b400), 3*b400_err/sqrt(-6*b400)))
print(sprintf("500\t%f\t%f\t%f\t%f\n",b500, b500_err, sqrt(-6*b500), 3*b500_err/sqrt(-6*b500)))
print(sprintf("600\t%f\t%f\t%f\t%f\n",b600, b600_err, sqrt(-6*b600), 3*b600_err/sqrt(-6*b600)))
print(sprintf("700\t%f\t%f\t%f\t%f\n",b700, b700_err, sqrt(-6*b700), 3*b700_err/sqrt(-6*b700)))
print(sprintf("800\t%f\t%f\t%f\t%f\n",b800, b800_err, sqrt(-6*b800), 3*b800_err/sqrt(-6*b800)))
print(sprintf("900\t%f\t%f\t%f\t%f\n",b900, b900_err, sqrt(-6*b900), 3*b900_err/sqrt(-6*b900)))
print(sprintf("1000\t%f\t%f\t%f\t%f\n",b1000, b1000_err, sqrt(-6*b1000), 3*b1000_err/sqrt(-6*b1000)))
print(sprintf("1100\t%f\t%f\t%f\t%f\n",b1100, b1100_err, sqrt(-6*b1100), 3*b1100_err/sqrt(-6*b1100)))
print(sprintf("1200\t%f\t%f\t%f\t%f\n",b1200, b1200_err, sqrt(-6*b1200), 3*b1200_err/sqrt(-6*b1200)))




set terminal pdfcairo size 4in, 4in
set out 'test.pdf'

set xlabel "q^2/fm^-2"
set ylabel "F"
#set yrange [-0.05:1.05]
plot for [lambda=300:1200:100] file using ((($1==lambda)&&($2==20))?$4:1/0):5 with linespoints title sprintf("lambda=%d", lambda)
set xrange [0:5]
plot for [lambda=300:1200:100] file using ((($1==lambda)&&($2==20))?$4:1/0):5 with linespoints title sprintf("lambda=%d", lambda)#, F1200(x)
set xrange [1:1.5]
plot for [lambda=300:1200:100] file using ((($1==lambda)&&($2==20))?$4:1/0):5 with linespoints title sprintf("lambda=%d", lambda)#, F1200(x)
set xrange [0:0.005]
plot for [lambda=300:1200:100] file using ((($1==lambda)&&($2==20))?$4:1/0):5 with linespoints title sprintf("lambda=%d", lambda)#, F1200(x)
set xrange [80:100]
unset yrange
plot for [lambda=300:1200:100] file using ((($1==lambda)&&($2==20))?$4:1/0):5 with linespoints title sprintf("lambda=%d", lambda)
unset xrange
set title 'accuracy'
set xlabel 'nx'
plot for [value=1:10:1] file using (($1==1200)&&($3==value)?$2:1/0):5 with linespoints title sprintf("q=%.1f", value)

set xlabel "lambda" 

plot radiusfile u 1:4:5 w yerrorbars, file using ($3==0&&$2==20?$1:1/0):5

set yrange [-1e-6:+1e-6]
zero(x)=0
plot file using ($3==0&&$2==20?$1:1/0):($5-1) ls 1 title '', zero(x) title ''
unset xrange
unset yrange


set ter epslatex size 15 cm, 10 cm color colortext
unset title
set out 'stability.tex'
set key top rmargin title '$q/\si{\per\femto\meter}$'
set xlabel 'nx'
set ylabel '$F(\vec{q}^2)/\si{\femto\meter^2}$'
plot file using (($1==1200)&&($3==1)?$2:1/0):5 with linespoints ls 1 title '1',\
	file using (($1==1200)&&($3==2)?$2:1/0):5 with linespoints ls 2 title '2',\
	file using (($1==1200)&&($3==3)?$2:1/0):5 with linespoints ls 3 title '3',\
	file using (($1==1200)&&($3==9)?$2:1/0):5 with linespoints ls 4 title '9'
	
set out 'radius.tex'
set key inside top right title ''
set xlabel '$\Lambda/\si{\mega\electronvolt}$'
set ylabel '$\sqrt{\langle r^2 \rangle}/\si{\femto\meter}$'
plot radiusfile u 1:4:5 w yerrorbars ls 1 title 'measured radius', '../data/rlecture.dat' u 1:2 ls 3 title 'radius from lecture'

set yrange [-1e-6:+1e-6]
set out 'fofzero.tex'
set key inside top right title ''
set xlabel '$\Lambda/\si{\mega\electronvolt}$'
set ylabel '$F(\vec{q}^2)/\si{\femto\meter^2}-1$'
plot file using ($3==0&&$2==20?$1:1/0):($5-1) ls 1 title '', zero(x) ls 1 dt 2 title ''
unset yrange

set out 'formfactor.tex'
set key top rmargin title '$\Lambda/\si{\mega\electronvolt}$'
set xlabel '$q^2/\si{\femto\meter^{-2}}$'
set ylabel '$F(\vec{q}^2)/\si{\femto\meter^2}$'

plot for [lambda=300:1200:100] file using ((($1==lambda)&&($2==20))?$4:1/0):5 with linespoints ls ((lambda-200)/100) title sprintf('%d', lambda)

set out 'formfactorsmall.tex'
set xrange [0:5]
set xlabel '$q^2/\si{\femto\meter^{-2}}$'
set ylabel '$F(\vec{q}^2)/\si{\femto\meter^2}$'

plot for [lambda=300:1200:100] file using ((($1==lambda)&&($2==20))?$4:1/0):5 with linespoints ls ((lambda-200)/100) title sprintf('%d', lambda)

set out 'formfactorbig.tex'
set xrange [80:100]
set xlabel '$q^2/\si{\femto\meter^{-2}}$'
set ylabel '$F(\vec{q}^2)/\si{\femto\meter^2}$'

plot for [lambda=300:1200:100] file using ((($1==lambda)&&($2==20))?$4:1/0):5 with linespoints ls ((lambda-200)/100) title sprintf('%d', lambda)

set output
