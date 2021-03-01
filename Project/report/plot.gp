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

set yrange[0:1]

plot '../data/plaquettethermhotstartepsilon0p4beta2p3size4.dat' u 0:2 ls 1 title 'during',\
	 '../data/plaquettethermhotstartepsilon0p4beta2p3size4.dat' u 0:3 ls 2 title 'after',\

set xrange [0:100]

set title 'acceptance rate early'
plot '../data/plaquettethermhotstartepsilon0point2beta2point3.dat' u 0:1 ls 1 title 'ar hot',\
	 '../data/plaquettethermcoldstartepsilon0point2beta2point3.dat' u 0:1 ls 4 title 'ar cold'

set title 'plaquette early'
plot '../data/plaquettethermhotstartepsilon0point2beta2point3.dat' u 0:2 w lines ls 2 title 'pd hot',\
	 '../data/plaquettethermhotstartepsilon0point2beta2point3.dat' u 0:3 w lines ls 3 title 'pa hot',\
	 '../data/plaquettethermcoldstartepsilon0point2beta2point3.dat' u 0:2 w lines ls 5 title 'pd cold',\
	 '../data/plaquettethermcoldstartepsilon0point2beta2point3.dat' u 0:3 w lines ls 6 title 'pa cold'

set xrange [900:1000]

set title 'acceptance rate late'
plot '../data/plaquettethermhotstartepsilon0point2beta2point3.dat' u 0:1 ls 1 title 'ar hot',\
	 '../data/plaquettethermcoldstartepsilon0point2beta2point3.dat' u 0:1 ls 4 title 'ar cold'

set title 'plaquette late'
plot '../data/plaquettethermhotstartepsilon0point2beta2point3.dat' u 0:2 ls 2 title 'pd hot',\
	 '../data/plaquettethermhotstartepsilon0point2beta2point3.dat' u 0:3 ls 3 title 'pa hot',\
	 '../data/plaquettethermcoldstartepsilon0point2beta2point3.dat' u 0:2 ls 5 title 'pd cold',\
	 '../data/plaquettethermcoldstartepsilon0point2beta2point3.dat' u 0:3 ls 6 title 'pa cold'


unset xrange

set mytics 4
set mxtics 2

set xrange [0:50]

set out 'comparisoncreutz.pdf'

set title "plaquette, measured during sweep"
plot "../data/beta12.dat" u 1:3 title "1.2",\
	 "../data/beta16.dat" u 1:3 title "1.6",\
	 "../data/beta20.dat" u 1:3 title "2.0",\
	 "../data/beta24.dat" u 1:3 title "2.4",\
	 "../data/beta28.dat" u 1:3 title "2.8",\
	 "../data/beta32.dat" u 1:3 title "3.2"	 

set title "1-plaquette, measured during sweep"
plot "../data/beta12.dat" u 1:(1-$3) title "1.2",\
	 "../data/beta16.dat" u 1:(1-$3) title "1.6",\
	 "../data/beta20.dat" u 1:(1-$3) title "2.0",\
	 "../data/beta24.dat" u 1:(1-$3) title "2.4",\
	 "../data/beta28.dat" u 1:(1-$3) title "2.8",\
	 "../data/beta32.dat" u 1:(1-$3) title "3.2"	 
	 
set title "plaquette, measured after sweep"
plot "../data/beta12.dat" u 1:4 title "1.2",\
	 "../data/beta16.dat" u 1:4 title "1.6",\
	 "../data/beta20.dat" u 1:4 title "2.0",\
	 "../data/beta24.dat" u 1:4 title "2.4",\
	 "../data/beta28.dat" u 1:4 title "2.8",\
	 "../data/beta32.dat" u 1:4 title "3.2"	 

set title "1-plaquette, measured after sweep"
plot "../data/beta12.dat" u 1:(1-$4) title "1.2",\
	 "../data/beta16.dat" u 1:(1-$4) title "1.6",\
	 "../data/beta20.dat" u 1:(1-$4) title "2.0",\
	 "../data/beta24.dat" u 1:(1-$4) title "2.4",\
	 "../data/beta28.dat" u 1:(1-$4) title "2.8",\
	 "../data/beta32.dat" u 1:(1-$4) title "3.2"	 

	 
set ter epslatex size 15 cm, 10 cm color colortext
unset title
set output
