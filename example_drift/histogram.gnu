set style data histogram
set ytics nomirror
set y2tics nomirror
set xlabel "Energy (eV)"
set ylabel "Occurance (-)"
set y2label "Distribution (-)"
set title "Exercise 8.4"
set xrange [0:]
set size ratio 1
set style fill solid 1 noborder
a = 1
f(x) = 2 * sqrt(x/a/pi) * exp(-x/a)
plot "eedf_2000.dat" w l t "Found distribution" axis x1y1, f(x) axis x1y2 title "Maxwellian Distribution"
pause 5
set term pdf
set output "../figs/Exercise 8.4_2k.pdf"
replot
unset output
unset terminal