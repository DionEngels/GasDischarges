set style data histogram
set ytics nomirror
set xlabel "Energy (eV)"
set ylabel "Occurance (-)"
set title "Exercise 8.4"
set xrange [0:]
set size ratio 1
set style fill solid 1 noborder
plot "eedf.dat" w l t "Found distribution" 
pause 5
set term pdf
set output "../figs/Exercise 8.4.pdf"
replot
unset output
unset terminal