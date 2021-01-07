set ytics nomirror
set xlabel "Energy (eV)"
set ylabel "Occurance (-)"
set title "Exercise 8.5"
set xrange [0:]
set size ratio 1
set style fill solid 1 noborder
plot "eedf.dat" w l t "8.5 distribution", "eedf_2000.dat" w l title "8.4 Distribution"
pause 5
set term pdf
set output "../figs/Exercise 8.5_2k.pdf"
replot
unset output
unset terminal