set logscale y
set xlabel "Energy (eV)"
set ylabel "Distribution (-)"
set title "Exercise 7.21"
set key left bottom
plot "eedf1.dat" w l lw 4, "eedf10.dat" w l lw 4, "eedf100.dat" w l lw 4
pause 2
set terminal pdf
set output "../figs/Exericse 7.21.pdf"
replot
unset output
unset terminal