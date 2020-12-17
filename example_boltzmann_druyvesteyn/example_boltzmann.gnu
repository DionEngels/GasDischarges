set logscale y
set xlabel "Energy (eV)"
set ylabel "Distribution (-)"
set title "Exercise 7.18"
plot "eedf.dat" every ::0::3333 w l lw 4, "maxwell.dat" every ::0::3333 w l lw 2
pause 2
set terminal pdf
set output "../figs/Exericse 7.18_zoom.pdf"
replot
unset output
unset terminal