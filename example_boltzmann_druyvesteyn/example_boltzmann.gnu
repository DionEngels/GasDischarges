set logscale y
set xlabel "Energy (eV)"
set ylabel "Distribution (-)"
set title "Exercise 7.20 - 5 eV, 200 * 10^{-22}"
plot "eedf.dat" w l lw 4, "maxwell.dat" w l lw 2
pause 2
set terminal pdf
set output "../figs/Exericse 7.20_5_200.pdf"
replot
unset output
unset terminal