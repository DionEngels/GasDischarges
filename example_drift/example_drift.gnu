set ytics nomirror
set y2tics nomirror
set xlabel "time (s)"
set ylabel "v (m/s)"
set y2label "number of particles (-)"
set title "Exercise 8.2"
set key left bottom
set key horizontal
set size ratio 0.15
plot "vel_drift.dat" w l lw 1 axis x1y1, "vel_average.dat" w l lw 1 axis x1y1, "n_particles.dat" w l lw 1 axis x1y2
pause 2
set term pdf size 6 in, 1.5 in font "Helvetica,8"
set output "../figs/Exercise 8.2 combi.pdf"
replot
unset output
unset terminal