set logscale y
plot "b.dat" using 1:(($1)**0.5*($2)) w l lw 4, "maxwell.dat" using 1:2 w l lw 2
pause -1
