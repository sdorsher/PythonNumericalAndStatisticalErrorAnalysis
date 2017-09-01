set xlabel "Time (M)"
set ylabel "Psir"
set title "Radial self force summed from l=0 to l=30"

set terminal eps enhanced color
set output "psirvtwfinfdgorders.eps"
plot 'genrawsummodes.csv' u 1:2 w l t "DG order 16", '' u 1:3 w l t "DG order 20", '' u 1:4 w l t "DG order 24", '' u 1:5 w l t "DG order 28", '' u 1:6 w l t "DG order 32", '' u 1:7 w l t "DG order 36", '' u 1:8 w l t "DG order 40", '' u 1:9 w l t "DG order 44", 'selfForceOverTime_mixed.csv' u 1:5 w l t "Finf"
set terminal x11
plot 'genrawsummodes.csv' u 1:2 w l t "DG order 16", '' u 1:3 w l t "DG order 20", '' u 1:4 w l t "DG order 24", '' u 1:5 w l t "DG order 28", '' u 1:6 w l t "DG order 32", '' u 1:7 w l t "DG order 36", '' u 1:8 w l t "DG order 40", '' u 1:9 w l t "DG order 44", 'selfForceOverTime_mixed.csv' u 1:5 w l t "Finf"