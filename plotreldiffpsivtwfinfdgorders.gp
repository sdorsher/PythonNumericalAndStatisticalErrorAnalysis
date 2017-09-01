set xlabel "Time (M)"
set ylabel "|(Psir-Finf)/Finf|"
set logscale y
set title "Relative difference of radial self force summed from l=0 to l=30"

set terminal eps enhanced color
set output "reldiffpsirvtwfinfdgorders.eps"
plot '< paste genrawsummodes2.csv selfForceOverTime_mixed.csv' u 1:(abs(($2-$14)/($14))) w l t "DG order 16 wrt Finf", '' u 1:(abs(($3-$14)/($14))) w l t "DG order 20 wrt Finf", ''  u 1:(abs(($4-$14)/($14))) w l t "DG order 24 wrt Finf", '' u 1:(abs(($5-$14)/($14))) w l t "DG order 28 wrt Finf", '' u 1:(abs(($6-$14)/($14))) w l t "DG order 32 wrt Finf", '' u 1:(abs(($7-$14)/($14))) w l t "DG order 36 wrt Finf", '' u 1:(abs(($8-$14)/($14))) w l t "DG order 40 wrt Finf", '' u 1:(abs(($9-$14)/($14))) w l t "DG order 44 wrt Finf" 
set terminal x11
plot '< paste genrawsummodes2.csv selfForceOverTime_mixed.csv' u 1:(abs(($2-$14)/($14))) w l t "DG order 16 wrt Finf", '' u 1:(abs(($3-$14)/($14))) w l t "DG order 20 wrt Finf", ''  u 1:(abs(($4-$14)/($14))) w l t "DG order 24 wrt Finf", '' u 1:(abs(($5-$14)/($14))) w l t "DG order 28 wrt Finf", '' u 1:(abs(($6-$14)/($14))) w l t "DG order 32 wrt Finf", '' u 1:(abs(($7-$14)/($14))) w l t "DG order 36 wrt Finf", '' u 1:(abs(($8-$14)/($14))) w l t "DG order 40 wrt Finf", '' u 1:(abs(($9-$14)/($14))) w l t "DG order 44 wrt Finf" 