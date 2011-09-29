# A simple gnuplot script to generate the tracer output file for the manual.

set key off
set term postscript eps color
set output "tracers.eps"
plot "./swarmOutput.00000.dat" using 3:4 w l lw 5 title "particle 0", "./swarmOutput.00001.dat" using 3:4 w l lw 5 title "particle 1", "./swarmOutput.00002.dat" using 3:4 w l lw 5 title "particle 2", "./swarmOutput.00003.dat" using 3:4 w l lw 5 title "particle 3", "./swarmOutput.00004.dat" using 3:4 w l lw 5 title "particle 4", "./swarmOutput.00005.dat" using 3:4 w l lw 5 title "particle 5", "./swarmOutput.00006.dat" using 3:4 w l lw 5 title "particle 6", "./swarmOutput.00007.dat" using 3:4 w l lw 5 title "particle 7"
