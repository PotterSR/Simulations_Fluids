# lattice.gp

set terminal pngcairo enhanced size 800,800
set output "lattice.png"

# Match your simulation box
L = 40.0
R = 1.0

set xrange [0:L]
set yrange [0:L]
set size ratio 1          # square plot — important for circles to look round

unset key
set title "Square Lattice — N=200, eta=0.393"
set xlabel "x"
set ylabel "y"

# Draw each particle as a filled circle of radius R
# 'with circles' needs: x  y  radius
plot "positions.dat" using 1:2:(R) with circles \
     lc rgb "#4A90D9"  \
     fs solid 0.5 border lc rgb "#1A5A99"