set yrange[0:0.11];
set terminal png
set output 'error3D.png'
set xlabel '\beta_{universal}'
set ylabel 'percent error'
set key left
plot "EHarselfscr14.dat" u 1:(100.0*(abs($2-$3)/$2)) title 'nr = 14'w lp, \
    0.05 title '' , 0.025 title '', 0.1 title '', 0.01 title '' 
