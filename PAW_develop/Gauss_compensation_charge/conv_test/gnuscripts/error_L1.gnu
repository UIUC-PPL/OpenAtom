set yrange[0:0.11];
set terminal png
set output 'error.png'
set xlabel '\beta_{universal}'
set ylabel 'percent error'
set key left
plot "EHarselfscr6.dat" u 1:(100.0*(abs($2-$3)/$2)) title 'nr = 6' w lp, \
	"EHarselfscr8.dat" u 1:(100.0*(abs($2-$3)/$2)) title 'nr = 8'w lp, \
	"EHarselfscr10.dat" u 1:(100.0*(abs($2-$3)/$2)) title 'nr = 10'w lp, \
	"EHarselfscr12.dat" u 1:(100.0*(abs($2-$3)/$2)) title 'nr = 12'w lp, \
	"EHarselfscr14.dat" u 1:(100.0*(abs($2-$3)/$2)) title 'nr = 14'w lp, \
    0.05 title '' , 0.025 title '', 0.1 title '', 0.01 title '' 
