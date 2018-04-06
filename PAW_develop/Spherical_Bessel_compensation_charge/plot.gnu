set terminal png
set output 'error.png'
set xlabel 'log_{10}n_r'
set ylabel 'log_{10}error'
set key right
plot "res2.0.dat" u (log10($2)):(log10($3-3.647597757)) title '\beta = 2.0' w lp, \
	 "res3.0.dat" u (log10($2)):(log10($3-3.816719151)) title '\beta = 3.0' w lp, \
	 "res4.0.dat" u (log10($2)):(log10($3-3.881322854)) title '\beta = 4.0' w lp
