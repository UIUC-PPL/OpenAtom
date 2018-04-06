set terminal png
set output 'errorsimple.png'
set xlabel 'log_{10}n_r'
set ylabel 'log_{10}error'
set key right
plot "test.dat" u (log10($2)):(log10($4-3.9690515)) title 'simple' w lp
