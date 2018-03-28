set terminal png
set output 'conv.png'
set xlabel 'n_r'
set ylabel 'log_{10}(Delta E/E)'
set key left
plot "conv.dat" u (log10($1)):(log10(abs(($2-$8)-($3-$9))/($3-$9))) w lp title 'e-N w/o self', \
        "conv.dat" u (log10($1)):(log10(abs(($4-$7)-($5-$6))/($5-$6))) w lp title 'Hartree w/o self'
