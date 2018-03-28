set terminal png enhanced
set output 'fconv.png'
set xlabel 'n_r'
set ylabel 'log_{10}(Delta f/f)'
set key left
plot "feN.dat" u (log10($1)):(log10(abs($2-$3)/$3)) w lp title 'e-N force', \
        "fHar.dat" u (log10($1)):(log10(abs($2-$3)/$3)) w lp title 'Har force'
