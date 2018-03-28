set terminal png
set output 'diffbeta.png'
set xlabel 'log_{10}(\beta_{universal})'
set ylabel 'log_{10}(percent error)'
plot "diff.dat" u (log10($1)):(log10(abs($2))) w lp title 'screening error', \
	"diff.dat" u (log10($1)):(log10(abs($3))) w lp title 'screening + grid error', \
	"diff.dat" u (log10($1)):(1.38-2*log10($1)) w lp title 'asymptotic scaling'	
