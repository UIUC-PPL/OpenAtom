set terminal png
#set output 'betaconvloglog.png'
set output 'betaconv.png'
#set xlabel 'log_{10}(n_{r})'
set xlabel 'n_{r}'
#set ylabel 'log_{10}(percent error)'
#plot "betaconv2.dat" u (log10($1)):(log10(abs($2))) w lp title 'beat = 2.0', \
#	"betaconv2.5.dat" u (log10($1)):(log10(abs($2))) w lp title 'beat = 2.5', \
#	"betaconv3.dat" u (log10($1)):(log10(abs($2))) w lp title 'beat = 3.0', \
#	"betaconv3.5.dat" u (log10($1)):(log10(abs($2))) w lp title 'beat = 3.5', \
#	"betaconv4.dat" u (log10($1)):(log10(abs($2))) w lp title 'beat = 4.0', \
#	"betaconv4.5.dat" u (log10($1)):(log10(abs($2))) w lp title 'beat = 4.5'
plot "betaconv2.dat" u 1:(log10(abs($2))) w lp title 'beat = 2.0', \
	"betaconv2.5.dat" u 1:(log10(abs($2))) w lp title 'beat = 2.5', \
	"betaconv3.dat" u 1:(log10(abs($2))) w lp title 'beat = 3.0', \
	"betaconv3.5.dat" u 1:(log10(abs($2))) w lp title 'beat = 3.5', \
	"betaconv4.dat" u 1:(log10(abs($2))) w lp title 'beat = 4.0', \
	"betaconv4.5.dat" u 1:(log10(abs($2))) w lp title 'beat = 4.5'
