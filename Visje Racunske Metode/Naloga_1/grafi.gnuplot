set terminal epslatex color solid
set style data lines
set datafile separator ","

set output "g_stabilnost.tex"
set xlabel "\"Cas"
set ylabel "$\\langle\\psi|\\psi\\rangle$"

set logscale xy

plot "g_eksplicitna_0.csv" u 1:2 t "Eksplicitna, $\\lambda=0$", \
"g_eksplicitna_0.1.csv" u 1:2 t "Eksplicitna, $\\lambda=0.1$", \
"g_implicitna_0.csv" u 1:2 t "Implicitna, $\\lambda=0$", \
"g_implicitna_0.1.csv" u 1:2 t "Implicitna, $\\lambda=0.1$"

set output "g_trajektorija.tex"
set xlabel "\"Cas"
set ylabel "$\\langle\\psi|x|\\psi\\rangle$"

unset logscale

plot[0:2000] "g_implicitna_0.csv" u 1:3 t "$\\lambda=0$", \
"g_implicitna_0.1.csv" u 1:3 t "$\\lambda=0.1$", \
"g_implicitna_1.csv" u 1:3 t "$\\lambda=1$"

