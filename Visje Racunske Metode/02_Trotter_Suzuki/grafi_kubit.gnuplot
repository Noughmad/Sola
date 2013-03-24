set terminal epslatex color solid


set output "g_casovna_zahtevnost.tex"
set xlabel "\"Stevilo delcev $n$"
set ylabel "\"Cas ra\"cunanja [s]"
set logscale y

A = 0.14
O(x) = A + B*x*2**x
fit O(x) "g_time.dat" u 1:2 via B

plot[5:19] "g_time.dat" u 1:2 t "Meritve", O(x) t "$\\mathcal{O}(N\\log N)$"


set output "g_energija.tex"

set xlabel "$\\beta$"
unset ylabel

plot "g_FE.dat" u 1:2:3 w errorbars t "Prosta energija $F(\\beta)$"
, '' u 1:4:4 w errorbars t "Povpre"cna energija $\\langle H \\rangle_\\beta$"