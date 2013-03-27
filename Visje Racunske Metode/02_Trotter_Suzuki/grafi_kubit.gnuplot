set terminal epslatex color solid
set datafile separator ","

set output "g_casovna_zahtevnost.tex"
set xlabel "\"Stevilo delcev $n$"
set ylabel "\"Cas ra\"cunanja [s]"
set logscale y

A = 0.14
O(x) = A + B*x*2**x
fit O(x) "g_time_brez.dat" u 1:2 via A,B

A2 = 0.1
O2(x) = A2 + B2*x*2**x
fit O2(x) "g_time_prem.dat" u 1:2 via B2

set key bottom right
plot[5:19] "g_time_brez.dat" u 1:2 t "Vsaka posebej", O(x) t "$\\mathcal{O}(N\\log N)$", \
"g_time_prem.dat" u 1:2 t "Vnaprej zmno\"zene", O2(x) t "$\\mathcal{O}(N\\log N)$"

set output "g_prosta_energija.tex"
unset logscale
set key top right

set xlabel "$\\beta$"
set ylabel "$F(\\beta)$"
set logscale x

E(x) = Ae + Be / (1 + Ce*x)
F(x) = Af + Bf / (1 + Cf*x)

fit E(x) "g_FE.dat" u 1:4:5 via Ae, Be, Ce
fit F(x) "g_FE.dat" u 1:2:3 via Af, Bf, Cf

plot "g_FE.dat" u 1:2:3 w errorbars t "Izra\"cun", F(x) t "Fit"

set output "g_energija.tex"
unset logscale
set key top right

set xlabel "$\\beta$"
set ylabel "$\\langle H \\rangle_\\beta$"
set logscale x

plot "g_FE.dat" u 1:4:5 w errorbars t "Izra\"cun", E(x) t "Fit"

set output "g_energija_obe.tex"

plot "g_FE.dat" u 1:2:3 w errorbars t "Prosta energija -- izra\"cun", \
F(x) t "Prosta energija -- fit", \
"g_FE.dat" u 1:(-$4):5 w errorbars t "Energija -- izra\"cun", \
(-E(x)) t "Energija -- fit"

unset logscale
# set logscale xy

set output "g_korelacija_mag.tex"

set xlabel "$t$"
set ylabel "$\\langle \\sigma_1^z(t)\\sigma_1^z(0)\\rangle$"

plot "g_korelacija.dat" u 1:2 w lines t "Izra\"cun"

set output "g_korelacija_tok.tex"

set xlabel "$t$"
set ylabel "$\\langle J(t)J(0)\\rangle$"

B = 1.5
k = 7.3
q = 0.0
p = 0.4
C(x) = B*cos(k*x)*exp(-p*x)

fit[1:] C(x) "g_tok.dat" u 1:2:($3/sqrt(100)) via B,k,p
set samples 300
plot[0:10] "g_tok.dat" u 1:2 w lines t "Izra\"cun", C(x) t "$C(t) \\propto \\cos(\\omega t)\\exp(-t/\\tau)$"

set output "g_tok_log.tex"
set logscale y
plot[][1e-5:] "g_tok_long.dat" u 1:(abs($2)) w lines notitle

set logscale y

set output "g_fft_mag.tex"
set xlabel "$\\omega$"
set ylabel "$|S(\\omega)|$"

plot[0:5] "g_fft_mag.dat" u 1:2 w lines notitle


set output "g_fft_tok.tex"
set xlabel "$\\omega$"
set ylabel "$|S(\\omega)|$"

plot[0:2.5] "g_fft_tok.dat" u 1:2 w lines notitle

