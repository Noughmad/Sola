set terminal epslatex color solid

set output "g_test_vodik.tex"
set xlabel "$r$"
set ylabel "$u(r)$"
plot[0:10] "data/g_test_vodik.dat" w lines notitle

p(x) = -(x + 1) * exp(-2*x) + 1

set output "g_test_vodik_pot.tex"
set ylabel "$U(r)$"
plot[0:5] "data/g_test_vodik_potencial.dat" t "Izra\"cun", p(x) t "Napoved"

set output "g_helij_e.tex"
set xlabel "$r$"
set ylabel "$u(r)$"
plot[0:10] "data/g_helij_e.dat" w lines notitle

p(x) = -(x + 1) * exp(-2*x) + 1

set output "g_helij_pot.tex"
set ylabel "$U(r)$"
plot[0:5] "data/g_helij_pot.dat" t "Izra\"cun", p(x) t "Napoved (vodik)"
