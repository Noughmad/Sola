set terminal epslatex color solid
set datafile separator ","

set output "g_entropija.tex"
set xlabel "$n_A$"
set ylabel "$E$"

plot "g_entent_random.dat" w errorbars t "Naklju\"cno stanje", \
"" w lines lt 1 notitle, \
"g_entent_ground.dat" w lines t "Osnovno stanje"

set output "g_velikost.tex"
set xlabel "$n$"

plot "g_entent_n.dat" u 1:2 w lines t "$n_A = 1$", \
"g_entent_n.dat" u 1:3 w lines t "$n_A = n/2$"

set output "g_nekompaktna.tex"
set xlabel "$k = n_A / n$"

set xrange [1.5:6.5]

f(x) = a/x**0.75
a = 2
fit f(x) "g_entent_nonc.dat" u 1:2 via a
plot "g_entent_nonc.dat" u 1:2 w points t "Podatki", f(x) t "$E \\propto k^{-3/4}$"