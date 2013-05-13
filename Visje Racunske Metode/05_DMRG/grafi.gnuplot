set terminal epslatex color solid
set datafile separator ","

set output "g_entropija.tex"
set xlabel "$n_A$"
set ylabel "$E$"

plot "g_entent_random.dat" w errorbars t "Naklju\"cno stanje", \
"" w lines lt 1 notitle, \
"g_entent_ground.dat" w lines t "Osnovno stanje"