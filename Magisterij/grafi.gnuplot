set terminal epslatex color solid
PI = 3.141592

set output "g_test_uniform.tex"
uni(x) = A * (sin(2*x))**2

# TODO: A se da izracunat, ne bi ga smel fittat
fit uni(x) "Podatki/test_uniform.dat" via A

set xlabel "$\\beta$"
set xtics ("0" 0, "$\\frac{\\pi}{2}$" PI/2, "$\\pi$" PI, "$\\frac{3\\pi}{2}$" 1.5*PI, "$2\\pi$" 2*PI)
set ylabel "Prepustnost $I/I_0$"
set yrange [0:0.4]
plot "Podatki/test_uniform.dat" title "Rezultat simulacije", uni(x) t "Napoved\\cite{kleman}"

reset
set output "g_test_periodic.tex"

set autoscale
set xtics
set xlabel "Frekvenca svetlobe $\\omega a/2\\pi$"
set ylabel "Prepustnost $I/I_0$"
set yrange [0:0.8]
plot "Podatki/test_bandgap_11.dat" w l title "$\\varepsilon_1 = 13,\\, \\varepsilon_2 = 11$", \
"Podatki/test_bandgap_1.dat" w l title "$\\varepsilon_1 = 13,\\, \\varepsilon_2 = 1$"

reset
set output "g_test_absorption.tex"
set style data lines
set logscale y

set autoscale
set xlabel "Povpre\"cne izgube v plasti"
set ylabel "Odbojnost plasti"

plot "Podatki/absorption.dat" u 1:2 t "$p=0$", \
"" u 1:3 t "$p=1$", \
"" u 1:4 t "$p=2$", \
"" u 1:5 t "$p=3$"