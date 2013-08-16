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
plot uni(x) w l lw 5 lc 2 t "Napoved\\cite{kleman}", "Podatki/test_uniform.dat" w p pt 7 ps 1.2 lc 1 title "Rezultat simulacije"

set xlabel "$\\beta$"
set xtics ("0" 0, "$\\frac{\\pi}{2}$" PI/2, "$\\pi$" PI, "$\\frac{3\\pi}{2}$" 1.5*PI, "$2\\pi$" 2*PI)
set ylabel "Transmittance $I/I_0$"
set yrange [0:0.4]
set output "g_test_uniform_en.tex"
plot uni(x) w l lw 5 lc 2 t "Prediction\\cite{kleman}", "Podatki/test_uniform.dat" w p pt 7 ps 1.2 lc 1 title "Simulation"

reset
set output "g_test_periodic.tex"

set xtics
set xlabel "Frekvenca svetlobe $\\omega a/2\\pi$"
set ylabel "Prepustnost $I/I_0$"
set yrange [0:0.4]

plot "Podatki/test_bandgap_12.dat" w l lw 5 title "$\\varepsilon_1 = 13,\\, \\varepsilon_2 = 12$", \
"Podatki/test_bandgap_1.dat" w l lw 5 title "$\\varepsilon_1 = 13,\\, \\varepsilon_2 = 1$"

set output "g_test_periodic_en.tex"

set xtics
set xlabel "Light frequency $\\omega a/2\\pi$"
set ylabel "Transmittance $I/I_0$"
set yrange [0:0.4]

plot "Podatki/test_bandgap_12.dat" w l lw 5 title "$\\varepsilon_1 = 13,\\, \\varepsilon_2 = 12$", \
"Podatki/test_bandgap_1.dat" w l lw 5 title "$\\varepsilon_1 = 13,\\, \\varepsilon_2 = 1$"

reset
set output "g_test_absorption.tex"
set style data lines
set logscale y

set autoscale
set xlabel "Povpre\"cne izgube v plasti"
set ylabel "Odbojnost plasti"

set logscale xy
set key at 0.002,0.002

plot "Podatki/absorption_0.dat" lw 5 t "$p=0$", \
"Podatki/absorption_1.dat" u 1:2 lw 5 t "$p=1$", \
"Podatki/absorption_2.dat" u 1:2 lw 5 t "$p=2$", \
"Podatki/absorption_3.dat" u 1:2 lw 5 t "$p=3$", \
"Podatki/absorption_4.dat" u 1:2 lw 5 t "$p=4$"

reset
set terminal png size 800,600 crop
set pm3d map
unset colorbox
unset tics

set lmargin 0
set rmargin 0
set bmargin 0
set tmargin 0

set output "g_test_plane.png"
splot "Podatki/cross_plane.dat" matrix notitle

set output "g_test_plane_profile.png"
splot[][0:150] "Podatki/profile.dat" matrix notitle

reset
set terminal epslatex color solid

set output "g_refraction_test.tex"
set pm3d map

set xlabel "$z$" offset 0,1.8
set ylabel "$x$" offset 0,0
unset tics
unset colorbox
unset border
splot "Podatki/brewster_refraction.dat" matrix notitle
