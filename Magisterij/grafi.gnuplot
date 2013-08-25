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
set xlabel "Frekvenca svetlobe $\\omega a/2\\pi c$"
set ylabel "Prepustnost $I/I_0$"
set yrange [0:0.4]

plot "Podatki/test_bandgap_12.dat" w l lw 5 title "$\\varepsilon_1 = 13,\\, \\varepsilon_2 = 12$", \
"Podatki/test_bandgap_1.dat" w l lw 5 title "$\\varepsilon_1 = 13,\\, \\varepsilon_2 = 1$"

set output "g_test_periodic_en.tex"

set xtics
set xlabel "Light frequency $\\omega a/2\\pi c$"
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

plot "Podatki/absorption_0.dat" u 1:($2*2) lw 5 t "$p=0$", \
"Podatki/absorption_1.dat" u 1:($2*2) lw 5 t "$p=1$", \
"Podatki/absorption_2.dat" u 1:($2*2) lw 5 t "$p=2$", \
"Podatki/absorption_3.dat" u 1:($2*2) lw 5 t "$p=3$", \
"Podatki/absorption_4.dat" u 1:($2*2) lw 5 t "$p=4$"

reset
set terminal png size 800,600 crop
set pm3d map interpolate 2,1
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
set terminal eps color

set output "g_refraction_test.eps"
set pm3d map interpolate 2,1

unset tics
unset border
unset colorbox

I = 170.0
angle = asin(6.0 * 24.0 / I)
EMF_K = 300.0 + 34.0 + 34.0
LINE_X = EMF_K * 0.5
LINE_Y = 80
ARROW_SIZE = 80
refracted = asin(cos(angle))

set size ratio I*1.0/EMF_K

set obj 5 circle arc [180:180+angle*180/3.141592] fs transparent solid 0 fc rgb "black" lw 3
set obj 5 circle at LINE_X,LINE_Y size 37.7 front
set obj 6 circle arc [0:refracted*180/3.141592] fs transparent solid 0 fc rgb "black" lw 3
set obj 6 circle at LINE_X,LINE_Y size 50 front


set arrow 3 from EMF_K/2-100,LINE_Y to EMF_K/2+100,LINE_Y lw 4 front nohead
set arrow 1 from (EMF_K/2-ARROW_SIZE*cos(angle)),(LINE_Y-ARROW_SIZE*sin(angle)) to EMF_K/2,LINE_Y lw 6 lc rgb "white" front
set arrow 2 from EMF_K/2,LINE_Y to (EMF_K/2+ARROW_SIZE*cos(refracted)),(LINE_Y+ARROW_SIZE*sin(refracted)) lw 6 lc rgb "white" front
show arrow

splot "Podatki/brewster_refraction.dat" matrix notitle
