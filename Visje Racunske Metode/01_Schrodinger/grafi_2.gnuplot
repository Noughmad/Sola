set terminal epslatex color solid
set style data lines
set datafile separator ","

set output "g_konvergenca_ho_1.tex"
plot for [i=3:22] "g_konvergenca_ho_1.csv" u 1:i notitle

set output "g_konvergenca_ho_01.tex"
plot for [i=3:22] "g_konvergenca_ho_0.1.csv" u 1:i notitle

set output "g_konvergenca_L_1.tex"
plot for [i=3:22] "g_konvergenca_L_1.csv" u 1:i notitle

set output "g_konvergenca_L_01.tex"
plot for [i=3:22] "g_konvergenca_L_0.1.csv" u 1:i notitle

set output "g_konvergenca.tex"
plot "g_konvergenca_ho_0.1.csv" u 1:2 t "HO, $\\lambda = 0.1$", \
"g_konvergenca_ho_1.csv" u 1:2 t "HO, $\\lambda = 1$", \
"g_konvergenca_L_0.1.csv" u 1:2 t "Lanczos, $\\lambda = 0.1$", \
"g_konvergenca_L_1.csv" u 1:2 t "Lanczos, $\\lambda = 1$"
