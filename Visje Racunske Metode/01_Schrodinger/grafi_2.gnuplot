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

h0(x) = a*x + A
h1(x) = b*x + B
l0(x) = c*x + C
l1(x) = d*x + D

fit[0:80] h0(x) "g_konvergenca_ho_0.1.csv" u 1:2 via a,A
fit[0:110] h1(x) "g_konvergenca_ho_1.csv" u 1:2 via b,B

plot[][0:20] "g_konvergenca_ho_0.1.csv" u 1:2 t "HO, $\\lambda = 0.1$", h0(x) t sprintf("$k = %.2g$", a), \
"g_konvergenca_ho_1.csv" u 1:2 t "HO, $\\lambda = 1$", h1(x) t sprintf("$k = %.2g$", b), \
"g_konvergenca_L_0.1.csv" u 1:2 t "Lanczos, $\\lambda = 0.1$", \
"g_konvergenca_L_1.csv" u 1:2 t "Lanczos, $\\lambda = 1$"
