set terminal epslatex color solid
set style data lines
set datafile separator ","

set xlabel "$N$"
set ylabel "$E_n^{(\\lambda)}$"

set output "g_konvergenca_ho_01.tex"
plot[20:100] for [i=3:22] "data/g_konvergenca_ho_0.1.csv" u 1:i notitle

set output "g_konvergenca_ho_1.tex"
plot[20:200][0:500] for [i=3:22] "data/g_konvergenca_ho_1.csv" u 1:i notitle

set output "g_konvergenca_L_1.tex"
plot[20:275][0:1000] for [i=3:22] "data/g_konvergenca_L_1.csv" u 1:i notitle

set output "g_konvergenca_L_01.tex"
plot[20:275][0:1000] for [i=3:22] "data/g_konvergenca_L_0.1.csv" u 1:i notitle

set output "g_konvergenca.tex"
set key bottom right

h0(x) = a*x + A
h1(x) = b*x + B
l0(x) = c*x + C
l1(x) = d*x + D

fit[0:60] h0(x) "data/g_konvergenca_ho_0.1.csv" u 1:2 via a,A
fit[0:110] h1(x) "data/g_konvergenca_ho_1.csv" u 1:2 via b,B

plot[20:120][0:20] "data/g_konvergenca_ho_0.1.csv" u 1:2 t "HO, $\\lambda = 0.1$", h0(x) t sprintf("$r = %.2g$", a), \
"data/g_konvergenca_ho_1.csv" u 1:2 t "HO, $\\lambda = 1$", h1(x) t sprintf("$r = %.2g$", b)

set output "g_energije.tex"
set key top left

A = 0.5
B = 1
p = 1.1

f(x) = A+B*(x**p)
fit[0:95] f(x) "data/g_energije.csv" u 1:11 via A,B,p
set xlabel "$n$"
set ylabel "$E_n^{(\\lambda)}$"
plot "data/g_energije.csv" u 1:2 t "$\\lambda = 0$", \
'' u 1:5 t "$\\lambda = 0.001$", \
'' u 1:6 t "$\\lambda = 0.003$", \
'' u 1:7 t "$\\lambda = 0.01$", \
'' u 1:8 t "$\\lambda = 0.03$", \
'' u 1:9 t "$\\lambda = 0.1$", \
'' u 1:10 t "$\\lambda = 0.3$", \
'' u 1:11 t "$\\lambda = 1$"
