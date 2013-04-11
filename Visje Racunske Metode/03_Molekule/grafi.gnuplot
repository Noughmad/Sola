set terminal epslatex color solid

set key top left
set ylabel "$T_i$"
set xlabel "$i$"

set output "g_hoover_lambda_10.tex"
plot for [lambda in "0 0.01 0.1 0.3 0.6 1"] sprintf("g_hoover_10_%s.dat", lambda) u 0:3 w lp t sprintf("$\\lambda = %s$", lambda)

set output "g_hoover_lambda_50.tex"
plot for [lambda in "0 0.01 0.1 0.3 0.6 1"] sprintf("g_hoover_50_%s.dat", lambda) u 0:3 w l t sprintf("$\\lambda = %s$", lambda)

set output "g_odv_lambda_30.tex"
set key top left
plot for [lambda in "0 0.01 0.1 0.3 0.6 1"] sprintf("g_maxwell_30_%s.dat", lambda) u 0:3 w lp t sprintf("$\\lambda = %s$", lambda)

set output "g_odv_lambda_100.tex"
set key top left
plot for [lambda in "0 0.01 0.1 0.3 0.6 1"] sprintf("g_maxwell_100_%s.dat", lambda) u 0:3 w l t sprintf("$\\lambda = %s$", lambda)

set output "g_odv_velikost.tex"
plot for [N in "10 20 30 50 100"] sprintf("g_maxwell_%s_0.3.dat", N) u (1+$0*10/(N-1)):3 w l t sprintf("$N = %s$", N)

set xlabel "$\\lambda$"
set ylabel "$T(1/5)$"
set key bottom left
set output "g_prehod.tex"
plot for [N in "10 20 30 50 100 200"] sprintf("g_prehod_%s.dat", N) w l t sprintf("$N = %s$", N)

set xlabel "$i$"
set ylabel "$J_i$"

set output "g_odv_lambda_30_tok.tex"
set key top left
plot for [lambda in "0 0.01 0.1 0.3 0.6 1"] sprintf("g_maxwell_30_%s.dat", lambda) u 0:1 w lp t sprintf("$\\lambda = %s$", lambda)

set output "g_odv_lambda_100_tok.tex"
set key top left
plot for [lambda in "0 0.01 0.1 0.3 0.6 1"] sprintf("g_maxwell_100_%s.dat", lambda) u 0:1 w l t sprintf("$\\lambda = %s$", lambda)

set output "g_odv_velikost_tok.tex"
plot for [N in "10 20 30 50 75 100"] sprintf("g_maxwell_%s_0.1.dat", N) u (1+$0*10/(N-1)):1 w l t sprintf("$N = %s$", N)

set ylabel "$\\overline{J}$"
set xlabel "$\\lambda$"

set key top right

set output "g_tok_lambda.tex"
plot for [N in "30 100"] sprintf("g_maxwell_tok_lambda_%s.dat",N) w lines t sprintf("$N = %s$", N)

set xlabel "$N$"
set key at 90,0.14

set output "g_tok_N.tex"
plot for [lambda in "0 0.1 1 2"] sprintf("g_maxwell_tok_N_%s.dat", lambda) w lines t sprintf("$\\lambda = %s$", lambda)

set terminal epslatex color solid size 18cm,8.5cm
set key outside
set output "g_tok_N_log.tex"
set logscale xy

A = 6
p = -1
f(x) = A * x**p
B = 0.1
q = -1
g(x) = B * x**q
fit f(x) "g_maxwell_tok_N_1.dat" via A,p
fit g(x) "g_maxwell_tok_N_2.dat" via B,q
plot for [lambda in "0 0.1"] sprintf("g_maxwell_tok_N_%s.dat", lambda) w lines t sprintf("$\\lambda = %s$", lambda), \
"g_maxwell_tok_N_1.dat" w lines t "$\\lambda = 1$", f(x) t sprintf("$p_1 = %.2g$", -p), \
"g_maxwell_tok_N_2.dat" w lines t "$\\lambda = 1$", g(x) t sprintf("$p_2 = %.2g$", -q)

