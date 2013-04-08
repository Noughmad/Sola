set terminal epslatex color solid

set output "g_odv_lambda_30.tex"
set key top left
plot for [lambda in "0 0.01 0.1 0.3 0.6 1"] sprintf("g_maxwell_30_%s.dat", lambda) u 0:3 w lp t sprintf("$\\lambda = %s$", lambda)

set output "g_odv_lambda_100.tex"
set key top left
plot for [lambda in "0 0.01 0.1 0.3 0.6 1"] sprintf("g_maxwell_100_%s.dat", lambda) u 0:3 w l t sprintf("$\\lambda = %s$", lambda)

set output "g_odv_velikost.tex"
plot for [N in "10 20 30 50 100"] sprintf("g_maxwell_%s_0.3.dat", N) u (1+$0*10/(N-1)):3 w l t sprintf("$N = %s$", N)

set output "g_prehod.tex"
plot for [N in "10 20 30 50 100 200"] sprintf("g_prehod_%s.dat", N) w l t sprintf("$N = %s$", N)
