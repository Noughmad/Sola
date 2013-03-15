set terminal epslatex color solid
set style data lines

# OHRANITEV ENERGIJE

set output "g_energija_0.tex"
plot "g_vse_0.dat" u 1:4 t "\\texttt{RK4}", '' u 1:7 t "$S_2$", '' u 1:10 t "$S_4$"

set output "g_energija_001.tex"
plot "g_vse_0.01.dat" u 1:4 t "\\texttt{RK4}", '' u 1:7 t "$S_2$", '' u 1:10 t "$S_4$"

set output "g_energija_01.tex"
plot "g_vse_0.1.dat" u 1:4 t "\\texttt{RK4}", '' u 1:7 t "$S_2$", '' u 1:10 t "$S_4$"

set output "g_energija_1.tex"
plot "g_vse_1.dat" u 1:4 t "\\texttt{RK4}", '' u 1:7 t "$S_2$", '' u 1:10 t "$S_4$"

set output "g_energija_5.tex"
plot "g_vse_5.dat" u 1:4 t "\\texttt{RK4}", '' u 1:7 t "$S_2$", '' u 1:10 t "$S_4$"

set output "g_energija_10.tex"
plot "g_vse_10.dat" u 1:4 t "\\texttt{RK4}", '' u 1:7 t "$S_2$", '' u 1:10 t "$S_4$"

# EKVIPARTICIJSKI IZREK

set output "g_ekviparticija_0.tex"
plot "g_eq_0.dat" u 1:2 t "$\\langle p_1^2 \\rangle$", '' u 1:3 t "$\\langle p_2^2 \\rangle$"

set output "g_ekviparticija_001.tex"
plot "g_eq_0.01.dat" u 1:2 t "$\\langle p_1^2 \\rangle$", '' u 1:3 t "$\\langle p_2^2 \\rangle$"

set output "g_ekviparticija_01.tex"
plot "g_eq_0.1.dat" u 1:2 t "$\\langle p_1^2 \\rangle$", '' u 1:3 t "$\\langle p_2^2 \\rangle$"

set output "g_ekviparticija_1.tex"
plot "g_eq_1.dat" u 1:2 t "$\\langle p_1^2 \\rangle$", '' u 1:3 t "$\\langle p_2^2 \\rangle$"

set output "g_ekviparticija_5.tex"
plot "g_eq_5.dat" u 1:2 t "$\\langle p_1^2 \\rangle$", '' u 1:3 t "$\\langle p_2^2 \\rangle$"

set output "g_ekviparticija_10.tex"
plot "g_eq_10.dat" u 1:2 t "$\\langle p_1^2 \\rangle$", '' u 1:3 t "$\\langle p_2^2 \\rangle$"

# MEJA ZA EKVIPARTICIJO
# TODO
