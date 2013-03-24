#!/bin/bash

for L in 0 1 2 3 1.16 1.18 1.2 4 5; do
 #   ./build/trotter vse ${L} 100000 2000 500 > g_vse_${L}.dat
    ./build/trotter eq ${L} 10000 2000 50 > g_eq_${L}.dat
done

gnuplot grafi.gnuplot
pdflatex miha_cancula_3.tex