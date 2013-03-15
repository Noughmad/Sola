#!/bin/bash

for L in 0 0.01 0.1 1 5 10; do
    ./build/trotter vse ${L} 100000 2000 500 > g_vse_${L}.dat
    ./build/trotter eq ${L} 10000 2000 50 > g_eq_${L}.dat
done