#!/bin/bash

I=1e8
A=1e7
M=100

for N in 10 20 30 50 100 200
do
    for lambda in 0 0.01 0.1 0.3 0.6 1
    do
        echo "Calculating for N = ${N}, lambda = ${lambda}"
        ./kopeli/build/kopeli maxwell $N $I $A $M $lambda > g_maxwell_${N}_${lambda}.dat
    #   ./kopeli/build/kopeli hoover $N $I $A $M $lambda > g_hoover_${N}_${lambda}.dat &
    done
done