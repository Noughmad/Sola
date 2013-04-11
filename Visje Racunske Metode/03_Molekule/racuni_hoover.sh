#!/bin/bash

I=1e6
A=1e6
M=10

for N in 10 50
do
    for lambda in 0 0.01 0.1 0.3 0.6 1
    do
        echo "Calculating for N = ${N}, lambda = ${lambda}"
        ./kopeli/build/kopeli hoover $N $I $A $M $lambda > g_hoover_${N}_${lambda}.dat
    done
done