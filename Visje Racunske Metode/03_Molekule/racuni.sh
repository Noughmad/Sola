#!/bin/bash

I=1e7
A=1e7
M=100

for N in ""
do
    for lambda in 0 0.01 0.1 0.3 0.6 1
    do
        echo "Calculating for N = ${N}, lambda = ${lambda}"
        ./kopeli/build/kopeli maxwell $N $I $A $M $lambda > g_maxwell_${N}_${lambda}.dat
    done
done

for N in 10 20 30 50 75 100
do
    for lambda in 2 # 0 0.1 1
    do
        echo "Calculating for N = ${N}, lambda = ${lambda}"
        ./kopeli/build/kopeli maxwell $N $I $A $M $lambda > g_maxwell_${N}_${lambda}.dat
    done
done
