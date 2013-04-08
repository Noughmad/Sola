#!/bin/bash

I=1e8
A=1e7
M=100

for N in 10 20 30 50
do
    P="g_prehod_${N}.dat"
    rm -f $P

    for lambda in 0 0.005 0.01 0.03 0.05 0.075 0.1 0.125 0.15
    do
        F="g_maxwell_${N}_${lambda}.dat"
        echo "Calculating for N = ${N}, lambda = ${lambda}"
        ./kopeli/build/kopeli maxwell $N $I $A $M $lambda > ${F}
        Value=`awk "NR==3*${N}/10" ${F} | awk '{print $3}'`
        echo "${lambda}, ${Value}" >> ${P}
    done

done

