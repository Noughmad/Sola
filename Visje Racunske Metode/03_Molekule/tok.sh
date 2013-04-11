#!/bin/bash

for N in 100 30
do
    F="g_maxwell_tok_lambda_${N}.dat"
    rm -f ${F}
    for lambda in 0 0.01 0.1 0.3 0.6 1
    do
        D="\$1"
        Avg=`awk "NR>1 && NR<${N} {sum += ${D}} END {printf sum / (${N}-2)}" g_maxwell_${N}_${lambda}.dat`
        echo "${lambda}, ${Avg}" >> ${F}
    done
done

    for lambda in 0 0.1 1 2
do
    F="g_maxwell_tok_N_${lambda}.dat"
    rm -f ${F}
    for N in 10 20 30 50 75 100
    do
        Avg=`awk "NR>1 && NR<${N} {sum += \\$1} END {printf sum / (${N}-2)}" g_maxwell_${N}_${lambda}.dat`
        echo "${N}, ${Avg}" >> ${F}
    done
done
