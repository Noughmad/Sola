#!/bin/sh

DELTA=$1

for MU in 0.0 0.01 0.02 0.05 0.1 0.2;
do
    OUTPUT="g_lambert_${MU}_${DELTA}.dat"
    rm -f ${OUTPUT}
    for T in {1..50}; 
    do
	echo "Starting ${MU}, ${T}"
        ../Planetki/build/planetki $MU $DELTA $T 500 0 2&>> ${OUTPUT};
    done
done
