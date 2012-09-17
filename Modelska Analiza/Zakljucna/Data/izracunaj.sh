#!/bin/sh

MU=$1
DELTA=$2

for TYPE in direct spiral spline orbit;
do
    OUTPUT="g_${TYPE}_${MU}_${DELTA}.dat"
    rm -f ${OUTPUT}
    for T in {1..50}; 
    do
	echo "Starting ${TYPE}, ${T}"
        ../Planetki/build/planetki $MU $DELTA $T 200 $TYPE 2&>> ${OUTPUT};
    done
done
