#!/bin/sh

MU=$1

for DELTA in -2.094 0 2.094;
do
    OUTPUT="g_zoom_${MU}_${DELTA}.dat"
    rm -f ${OUTPUT}
    for T in {1..50};
    do
	echo "Starting ${MU}, ${DELTA}, ${T}"
        ../Planetki/build/planetki $MU $DELTA `echo "scale=1; $T/5" | bc` 500 0 2&>> ${OUTPUT};
    done
done
