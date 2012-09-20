#!/bin/sh

OUTPUT="g_min_hitrost.dat"
rm -f ${OUTPUT}
for MU in 0.0 0.001 0.003 0.01 0.02 0.05 0.1 0.2;
do
	echo "${MU} `cat g_{lambert,zoom}_${MU}_*.dat | awk 'BEGIN { min=1000 } $5 < min { min=$5 }; END {print min}'`" >> ${OUTPUT};
done
