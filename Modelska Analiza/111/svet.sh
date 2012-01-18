#!/bin/sh
SIZES="2000 5000"
for i in $SIZES
do
    NAME=g_svet_${i}.dat
    rm -rf $NAME
    python svet.py $i 3 > $NAME
    echo "Koncal za N=${i}"
done