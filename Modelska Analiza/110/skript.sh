#!/bin/sh

S=1000000

for t in 1 0.5 0.1 0.01
do
	for n in 25 250 
	do
		./Populacije/build/populacije g_iz_exp_${t}_${n}.dat $t $S izumrtje exp $n
                ./Populacije/build/populacije g_iz_rs_${t}_${n}.dat $t $S izumrtje rs $n
	done

#	for l in 40 50 60
#	do
#		./Populacije/build/populacije g_zl_${t}_200_${l}.dat $t $S zajlis 200 $l
#	done
done

