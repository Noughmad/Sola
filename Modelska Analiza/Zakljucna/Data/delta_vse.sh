#!/bin/sh

for MU in 0.0 0.001 0.01 0.05 0.1 0.2;
do
    ./delta.sh ${MU};
done

