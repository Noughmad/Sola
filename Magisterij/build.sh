#!/bin/bash

# First generated the plots
gnuplot grafi.gnuplot

# Then fix EPS files
TEMP="temp.eps"
for f in g_test_plane.eps g_test_plane_profile.eps g_refraction_test.eps
do
    mv $f $TEMP
    eps2eps $TEMP $f
done

rm -f $TEMP

# Generate line defect illustrations
# python Defekti/defects.py

pdflatex magisterij.tex
