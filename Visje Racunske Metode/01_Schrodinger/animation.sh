#!/bin/sh

rm -f *.avi *.mp4
for L in 0 0.1 1; do
  ffmpeg -r 10 -i "g_animation_2D_${L}_%03d.png" -q:v 0 "miha_cancula_1_2D_$L.mp4"
  ffmpeg -r 10 -i "g_animation_1D_${L}_%03d.png" -q:v 0 "miha_cancula_1_1D_$L.mp4"
done