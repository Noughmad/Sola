#!/bin/sh

for L in 0 0.1 1; do
  ffmpeg -r 10 -i "g_animation_${L}_%03d.png" -q:v 0 "g_animation_$L.avi"
done