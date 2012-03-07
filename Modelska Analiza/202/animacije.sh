#!/bin/sh

mkdir temp &> /dev/null
cd temp
../Integrator/build/integrator
echo "Integrator koncal"
ls
for i in {0..62}
do
  cd zvezda_${i}
  rm -f g_zvezda-${i}.avi
  ffmpeg -r 16 -i g_frame_%03d.png g_zvezda-${i}.avi &> /dev/null
  cp g_zvezda-${i}.avi ../..
  cd ..
  echo "Koncal kodiranje za zvezda-${i}"
done
