#!/bin/bash

for FILE in $* ; do

#  BAK=${FILE}~
#  if [ "$TMP" = "" ]; then TMP="."; fi; BAK=$TMP\fixbb__.ps
  BAK='fixbb.eps'

  mv $FILE $BAK
  eps2eps --ignoreBB $BAK $FILE

  # delete the backup file
  rm -f $BAK

done

