#!/bin/sh

cd Jarzynski/build
make

# NOTE: Graf za 0.2 / 0.3 ze zgleda v redu
# ./jarzynski 0.2 0.3 &
./jarzynski c 0.1 &
./jarzynski c 0.3 &
./jarzynski 0.7 0.1