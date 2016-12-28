#!/bin/bash

./bin/makefield.out
./bin/main.out

cat log | perl -lne 'if(/(.*\/)psi(.*?\.bin)/) { print "$_=>$1pres$2"; }' | ./bin/invert_pres.out
