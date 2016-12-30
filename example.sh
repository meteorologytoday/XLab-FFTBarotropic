#!/bin/bash

#./bin/makefield-elliptic-vortex.out
./bin/makefield-const-vortex.out
#./bin/makefield-gaussian.out
./bin/main.out

cat log | perl -lne 'if(/(.*\/)psi(.*?\.bin)/) { print "$_=>$1pres$2"; }' | ./bin/invert_pres.out
ls output | grep pres_step | awk '{ print "output/" $1 }' | ./bin/find_min.out 1> output/pres_timeseries.txt

