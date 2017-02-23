#!/bin/bash

fifo=vort_src_fifo
#./bin/makefield-elliptic-vortex.out
#./bin/makefield-const-vortex.out
#./bin/makefield-gaussian.out
makefield-Kuo2004.out

rm -f $fifo
mkfifo $fifo

vort_src_input.out > $fifo &
main.out -f${fifo}

cat log | perl -lne 'if(/(.*\/)psi(.*?\.bin)/) { print "$_=>$1pres$2"; }' | invert_pres.out
ls -v output | grep pres_step | awk '{ print "output/" $1 }' | ./bin/find_min.out 1> output/pres_timeseries.txt
