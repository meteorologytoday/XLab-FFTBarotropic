#!/bin/bash

makefield-Kuo2004.out
main.out

cat log | perl -lne 'if(/(.*\/)psi(.*?\.bin)/) { print "$_=>$1pres$2"; }' | invert_pres.out
ls -v output | grep pres_step | awk '{ print "output/" $1 }' | find_min.out 1> output/pres_timeseries.txt
