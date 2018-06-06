#!/bin/bash

for i in {1..10} ; do ./test ; done | grep Elapsed | gawk 'BEGIN {a=0.0;n=0} {a=a+$3; n++; print $3, a, n} END {a=a/n; print "Average ", a}'

