#!/bin/bash

echo "O2 version"
for i in {1..10} ; do ./test_O2 ; done | grep Elapsed | gawk 'BEGIN {a=0.0;n=0} {a=a+$3; n++; print $3, a, n} END {a=a/n; print "Average ", a}'

echo "O3 version"
for i in {1..10} ; do ./test_O3 ; done | grep Elapsed | gawk 'BEGIN {a=0.0;n=0} {a=a+$3; n++; print $3, a, n} END {a=a/n; print "Average ", a}'

echo "Ofast version"
for i in {1..10} ; do ./test_Ofast ; done | grep Elapsed | gawk 'BEGIN {a=0.0;n=0} {a=a+$3; n++; print $3, a, n} END {a=a/n; print "Average ", a}'

echo "ndebug version"
for i in {1..10} ; do ./test_ndebug ; done | grep Elapsed | gawk 'BEGIN {a=0.0;n=0} {a=a+$3; n++; print $3, a, n} END {a=a/n; print "Average ", a}'

echo "unsafe version"
for i in {1..10} ; do ./test_unsafe ; done | grep Elapsed | gawk 'BEGIN {a=0.0;n=0} {a=a+$3; n++; print $3, a, n} END {a=a/n; print "Average ", a}'

