#!/usr/bin/env gnuplot 

# This script creates a figure with two plots, both showing electronic density, 
#   but one is in log scale

set term pdfcairo enhanced size 24cm,12cm enhanced font "Sans Serif, 18" mono lw 2
set encoding utf8

set grid
set key off
set xtics
set ytics

r_cut = 5.5

set xlabel "distance [Å]"
set ylabel "electronic density [e/Å^3]"

set output "xx.pdf"
set multiplot title "XX density profile" layout 1,2
set logscale y
plot [0:r_cut] \
  "xx_rho.data" u 1:2 w l lw 2 lc -1

unset logscale y
unset ylabel
plot [0:r_cut] \
  "xx_rho.data" u 1:2 w l lw 2 lc -1
unset multiplot
