#!/bin/bash

IN="xx" 
OU="xx"  

OPIUM=opium_pclinux_static38.exe 
$OPIUM ${IN}.param ${IN}.log all

L1=$(cat -n ${IN}.pcc_plt | grep '@' | head -n 1 | gawk '{print $1}')
L2=$(cat -n ${IN}.pcc_plt | grep '@' | tail -n 1 | gawk '{print $1}')

head -n $[L1-1] ${IN}.pcc_plt > ${OU}_core_den.data
cat ${IN}.pcc_plt | tail -n $[L2-L1] | head -n $[L2-L1-1] > ${OU}_valence_den.data

