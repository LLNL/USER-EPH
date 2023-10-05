#!/bin/bash

XX=H # Element name with capital
xx=h # Element name non-cap
ZZ=1 # Z of the element

sed -i -e "s/xx/${xx}/g" 01_extract.sh
sed -i -e "s/XX/${XX}/g" \
       -e "s/xx/${xx}/g" \
       -e "s/ZZ/${ZZ}/g" 02_main.m
sed -i -e "s/XX/${XX}/g" \
       -e "s/xx/${xx}/g" 03_create.sh
sed -i -e "s/xx/${xx}/g" \
       -e "s/XX/${XX}/g" 04_figure.plot
sed -i -e "s/XX/${XX}/g" \
       -e "s/ZZ/${ZZ}/g" xx_header.data
sed -i -e "s/XX/${XX}/g" \
       -e "s/ZZ/${ZZ}/g" xx.param

mv xx_beta.data ${xx}_beta.data
mv xx_header.data ${xx}_header.data
mv xx.param ${xx}.param


