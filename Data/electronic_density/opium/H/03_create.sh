#!/bin/bash

cat h_header.data > H.beta
cat h_rho.data | gawk '{printf "%.12e\n", $2}' >> H.beta
cat h_beta.data | gawk '{printf "%.12e\n", $2}' >> H.beta


