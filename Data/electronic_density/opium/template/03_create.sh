#!/bin/bash

cat xx_header.data > XX.beta
cat xx_rho.data | gawk '{printf "%.12e\n", $2}' >> XX.beta
cat xx_beta.data | gawk '{printf "%.12e\n", $2}' >> XX.beta


