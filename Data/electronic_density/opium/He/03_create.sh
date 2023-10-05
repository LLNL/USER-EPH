#!/bin/bash

cat he_header.data > He.beta
cat he_rho.data | gawk '{printf "%.12e\n", $2}' >> He.beta
cat he_beta.data | gawk '{printf "%.12e\n", $2}' >> He.beta


