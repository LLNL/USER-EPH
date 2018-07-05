#!/bin/bash
echo "Updating USER-EPH private repo"
IN=~/Work/01_Electronic_Stopping/Software/USER-EPH
cd $IN
git pull

OUT=~/Build/lammps_LLNL/src/USER-EPH
cd $OUT
git pull

echo "Syncing files in private repo with public repo"
SRC=("README.md" "LICENSE" "eph_spline.h" "eph_spline.cpp" "eph_beta.h" "eph_beta.cpp" "eph_fdm.h" "eph_fdm.cpp" "fix_eph.h" "fix_eph.cpp")
SRC=(${SRC[*]} "Doc/Beta/input.beta" "Doc/FDM/T_input.fdm")

# documentation on input files
SRC=(${SRC[*]} "Examples/Beta/Ni.beta" "Examples/Beta/Ni_model_4.beta")
# Example 1
SRC=(${SRC[*]} "Examples/Example_1/Ni.eam" "Examples/Example_1/Ni_model_4.beta" "Examples/Example_1/run.lmp" "Examples/Example_1/README" "Examples/Example_1/Tout.pdf" "Examples/Example_1/Tout.png")
# Example 2
SRC=(${SRC[*]} "Examples/Example_2/Ni.eam" "Examples/Example_2/Ni_model_4.beta" "Examples/Example_2/run.lmp" "Examples/Example_2/README" "Examples/Example_2/Tout.pdf" "Examples/Example_2/Tout.png")
# Example 3
SRC=(${SRC[*]} "Examples/Example_3/Ni.eam" "Examples/Example_3/Ni_model_4.beta" "Examples/Example_3/run.lmp" "Examples/Example_3/README" "Examples/Example_3/Tout.pdf" "Examples/Example_3/Tout.png")
# Example 4
SRC=(${SRC[*]} "Examples/Example_4/Ni.eam" "Examples/Example_4/Ni_model_4.beta" "Examples/Example_4/run.lmp" "Examples/Example_4/T.in" "Examples/Example_4/README" "Examples/Example_4/Tfieldout.pdf" "Examples/Example_4/Tfieldout.png")

N=${#SRC[*]}
echo $N

for i in $(seq 0 $[N-1]) ; do
  file=${SRC[$i]}
  ifile=$(md5sum $IN/$file | gawk '{print $1}')
  ofile=$(md5sum $OUT/$file | gawk '{print $1}')

  echo $file $ifile $ofile
  if [ "$ifile" != "$ofile" ] ; then
    cp -v $IN/$file $OUT/$file
    git add $file
  fi

  done

cd $OUT
git commit
git push

