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
SRC=(${SRC[*]} "Examples/Beta/Ni.beta" "Examples/Beta/Ni_model_4.beta")
SRC=(${SRC[*]} "Examples/Example_1/Ni.eam" "Examples/Example_1/Ni_model_4.beta" "Examples/Example_1/run.lmp" "Examples/Example_1/README" "Examples/Example_1/Tout.pdf" "Examples/Example_1/Tout.png")

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

