# USER-EPH
LAMMPS extension to capture electron-ion interaction

Artur Tamm and Alfredo A. Correa (LLNL)

## Installation Instructions

Get LAMMPS (source code)

`git clone https://github.com/lammps/lammps.git lammps
Get USER-EPH (this plugin)

```
git clone https://github.com/artuuuro/USER-EPH.git lammps/src/USER-EPH`
cd lammps/src
```

Edit `Makefile` add string ` user-eph ` to the end of `PACKUSER` variable, near line 66.
Edit `MAKE/Makefile.mpi` and `MAKE/Makefile.serial` to read `CCFLAGS = -g -O3 -std=c++11` near line 10.

Execute
```
make yes-manybody
make yes-user-eph
make -j 8 serial
make -j 8 mpi
```

The executable are `lmp_mpi lmp_serial`, you can copy them elsewhere.

## Usage

* Take you MD input file
* Add a line at the correct place, 

```
fix ID group-ID eph seed ...
```
