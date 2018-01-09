# USER-EPH
LAMMPS extension to capture electron-ion interaction

Artur Tamm and Alfredo A. Correa (LLNL)

## Installation Instructions

Get LAMMPS (source code)
```
git clone https://github.com/lammps/lammps.git lammps
```

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
fix ID group-ID eph seed flags model C_e rho_e kappa_e T_e NX NY NZ T_infile N T_map beta_infile A B C
```
Where:
* seed -> seed for random number generator [integer]
* flags -> control of different terms or'd together [integer]
** 1 -> enable friction
** 2 -> enable random force
** 4 -> enable FDM
* model -> select model for friction and random force [integer]
** 1 -> Standard Langevin
** 2 -> Tamm2016 PRB paper version
** 3 -> New version with CM correction only
** 4 -> New version full model
* rho_e -> electronic density for the FDM grid [float]
* C_e -> electronic heat capacity per density [float]
* kappa_e -> electronic heat conduction [float]
* T_e -> electronic temperature [float]
* NX, NY, NZ -> grid size in x, y, and z direction [integer]
* T_infile -> Input file for the FDM grid parameters and initial values [string]
* N -> heat map output frequency [integer]
* T_map -> output heat map file name [string]
* beta_infile -> beta(rho) input file [string]
* A, B, C -> element type mapping [string]

This fix produces two types of output:
* vector with the energy and temperature of the electronic system
** f_ID[1] -> 
** f_ID[2]
* per atom values for 

## Notes
* There is an incompatibility between model 4 and models 1, 2, 3, where the same beta(rho) function cannot be used. This is due to angular correction.
* C_e, rho_e, and kappa_e are constants.
* If T_infile is not NULL then C_e, rho_e, kappa_e, T_e, NX, NY, NZ are ignored and are read from the filename supplied. If NULL is provided as the file name then the FDM grid is initialised with the parameters provided in the command.
* In order to disable output of the electronic temperature zero can be provided as the value of N.



