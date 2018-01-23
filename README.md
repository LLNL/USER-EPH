# USER-EPH

LAMMPS extension ("fix") to capture electron-ion interaction

Artur Tamm and Alfredo A. Correa (LLNL)

## Installation Instructions

Get LAMMPS (source code)
```
$ cd mywork
$ git clone https://github.com/lammps/lammps.git
```

Get USER-EPH (this plugin, you have to have access to the repository)
```
$ cd lammps/src
$ git clone https://github.com/artuuuro/USER-EPH.git
```

Edit `Makefile` add string ` user-eph ` to the end of `PACKUSER` variable (near line 66).

Edit `MAKE/Makefile.mpi` and `MAKE/Makefile.serial` to read `CCFLAGS = -g -O3 -std=c++11` (near line 10).

Execute:
```
$ make yes-manybody yes-user-eph
$ make -j 8 serial
$ make -j 8 mpi
```

The executables are `./lmp_mpi` (for parallel runs) `./lmp_serial` (for serial runs, testing), you can copy them elsewhere.

## Usage

* Take you MD input file
* Add a line at the correct place, 
```
fix ID group-ID eph seed flags model C_e rho_e kappa_e T_e NX NY NZ T_infile N T_map beta_infile A B C
```
Where:
* `seed` -> seed for random number generator [integer, for example 123]
* `flags`: control of different terms or'd together [integer]
  * `1` -> enable friction
  * `2` -> enable random force
  * `4` -> enable FDM
* `model`: select model for friction and random force [integer]
  * `1` -> Standard Langevin
  * `2` -> Tamm2016 PRB paper version
  * `3` -> New version with CM correction only
  * `4` -> New version full model (Tamm2018 paper)
* `rho_e` -> scaling parameter for the FDM grid [float] [unitless]
* `C_e` -> electronic heat capacity per density [float] [in eV/K/Ang^3]
* `kappa_e` -> electronic heat conduction [float] [in eV/K/Ang/ps]
* `T_e` -> electronic temperature [float] [in K]
* `NX`, `NY`, `NZ` -> grid size in x, y, and z direction [integer, for example 1 1 1]
* `T_infile` -> Input file for the FDM grid parameters and initial values [string]
* `N` -> heat map output frequency [integer, for example 10]
* `T_map` -> output heat map file name [string]
* `beta_infile` -> beta(rho) input file [string]
* `A`, `B`, `C` -> element type mapping [string]

For example the following line will run MD including the coupling for electrons, within the spatially correlated Langevin bath.
The electronic specific heat is assumed to be 2.5e-6 eV/K/Ang^3 (400000 J/m³/K) (see LinPRB772008) which is a good approximation for a renge of electronic temperatures in the range of 500 to 1500K. 
Initial electron temperature is set to 300K (and not from a file).
We use uniform tempetures (one grid element), therefore the heat conductivity is not relevant in this case.

```
fix ID group-ID eph 123 4 4 2.5e-6 rho_e kappa_e 300 1 1 1 NULL 10 T_map.heat.out.dat Ni.betarho.in.dat A B C
```

This fix produces two types of output:
* vector with the energy and temperature of the electronic system
  * `f_ID[1]` -> Net energy transfer between electronic and ionic system
  * `f_ID[2]` -> Average electronic temperature
* per atom values:
  * `f_ID[i][1]` -> site density
  * `f_ID[i][2]` -> coupling parameter

### Beta(rho) input file

This file provides the electronic densities and beta(rho) functions for individual species.
The format is described in `Doc/Beta/input.beta`. 
The file is similar to eam setfl format. 
The beta(rho) function has units of [eV ps/Ang^2]. An example is provided in `Examples/Beta/Ni_model_4.beta`.

### FDM grid input file

This file is used to initialise FDM grid for the electronic system. 
The format is described in `Doc/FDM/T_input.fdm`. 
This allows fine control of the properties at various grid points. 
Also, the grid can be larger than ionic system if energy reservoir far away is needed. 
Grid points can act as energy sources, sinks or walls (wall feature is not implemented correctly). 
An example of input file is provided in `Examples/FDM/T.in`.

## Notes and limitations

* There is an incompatibility between model 4 and models 1, 2, 3, where the same beta(rho) function cannot be used. 
This is due to angular correction.
* `C_e`, `rho_e`, and `kappa_e` are constants at the moment.
* If `T_infile` is not `NULL` then `C_e`, `rho_e`, `kappa_e`, `T_e`, `NX`, `NY`, `NZ` are ignored and are read from the filename supplied. 
If `NULL` is provided as the filename then the FDM grid is initialised with the parameters provided in the command.
* In order to disable output of the electronic temperature zero can be provided as the value of N.

# Tutorial
Tutorials can be found in `Examples/` directory. To run them type `lmp_serial -i run.lmp` in the appropriate example directory and assuming executable is in `PATH`. Some of the examples may take long on older machines, so tweak the input file (`run.lmp`) accordingly. Every example contains a README file that describes what the runscript does.

## Example 1
`Examples/Example_1/`
In this example a crystal structure is created and the model is applied with both the friction and random force terms. The electrons are kept at constant temperature (300K). This example illustrates the thermalisation process from 0K through electron-ion interaction only.

## Example 2
`Examples/Example_2/`
This example illustrates the use cooling of the ionic systems due to electrons only. This means that only the friction term acts on atoms and removes energy. This is equivalent to having electrons at 0K.

## Example 3
`Examples/Example_3/`
In this example full model with electronic FDM grid is used. The crystal is created with 0K motion and is heated by electrons. During the simulation the electronic system will cool and ionic system heat. At equilibrium both systems have same temperature on average. Also, this example illustrates the automatic initialisation of the FDM grid with constant parameters. The electronic temperature
at various grid points is written to files. Final state of the grid is stored and can be reused in later simulations (`T.restart`).

## Example 4
`Examples/Example_4/`
This example uses reads the FDM grid parameters from a file (`T.in`). In this file a source term is added at grid point `(0 0 0)` representing for example a laser. During the simulation the whole system will heat and due to gradient in the electronic system forces acting on atoms at different grid points will 'feel' different random forces in magnitude.
