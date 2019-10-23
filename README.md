# USER-EPH

LAMMPS extension (LAMMPS "fix") to capture electron-ion interaction and electronic stopping. LLNL-CODE-750832

Artur Tamm and Alfredo A. Correa (LLNL)

## Introduction

In LAMMPS, a "fix" is a plugin or extension to the main code that performs a specific operation to the atomistic system during
timestepping or minimization. 
We use this extension mechanism to generalize the two-temperature model to include electron-phonon coupling.
The theory behind this extension is designed to represent cascades, laser heating and equilibration and study energy transport with realistic electronic stopping power and electron-phonon coupling.
The theory is developed in the papers "Langevin dynamics with spatial correlations as a model for electron-phonon coupling" (https://dx.doi.org/10.1103/PhysRevLett.120.185501) and "Electron-phonon interaction within classical molecular dynamics" (https://link.aps.org/doi/10.1103/PhysRevB.94.024305).

## Installation Instructions

First you need a version of LAMMPS (source code)
```
$ mkdir mywork
$ cd mywork
$ git clone https://github.com/lammps/lammps.git
```

### Compile the USER-EPH fix for CPUs
 
Get USER-EPH (this plugin, you have to have access to the repository)
```
$ cd lammps/src
$ git clone https://github.com/LLNL/USER-EPH.git
```

Edit LAMMPS' `Makefile` add the string ` user-eph ` to the `PACKUSER` variable (near line 68), for example:

```Makefile
PACKUSER = user-adios user-atc user-awpmd user-bocs user-cgdna user-cgsdk user-colvars \
    ... 
    user-sdpd user-sph user-tally user-uef user-vtk user-yaff \
    user-eph
```

Edit `MAKE/Makefile.mpi` and `MAKE/Makefile.serial` and add `-std=c++11` to the `CCFLAGs` varialbe to read `CCFLAGS = -g -O3 -std=c++11` (near line 10).

Execute:
```bash
$ make yes-manybody yes-user-eph
$ make -j 8 serial
```

(You can also enable other packages as needed)

Make sure your MPI enviroment is setup already (`mpicxx` compiler wrapper works), this may require for example running `$ module load mpi/mpich-x86_64`

```
$ make -j 8 mpi
```

The executables are `./lmp_mpi` (for parallel runs) `./lmp_serial` (for serial runs, testing), you can copy them elsewhere.

### Compile for CUDA-enabled GPUs (optional)

The code is ported to GPUs, a CUDA toolkit is required to compile this version and a CUDA card(s) supporting architecture *at least* 6.0 (`sm_60`, like
[Pascal, Volta, Turing, etc](https://en.wikipedia.org/wiki/CUDA#GPUs_supported). 
The command `nvidia-smi` will give you details.

If not defined, set the CUDA environment variable (e.g. `/usr/local/cuda` or `/usr`)
```bash
$ echo $CUDA_HOME
$ export CUDA_HOME=/usr/local/cuda 
```

Go back to the LAMMPS GPU library directory (`cd ../../lammps/lib/gpu`) and modify the file `Makefile.linux.double` add and activate your CUDA architecture and if needed `CUDA_HOME` and `CUDA_INCLUDE`. For example (after line 10),

```Makefile
...
#CUDA_HOME = /usr/local/cuda
NVCC = nvcc -ccbin=cuda-c++   ### or, for example -ccbin=mpicxx

# Kepler CUDA
#CUDA_ARCH = -arch=sm_35
# Tesla CUDA
# CUDA_ARCH = -arch=sm_21
# newer CUDA
#CUDA_ARCH = -arch=sm_13
# older CUDA
#CUDA_ARCH = -arch=sm_10 -DCUDA_PRE_THREE
# Pascal (your architecture)
CUDA_ARCH = -arch=sm_60
...
```

and compile the library and return to the USER-EPH

```bash
$ make -f Makefile.linux.double
```

Go to the USER-EPH directory
```bash
$ cd ../../../lammps/src/USER-EPH/lib
```

Modify `Makefile` if needed (`CUDA_ARCH`, `CUDA_CODE`, `NVCCFLAGS`, `INCFLAGS`) and load necessary modules (e.g. `module load mpi/openmpi-x86_64`) buid

```bash
$ make
```

Go back to the source directory `lammps/src` and create a new file `MAKE/Makefile.mpi_gpu` 

```bash
$ cd ../..
$ cp MAKE/Makefile.mpi MAKE/Makefile.mpi_gpu
```

Modify the `CCFLAGS` and `LIB` variables in the new `MAKE/Makefile.mpi_gpu` to read

```Makefile
...
CCFLAGS =	-g -O3 -std=c++11 -DFIX_EPH_GPU
...
LIB = -L../USER-EPH/lib -leph_gpu -lcuda -lcudart
...
```

```bash
$ make yes-user-eph
$ make yes-gpu
$ make -j mpi_gpu
```

The GPU executable will be in `lmp_mpi_gpu`.

The use the GPU accelerated potentials you enable the GPU package when running LAMMPS
either by supplying it on the command line or through run scripts.

```bash
$ cd USER-EPH/Examples/Example_5
$ lmp_mpi_gpu -pk gpu 1 -sf gpu -i run.lmp
```

or

```
# LAMMPS input file run.lmp
package gpu 1
...
pair eam/alloy/gpu 
...
fix friction all eph/gpu
...
```
The `-pk gpu 1` and `-sf gpu` flags allow the use of a single lammps run script where
the last will try to substitute all fixes and pairs supporting gpu extension.

See also Example-5 for a possible input script using the GPU version. 

#### Note

The current GPU version of the fix is not multi-device aware, so in order to utilise all GPUs on a node, 
an equal number of tasks has to be used. Thus, `mpirun` has to be aware of gpus in order to assign correct GPUs 
to each task. As a workaround, GPUs can be set into a special mode to block multiple tasks running on one GPU card.

## Usage

* Take your MD input file
* Add a line at the correct place, 
```
fix [ID] [group-ID] eph [seed] [flags] [model] [rho_e] [C_e] [kappa_e] [T_e] [NX] [NY] [NZ] [T_infile] [freq] [Te_outfile] [beta_infile] [A] [B] [C...]
```
Where:

* `ID` -> user-assigned name for the fix, [string, e.g. `ephttm` or `friction`]
* `group-ID` -> group of atoms to which this fix will be applied, [string, e.g. `all`]
* `seed` -> seed for random number generator [integer, e.g. 123]
* `flags`: control of different terms or'd together [integer or bitmask]
  * `1` -> enable friction (only) (pure damping)
  * `2` -> enable random force (only) (not recommended)
  * `3` -> enable friction and random force (fixed e-temperature) 
  * `4` -> heat equation (by FDM finite difference method) (decoupled from ions)
  * `5` -> enable friction and heat equation (no feedback from e-)
  * `7` -> enable friction, random force, and heat equation (coupled e-ions)
* `model`: select model for friction and random force [integer]
  * `1` -> standard Langevin (for vanilla TTM with beta(rho))
  * `2` -> simple e-ph model (https://link.aps.org/doi/10.1103/PhysRevB.94.024305) (not recommended)
  * `3` -> e-ph with spatial correlations, with CM-correction only (https://arxiv.org/abs/1801.06610)
  * `4` -> e-ph with spatial correlations, full model (https://arxiv.org/abs/1801.06610)
* `rho_e` -> scaling parameter for the FDM grid [float, recommended `1.0`] [unitless]
* `C_e` -> electronic heat capacity per volume [float, e.g. `2.5e-6`] [in eV/K/Ang^3]
* `kappa_e` -> electronic thermal conductivity [float,  ignored for single grid point] [in eV/K/Ang/ps]
* `T_e` -> electronic temperature [float, e.g. `300`] [in K]
* `NX`, `NY`, `NZ` -> grid size in x, y, and z direction [integer, e.g. `1` `1` `1` sets single grid point]
* `T_infile` -> input filename for the FDM grid parameters and initial values [string or NULL]
* `freq` -> heat map output (`T_output`) frequency, `0` to disable [integer, e.g. `10`]
* `Te_output` -> output heat map filename (CUBE format) [string, e.g. `Te_output.cub`]
* `beta_infile` -> beta(rho) input filename [string, e.g. `NiFe.beta`]
* `A`, `B`, `C...` -> element type mapping [1 or more strings, `Ni Ni Fe`]
~/Lassen/09_Techbase/TB_Bench/02_Bench/02_GPU/TB_12/02_2_nodes
For example the following line in LAMMPS input script, 
will run the MD including the coupling to electrons, 
within the spatially correlated Langevin bath.
The electronic specific heat is assumed to be 2.5e-6 eV/K/Ang^3 (400000 J/m³/K) (see LinPRB772008) which is a good approximation for a range of electronic temperatures from 500 to 1500K. 
Initial electron temperature is set to 300K (and not from a file).
We use uniform tempetures (one grid element), therefore the heat conductivity is not relevant in this case.

```
fix ephttm all eph 123 7 4 1.0 2.5e-6 1.0 300.0 1 1 1 NULL 10 Te_output.cub Ni.beta Ni 
```

This fix produces two types of Lammps-internal results in addition to the normal MD:

* vector with the energy and temperature of the electronic system
  * `f_ID[1]` -> Net energy transfer between electronic and ionic system
  * `f_ID[2]` -> Average electronic temperature
 
* per atom values:
  * `f_ID[i][1]` -> site density
  * `f_ID[i][2]` -> coupling parameter

To access them in the output file add this to the LAMMPS input script:
```
fix out all print 1000 "$(step) $(time) $(temp) $(f_ephttm[1]) $(f_ephtmm[2])" file out.data screen no
dump out all custom 10 strucs_out.dump.gz type x y z f_ephttm[1] f_ephttm[2]
```

### Beta(rho) input file

This file provides the electronic densities and beta(rho) functions for individual species (see https://dx.doi.org/10.1103/PhysRevLett.120.185501).
The format is described in `Doc/Beta/input.beta`. 
The file is similar to EAM setfl format. 
The beta(rho) function has units of [eV ps/Ang^2]. 
An example is provided in `Examples/Beta/Ni_model_4.beta`.

### heat equation FDM grid input file

This file is used to initialise FDM grid for the electronic system. 
The format is described in `Doc/FDM/T_input.fdm`. 
This allows fine control of the properties at various grid points. 
Also, the grid can be larger than ionic system if energy reservoir far away is needed. 
Grid points can act as energy sources, sinks or walls (individual grid points updates can be deactivated).
An example of input file is provided in `Examples/FDM/T.in`.
Units are [Kelvin], [eV/Ang^3/ps], [unitless], [in eV/K/Ang^3] [eV/K/Ang/ps] for T, local source term, rho_e, Ce, kappa_e respectively.

## Notes and limitations

* The exact physical interpretation of beta(rho) changes with the precise model. 
* `C_e`, `rho_e`, and `kappa_e` are constants at the moment (not temperature dependent).
* If `T_infile` is not `NULL` then `C_e`, `rho_e`, `kappa_e`, `T_e`, `NX`, `NY`, `NZ` are ignored and are read from the filename supplied. 
If `NULL` is provided as the filename then the FDM grid is initialised with the parameters provided in the command.
* The implementation of the model is applicable to alloys, but this has not been tested thoroughly yet.

# Tutorial

Examples can be found in `Examples/` directory. 
To run them type `lmp_serial -i run.lmp` in the appropriate example directory and assuming executable is in `PATH`. 
Some of the examples may take long on older machines, so tweak the input file (`run.lmp`) accordingly. Every example contains a `README` file that describes what the runscript does.

## Example 1

`Examples/Example_1/`: 
In this example a crystal structure is created and the model is applied with both the friction and random force terms. 
The electrons are *kept* at constant temperature (300K). 
This example illustrates the thermalisation process from 0K to the target temperature through electron-ion interaction only.

```bash
$ cd Examples/Example_1
$ mypath/lmp_serial -i run.lmp
```

The run will write to the file `out.data`, column 2 has the time (in ps), columns 3 and 5 have the ionic and the electronic temperature respectively.
```
$ gnuplot
> plot "out.data" u 2:3 w lp lw 2 t "Tion", "out.data" u 2:5 w lp lw 2 t "Te"
```
You will see that the ionic temperature increases and approaches the (fixed) electronic temperature (300K).

![Alt text](Examples/Example_1/Tout.png?raw=true "Temperature Example 1")

## Example 2

`Examples/Example_2/`: 
This example illustrates the use cooling of the ionic systems due to electrons only. 
This means that only the friction term acts on atoms and removes energy. 
This is equivalent to having electrons at 0K.

![Alt text](Examples/Example_2/Tout.png?raw=true "Temperature Example 2")

## Example 3
`Examples/Example_3/`: 
In this example the full model with electronic heat equation FDM grid is used. 
The crystal is created at equilibrium positions (0K) and it is heated by electrons. 
During the simulation the electronic system will cool and the ionic system heat. 
At equilibrium both systems end up at the same temperature on average. 
Also, this example illustrates the automatic initialisation of the FDM grid with constant parameters. 
The electronic temperature at various grid points is written to files (one per step) (`T_out_XXXXXX`). 
Final state of the grid is stored and can be reused in later simulations (`T.restart`).

Although you can run in serial mode (like above), 
you can try to run this example with in parallel, for example in 4 processes.
Make sure your MPI environment works, for example you may need `$ module load mpi/mpich-x86_64`
```
$ mpirun -np 4 mypath/lmp_mpi -i run.lmp
$ gnuplot
> plot "out.data" u 2:3 w lp lw 2 t "Tion", "out.data" u 2:5 w lp lw 2 t "Te"
```

![Alt text](Examples/Example_3/Tout.png?raw=true "Temperature Example 2")

## Example 4
`Examples/Example_4/`: 
This example reads the FDM grid parameters from a file (`T.in`). 
In this file a source term is added at line in the grid representing the energy dumped by swift ion. 
During the simulation the ionic system will heat while electron temperature will diffuse and due to gradient in the electronic system forces acting on atoms at different grid points will 'feel' different random forces in magnitude.

After a few MD-TTM steps the electronic temperature field will look like this:

![Alt text](Examples/Example_4/Tfieldout.png?raw=true "Temperature Example 4")

## Example 5
`Examples/Example_5/`:

This the same as Example 1, but with GPU-enabled acceleration.

```bash
$ cd Examples/Example_5
$ mpirun -n 2 mypath/lmp_mpi_gpu -i run.lmp
```

On a specialized cluster like LLNL/lassen you need a submission script

```bash
#!/bin/bash
## job.bsub
#BSUB -nnodes 2    # nodes
#BSUB -W 120       # walltime in min
#BSUB -G Bank      # your bank
#BSUB -J Jobname   # name of job
#BSUB -q pbatch    # queue to use (or pdebug)

# Run info and srun job launch
jsrun -r4 -a1 -c1 -g1 -E OMP_NUM_THREADS=1 lmp_mpi_gpu -i run.lmp
```
```
$ bsub < job.bsub
```

To run interactively:

```bash
$ bsub -Is jsrun -r4 -a1 -c1 -g1 ~/mywork/lammps/src/lmp_mpi_gpu -i run.lmp
```
(the arguments are: `-r`:num_resources (cpus), `-a`:mpi_tasks_per_resource, `-c`:cpus_per_resources, `-g`:gpus_per_resource)

# Benchmark, CPUs vs GPU


![Alt text](Doc/Benchmark/CPU_Timing_Quartz.png?raw=true "Code benchmark on Quartz CPU")
![Alt text](Doc/Benchmark/CPU_Timing_Lassen.png?raw=true "Code benchmark on Lassen CPU")
![Alt text](Doc/Benchmark/GPU_Timing_Lassen.png?raw=true "Code benchmark on Lassen GPU")

# Release

## History

- 2018/05/10 Initial Release
- 2019/09/15 GPU port Release

## TODO

- Implement CUBE format output for temperature field

## License and Copying

USER-EPH is licensed under the terms of the [GPL v3 License](/COPYING).

USER-EPH is not part of the LAMMPS code https://github.com/lammps/lammps

If you have any questions contact Artur Tamm <tamm3@llnl.gov> or Alfredo Correa <correaa@llnl.gov>

``LLNL-CODE-750832``
