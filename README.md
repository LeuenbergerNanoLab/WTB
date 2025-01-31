# WTB

This is the improved version of [WanTiBexos](https://github.com/ac-dias/wantibexos) code,

see A.C. Dias, J.F.R.V. Silveira, F. Qu "WanTiBEXOS: A Wannier based Tight Binding code for electronic band structure, excitonic and optoelectronic properties of solids" Comp. Phys. Comm., 285, 108636 (2022) [DOI](https://doi.org/10.1016/j.cpc.2022.108636)

## Version 1.0
What's new:

- the MKL solver for Bethe-Salpeter equation was replaced by MPI+GPU paralel solver from the [Scalable Library for Eigenvalue Problem Computations (SLEPc)](https://slepc.upv.es)

The calculation time for [Example-3](examples/example-3) is reduced from 680 s using MKL solver to 3 s using SLEPc solver on 64 cores.

## How to install the code:

# SLEPc installation and examples

Load the appropriate compilers, for example:

> module load openmpi/4.1.5:intel-2020

Download the code, untar and go to *src* folder.

Download PETSc library archive:

> wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.22.2.tar.gz

untar the archive

> tar -xvf petsc-3.22.2.tar.gz

> cd petsc-3.22.2

Configure the package with complex scalars

> ./configure --with-scalar-type=complex --with-debugging=no

Make installation:

> make PETSC_DIR=/home/3811/petsc-3.22.2 PETSC_ARCH=arch-linux-c-debug all

Return to *src* directory

> cd ..

Download SLEPc library archive:

> wget https://slepc.upv.es/download/distrib/slepc-3.22.2.tar.gz

and untar

> tar -xvf slepc-3.22.2.tar.gz

> cd slepc-3.22.2/

> export SLEPC_DIR=$PWD

> export PETSC_DIR={full path to petsc folder}

> export PETSC_ARCH=arch-linux-c-debug

> ./configure

> make SLEPC_DIR=/home/3811/slepc-3.22.2 PETSC_DIR=/home/3811/petsc-3.22.2 PETSC_ARCH=arch-linux-c-debug




```
module load openmpi/4.0.2-intel20.4
```

The remaining installation will be the same


## Examples on Fortran

[Example 1](Example) solve large system for 3000 x 3000 and calculate 1000 first lowest eigenvalues.

[Example 2](Example_2) calculate eigenvalues and eigenvectors, store them into arrays, and check the solutions. 

[Example 3](Example_3) use slepc_solver subroutine to solve 8000x8000 Hamiltonian and print convergence process.


  
