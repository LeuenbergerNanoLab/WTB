

## How to install the code:

Load the appropriate compiler and MPI package, for example:

> module load openmpi/4.1.5:intel-2020

Download the code, untar and go to *src* folder. 

Before installation of WTB code, it is necessary to install two scientific libraries, [PETSc](https://petsc.org/release/) and [SLEPc](https://slepc.upv.es).

Download PETSc library archive:

> wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.22.2.tar.gz

untar the archive

> tar -xvf petsc-3.22.2.tar.gz

> cd petsc-3.22.2

Configure the package with complex scalars

> ./configure --with-scalar-type=complex --with-debugging=no

Make installation (the actual installation command will be provided by *configure*):

> make PETSC_DIR=/lustre/fs1/home/dm606074/WTB/src/petsc-3.22.2 PETSC_ARCH=arch-linux-c-opt all

Return to *src* directory

> cd ..

Download SLEPc library archive:

> wget https://slepc.upv.es/download/distrib/slepc-3.22.2.tar.gz

and untar

> tar -xvf slepc-3.22.2.tar.gz

> cd slepc-3.22.2/

Configure few variables for SLEPc installation:

> export SLEPC_DIR=$PWD

> export PETSC_DIR={full path to petsc folder}

> export PETSC_ARCH=arch-linux-c-opt

Configure the package:

> ./configure

Make installation of SLEPc library (the actual installation command will be provided by *configure*)

> make SLEPC_DIR=/home/3811/slepc-3.22.2 PETSC_DIR=/home/3811/petsc-3.22.2 PETSC_ARCH=arch-linux-c-opt

Go to *src* folder and make installation of WTB package

> cd ..

> mkdir bin

> make

## Installation on STOKES cluster of UCF, 1/16/2026

Load the modules:

> module load oneapi/oneapi-2023.1.0/mpi/mpi-2021.9.0 

> module load openblas/openblas-0.3.25-oneapi-2023.1.0

Deactivate Conda:

> conda deactivate

Download PETSc library archive in *src* folder:

> wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.22.2.tar.gz

> tar -xvf petsc-3.22.2.tar.gz

> cd petsc-3.22.2

Configure the package with SLEPc library:
```
./configure --with-scalar-type=complex \
--with-debugging=no \
--download-slepc=https://slepc.upv.es/download/distrib/slepc-3.22.2.tar.gz \
--with-clean=1
```

Then you need to request interactive job because installation of big packages on STOKES is prohibited:

> salloc --time=1:00:00 --cpus-per-task=4 --mem=8G

> conda deactivate

Make installation:

> make PETSC_DIR=/lustre/fs1/home/dm606074/WTB/src/petsc-3.22.2 PETSC_ARCH=arch-linux-c-opt all

This will install the PETSc and SLEPs libraries at petsc-3.22.2 and petsc-3.22.2/arch-linux-c-opt/externalpackages/slepc-3.22.2 folders. 

Then you need to modify *makefile* accordingly:
```
export DIR=${PWD}
export PETSC_DIR=${DIR}/petsc-3.22.2
export SLEPC_DIR=${DIR}/petsc-3.22.2/arch-linux-c-opt/externalpackages/slepc-3.22.2
export PETSC_ARCH=arch-linux-c-opt
```
Go to WTB directory and make *bin* folder:

> cd WTB

> mkdir bin

> cd src

Make installation of WTB package:

> make
