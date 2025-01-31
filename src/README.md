

## How to install the code:

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

Make installation (the actual installation command will be provided by *configure*):

> make PETSC_DIR=/home/3811/petsc-3.22.2 PETSC_ARCH=arch-linux-c-opt all

Return to *src* directory

> cd ..

Download SLEPc library archive:

> wget https://slepc.upv.es/download/distrib/slepc-3.22.2.tar.gz

and untar

> tar -xvf slepc-3.22.2.tar.gz

> cd slepc-3.22.2/

> export SLEPC_DIR=$PWD

> export PETSC_DIR={full path to petsc folder}

> export PETSC_ARCH=arch-linux-c-opt

> ./configure

Make installation of SLEPc library (the actual installation command will be provided by *configure*)

> make SLEPC_DIR=/home/3811/slepc-3.22.2 PETSC_DIR=/home/3811/petsc-3.22.2 PETSC_ARCH=arch-linux-c-opt

Go to *src* folder and make installation of WTB package

> mkdir bin

> make




  
