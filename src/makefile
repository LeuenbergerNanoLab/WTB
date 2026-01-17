default: wtb

export DIR=${PWD}
export PETSC_DIR=${DIR}/../src/petsc-3.22.2
export SLEPC_DIR=${DIR}/../src/petsc-3.22.2/arch-linux-c-opt/externalpackages/slepc-3.22.2
export PETSC_ARCH=arch-linux-c-opt

OBJS = mpi_p.o \
       module_slepc_solver.o \
       module_input.o   \
       module_hamiltonian.o \
       module_bse.o \
       module_pce.o \
       module_coulomb.o \
       module_math.o \
       module_berry.o \
       module_optics.o \
       module_tools.o \
       module_kmesh.o \
       main.o

wtb: $(OBJS)
	-${CLINKER} -o $(DIR)/../bin/wtb.x $(OBJS) ${SLEPC_SYS_LIB}
	${RM} main.o

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common



mpi_p.f90:
module_math.o:
module_optics.o:      module_math.o
module_hamiltonian.o: module_math.o
module_pce.o:         module_math.o
module_kmesh.o:       module_math.o
module_coulomb.o:     module_math.o
module_input.o:       module_hamiltonian.o
module_tools.o:       module_hamiltonian.o module_kmesh.o
module_berry.o:       module_hamiltonian.o module_kmesh.o module_optics.o
module_bse.o:         module_hamiltonian.o module_kmesh.o module_optics.o module_coulomb.o   
main.o:               module_hamiltonian.o \
                      module_input.o \
                      module_bse.o \
                      module_pce.o \
                      module_tools.o \
                      module_berry.o
