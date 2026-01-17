Module module_slepc_solver
#include <slepc/finclude/slepceps.h>
      use slepceps
      implicit none
      PetscLogDouble           :: ts0,tsf
      PetscMPIInt              :: rank
      PetscErrorCode           :: ierr
      PetscReal                :: tols
contains


     subroutine slepc_solver(HX,EX,HE,NX,NE)
      EPS                      :: eps
      Mat                      :: HH
      Vec                      :: xx,vtr
      VecScatter               :: ctr
      PetscInt                 :: ns, is, js, Istart, Iend
      PetscInt                 :: nev,ncv,mpd,max_it
      PetscInt                 :: row(1)
      PetscInt, allocatable    :: col(:)
      PetscScalar              :: kr 
      PetscScalar, allocatable :: H1(:) 
      PetscScalar, pointer     :: xxs(:)
      integer                  :: NX,NE
      complex(8)               :: HX(NX,NX)
      complex(8)               :: HE(NX,NE)
      real(8)                  :: EX(NE)

      ns   = NX                                                           ! size of the matrix
      nev  = NE                                                           ! number of the solutions
      allocate(col(ns))
      allocate(H1(ns))
      if(nev<1000) then
       mpd  = nev                                                         ! working projections (for small nev mpd=nev, for large nev mpd<<nev) 
      else
       mpd = nev/5
      endif 
      ncv  = nev + mpd                                                    ! working subspace
      tols = 1.d-8                                                        ! accuracy of the solution
      max_it = 200                                                        ! max number of iterations
      PetscCallA(SlepcInitialize(PETSC_NULL_CHARACTER,ierr))              ! initialize SLEPc
      PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))            ! Initialize MPI processes
      PetscCallA(PetscTime(ts0,ierr))
      PetscCallA(MatCreate(PETSC_COMM_WORLD,HH,ierr))                     ! create matrix HH
      PetscCallA(MatSetSizes(HH,PETSC_DECIDE,PETSC_DECIDE,ns,ns,ierr))    ! set size of the matrix HH (ns x ns)
      PetscCallA(MatSetFromOptions(HH,ierr))                              ! type of matrix HH (default MATAIJ)
      PetscCallA(MatGetOwnershipRange(HH,Istart,Iend,ierr))               ! start to fill matrix HH in parallel
      PetscCallA(MatCreateVecs(HH,PETSC_NULL_VEC,xx,ierr))                ! create vector xx for eigenvectors

      do is=Istart,Iend-1                                                 ! fill Hamiltonian in parallel regime
       row(1) = is
       do js=1,ns
        col(js) = js-1
       enddo 
       do js=1,ns
        H1(js) = HX(is+1,js)
       enddo
       PetscCallA(MatSetValues(HH,1,row,ns,col,H1,INSERT_VALUES,ierr))     ! insert one row in HH (it may be a part of the row)
      enddo
      PetscCallA(MatAssemblyBegin(HH,MAT_FINAL_ASSEMBLY,ierr))             ! MPI assembling of matrix HH
      PetscCallA(MatAssemblyEnd(HH,MAT_FINAL_ASSEMBLY,ierr))               !
      call SLEPc_print('SLEPc matrix assembling')

      PetscCallA(EPSCreate(PETSC_COMM_WORLD,eps,ierr))                    ! create object EPS
      PetscCallA(EPSSetOperators(eps,HH,PETSC_NULL_MAT,ierr))             ! set matrix HH to EPS problem AX=lX
      PetscCallA(EPSSetProblemType(eps,EPS_HEP,ierr))                     ! Hermitian Hamiltonian
      PetscCallA(EPSSetWhichEigenpairs(eps,EPS_SMALLEST_MAGNITUDE,ierr))  ! calculate smallest eigenvalues
      PetscCallA(EPSSetType(eps,EPSKRYLOVSCHUR,ierr))                     ! Krylov-Schur method (default)
      PetscCallA(EPSSetFromOptions(eps,ierr))                             ! set all parameters for EPS
      PetscCallA(EPSSetDimensions(eps,nev,ncv,mpd,ierr))                  ! set number of eigenvalues to find (nev)
      PetscCallA(EPSSetTolerances(eps,tols,max_it,ierr))                  ! set max number of iterations and tolerance 
      PetscCallA(EPSMonitorSet(eps,MyEPSMonitor,0,PETSC_NULL_FUNCTION,ierr))    ! print iterations

      PetscCallA(EPSSolve(eps,ierr))                                      ! solve AX = lX 
      call SLEPc_print('SLEPc solve AX=lX')

      do is=0,nev-1                                                       ! in SLEPc indexes start from zero
       PetscCallA(EPSGetEigenpair(eps,is,kr,PETSC_NULL_SCALAR,xx,PETSC_NULL_VEC,ierr)) ! read solution is
       PetscCallA(VecScatterCreateToAll(xx,ctr,vtr,ierr))                 ! request parts of xx from all processes 
       PetscCallA(VecScatterBegin(ctr,xx,vtr,INSERT_VALUES,SCATTER_FORWARD,ierr)) ! combine parts of xx from all processes
       PetscCallA(VecScatterEnd(ctr,xx,vtr,INSERT_VALUES,SCATTER_FORWARD,ierr))   ! into vtr
       PetscCallA(VecGetArrayReadF90(vtr,xxs,ierr))                       ! get access to Vec vtr
        EX(is+1) = PetscRealPart(kr)                                      ! store eigenvalues
        HE(:,is+1) = xxs                                                  ! store eigenvectors
       PetscCallA(VecRestoreArrayReadF90(xx,xxs,ierr))                    ! restore Vec xx
      enddo
      call SLEPc_print('SLEPc extract solution')

      deallocate(col)
      deallocate(H1)
      PetscCallA(EPSDestroy(eps,ierr))                                    ! destroy EPS object
      PetscCallA(VecDestroy(xx,ierr))                                     ! destroy vector
      PetscCallA(VecScatterDestroy(ctr,ierr))
      PetscCallA(VecDestroy(vtr,ierr))
      PetscCallA(MatDestroy(HH,ierr))                                     ! destroy matrix
      PetscCallA(SlepcFinalize(ierr))                                     ! finalize SLEPc MPI calculation
     end subroutine slepc_solver



     subroutine MyEPSMonitor(eps,its,nconv,eigr,eigi,errest,nest,dummy,ierr)    ! user-defined routine for monitoring
      EPS            eps                          ! eigensolver context
      PetscInt       its                          ! iteration number
      PetscInt       nconv                        ! number of converged eigenpairs
      PetscScalar    eigr(*)                      ! real part of the eigenvalues
      PetscScalar    eigi(*)                      ! imaginary part of the eigenvalues
      PetscReal      errest(*)                    ! relative error estimates for each eigenpair
      PetscInt       nest                         ! number of error estimates
      PetscInt       dummy                        ! optional user-defined monitor context
      PetscErrorCode ierr
      ierr = 0
      if (its .gt. 0 .and. rank .eq. 0) then
       print 140,its,nconv+1,PetscRealPart(eigr(nconv+1)),errest(nconv+1)
      endif
 140  format('SLEPc ',i5,' EPS pair=',i4,' value (error) ',E15.4,' (',g10.3,')')
     end subroutine MyEPSMonitor



   subroutine slepc_check_results(HX,EX,HE,NX,NE)          ! check NE solutions AX=lX
    integer                 :: NX                          ! size of Hamiltonian matrix
    integer                 :: NE                          ! number of solutions
    complex(8)              :: HX(NX,NX)                   ! Hamiltonian matrix
    complex(8)              :: HE(NX,NE)                   ! NE eigenvectors
    real(8)                 :: EX(NE)                      ! NE eigenvalues
    integer                 :: q,l,i
    complex(8)              :: s
    logical                 :: l_wrong,l_one
    l_one = .false.
    do q = 1,NE
     if(rank==0) print *,'SLEPc check solution q=',q
     l_wrong = .false.
     do i=1,NX
      s = (0.d0,0.d0)
      do l=1,NX
       s = s + HX(i,l)*HE(l,q)                             ! calculate sum A(i,l)*X(l,q)
      enddo
      if(zabs(s-EX(q)*HE(i,q)) > tols*10) then             ! check i-component of AX-lX is zero or not, for i=1,NE
       l_wrong = .true.
      endif
     enddo
     if(l_wrong) then
      if(rank==0) print *,'SLEPc solution q=',q,' is wrong'
      l_one = .true.
     endif 
    enddo 
    if(.not.l_one .and. rank==0) print *,'SLEPc all solutions are correct'
   end subroutine slepc_check_results



   subroutine SLEPc_print(comment)
    character(*)          :: comment
    integer               :: l1     !,l2
    character(len=60)     :: cp
    PetscCallA(PetscTime(tsf,ierr))
    cp = ' '
    l1 = len(trim(adjustl(comment)))
    if(l1>=60) then
     if(rank==0) print 1,trim(adjustl(comment)),(tsf-ts0)
    else
  !   l2 = 60 - l1
     cp(1:l1) = trim(adjustl(comment))
     if(rank==0) print 2,cp,(tsf-ts0)
    endif
    ts0 = tsf
 1  format(A,' (',F10.3,' s)')
 2  format(A60,' (',F10.3,' s)')
   end subroutine SLEPc_print



end Module module_slepc_solver




