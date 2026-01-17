


module module_hamiltonian
    use MPI_P                                                             ! MPI module 
    use module_math
    implicit none
    integer                               :: w90basis,ntype,nocp
    double precision                      :: efermi,scs
    integer                               :: nvec
    double precision,dimension(3,3)       :: rlat
    integer,allocatable,dimension(:,:)    :: rvec
    real(8),allocatable,dimension(:,:,:)  :: hopmatrices,ihopmatrices
    integer,allocatable,dimension(:)      :: ffactor
    character(len=4)                      :: systype
    real(8),allocatable,dimension(:,:,:)  :: hopmatricesu,ihopmatricesu   !variables for spin polarized hamiltonian
    real(8),allocatable,dimension(:,:,:)  :: hopmatricesd,ihopmatricesd
    integer,allocatable,dimension(:)      :: ffactoru,ffactord
    integer,allocatable,dimension(:,:)    :: rvecu,rvecd
    integer                               :: w90basisu,w90basisd
    integer                               :: nvecu,nvecd
    complex(8), allocatable               :: H_k(:,:,:)                                ! full Hamiltonian H(k) from Wannier90
    complex(8), allocatable               :: Psi_ik(:,:,:)                             ! |ik>     eigenvectors of H(k)
    real(8), allocatable                  :: E_ik(:,:)                                 ! E_ik     eigenvalues  of H(k)
    integer                               :: ik
    logical                               :: lkbe
    complex(8), allocatable               :: htb(:,:)
    complex(8), allocatable               :: fermimt(:,:)
    complex(8), allocatable               :: scor(:,:)
    complex(8), allocatable               :: hexc(:,:)
    complex(8), allocatable               :: kvecs(:)
Contains


subroutine hamiltonian_input_read(unidade,inputfile)                 ! read Hamiltonian in parallel regime
    integer                     :: i,j,k
    integer                     :: unidade,erro
    integer                     :: flagaux
    character(len=70)           :: inputfile,charflag

!    if(Process==0) print *,'  hamiltonian input read'
    call MSG('MODULE_HAMILTONIAN: open input file')
    if(Process==0) then
!     print *,'hamiltonian_input_read: lkbe=',lkbe
     OPEN(UNIT=unidade, FILE= inputfile,STATUS='old', IOSTAT=erro)
        if (erro/=0) stop "Error opening wannier90 hr input file"
    
     read(unidade,*) systype
     read(unidade,*) scs
     read(unidade,*) efermi
     read(unidade,*) rlat(1,1),rlat(1,2),rlat(1,3)
     read(unidade,*) rlat(2,1),rlat(2,2),rlat(2,3)
     read(unidade,*) rlat(3,1),rlat(3,2),rlat(3,3)
    endif 
    call P_sendA(systype)                                 ! send array to all other processes
    !if(Process==0) print *,' systype=',systype
    call P_sendR(scs)                                     ! send variable to all other processes
    call P_sendR(efermi)                                  ! send variable to all other processes
    call P_sendR2(rlat,3,3)                               ! send array to all other processes
    
    select case (systype)
    
    case("SP")
    
    !lendo parte spin up
    
    if(Process==0) read(unidade,*) charflag
    call P_sendA(charflag)                               ! send variable to all other processes
    if(Process==0) print *,' charflag=',charflag
    if(Process==0) read(unidade,*) w90basisu
    call P_sendI(w90basisu)                              ! send variable to all other processes
    if(Process==0) read(unidade,*) nvecu
    call P_sendI(nvecu)                                  ! send variable to all other processes

    allocate(ffactoru(nvecu))

    if(Process==0) read(unidade,*)(ffactoru(i), i=1,nvecu)
    call P_sendI1(ffactoru,nvecu)                        ! send array to all other processes

    allocate(rvecu(nvecu,3))
    allocate(hopmatricesu(nvecu,w90basisu,w90basisu))
    allocate(ihopmatricesu(nvecu,w90basisu,w90basisu))

    if(Process==0) then
     do i=1,nvecu
        do j=1,w90basisu
            do k=1,w90basisu
            read(unidade,*) rvecu(i,1),rvecu(i,2),rvecu(i,3),flagaux,flagaux,hopmatricesu(i,k,j),ihopmatricesu(i,k,j)
            end do
        end do
     end do
    endif 
    call P_sendI2(rvecu,nvecu,3)                           ! send array to all other processes
    call P_sendI(flagaux)                                  ! send integer variable to all other processes
    call P_sendR3(hopmatricesu,nvecu,w90basisu,w90basisu)    ! send array to all other processes
    call P_sendR3(ihopmatricesu,nvecu,w90basisu,w90basisu)   ! send array to all other processes



    !lendo parte spin down  
    
    if(Process==0) read(unidade,*) charflag
    call P_sendA(charflag)                               ! send variable to all other processes

    if(Process==0) read(unidade,*) w90basisd
    call P_sendI(w90basisd)                                  ! send integer variable to all other processes
    if(Process==0) read(unidade,*) nvecd
    call P_sendI(nvecd)                                  ! send integer variable to all other processes

    allocate(ffactord(nvecd))

    if(Process==0) read(unidade,*)(ffactord(i), i=1, nvecd)
    call P_sendI1(ffactord,nvecd)                                  ! send integer variable to all other processes

    allocate(rvecd(nvecd,3))
    allocate(hopmatricesd(nvecd,w90basisd,w90basisd))
    allocate(ihopmatricesd(nvecd,w90basisd,w90basisd))

    if(Process==0) then
     do i=1,nvecd
        do j=1,w90basisd
            do k=1,w90basisd
            read(unidade,*) rvecd(i,1),rvecd(i,2),rvecd(i,3),flagaux,flagaux,hopmatricesd(i,k,j),ihopmatricesd(i,k,j)
            end do
        end do
     enddo
    endif 
    call P_sendI2(rvecd,nvecd,3)                           ! send array to all other processes
    call P_sendI(flagaux)                                  ! send integer variable to all other processes
    call P_sendR3(hopmatricesd,nvecd,w90basisd,w90basisd)    ! send array to all other processes
    call P_sendR3(ihopmatricesd,nvecd,w90basisd,w90basisd)   ! send array to all other processes


    !juntando as partes up e down
    
    w90basis = w90basisd+w90basisu
    nvec = nvecd+nvecu
    
    allocate(ffactor(nvec)) 
    allocate(rvec(nvec,3))  
    allocate(hopmatrices(nvec,w90basis,w90basis))
    allocate(ihopmatrices(nvec,w90basis,w90basis))

    ffactor(1:nvecu) = ffactoru
    ffactor(nvecu+1:nvec) = ffactord
    
    rvec(1:nvecu,:) = rvecu
    rvec(nvecu+1:nvec,:) = rvecd
    
    hopmatrices = 0.0
    ihopmatrices = 0.0 
    
    hopmatrices(1:nvecu,1:w90basisu,1:w90basisu) = hopmatricesu
    ihopmatrices(1:nvecu,1:w90basisu,1:w90basisu) = ihopmatricesu
    
    hopmatrices(nvecu+1:nvec,w90basisu+1:w90basis,w90basisu+1:w90basis) = hopmatricesd
    ihopmatrices(nvecu+1:nvec,w90basisu+1:w90basis,w90basisu+1:w90basis) = ihopmatricesd
    
    deallocate(ffactord)
    deallocate(rvecd,hopmatricesd,ihopmatricesd)
    deallocate(ffactoru)
    deallocate(rvecu,hopmatricesu,ihopmatricesu)
    
    case default

    if(Process==0) read(unidade,*) charflag
    call P_sendA(charflag)                                  ! send character variable to all other processes
    if(Process==0) read(unidade,*) w90basis
    call P_sendI(w90basis)                                  ! send character variable to all other processes
    if(Process==0) read(unidade,*) nvec
    call P_sendI(nvec)                                  ! send character variable to all other processes
    allocate(ffactor(nvec))
    if(Process==0) read(unidade,*)(ffactor(i), i=1, nvec)
    call P_sendI1(ffactor,nvec)                                  ! send integer array to all other processes

    allocate(rvec(nvec,3))
    allocate(hopmatrices(nvec,w90basis,w90basis))
    allocate(ihopmatrices(nvec,w90basis,w90basis))
    if(Process==0) then
     do i=1,nvec
        do j=1,w90basis
            do k=1,w90basis
            read(unidade,*) rvec(i,1),rvec(i,2),rvec(i,3),flagaux,flagaux,hopmatrices(i,k,j),ihopmatrices(i,k,j)
            end do
        end do
     end do 
    endif 
    call P_sendI2(rvec,nvec,3)                           ! send array to all other processes
    call P_sendI(flagaux)                                  ! send integer variable to all other processes
    call P_sendR3(hopmatrices,nvec,w90basis,w90basis)    ! send array to all other processes
    call P_sendR3(ihopmatrices,nvec,w90basis,w90basis)   ! send array to all other processes
   end select
!   call P_stop                                          ! stop parallel calculation
end subroutine hamiltonian_input_read





subroutine eigsys(nthreads,scs,exc,nocp,ffactor,kx,ky,kz,w90basis,nvec,rlat,rvec,hopmatrices,&
           ihopmatrices,efermi,energias,autovetores)
    integer :: i,j,k,nthreads
    integer :: w90basis,nvec
    double precision :: kx,ky,kz
    double precision :: scs,exc
    integer :: nocp
    double precision,dimension(3,3) :: rlat
    integer,dimension(nvec,3) :: rvec
    integer,dimension(w90basis) :: ffactor
    double precision,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices
    double precision :: efermi
    double precision,dimension(w90basis) :: energias
    double complex,dimension(w90basis,w90basis) :: autovetores
    double complex,parameter :: imag=cmplx(0.0,1.0)
    INTEGER   ::       ifail
    double precision,parameter :: ABSTOL=1.0e-6
    INTEGER          LWMAX
    INTEGER ::         LIWORK, LRWORK
    INTEGER  ::        INFO, LWORK
    real(8), allocatable               :: RWORK(:)
    complex(8), allocatable            :: WORK(:)
    INTEGER, allocatable               :: IWORK(:)

    allocate(RWORK(1 + 5*w90basis + 2*w90basis**2))
    allocate(WORK(2*w90basis+w90basis**2))
    allocate(IWORK(3 + 5*w90basis))
    allocate(htb(w90basis,w90basis))
    allocate(fermimt(w90basis,w90basis))
    allocate(scor(w90basis,w90basis))
    allocate(hexc(w90basis,w90basis))
    allocate(kvecs(nvec))

    LWORK  = 2*w90basis+w90basis**2
    LIWORK = 3 + 5*w90basis
    LRWORK = 1 + 5*w90basis + 2*w90basis**2
    LWMAX  = 2*w90basis-1

!   call OMP_SET_NUM_THREADS(nthreads)

        do i=1,nvec
        kvecs(i) = cmplx(0.0,kx*(rvec(i,1)*rlat(1,1)+rvec(i,2)*rlat(2,1)+rvec(i,3)*rlat(3,1))+&
               ky*(rvec(i,1)*rlat(1,2)+rvec(i,2)*rlat(2,2)+rvec(i,3)*rlat(3,2))+&
               kz*(rvec(i,1)*rlat(1,3)+rvec(i,2)*rlat(2,3)+rvec(i,3)*rlat(3,3)))
        end do

        htb = 0.0
        do i=1,nvec
            htb=htb+cdexp(kvecs(i))*(hopmatrices(i,:,:)+imag*ihopmatrices(i,:,:))*(1.0/dble(ffactor(i)))
        end do

        call fermilvl(w90basis,efermi,fermimt)

        !call syscor(w90basis,scor)

        !call exch(w90basis,exc,exc2,mag,hexc)

        !call spincor(exc,w90basis,scor)

        htb = htb+fermimt!+scor+hexc
        
        if(lkbe) then                    ! store full Hamiltonian H(k)
        !  print *,'store H in H_k'
         if(ik==0) then
          print *,'eigsys ERROR ik= ',ik
         else 
          do i=1,w90basis
           do j=1,w90basis
            H_k(j,i,ik) = htb(j,i)
           enddo
          enddo 
!    print *
!    print *,'Hamiltonian H_k for k-point ik=',ik
!    do j=1,w90basis
!     print 111,(H_k(j,i,ik),i=1,w90basis)
!    enddo 
!111 format(22(E12.4,',',4x))
!    stop   
         endif
        endif

            !call LA_HEEVR( htb, energias, JOBZ='V', UPLO='U', ABSTOL=ABSTOL, INFO=INFO)

            !WRITE(*,*)'ZHEEV Example Program Results'
            LWORK = -1
            CALL ZHEEV( 'Vectors', 'Lower', w90basis, htb, w90basis, energias, WORK, LWORK, RWORK, INFO )
            LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
            CALL ZHEEV( 'Vectors', 'Lower', w90basis, htb, w90basis, energias, WORK, LWORK, RWORK,INFO )
        !CALL ZHEEVD('V','U',w90basis,htb,w90basis,energias,WORK, LWORK,&
        !       RWORK,LRWORK,IWORK,INFO)
            IF( INFO.GT. 0 ) THEN
            WRITE(*,*)'The algorithm failed to compute eigenvalues.'
            STOP
            END IF

        do i=1,w90basis
            do j=1,w90basis
                autovetores(i,j) = htb(j,i)
            end do
        end do

        do i=1,w90basis
            if (energias(i) .gt. 0.0 ) then
                nocp = i-1
                EXIT
            end if
        end do

        do i=nocp+1,w90basis
            energias(i)=energias(i)+scs
        end do

    deallocate(RWORK)
    deallocate(WORK)
    deallocate(IWORK)
    deallocate(htb)
    deallocate(fermimt)
    deallocate(scor)
    deallocate(hexc)
    deallocate(kvecs)
end subroutine eigsys



subroutine eigsys2(nthreads,scs,exc,nocp,ffactor,kx,ky,kz,w90basis,nvec,rlat,rvec,hopmatrices,&
           ihopmatrices,efermi,energias,autovetores)
    integer :: i,j,k,nthreads
    integer :: w90basis,nvec
    double precision :: kx,ky,kz
    !double precision :: exc,exc2
    double precision :: scs,exc
    integer :: nocp
    double precision,dimension(3,3) :: rlat
    integer,dimension(nvec,3) :: rvec
    integer,dimension(w90basis) :: ffactor
    double precision,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices
    double precision :: efermi
    double precision,dimension(w90basis) :: energias
    double complex,dimension(w90basis,w90basis) :: autovetores
    double complex,parameter :: imag=cmplx(0.0,1.0)
    INTEGER   ::       ifail
    double precision,parameter :: ABSTOL=1.0e-6
    INTEGER          LWMAX
    INTEGER ::         LIWORK, LRWORK
    INTEGER  ::        INFO, LWORK
!        DOUBLE PRECISION RWORK( 1 + 5*w90basis + 2*w90basis**2 )
!        double complex, dimension (2*w90basis+w90basis**2) :: WORK
!   INTEGER,dimension(3 + 5*w90basis) :: IWORK
    real(8), allocatable               :: RWORK(:)
    complex(8), allocatable            :: WORK(:)
    INTEGER, allocatable               :: IWORK(:)

    allocate(RWORK(1 + 5*w90basis + 2*w90basis**2))
    allocate(WORK(2*w90basis+w90basis**2))
    allocate(IWORK(3 + 5*w90basis))

        LWORK = 2*w90basis+w90basis**2
        LIWORK = 3 + 5*w90basis
        LRWORK = 1 + 5*w90basis + 2*w90basis**2
    LWMAX= 2*w90basis-1

    !call OMP_SET_NUM_THREADS(nthreads)


        do i=1,nvec

        kvecs(i) = cmplx(0.0,kx*(rvec(i,1)*rlat(1,1)+rvec(i,2)*rlat(2,1)+rvec(i,3)*rlat(3,1))+&
               ky*(rvec(i,1)*rlat(1,2)+rvec(i,2)*rlat(2,2)+rvec(i,3)*rlat(3,2))+&
               kz*(rvec(i,1)*rlat(1,3)+rvec(i,2)*rlat(2,3)+rvec(i,3)*rlat(3,3)))

        end do


        htb = 0.0

        do i=1,nvec

            htb=htb+cdexp(kvecs(i))*(hopmatrices(i,:,:)+imag*ihopmatrices(i,:,:))*(1.0/dble(ffactor(i)))

        end do

        call fermilvl(w90basis,efermi,fermimt)

        !call syscor(w90basis,scor)

        !call exch(w90basis,exc,exc2,mag,hexc)

        !call spincor(exc,w90basis,scor)

        htb = htb+fermimt!+scor+hexc

            !call LA_HEEVR( htb, energias, JOBZ='V', UPLO='U', ABSTOL=ABSTOL, INFO=INFO)

            !WRITE(*,*)'ZHEEV Example Program Results'
            LWORK = -1
            CALL ZHEEV( 'Vectors', 'Lower', w90basis, htb, w90basis, energias, WORK, LWORK, RWORK, INFO )
            LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
            CALL ZHEEV( 'Vectors', 'Lower', w90basis, htb, w90basis, energias, WORK, LWORK, RWORK,INFO )
        !CALL ZHEEVD('V','U',w90basis,htb,w90basis,energias,WORK, LWORK,&
        !       RWORK,LRWORK,IWORK,INFO)
            IF( INFO.GT. 0 ) THEN
            WRITE(*,*)'The algorithm failed to compute eigenvalues.', INFO
            !STOP
            END IF

        do i=1,w90basis
            do j=1,w90basis
                autovetores(i,j) = htb(j,i)
            end do
        end do

        do i=1,w90basis

            if (energias(i) .gt. 0.0 ) then
                nocp = i-1
                EXIT
            end if
        end do

        do i=nocp+1,w90basis
            energias(i)=energias(i)+scs
        end do

    deallocate(RWORK)
    deallocate(WORK)
    deallocate(IWORK)
    deallocate(htb)
    deallocate(fermimt)
    deallocate(scor)
    deallocate(hexc)
    deallocate(kvecs)
end subroutine eigsys2






subroutine spinvl(w90basis,vec,dft,systype,spz)
    integer :: w90basis,w2
    double complex,dimension(w90basis) :: vec
    double complex,dimension(w90basis) :: vcconj
    double complex,dimension(w90basis) :: zvvaux
    double complex,dimension(w90basis,w90basis) :: spmz
    double complex :: spz
    integer :: i,j
    character(len=1) :: dft
    character(len=4) :: systype 

    spmz=0.0
    w2=w90basis/2
    
    if (systype .eq. "NP") then 
     do i=1,w90basis
        spmz(i,i) = cmplx(1.,0.)
     end do
    else if (systype .eq. "SP") then
     do i=1,w2
        spmz(i,i) = cmplx(1.,0.)
        spmz(w2+i,w2+i) = cmplx(-1.,0.)
     end do
    else if ((systype .eq. "SOC") .and. (dft .eq. "V") ) then
     do i=1,w2
        spmz(i,i) = cmplx(1.,0.)
        spmz(w2+i,w2+i) = cmplx(-1.,0.)
     end do
    else
     do i=1,w2
        spmz((2*i)-1,(2*i)-1) = cmplx(1.,0.)
        spmz(2*i,2*i) = cmplx(-1.,0.)
     end do
    end if

    call vecconjg(vec,w90basis,vcconj)
    call matvec(spmz,vec,w90basis,zvvaux)
    call prodintsq(vcconj,zvvaux,w90basis,spz)
end subroutine spinvl



subroutine fermilvl(w90basis,efermi,fermimt)
    integer :: w90basis,i
    double precision :: efermi
    double complex,dimension(w90basis,w90basis) :: fermimt
    fermimt = (0.0,0.0)
    do i=1,w90basis
        fermimt(i,i) = cmplx(-efermi,0.0)
    end do
end subroutine fermilvl







end module module_hamiltonian




