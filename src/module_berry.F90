


   module module_berry
    use MPI_P
!   use omp_lib
    use module_hamiltonian  
    use module_kmesh
    use module_optics
    implicit none
   
   
   contains

!INCLUDE "./subroutines/berry_curvature_subs.f90"   
subroutine berryct2(nocp,kx,ky,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
          ihopmatrices,efermi,energias,autovetores,nthread,gammas,bxx,bxy,bxz,byy,byz,bzz)
    integer :: nthread
    double precision :: kx,ky,kz,ezz
    integer :: nocp,w90basis
    integer :: nvec
    integer :: n,i,j
    double precision,parameter :: ktol = 0.001
    integer,dimension(w90basis) :: ffactor
    double precision,dimension(3,3) :: rlat
    integer,dimension(nvec,3) :: rvec
    integer :: nsite
    double complex, dimension(w90basis,w90basis) :: hx,hy,hz
    double complex, dimension(w90basis,w90basis) :: autovetores,autovetores2
    double precision,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices
    double precision, dimension(w90basis) :: energias,energias2
    !double precision,dimension(3) :: mag
    double precision :: efermi
    double complex, dimension(w90basis) :: vflagx,vflagy,vflagz
    double complex, dimension(w90basis) :: vn,vm,vl
    double complex,dimension(w90basis) :: vx,vy,vz
    double complex,dimension(w90basis) :: vxa,vya
    double complex, dimension(w90basis) :: vflagxa,vflagya
    double complex, dimension(w90basis) :: vna,vma

    double complex :: bxx,bxy,bxz,byy,byz,bzz
    double complex :: bfxx,bfxy,bfxz,bfyy,bfyz,bfzz

    double precision:: gammas

    double precision :: berryaux


    
    call hlm(kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
         ihopmatrices,hx,hy,hz)
    

    bxx=cmplx(0.0,0.0)
    bxy=cmplx(0.0,0.0)
    bxz=cmplx(0.0,0.0)
    byy=cmplx(0.0,0.0)
    byz=cmplx(0.0,0.0)
    bzz=cmplx(0.0,0.0)      

    do n=1,nocp,1
    

    do i=nocp+1,w90basis,1

    call vecconjg(autovetores(n,:),w90basis,vn)

    call matvec(hx,autovetores(i,:),w90basis,vflagx)

    call prodintsq(vn,vflagx,w90basis,vx(i))

    !write(*,*) vx(i)


    call vecconjg(autovetores(n,:),w90basis,vm)

    call matvec(hy,autovetores(i,:),w90basis,vflagy)

    call prodintsq(vm,vflagy,w90basis,vy(i))
    

    !write(*,*) vy(i)
    
    call vecconjg(autovetores(n,:),w90basis,vl)

    call matvec(hz,autovetores(i,:),w90basis,vflagz)

    call prodintsq(vl,vflagz,w90basis,vz(i))
    

    if (i == n) then

    bfxx=0.
    bfxy=0.
    bfxz=0.
    bfyy=0.
    bfyz=0.
    bfzz=0.         

    else if (abs(energias(i)-energias(n)) .lt. 0.0001 ) then

    bfxx=0.
    bfxy=0.
    bfxz=0.
    bfyy=0.
    bfyz=0.
    bfzz=0. 

    else
    
    berryaux= (energias(n)-energias(i))**2+(gammas)**2
    
    bfxx=-2.*(((vx(i))*conjg(vx(i)))/(berryaux))
    bfxy=-2.*(((vx(i))*conjg(vy(i)))/(berryaux))
    bfxz=-2.*(((vx(i))*conjg(vz(i)))/(berryaux))
    bfyy=-2.*(((vy(i))*conjg(vy(i)))/(berryaux))
    bfyz=-2.*(((vy(i))*conjg(vz(i)))/(berryaux))
    bfzz=-2.*(((vz(i))*conjg(vz(i)))/(berryaux))

    end if

    bxx=bxx+bfxx
    bxy=bxy+bfxy
    bxz=bxz+bfxz
    byy=byy+bfyy
    byz=byz+bfyz
    bzz=bzz+bfzz

    end do

    end do



end subroutine berryct2


subroutine berryct(nocp,kx,ky,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
          ihopmatrices,efermi,energias,autovetores,nthread,gammas,berry)
    integer :: nthread
    double precision :: kx,ky,kz,ezz
    integer :: nocp,w90basis
    integer :: nvec
    integer :: n,i,j
    double precision,parameter :: ktol = 0.001
    integer,dimension(w90basis) :: ffactor
    double precision,dimension(3,3) :: rlat
    integer,dimension(nvec,3) :: rvec
    integer :: nsite
    double complex, dimension(w90basis,w90basis) :: hx,hy,hz
    double complex, dimension(w90basis,w90basis) :: autovetores,autovetores2
    double precision,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices
    double precision, dimension(w90basis) :: energias,energias2
    !double precision,dimension(3) :: mag
    double precision :: efermi
    double complex, dimension(w90basis) :: vflagx,vflagy
    double complex, dimension(w90basis) :: vn,vm

    double complex,dimension(w90basis) :: vx,vy

    double complex,dimension(w90basis) :: vxa,vya
    double complex, dimension(w90basis) :: vflagxa,vflagya
    double complex, dimension(w90basis) :: vna,vma

    double complex :: berry,berryf

    double precision:: gammas

    double precision :: berryaux


    
    call hlm(kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
         ihopmatrices,hx,hy,hz)
    

    berry=cmplx(0.,0.)

    do n=1,nocp,1
    

    do i=nocp+1,w90basis,1

    call vecconjg(autovetores(n,:),w90basis,vn)

    call matvec(hx,autovetores(i,:),w90basis,vflagx)

    call prodintsq(vn,vflagx,w90basis,vx(i))

    !write(*,*) vx(i)


    call vecconjg(autovetores(i,:),w90basis,vm)

    call matvec(hy,autovetores(n,:),w90basis,vflagy)

    call prodintsq(vm,vflagy,w90basis,vy(i))

    !write(*,*) vy(i)

    if (i == n) then

    berryf=0.


    else if (abs(energias(i)-energias(n)) .lt. 0.0001 ) then

    
    !call eigsys(nthread,kx+ktol,ky+ktol,w90basis,nvec,rlat1,rlat2,rlat3,rvec,hopmatrices,&
    !      ihopmatrices,efermi,energias2,autovetores2)

    !call vecconjg(autovetores2(n,:),w90basis,vna)

    !call matvec(hx,autovetores2(i,:),w90basis,vflagxa)

    !call prodintsq(vna,vflagxa,w90basis,vxa(i))


    !call vecconjg(autovetores2(i,:),w90basis,vma)

    !call matvec(hy,autovetores2(n,:),w90basis,vflagya)

    !call prodintsq(vma,vflagya,w90basis,vya(i))


    !berryf=-2.*(((vxa(i))*vya(i))/((energias2(n)-energias2(i))**2+(gammas)**2))

    berryf=0.

    else
    
    berryaux= (energias(n)-energias(i))**2+(gammas)**2
    berryf=-2.*(((vx(i))*vy(i))/(berryaux))

    end if

    berry=berry+berryf

    end do

    end do



end subroutine berryct


subroutine berryc(n,nocp,kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
          ihopmatrices,efermi,energias,autovetores,nthread,gammas,berry)
           

!   implicit none

    integer :: nthread

    double precision :: kx,ky,kz,ezz
    integer :: nocp,w90basis
    integer :: nvec
    integer :: n,i,j

    double precision,parameter :: ktol = 0.001


    double precision,dimension(3,3) :: rlat
    integer,dimension(nvec,3) :: rvec
    integer,dimension(w90basis) :: ffactor
    integer :: nsite

    double complex, dimension(w90basis,w90basis) :: hx,hy,hz
    double complex, dimension(w90basis,w90basis) :: autovetores,autovetores2

    double precision,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices

    double precision, dimension(w90basis) :: energias,energias2

    !double precision,dimension(3) :: mag

    double precision :: efermi

    double complex, dimension(w90basis) :: vflagx,vflagy
    double complex, dimension(w90basis) :: vn,vm

    double complex,dimension(w90basis) :: vx,vy

    double complex,dimension(w90basis) :: vxa,vya
    double complex, dimension(w90basis) :: vflagxa,vflagya
    double complex, dimension(w90basis) :: vna,vma

    double complex :: berry,berryf

    double precision:: gammas


    
    call hlm(kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
         ihopmatrices,hx,hy,hz)
    

    berry=cmplx(0.,0.)
    

    do i=nocp+1,w90basis,1

    call vecconjg(autovetores(n,:),w90basis,vn)

    call matvec(hx,autovetores(i,:),w90basis,vflagx)

    call prodintsq(vn,vflagx,w90basis,vx(i))

    !write(*,*) vx(i)


    call vecconjg(autovetores(i,:),w90basis,vm)

    call matvec(hy,autovetores(n,:),w90basis,vflagy)

    call prodintsq(vm,vflagy,w90basis,vy(i))

    !write(*,*) vy(i)

    if (i == n) then

    berryf=0.


    else if (abs(energias(i)-energias(n)) .lt. 0.0001 ) then

    
    !call eigsys(nthread,kx+ktol,ky+ktol,w90basis,nvec,rlat1,rlat2,rlat3,rvec,hopmatrices,&
    !      ihopmatrices,efermi,energias2,autovetores2)

    !call vecconjg(autovetores2(n,:),w90basis,vna)

    !call matvec(hx,autovetores2(i,:),w90basis,vflagxa)

    !call prodintsq(vna,vflagxa,w90basis,vxa(i))


    !call vecconjg(autovetores2(i,:),w90basis,vma)

    !call matvec(hy,autovetores2(n,:),w90basis,vflagya)

    !call prodintsq(vma,vflagya,w90basis,vya(i))


    !berryf=-2.*(((vxa(i))*vya(i))/((energias2(n)-energias2(i))**2+(gammas)**2))

    berryf=0.

    else

    berryf=-2.*(((vx(i))*vy(i))/((energias(n)-energias(i))**2+(gammas)**2))

    end if

    berry=berry+berryf

    end do
end subroutine berryc





!INCLUDE "./subprograms/berry_curvature_bz-tool.f90"

!gfortran berry_curvature_kpath-tool.f90 -o berry_ct-tool.x -llapack95 -lopenblas -fopenmp

subroutine berrycurvbz(nthreads,outputfolder,params,sme,ngrid,mshift)

!   use omp_lib
!   use module_hamiltonian   !_input_variables

!   implicit none

    integer :: i,j,k,erro,j2,m
    double precision,parameter :: pi=acos(-1.)


    double complex,allocatable,dimension(:,:) :: htb
    double precision,allocatable,dimension(:) :: eigv

    double complex,allocatable,dimension(:,:) :: autovetores

    double precision,allocatable,dimension(:,:) :: kpts

    double precision :: kx,ky,kz,kp

    double precision,allocatable,dimension(:,:) :: eigvf

    double complex,allocatable,dimension(:,:,:) :: vector

    double complex :: bxx,bxy,bxz,byy,byz,bzz
    double precision,parameter :: gammas= 0.001

    double precision,allocatable,dimension(:,:) :: output

    !double precision,dimension(3) :: shift !shift no mhkpack

    !variaveis kpoints
    integer :: nk
    integer :: nks
    integer :: nkpts
    double precision,allocatable,dimension(:,:) :: ks

    integer,allocatable,dimension(:) :: nocpk 

    !modificacoes versao 2.1

    integer :: nthreads
    integer,dimension(3) :: ngrid
    integer :: nc,nv
    !integer :: ncrpa,nvrpa
    !integer :: ncbz,nvbz
    double precision :: edos0,edosf,numdos
    double precision :: ebse0,ebsef,numbse
    double precision :: sme,exc,rk
    double precision,dimension(3) :: mshift
    double precision :: ktol
    character(len=70) :: params   !parametros TB
    character(len=70) :: orbw     !peso orbitais lcount
    character(len=70) :: kpaths    !kpath
    character(len=70) :: kpathsbse    !kpath
    !character(len=70) :: diein    !ambiente dieletrico
    character(len=70) :: outputfolder    !pasta saida
    character(len=70) :: calcparms
    character(len=70) :: meshtype
    character(len=5) :: coultype
    double precision,dimension(3) :: ediel

    !fim modificacoes versao 2.1

    !call input_read


    !OPEN(UNIT=202, FILE= kpaths,STATUS='old', IOSTAT=erro)
        !if (erro/=0) stop "Erro na abertura do arquivo de entrada kpath"

    !OUTPUT : criando arquivos de saida
    OPEN(UNIT=300, FILE=trim(outputfolder)//"berry_curv_bz.dat",STATUS='unknown', IOSTAT=erro)
        if (erro/=0) stop "Error opening berry_curv_bz output file"

!   call OMP_SET_NUM_THREADS(nthreads)
    call hamiltonian_input_read(200,params)

    !ediel(2) = edielh

    !read(202,*) nks
    !read(202,*) nkpts

    !allocate(ks(nks,3))

    !do i=1,nks

    !   read(202,*) ks(i,1),ks(i,2),ks(i,3)
    
    !end do


    !termino leitura parametros

    !parametros do calculo


    !termino parametros calculo

 
    allocate(autovetores(w90basis,w90basis),eigv(w90basis))


    allocate(kpts(ngrid(1)*ngrid(2)*ngrid(3),3))

    !shift = 0.0
    call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpts)

    !definindo kpath


    allocate(eigvf(ngrid(1)*ngrid(2)*ngrid(3),w90basis),vector(ngrid(1)*ngrid(2)*ngrid(3),w90basis,w90basis))

    allocate(nocpk(ngrid(1)*ngrid(2)*ngrid(3)))

    ! $omp parallel do
    do i=1,ngrid(1)*ngrid(2)*ngrid(3)

        call eigsys(nthreads,scs,exc,nocpk(i),ffactor,kpts(i,1),kpts(i,2),kpts(i,3),w90basis,nvec,rlat,rvec,hopmatrices,&
           ihopmatrices,efermi,eigv,autovetores)

            do k=1,w90basis


                eigvf(i,k)=eigv(k)

                do m=1,w90basis

                vector(i,k,m)=autovetores(k,m)

                end do

                
            end do

            

    end do
    ! $omp end parallel do

    !termino definicao kpath

    allocate(output(ngrid(1)*ngrid(2)*ngrid(3),9))
    
    if(Process==0) write(300,*) "kx ky kz yz xz xy" 

    !$omp parallel do private(j,bxx,bxy,bxz,byy,byz,bzz)    

    do j=1,ngrid(1)*ngrid(2)*ngrid(3)




        call berryct2(nocpk(j),kpts(j,1),kpts(j,2),ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
            ihopmatrices,efermi,eigvf(j,:),vector(j,:,:),nthreads,gammas,bxx,bxy,bxz,byy,byz,bzz)
        ! $omp ordered
             !write(300,*) kpts(j,1),kpts(j,2),aimag(berry)
        ! $omp end ordered
            output(j,1) = kpts(j,1)
            output(j,2) = kpts(j,2)
            output(j,3) = kpts(j,3)
            output(j,4) = aimag(bxx)
            output(j,5) = aimag(bxy)
            output(j,6) = aimag(bxz)
            output(j,7) = aimag(byy)
            output(j,8) = aimag(byz)
            output(j,9) = aimag(bzz)
    

    end do
    !$omp end parallel do

 if(Process==0) then
    do j=1,ngrid(1)*ngrid(2)*ngrid(3)
        write(300,"(6F15.4)") real(output(j,1)),real(output(j,2)),real(output(j,3)),output(j,8),output(j,6),output(j,5)
    end do
 endif

    deallocate(rvec,hopmatrices,ihopmatrices,eigv)
    deallocate(ffactor)
    deallocate(kpts)
    deallocate(autovetores)
    deallocate(eigvf,vector)
    deallocate(output)


 if(Process==0) then
    close(200)
    close(201)
    close(300)
 endif
end subroutine berrycurvbz


!INCLUDE "./subprograms/berry_curvature_kpath-tool.f90"
!gfortran berry_curvature_kpath-tool.f90 -o berry_ct-tool.x -llapack95 -lopenblas -fopenmp

subroutine berrycurv(nthreads,outputfolder,params,kpaths,sme)

!   use omp_lib
!   use module_hamiltonian   !_input_variables

!   implicit none

    integer :: i,j,k,erro,j2,m
    double precision,parameter :: pi=acos(-1.)


    double complex,allocatable,dimension(:,:) :: htb
    double precision,allocatable,dimension(:) :: eigv

    double complex,allocatable,dimension(:,:) :: autovetores

    double precision,allocatable,dimension(:,:) :: kpts

    double precision :: kx,ky,kz,kp

    double precision,allocatable,dimension(:,:) :: eigvf

    double complex,allocatable,dimension(:,:,:) :: vector

    double complex :: bxx,bxy,bxz,byy,byz,bzz
    double precision,parameter :: gammas= 0.001

    integer,allocatable,dimension(:) :: nocpk 

    !variaveis kpoints

    integer :: nks,nocp1
    integer :: nkpts
    double precision,allocatable,dimension(:,:) :: ks

    !modificacoes versao 2.1

    integer :: nthreads
    integer,dimension(3) :: ngrid
    integer :: nc,nv
    !integer :: ncrpa,nvrpa
    !integer :: ncbz,nvbz
    double precision :: edos0,edosf,numdos
    double precision :: ebse0,ebsef,numbse
    double precision :: sme,exc,rk
    double precision,dimension(3) :: mshift
    double precision :: ktol
    character(len=70) :: params   !parametros TB
    character(len=70) :: orbw     !peso orbitais lcount
    character(len=70) :: kpaths    !kpath
    character(len=70) :: kpathsbse    !kpath
    !character(len=70) :: diein    !ambiente dieletrico
    character(len=70) :: outputfolder    !pasta saida
    character(len=70) :: calcparms
    character(len=70) :: meshtype
    character(len=5) :: coultype
    double precision,dimension(3) :: ediel

    !fim modificacoes versao 2.1

    !call input_read


    OPEN(UNIT=202, FILE= kpaths,STATUS='old', IOSTAT=erro)
        if (erro/=0) stop "Error opening kpath input file"

    !OUTPUT : criando arquivos de saida
    OPEN(UNIT=300, FILE=trim(outputfolder)//"berry_curv_kpath.dat",STATUS='unknown', IOSTAT=erro)
        if (erro/=0) stop "Error opening berry_curv_kpath output file"

    call hamiltonian_input_read(200,params)

    !ediel(2) = edielh


    read(202,*) nks
    read(202,*) nkpts

    allocate(ks(nks,3))

    do i=1,nks

        read(202,*) ks(i,1),ks(i,2),ks(i,3)
    
    end do



    !termino leitura parametros

    !parametros do calculo


    !termino parametros calculo 

    allocate(autovetores(w90basis,w90basis),eigv(w90basis))
    
    !allocate(kpts(nkpts*(nks-1),4))
    !allocate(nocpk(nkpts*(nks-1)))
    
    allocate(kpts((nks/2)*nkpts,4))
    allocate(nocpk((nks/2)*nkpts))

    !definindo kpath


    call kpath2(outputfolder,rlat(1,:),rlat(2,:),rlat(3,:),nks,ks,nkpts,kpts)

    allocate(eigvf(nkpts*(nks-1),w90basis),vector(nkpts*(nks-1),w90basis,w90basis))


    ! $omp parallel do
    !do i=1,nkpts*(nks-1)
    do i=1,(nks/2)*nkpts

        call eigsys(nthreads,scs,exc,nocpk(i),ffactor,kpts(i,2),kpts(i,3),kpts(i,4),w90basis,nvec,rlat,rvec,hopmatrices,&
           ihopmatrices,efermi,eigv,autovetores)

            do k=1,w90basis


                eigvf(i,k)=eigv(k)

                do m=1,w90basis

                vector(i,k,m)=autovetores(k,m)

                end do

                
            end do

            

    end do
    ! $omp end parallel do
    !termino definicao kpath
    
    if(Process==0) write(300,*) "kp yz xz xy"

    !$omp do ordered
    !do j=1,nkpts*(nks-1)
    do j=1,(nks/2)*nkpts
            kp= kpts(j,1)
        kx= kpts(j,2)
        ky= kpts(j,3)
        kz= kpts(j,4)
        call berryct2(nocpk(j),kx,ky,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
            ihopmatrices,efermi,eigvf(j,:),vector(j,:,:),nthreads,gammas,bxx,bxy,bxz,byy,byz,bzz)
        !$omp ordered
             if(Process==0) write(300,"(4F15.4)") kp,aimag(byz),aimag(bxz),aimag(bxy)
        !$omp end ordered
    end do
    !$omp end do

    deallocate(rvec,hopmatrices,ihopmatrices,eigv)
    deallocate(ffactor)
    deallocate(kpts)
    deallocate(autovetores)
    deallocate(eigvf,vector)
    deallocate(nocpk)

 if(Process==0) then
    close(200)
    close(202)
    close(300)
 endif
end subroutine berrycurv


  end module module_berry
  
  


