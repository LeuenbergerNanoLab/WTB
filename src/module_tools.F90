


   module module_tools
    use MPI_P
	use module_hamiltonian
	use module_kmesh
	implicit none
   
   contains

!INCLUDE "./subroutines/dos_subs.f90"
subroutine lcount(w90basis,autovetor,dcont,pcont)
	integer :: i
	integer :: w90basis,wbasisaux
	double complex,dimension(w90basis) :: autovetor
	double precision :: dcont,pcont

	dcont = 0.0
	pcont = 0.0

	wbasisaux=w90basis/2

	do i=1,5

		dcont=dcont+REAL(autovetor(i)*conjg(autovetor(i)))
		dcont=dcont+REAL(autovetor(wbasisaux+i)*conjg(autovetor(wbasisaux+i)))

	end do

	do i=1,6

		pcont=pcont+REAL(autovetor(5+i)*conjg(autovetor(5+i)))
		pcont=pcont+REAL(autovetor(wbasisaux+5+i)*conjg(autovetor(wbasisaux+5+i)))

	end do
end subroutine lcount



subroutine spincont(w90basis,autovetor,cup,cdown,dft,systype)
	integer :: i,j
	integer :: w90basis,wbasisaux
	double complex,dimension(w90basis) :: autovetor
	double precision :: cup,cdown
	character(len=1) :: dft
	character(len=4) :: systype	

	cup=0d0
	cdown=0d0
	
	wbasisaux=w90basis/2
	
	if (systype .eq. "NP") then 
	
		cup=0d0
		cdown=0d0
	
	else if (systype .eq. "SP") then
	
	 do i=1,wbasisaux

		cup=cup+REAL(autovetor(i)*conjg(autovetor(i)))

		cdown=cdown+REAL(autovetor(wbasisaux+i)*conjg(autovetor(wbasisaux+i)))


	 end do	
	
	
	else if ((systype .eq. "SOC") .and. (dft .eq. "V") ) then
	
	 do i=1,wbasisaux

		cup=cup+REAL(autovetor(i)*conjg(autovetor(i)))

		cdown=cdown+REAL(autovetor(wbasisaux+i)*conjg(autovetor(wbasisaux+i)))


	 end do
	
	else
	
	 do i=1,wbasisaux

		cup=cup+REAL(autovetor((2*i)-1)*conjg(autovetor((2*i)-1)))

		cdown=cdown+REAL(autovetor(2*i)*conjg(autovetor(2*i)))


	 end do
	
	
	end if
	
	

	
end subroutine spincont



subroutine endos(rlat,dimse,n1,n2,n3,en,etb,fosc,sme,ints)
	character(len=2) :: systype

	double precision :: en
	integer :: dimse !dimensão do sistema => nkpts*nbands
	double precision :: sme
	double precision ::  ints
	integer :: n1,n2,n3,i,j

	double precision :: vc

	double complex :: intsaux

	double precision,dimension(3,3) :: rlat

	double precision,dimension(dimse) :: etb
	double precision,dimension(dimse) :: fosc
	
	double precision :: caux,dimk    !,gaussian

	double precision,parameter :: pi=acos(-1.)
	

	dimk = dble(n1*n2*n3)

	caux=(1.0/dimk)


	ints=0.


	do i=1,dimse

		
		intsaux = fosc(i)*gaussian(en-etb(i),sme)


		ints=ints+intsaux


	end do


	ints=ints*caux!*(1./elux)


end subroutine endos



subroutine layercont(w90basis,wf,orbweight,lc)
	integer :: w90basis
	double complex,dimension(w90basis) :: wf
	double precision,dimension(w90basis) :: orbweight 
	integer :: i,j

	double complex :: lcaux
	double precision :: lc

	lcaux = 0.0


	do i=1,w90basis

		lcaux = wf(i)*orbweight(i)*conjg(wf(i))+lcaux

	end do	


	lc  = real(lcaux)


end subroutine layercont




!INCLUDE "./subprograms/bands-kpath-tool.f90"


!gfortran bands-kpath-tool.f90 -o bands_tool.x -llapack95 -lopenblas -fopenmp

subroutine bandstool(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,rk,dft)

!	use omp_lib
!	use module_hamiltonian    !_input_variables
!	implicit none

	integer :: i,j,k,erro,j2
	double precision,parameter :: pi=acos(-1.)

	!double complex,allocatable,dimension(:,:) :: htb
	double precision,allocatable,dimension(:) :: eigv

	double precision,allocatable,dimension(:,:,:) :: ebands

	double complex,allocatable,dimension(:,:) :: autovetores

	double precision,allocatable,dimension(:,:) :: kpts

	double precision,allocatable,dimension(:) :: orbweight

	double precision :: kx,ky,kz,kp

	double precision :: lco

	double complex :: spz

	!variaveis kpoints

	integer :: nks
	integer :: nkpts
	double precision,allocatable,dimension(:,:) :: ks

	!call input_read

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
	character(len=1) :: dft	

	!fim modificacoes versao 2.1

	! INPUT : lendo os parametros do modelo de tight-binding
	
	OPEN(UNIT=202, FILE= kpaths,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening kpath input file"


	!OUTPUT : criando arquivos de saida
	OPEN(UNIT=300, FILE=trim(outputfolder)//"bands.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bands output file"

	!OPEN(UNIT=500, FILE=trim(outputfolder)//"nocp.dat",STATUS='unknown', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de saida nocp"

!	call OMP_SET_NUM_THREADS(nthreads)

	call hamiltonian_input_read(200,params)

	!ediel(2) = edielh


	

	read(202,*) nks
	read(202,*) nkpts

	allocate(ks(nks,3))

	do i=1,nks

		read(202,*) ks(i,1),ks(i,2),ks(i,3)
	
	end do

	allocate(orbweight(w90basis))

	OPEN(UNIT=203, FILE= orbw,STATUS='old', IOSTAT=erro)
    	  	if (erro/=0) then
	
	orbweight= 1.0

	else

	do i=1,w90basis
		read(203,*) orbweight(i)
	end do

	end if



	!termino leitura parametros

	!parametros do calculo

	!OPEN(UNIT=2077, FILE= trim(calcparms)//"bands_calc.dat",STATUS='unknown', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de saida parametros calculo"

	!call param_out(2077,nthreads,outputfolder,calcparms,ngrid,nc,nv,edos0,edosf,numdos, &
	!	     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
	!	     exc,mshift,coultype)

	!close(2077)

	!termino parametros calculo
 

	allocate(autovetores(w90basis,w90basis),eigv(w90basis))
	allocate(ebands(w90basis,nkpts*(nks-1),4))
	!allocate(kpts(nkpts*(nks-1),4))
	allocate(kpts((nks/2)*nkpts,4))

	!definindo kpath

	call kpath2(outputfolder,rlat(1,:),rlat(2,:),rlat(3,:),nks,ks,nkpts,kpts)

	!termino definicao kpath

	ebands=0.0

	!do j=1,nkpts*(nks-1)
	do j=1,(nks/2)*nkpts

	        kp= kpts(j,1)
		kx= kpts(j,2)
		ky= kpts(j,3)
		kz= kpts(j,4)

		call eigsys(nthreads,scs,exc,nocp,ffactor,kx,ky,kz,w90basis,nvec,rlat,rvec,hopmatrices,&
		    ihopmatrices,efermi,eigv,autovetores)

		!write(500,*) nocp,kx,ky,kz 


		do i=1,w90basis

			call layercont(w90basis,autovetores(i,:),orbweight,lco)
			call spinvl(w90basis,autovetores(i,:),dft,systype,spz)

			ebands(i,j,1) = kp
			ebands(i,j,2) = eigv(i)
			ebands(i,j,3) = real(spz)
			ebands(i,j,4) = lco
			


		end do





	end do


 if(Process==0) then
	do i=1,w90basis
			write(300,*) "#band",i
		do j=1,(nks/2)*nkpts
			write(300,*) real(ebands(i,j,1)),real(ebands(i,j,2)),real(ebands(i,j,3)),real(ebands(i,j,4))
		end do
			write(300,*)
	end do
 endif
 
	deallocate(rvec,hopmatrices,ihopmatrices,ffactor)
	deallocate(eigv)
	deallocate(kpts)
	deallocate(orbweight)
	deallocate(autovetores)
	deallocate(ebands)
	deallocate(ks)

 if(Process==0) then
	close(200)
	close(300)
	close(202)
	close(203)
 endif
end subroutine bandstool


!INCLUDE "./subprograms/tdos-tool.f90"
!gfortran tdos-tool.f90 -o tdos-tool.x -llapack95 -lopenblas -fopenmp
subroutine dostool(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,rk,meshtype,dft)
	integer :: i,j,k,erro
	double precision,parameter :: pi=acos(-1.)
	double complex,parameter :: imag=cmplx(0.0,1.0)
	integer :: ngkpt
	!double precision,dimension(3) :: shift !shift no mhkpack
	double complex,allocatable,dimension(:,:) :: autovetores
	double precision,allocatable,dimension(:) :: eigv
	double precision,allocatable,dimension(:,:) :: kpts
	double precision,allocatable,dimension(:,:) :: res
	double precision,allocatable,dimension(:) :: orbweight
	double precision :: kx,ky,kz,kp
	integer  :: ndim,counter
	!variaveis para dos
	double precision :: tdos,wdos
	double precision,allocatable,dimension(:)   :: atdos,awdos
	double precision  :: cup,cdown
	double precision,allocatable,dimension(:)  :: updos,dndos,en
	double precision :: deltaen
	double precision :: sumdos
	double precision :: spinf
	!variaveis para grafico dos
	double precision :: at0
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
	character(len=1) :: dft
	double precision,dimension(3) :: ediel

	!fim modificacoes versao 2.1

	!call input_read


	! INPUT : lendo os parametros do modelo de tight-binding

	!OPEN(UNIT=201, FILE= diein,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada ambiente dieletrico"



	!OUTPUT : criando arquivos de saida
	OPEN(UNIT=300, FILE=trim(outputfolder)//"dos.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening dos output file"


!	call OMP_SET_NUM_THREADS(nthreads)
	!call OPENBLAS_SET_NUM_THREADS(nthread)

	call hamiltonian_input_read(200,params)

	!ediel(2) = edielh
	
	
	if ( systype .eq. "NP" ) then
	 spinf = 2.0
	 if(Process==0) write(300,*) "#energy tdos sumdos awdos"
	else
	 spinf = 1.0
	 if(Process==0) write(300,*) "#energy tdos sumdos updos dndos awdos"
	end if

	allocate(orbweight(w90basis))

 if(Process==0) then
	OPEN(UNIT=203, FILE= orbw,STATUS='old', IOSTAT=erro)
    if (erro/=0) then
	 orbweight= 1.0
	else
	 do i=1,w90basis
		read(203,*) orbweight(i)
	 end do
	end if
 endif
    call P_sendR1(orbweight,w90basis)

	!if (meshtype .eq. "RK3D") then

	!	call rkmesh(rk,rlat,ngrid)

	!else if (meshtype .eq. "RK2D") then

	!	call rkmesh2D(rk,rlat,ngrid)
	!else

	!	continue
	!end if


	allocate(eigv(w90basis))
	allocate(autovetores(w90basis,w90basis))

	ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
	allocate(kpts(ngkpt,3))


	!parametros do calculo

	!OPEN(UNIT=2077, FILE= trim(calcparms)//"dos_calc.dat",STATUS='unknown', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de saida parametros calculo"

	!call param_out(2077,nthreads,outputfolder,calcparms,ngrid,nc,nv,edos0,edosf,numdos, &
	!	     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
	!	     exc,mshift,coultype)

	!close(2077)

	!termino parametros calculo

	!shift = 0.0
	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpts)


	counter=0

	ndim=ngkpt*w90basis

	allocate(res(ndim,5))

	edos0= 0.0
	edosf= 0.0

	do j=1,ngkpt

		kx= kpts(j,1)
		ky= kpts(j,2)
		kz= kpts(j,3)

		call eigsys(nthreads,scs,exc,nocp,ffactor,kx,ky,kz,w90basis,nvec,rlat,rvec,hopmatrices,&
		    ihopmatrices,efermi,eigv,autovetores)

		if (eigv(1) .lt. edos0 ) then
		
			edos0 = eigv(1)
			
		end if
		
		
		if (eigv(w90basis) .gt. edosf ) then
		
			edosf = eigv(w90basis)
			
		end if		


		do i=1,w90basis

			
			tdos = 1.0
			call layercont(w90basis,autovetores(i,:),orbweight,wdos)
			call spincont(w90basis,autovetores(i,:),cup,cdown,dft,systype)

			counter=counter+1

			res(counter,1) = eigv(i)
			res(counter,2) = tdos
			res(counter,3) = cup
			res(counter,4) = cdown
			res(counter,5) = wdos
			!res(counter,4) = dcont
			!res(counter,5) = pcont
			!res(counter,6) = 1.0


		end do


	end do
	
	allocate(atdos(int(numdos)+1),updos(int(numdos)),dndos(int(numdos)),awdos(int(numdos)))
	allocate(en(int(numdos)))
	
	deltaen= dble(((edosf-edos0)))/(numdos-1.)
	
	sumdos = 0.0
	atdos = 0.0
	
	! $omp do ordered
	do i=1,int(numdos)

		en(i) = edos0 + dble(((edosf-edos0)*(i-1))/(numdos-1.))


		call endos(rlat,ndim,ngrid(1),ngrid(2),ngrid(3),en(i),res(:,1),res(:,2),sme,atdos(i))
		call endos(rlat,ndim,ngrid(1),ngrid(2),ngrid(3),en(i),res(:,1),res(:,3),sme,updos(i))
		call endos(rlat,ndim,ngrid(1),ngrid(2),ngrid(3),en(i),res(:,1),res(:,4),sme,dndos(i))
		call endos(rlat,ndim,ngrid(1),ngrid(2),ngrid(3),en(i),res(:,1),res(:,5),sme,awdos(i))
		!call endos(a0,ndim,ngrid,ngrid,1,en,res(:,1),res(:,6),sme,tcont)

		!sumdos = sumdos+ deltaen*(at0)

		! $omp ordered
		!write(300,"(6F15.4)") en,atdos,sumdos,updos,-dndos,awdos
		!write(301,*) en,d2cont,p2cont,tcont
		! $omp end ordered

	end do
	! $omp end do
	

	
	atdos = spinf*atdos
	awdos = spinf*awdos
	
	do i=1,int(numdos)
	
		if ( systype .eq. "NP" ) then
	
		if(Process==0) write(300,"(4F15.4)") en(i),atdos(i),sumdos,awdos(i)
		
		else
		
		if(Process==0) write(300,"(6F15.4)") en(i),atdos(i),sumdos,updos(i),-dndos(i),awdos(i)	

		end if
		
		if (i .eq. 1) then
		
		sumdos= sumdos 
		
		else
		
		sumdos= sumdos + (deltaen*(atdos(i-1)+atdos(i))*0.5)
		
		end if
	
	end do

	deallocate(rvec,hopmatrices,ihopmatrices,ffactor)

	deallocate(eigv)
	deallocate(autovetores)
	deallocate(kpts)
	deallocate(res)
	deallocate(orbweight)
	deallocate(atdos,updos,dndos,awdos)


 if(Process==0) then
	close(200)
	close(300)
 endif
end subroutine dostool



   end module module_tools



