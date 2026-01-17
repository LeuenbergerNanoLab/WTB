



  module module_bse
   use MPI_P                                                             ! MPI module 
   use module_slepc_solver
   use module_hamiltonian
   use module_coulomb
   use module_kmesh
   use module_optics
   implicit none
   integer                            :: dimbse                        ! size of HBSE
   COMPLEX(8),allocatable             :: hbse(:,:)                     ! BSE Hamiltonian
   real(8),allocatable                :: Ws(:)                          ! eigenvalues
   complex(8), allocatable            :: HE(:,:)           ! for slepc test
   logical                            :: ltest
    integer                           :: ngkpt                  ! number of k-points
    real(8),allocatable               :: kpt(:,:)               ! k-points
  contains


!INCLUDE "./subroutines/bse_subs.f90"
complex(8) function matrizelbse(coultype,tolr,w90basis,ediel,lc,ez,w,r0,ngrid,rlat,est1,ec1,ev1,vbc1 &
         ,vbv1,kpt1,est2,ec2,ev2,vbc2,vbv2,kpt2) !funcao para calcular o elemento de matriz da matriz bse
	implicit none
	character(len=5) :: coultype
	integer,dimension(3) :: ngrid
	integer :: w90basis
	double precision :: a,vcell1
	double precision :: ez,w
	integer, dimension(4) :: est1,est2
	double precision :: ec1,ec2,ev1,ev2
	double precision,dimension(3) :: kpt1,kpt2
	double complex, dimension(w90basis) :: vbc1,vbc2,vbv1,vbv2
	double complex, dimension(w90basis) :: vbc,vbv
	double precision,dimension(3,3) :: rlat
	double precision :: tolr
	integer :: ktol
	double precision,parameter:: pi=acos(-1.)
	double precision :: auxi
	double precision :: modk
	double precision,dimension(3) :: ediel
	double precision :: lc
	double complex :: vc,vv
	double precision :: vcoul1
	double precision :: r0

	select case (coultype)

	case("V2DK")
		vcoul1= v2dk(kpt1,kpt2,ediel,rlat,ngrid,lc,tolr)
	case("V3D")
		vcoul1= vcoul(kpt1,kpt2,rlat,ngrid,tolr)
	case("V3DL")
		vcoul1= v3diel(kpt1,kpt2,ediel,rlat,ngrid,tolr)
	case("V2DT")
		vcoul1= v2dt(kpt1,kpt2,ngrid,rlat,tolr)
	case("V2DT2")
		vcoul1= v2dt2(kpt1,kpt2,ngrid,rlat,lc,tolr)
	case("V2DOH")
		vcoul1= v2dohono(kpt1,kpt2,ngrid,rlat,ediel,w,ez,tolr)
	case("V2DRK")
		vcoul1= v2drk(kpt1,kpt2,ngrid,rlat,ediel,lc,ez,w,r0,tolr)
	case("V1DT")
		vcoul1=	v1dt(kpt1,kpt2,ngrid,rlat,tolr,lc)	
	case("V0DT")
		vcoul1= v0dt(kpt1,kpt2,ngrid,rlat,tolr)
	case default
		if(Process==0) write(*,*) "Wrong Coulomb Potential"
		STOP
	end select

	if (est1(1) .eq. est2(1)) then
		matrizelbse= (ec1-ev1) + vcoul1
	else
		call vecconjg(vbc1,w90basis,vbc)
		call vecconjg(vbv1,w90basis,vbv)
		call prodintsq(vbc,vbc2,w90basis,vc)
		call prodintsq(vbv,vbv2,w90basis,vv)
		matrizelbse=  vcoul1*vc*vv
	end if
end function matrizelbse




subroutine opticalactivity(dimse,excitonvec,hopt,activity,description)
	integer :: dimse
	double precision,dimension(dimse) :: activity
	double complex :: actaux
	double precision :: description(dimse)
	double complex,dimension(dimse) :: hopt
	double complex,dimension(dimse,dimse) :: excitonvec
	integer :: i,j

	activity=0.
	actaux=cmplx(0.,0.)

	do i=1,dimse
		do j=1,dimse
        actaux=actaux+(excitonvec(j,i)*hopt(j))
		end do
		activity(i)=actaux*conjg(actaux)
		if (activity(i) .gt. 0.1) then
			description(i)= 1.
		else if ( (activity(i) .lt. 0.1) .and. (activity(i) .gt. 1.0e-8)) then
			description(i)= 0.
		else
			description(i)= -1.
		end if
		actaux=cmplx(0.,0.)
	end do
end subroutine opticalactivity



subroutine excitonil(i,w90basis,nc,nv,ncount,stt,ndim,wf,lc,pinter,pintra)
	integer :: i !numero do estado excitonico
	integer :: j,k,nc,nv
	integer :: w90basis
	integer :: ncount,ndim
	integer,dimension(ncount*nc*nv,4) :: stt
	double complex,dimension(ndim) :: wf
	double precision,dimension(ncount,2*w90basis) :: lc
	double precision :: pinter,pintra
	double precision :: inter,intra
	double precision :: p1c,p1v,p2c,p2v

	pinter = 0d0
	pintra = 0d0

	do j=1,ndim
		p1c = 1d0-lc(stt(j,4),stt(j,3))
		p1v = 1d0-lc(stt(j,4),stt(j,2))
		p2c = lc(stt(j,4),stt(j,3))
		p2v = lc(stt(j,4),stt(j,2))
		intra = (p1c*p1v)+(p2c*p2v)
		inter = (p1c*p2v)+(p2c*p1v)
		pinter = pinter + Real(conjg(wf(j))*wf(j)*inter,8)
		pintra = pintra + Real(conjg(wf(j))*wf(j)*intra,8)
	end do
end subroutine excitonil



subroutine dielbse(dimse,dimse1,excitonvec,hopt1,hopt2,activity)
	integer                       :: dimse,dimse1
	double precision,dimension(:) :: activity           !, exciton
	double precision              :: actaux
	double complex,dimension(:)   :: hopt1,hopt2
	double complex,dimension(:,:) :: excitonvec
	integer                       :: i,j,k

	activity=0.0
	actaux=0.0

	do i=1,dimse1
		do j=1,dimse
			do k=1,dimse
		actaux=actaux+(excitonvec(j,i)*hopt1(j)*conjg(excitonvec(k,i))*conjg(hopt2(k)))
		end do
			end do
		activity(i)=actaux
		actaux=0.0
	end do
end subroutine dielbse



subroutine dielbsep(nthread,dimse,dimse1,excitonvec,hopt1,hopt2,activity)
	integer                       :: dimse,dimse1,nthread
	double precision,dimension(:) :: activity    !, exciton
	double complex,dimension(:)   :: hopt1,hopt2
	double complex,dimension(:,:) :: excitonvec
!	integer,dimension(dimse,4)    :: stt
	double complex                :: actaux,actaux2
	integer                       :: i,j   !,k   !,no

!	call OMP_SET_NUM_THREADS(nthread)

	activity=0.0

	! $OMP DO PRIVATE(actaux)
	!$OMP PARALLEL DO PRIVATE(actaux,actaux2)
	do i=1,dimse1
		actaux=0.0
		actaux2=0.0
		do j=1,dimse
			actaux=actaux+(excitonvec(j,i)*hopt1(j))
			actaux2=actaux2+(excitonvec(j,i)*hopt2(j))
		end do
		activity(i)=actaux*conjg(actaux2)
	end do
	!$OMP END PARALLEL DO
end subroutine dielbsep



!	call   dielbsev(nthreads,dimbse,dimbse1,hbse,hrx,actxx)                     ! change dimbse to dimbse1
subroutine dielbsev(nthread,dimse,dimse1,excitonvec,hopt,activity)
	integer                   :: dimse,dimse1,nthread
	real(8),dimension(:)      :: activity         !, exciton
	complex(8)                :: actaux
	complex(8),dimension(:)   :: hopt
	complex(8),dimension(:,:) :: excitonvec
	integer                   :: i,j,k

!	call OMP_SET_NUM_THREADS(nthread)	

	activity=0.0
	actaux=0.0
	
!	if(Process==0) print *,'dielbsev: size(excitonvec,1)=',size(excitonvec,1)
!	if(Process==0) print *,'dielbsev: size(excitonvec,2)=',size(excitonvec,2)
!	if(Process==0) print *,'dielbsev: size(hopt,1)=',size(hopt,1)
!	if(Process==0) print *,'dielbsev: size(activity,1)=',size(activity,1)
!	if(Process==0) print *,'dielbsev: dimse=',dimse
!	if(Process==0) print *,'dielbsev: dimse1=',dimse1
! if(process==0.and.ltest) then
!  print *,'dielbsev: excitonvec'
!  do i=1,1
!   print 1,excitonvec(1:dimse,i)
!  enddo
!  print *,'dielbsev: hopt'
!  print 1,(hopt(1:dimse))
! endif
!1 format(100(2E12.3,3x))
!	hopt(:) = (0.d0,0.d0)           ! TEST

	!$OMP PARALLEL DO PRIVATE(actaux)
	do i=1,dimse1
		actaux=(0.d0,0.d0)
		do j=1,dimse
		 actaux=actaux+(excitonvec(j,i)*hopt(j))
		end do
		activity(i)=actaux*conjg(actaux)
	end do
	!$OMP END PARALLEL DO
! if(process==0.and.ltest) then
!  print *,'dielbsev: activity'
!  print 1,(activity(1:dimse))
!  ltest = .false.
! endif
end subroutine dielbsev



subroutine excwf(outputfolder,ngkpt,kpt,nc,nv,nocp,stt,excenergy,excnum,ewf)
	integer :: i,erro
	character(len=70) :: outputfolder    !pasta saida
	character(len=200) :: file1
	CHARACTER(LEN=30) :: Format
	integer :: ngkpt,nc,nv,excnum
	integer,dimension(ngkpt) :: nocp
	integer,dimension(ngkpt*nc*nv,4) :: stt
	double precision :: excenergy
	double precision,dimension(ngkpt,3) :: kpt
	double complex,dimension(ngkpt*nc*nv) :: ewf

 if(Process==0) then
	WRITE (file1,'(a7,I0,a4)') ,'exc_wf_',excnum,'.dat'
	Format = "(I,3F15.4,2I,2E15.4)"
	
	OPEN(UNIT=700+excnum, FILE=trim(outputfolder)//trim(file1),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening wf BSE output file"
    	
    	write(700+excnum,*) "#"," ","excitonic state",excnum
    	write(700+excnum,*) "#"," ","exciton energy:",excenergy
    	write(700+excnum,*) "#"," ","number occupied states",nocp(1),"(not valid for metallic systems)"
    	write(700+excnum,*) "#"," ","Number of conduction states","",nc
    	write(700+excnum,*) "#"," ","Number of valence states","",nv
    	write(700+excnum,*) "#"," ","Number of kpoints","",ngkpt    	
    	write(700+excnum,*) "#"," ","nocpk"," ","kx"," ","ky"," ","kz"," ","nc"," ","nv"," ","re_wf"," ","imag_wf"
    	

	do i=1,ngkpt*nc*nv
 	write(700+excnum,Format) nocp(stt(i,4)),kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),&
 			   	  stt(i,3),stt(i,2),real(ewf(i),8),aimag(ewf(i))
	end do
	close(700+excnum)
 endif
end subroutine excwf



subroutine excwfi(outputfolder,ngkpt,kpt,qpt,nc,nv,nocp,stt,excenergy,excnum,qptnum,ewf)
	integer :: i,erro
	character(len=70) :: outputfolder    !pasta saida
	character(len=200) :: file1
	CHARACTER(LEN=30) :: Format
	integer :: ngkpt,nc,nv,excnum,qptnum
	integer,dimension(ngkpt) :: nocp
	integer,dimension(ngkpt*nc*nv,4) :: stt
	double precision :: excenergy
	double precision,dimension(ngkpt,3) :: kpt
	double precision,dimension(4) :: qpt
	double complex,dimension(ngkpt*nc*nv) :: ewf

 if(Process==0) then
	WRITE (file1,'(a7,I0,a1,I0,a4)') ,'exc_wf_',excnum,"_",qptnum,'.dat'
	Format = "(I,3F15.4,2I,2E15.4)"
	
	OPEN(UNIT=700+excnum*qptnum, FILE=trim(outputfolder)//trim(file1),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening wf BSE output file"
    	
    	write(700+excnum*qptnum,*) "#"," ","excitonic momentum",real(qpt(2)),real(qpt(3)),real(qpt(4))
    	write(700+excnum*qptnum,*) "#"," ","excitonic state",excnum
    	write(700+excnum*qptnum,*) "#"," ","exciton energy:",excenergy
    	write(700+excnum*qptnum,*) "#"," ","number occupied states",nocp(1),"(not valid for metallic systems)"
    	write(700+excnum*qptnum,*) "#"," ","Number of conduction states","",nc
    	write(700+excnum*qptnum,*) "#"," ","Number of valence states","",nv
    	write(700+excnum*qptnum,*) "#"," ","Number of kpoints","",ngkpt    	
    	write(700+excnum*qptnum,*) "#"," ","nocpk"," ","kx"," ","ky"," ","kz"," ","nc"," ","nv"," ","re_wf"," ","imag_wf"

	do i=1,ngkpt*nc*nv
 	write(700+excnum*qptnum,Format) nocp(stt(i,4)),kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),&
 			   	  stt(i,3),stt(i,2),real(ewf(i)),aimag(ewf(i))
	end do
	close(700+excnum*qptnum)
 endif
end subroutine excwfi



subroutine exclft(sysdim,ngrid,rlat,fosc,enexc,lft) 
	integer,dimension(3) :: ngrid
	double precision :: lft,fosc,refr,enexc
	double precision,dimension(3,3) :: rlat
	double precision,parameter :: cfe= 7.297E-3 !constante estrutura fina
	double precision,parameter :: lspd= 299.792E16 !vel luz Ang/s^2
	double precision,parameter :: hbar= 6.582E-16 !const planck eV.s
	double precision,parameter :: cic= -(0.0904756)*10**(3)
	double precision,parameter:: pi=acos(-1.)	
	double precision :: aux1,aux2,aux3,vc
	character(len=5) :: sysdim

	select case (sysdim)
	
	 case("3D")
	 
	 	!for 3D systems- DOI: 10.1103/PhysRevLett.95.247402
		aux1 = (4.0*cfe)*(enexc*enexc*enexc)
		aux2 = fosc/dble(ngrid(1)*ngrid(2)*ngrid(3))
		aux3 = 3.0*(lspd*lspd)*(hbar*hbar*hbar)
		lft = (aux1*aux2)/aux3
		lft = (1.0/lft)
	 
	 case("2D")
	 
	 	 !for 2D systems- DOI:10.1021/nl503799t
	 	call vcell2D(rlat,vc)
		aux1 = dble(ngrid(1)*ngrid(2)*ngrid(3))*vc*hbar
		aux2 = (8.0*pi*cfe*enexc)*(fosc)	
		lft = aux1/aux2
	 
	 case("1D")
	 
	 	!for 1D systems- DOI: 10.1103/PhysRevLett.95.247402
		aux2 = (2.0*pi*cfe*enexc*enexc)*fosc
		aux1 = dble(ngrid(1)*ngrid(2)*ngrid(3))*rlat(1,1)*hbar*hbar*lspd
		lft = aux1/aux2
	 
	 case default

		if(Process==0) write(*,*) "Wrong System dimension"
		STOP

	end select 
end subroutine exclft



subroutine excwf2(outputfolder,ngkpt,kpt,nc,nv,nocp,stt,excenergy,excnum,ewf)
	integer :: i,erro
	character(len=70) :: outputfolder    !pasta saida
	character(len=200) :: file1
	CHARACTER(LEN=30) :: Format
	integer :: ngkpt,nc,nv,excnum
	integer,dimension(ngkpt) :: nocp
	integer,dimension(ngkpt*nc*nv,4) :: stt
	double precision :: excenergy
	double precision,dimension(ngkpt,3) :: kpt
	double complex,dimension(ngkpt*nc*nv) :: ewf

 if(Process==0) then
	WRITE (file1,'(a7,I0,a4)') ,'exc_wf_',excnum,'.dat'
	Format = "(I,3F15.4,2I,2E15.4)"
	
	OPEN(UNIT=700+excnum, FILE=trim(outputfolder)//trim(file1),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening wf BSE output file"
    	
    	write(700+excnum,*) "#"," ","excitonic state",excnum
    	write(700+excnum,*) "#"," ","exciton energy:",excenergy
    	write(700+excnum,*) "#"," ","number occupied states",nocp(1),"(not valid for metallic systems)"
    	write(700+excnum,*) "#"," ","Number of conduction states","",nc
    	write(700+excnum,*) "#"," ","Number of valence states","",nv
    	write(700+excnum,*) "#"," ","Number of kpoints","",ngkpt    	
    	write(700+excnum,*) "#"," ","nocpk"," ","kx"," ","ky"," ","kz"," ","nc"," ","nv"," ","re_wf"," ","imag_wf"
    	

	do i=1,ngkpt*nc*nv
 	write(700+excnum,Format) nocp(stt(i,4)),kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),&
 			   	  nocp(stt(i,4))+stt(i,3)-nv,nocp(stt(i,4))-nv+stt(i,2),real(ewf(i)),aimag(ewf(i))
	end do
	close(700+excnum)
 endif
end subroutine excwf2



subroutine excwfi2(outputfolder,ngkpt,kpt,qpt,nc,nv,nocp,stt,excenergy,excnum,qptnum,ewf)
	integer :: i,erro
	character(len=70) :: outputfolder    !pasta saida
	character(len=200) :: file1
	CHARACTER(LEN=30) :: Format
	integer :: ngkpt,nc,nv,excnum,qptnum
	integer,dimension(ngkpt) :: nocp
	integer,dimension(ngkpt*nc*nv,4) :: stt
	double precision :: excenergy
	double precision,dimension(ngkpt,3) :: kpt
	double precision,dimension(4) :: qpt
	double complex,dimension(ngkpt*nc*nv) :: ewf

 if(Process==0) then
	WRITE (file1,'(a7,I0,a1,I0,a4)') ,'exc_wf_',excnum,"_",qptnum,'.dat'
	Format = "(I,3F15.4,2I,2E15.4)"
	
	OPEN(UNIT=700+excnum*qptnum, FILE=trim(outputfolder)//trim(file1),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening wf BSE output file"
    	
    	write(700+excnum*qptnum,*) "#"," ","excitonic momentum",real(qpt(2)),real(qpt(3)),real(qpt(4))
    	write(700+excnum*qptnum,*) "#"," ","excitonic state",excnum
    	write(700+excnum*qptnum,*) "#"," ","exciton energy:",excenergy
    	write(700+excnum*qptnum,*) "#"," ","number occupied states",nocp(1),"(not valid for metallic systems)"
    	write(700+excnum*qptnum,*) "#"," ","Number of conduction states","",nc
    	write(700+excnum*qptnum,*) "#"," ","Number of valence states","",nv
    	write(700+excnum*qptnum,*) "#"," ","Number of kpoints","",ngkpt    	
    	write(700+excnum*qptnum,*) "#"," ","nocpk"," ","kx"," ","ky"," ","kz"," ","nc"," ","nv"," ","re_wf"," ","imag_wf"

	do i=1,ngkpt*nc*nv
 	write(700+excnum*qptnum,Format) nocp(stt(i,4)),kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),&
 			   	  nocp(stt(i,4))+stt(i,3)-nv,nocp(stt(i,4))-nv+stt(i,2),real(ewf(i)),aimag(ewf(i))
	end do
	close(700+excnum*qptnum)
 endif
end subroutine excwfi2



!INCLUDE "./subroutines/bse_subs_kpath.f90"
complex(8) function matrizelbsekq(coultype,tolr,w90basis,ediel,lc,ez,w,r0,ngrid,q,rlat,est1,ec1,ev1,vbc1 &
                       ,vbv1,kpt1,est2,ec2,ev2,vbc2,vbv2,kpt2) !funcao para calcular o elemento de matriz da matriz bse
	character(len=5) :: coultype
	integer,dimension(3) :: ngrid
	double precision,dimension(3,3) :: rlat
	double precision :: ez,w	
	integer :: w90basis
	double precision :: a,vcell1
	integer, dimension(4) :: est1,est2
	double precision :: ec1,ec2,ev1,ev2
	double precision,dimension(3) :: kpt1,kpt2
	double precision,dimension(4) :: q
	double precision,dimension(3) :: vq,v0
	double complex, dimension(w90basis) :: vbc1,vbc2,vbv1,vbv2
	double complex, dimension(w90basis) :: vbc,vbv,vbvkp
	double precision :: tolr
	integer :: ktol
	double precision :: modk,modq
	double precision,dimension(3) :: ediel
	double precision :: lc
	double complex :: vc,vv
	double complex :: vcv,vvc
	double precision :: vcoulk,vcoulq
	double precision :: r0

	v0 = 0.0
	vq = q(2:4)
	call modvec(kpt1,kpt2,modk)
	call modvecq(q,modq)

	select case (coultype)

	case("V2DK")

		vcoulk= v2dk(kpt1,kpt2,ediel,rlat,ngrid,lc,tolr)
		vcoulq= v2dk(vq,v0,ediel,rlat,ngrid,lc,tolr)

	case("V3D")

		vcoulk= vcoul(kpt1,kpt2,rlat,ngrid,tolr)
		vcoulq= vcoul(vq,v0,rlat,ngrid,tolr)

	case("V3DL")

		vcoulk= v3diel(kpt1,kpt2,ediel,rlat,ngrid,tolr)
		vcoulq= v3diel(vq,v0,ediel,rlat,ngrid,tolr)

	case("V2DT")

		vcoulk= v2dt(kpt1,kpt2,ngrid,rlat,tolr)
		vcoulq= v2dt(vq,v0,ngrid,rlat,tolr)

	case("V2DT2")

		vcoulk= v2dt2(kpt1,kpt2,ngrid,rlat,lc,tolr)
		vcoulq= v2dt2(vq,v0,ngrid,rlat,lc,tolr)
		
	case("V2DOH")

		vcoulk= v2dohono(kpt1,kpt2,ngrid,rlat,ediel,w,ez,tolr)
		vcoulq= v2dohono(vq,v0,ngrid,rlat,ediel,w,ez,tolr)
		
	case("V2DRK")

		vcoulk= v2drk(kpt1,kpt2,ngrid,rlat,ediel,lc,ez,w,r0,tolr)
		vcoulq= v2drk(vq,v0,ngrid,rlat,ediel,lc,ez,w,r0,tolr)
		
	case("V1DT")
	
		vcoulk=	v1dt(kpt1,kpt2,ngrid,rlat,tolr,lc)	
		vcoulq= v1dt(vq,v0,ngrid,rlat,tolr,lc)		

	case("V0DT")

		vcoulk= v0dt(kpt1,kpt2,ngrid,rlat,tolr)
		vcoulq= v0dt(vq,v0,ngrid,rlat,tolr)

	case default

		if(Process==0) write(*,*) "Wrong Coulomb Potential"
		STOP

	end select

if (modq .eq. 0.) then
	if (est1(1) .eq. est2(1)) then
		matrizelbsekq= (ec1-ev1) + vcoulk
	else
		call vecconjg(vbc1,w90basis,vbc)
		call vecconjg(vbv1,w90basis,vbv)
		call prodintsq(vbc,vbc2,w90basis,vc)
		call prodintsq(vbv,vbv2,w90basis,vv)
		matrizelbsekq= vcoulk*vc*vv
	end if
else 
	if (est1(1) .eq. est2(1)) then
		call vecconjg(vbc1,w90basis,vbc)
		call vecconjg(vbv1,w90basis,vbv)
		call vecconjg(vbv2,w90basis,vbvkp)
		call prodintsq(vbc,vbc2,w90basis,vc)
		call prodintsq(vbv,vbv2,w90basis,vv)
		call prodintsq(vbc,vbv1,w90basis,vcv)
		call prodintsq(vbvkp,vbc2,w90basis,vvc)
		matrizelbsekq= (ec1-ev1) + vcoulk &
			     - vcoulq*vcv*vvc
	else
		call vecconjg(vbc1,w90basis,vbc)
		call vecconjg(vbv1,w90basis,vbv)
		call vecconjg(vbv2,w90basis,vbvkp)
		call prodintsq(vbc,vbc2,w90basis,vc)
		call prodintsq(vbv,vbv2,w90basis,vv)
		call prodintsq(vbc,vbv1,w90basis,vcv)
		call prodintsq(vbvkp,vbc2,w90basis,vvc)
		matrizelbsekq= vcoulk*vc*vv&
				- vcoulq*vcv*vvc
	end if
end if
end function matrizelbsekq








!INCLUDE "./subroutines/bse_subs_temp.f90"

complex(8) function matrizelbsetemp(coultype,tolr,w90basis,ediel,lc,ez,w,r0,ngrid,rlat,est1,ec1,ev1,vbc1 &
         ,vbv1,kpt1,est2,ec2,ev2,vbc2,vbv2,kpt2,temp) !funcao para calcular o elemento de matriz da matriz bse
	character(len=5) :: coultype
	integer,dimension(3) :: ngrid
	integer :: w90basis
	double precision :: a,vcell1
	double precision :: ez,w
	integer, dimension(4) :: est1,est2
	double precision :: ec1,ec2,ev1,ev2
	double precision,dimension(3) :: kpt1,kpt2
	double complex, dimension(w90basis) :: vbc1,vbc2,vbv1,vbv2
	double complex, dimension(w90basis) :: vbc,vbv
	double precision,dimension(3,3) :: rlat
	double precision :: tolr
	integer :: ktol
	double precision,parameter:: pi=acos(-1.)
	double precision :: auxi
	double precision :: modk
	double precision,dimension(3) :: ediel
	double precision :: lc
	double complex :: vc,vv
	double precision :: vcoul1
	double precision :: r0
	double precision :: temp      !, fermidisteh

	select case (coultype)

	case("V2DK")

		vcoul1= v2dk(kpt1,kpt2,ediel,rlat,ngrid,lc,tolr)

	case("V3D")

		vcoul1= vcoul(kpt1,kpt2,rlat,ngrid,tolr)

	case("V3DL")

		vcoul1= v3diel(kpt1,kpt2,ediel,rlat,ngrid,tolr)

	case("V2DT")

		vcoul1= v2dt(kpt1,kpt2,ngrid,rlat,tolr)

	case("V2DT2")

		vcoul1= v2dt2(kpt1,kpt2,ngrid,rlat,lc,tolr)
		
	case("V2DOH")

		vcoul1= v2dohono(kpt1,kpt2,ngrid,rlat,ediel,w,ez,tolr)
		
	case("V2DRK")

		vcoul1= v2drk(kpt1,kpt2,ngrid,rlat,ediel,lc,ez,w,r0,tolr)
		
	case("V1DT")
	
		vcoul1=	v1dt(kpt1,kpt2,ngrid,rlat,tolr,lc)	
		
		
	case("V0DT")

		vcoul1= v0dt(kpt1,kpt2,ngrid,rlat,tolr)

	case default

		if(Process==0) write(*,*) "Wrong Coulomb Potential"
		STOP

	end select

	if (est1(1) .eq. est2(1)) then
		matrizelbsetemp= (ec1-ev1) + vcoul1*fermidisteh(ec1,ev1,temp)
	else
		call vecconjg(vbc1,w90basis,vbc)
		call vecconjg(vbv1,w90basis,vbv)
		call prodintsq(vbc,vbc2,w90basis,vc)
		call prodintsq(vbv,vbv2,w90basis,vv)
		matrizelbsetemp=  vcoul1*vc*vv*fermidisteh(ec1,ev1,temp)
	end if
end function matrizelbsetemp



complex(8) function matrizelbsekqtemp(coultype,tolr,w90basis,ediel,lc,ez,w,r0,ngrid,q,rlat,est1,ec1,ev1,vbc1 &
                       ,vbv1,kpt1,est2,ec2,ev2,vbc2,vbv2,kpt2,temp) !funcao para calcular o elemento de matriz da matriz bse
	character(len=5) :: coultype
	integer,dimension(3) :: ngrid
	double precision,dimension(3,3) :: rlat
	double precision :: ez,w	
	integer :: w90basis
	double precision :: a,vcell1
	integer, dimension(4) :: est1,est2
	double precision :: ec1,ec2,ev1,ev2
	double precision,dimension(3) :: kpt1,kpt2
	double precision,dimension(4) :: q
	double precision,dimension(3) :: vq,v0
	double complex, dimension(w90basis) :: vbc1,vbc2,vbv1,vbv2
	double complex, dimension(w90basis) :: vbc,vbv,vbvkp
	double precision :: tolr
	integer :: ktol
	double precision :: modk,modq
	double precision,dimension(3) :: ediel
	double precision :: lc
	double complex :: vc,vv
	double complex :: vcv,vvc
	double precision :: vcoulk,vcoulq
	double precision :: r0
	double precision :: temp          !, fermidisteh

	v0 = 0.0
	vq = q(2:4)
	call modvec(kpt1,kpt2,modk)
	call modvecq(q,modq)

	select case (coultype)

	case("V2DK")

		vcoulk= v2dk(kpt1,kpt2,ediel,rlat,ngrid,lc,tolr)
		vcoulq= v2dk(vq,v0,ediel,rlat,ngrid,lc,tolr)

	case("V3D")

		vcoulk= vcoul(kpt1,kpt2,rlat,ngrid,tolr)
		vcoulq= vcoul(vq,v0,rlat,ngrid,tolr)

	case("V3DL")

		vcoulk= v3diel(kpt1,kpt2,ediel,rlat,ngrid,tolr)
		vcoulq= v3diel(vq,v0,ediel,rlat,ngrid,tolr)

	case("V2DT")

		vcoulk= v2dt(kpt1,kpt2,ngrid,rlat,tolr)
		vcoulq= v2dt(vq,v0,ngrid,rlat,tolr)

	case("V2DT2")

		vcoulk= v2dt2(kpt1,kpt2,ngrid,rlat,lc,tolr)
		vcoulq= v2dt2(vq,v0,ngrid,rlat,lc,tolr)
		
	case("V2DOH")

		vcoulk= v2dohono(kpt1,kpt2,ngrid,rlat,ediel,w,ez,tolr)
		vcoulq= v2dohono(vq,v0,ngrid,rlat,ediel,w,ez,tolr)
		
	case("V2DRK")

		vcoulk= v2drk(kpt1,kpt2,ngrid,rlat,ediel,lc,ez,w,r0,tolr)
		vcoulq= v2drk(vq,v0,ngrid,rlat,ediel,lc,ez,w,r0,tolr)
		
	case("V1DT")
	
		vcoulk=	v1dt(kpt1,kpt2,ngrid,rlat,tolr,lc)	
		vcoulq= v1dt(vq,v0,ngrid,rlat,tolr,lc)
		
	case("V0DT")

		vcoulk= v0dt(kpt1,kpt2,ngrid,rlat,tolr)
		vcoulq= v0dt(vq,v0,ngrid,rlat,tolr)

	case default

		if(Process==0) write(*,*) "Wrong Coulomb Potential"
		STOP

	end select

if (modq .eq. 0.) then
	if (est1(1) .eq. est2(1)) then
		matrizelbsekqtemp= (ec1-ev1) + vcoulk*fermidisteh(ec1,ev1,temp)
	else
		call vecconjg(vbc1,w90basis,vbc)
		call vecconjg(vbv1,w90basis,vbv)
		call prodintsq(vbc,vbc2,w90basis,vc)
		call prodintsq(vbv,vbv2,w90basis,vv)
		matrizelbsekqtemp= vcoulk*vc*vv*fermidisteh(ec1,ev1,temp)
	end if
else 
	if (est1(1) .eq. est2(1)) then
		call vecconjg(vbc1,w90basis,vbc)
		call vecconjg(vbv1,w90basis,vbv)
		call vecconjg(vbv2,w90basis,vbvkp)
		call prodintsq(vbc,vbc2,w90basis,vc)
		call prodintsq(vbv,vbv2,w90basis,vv)
		call prodintsq(vbc,vbv1,w90basis,vcv)
		call prodintsq(vbvkp,vbc2,w90basis,vvc)
		matrizelbsekqtemp= (ec1-ev1) + (vcoulk - vcoulq*vcv*vvc)*fermidisteh(ec1,ev1,temp)
	else
		call vecconjg(vbc1,w90basis,vbc)
		call vecconjg(vbv1,w90basis,vbv)
		call vecconjg(vbv2,w90basis,vbvkp)
		call prodintsq(vbc,vbc2,w90basis,vc)
		call prodintsq(vbv,vbv2,w90basis,vv)
		call prodintsq(vbc,vbv1,w90basis,vcv)
		call prodintsq(vbvkp,vbc2,w90basis,vvc)
		matrizelbsekqtemp= (vcoulk*vc*vv- vcoulq*vcv*vvc)*fermidisteh(ec1,ev1,temp)
	end if
end if
end function matrizelbsekqtemp



subroutine dielbseptemp(nthread,dimse,dimse1,excitonvec,hopt1,hopt2,fdeh,activity)
	integer                         :: dimse,dimse1,nthread
	double precision,dimension(:)   :: activity           !, exciton
	double complex                  :: actaux,actaux2
	double complex,dimension(:)     :: hopt1,hopt2
	double precision,dimension(:)   :: fdeh
	double complex,dimension(:,:)   :: excitonvec
	integer                         :: i,j             !,k,no

!	call OMP_SET_NUM_THREADS(nthread)

	activity=0.0
	!actaux=0.0

	! $OMP DO PRIVATE(actaux)
	!$OMP PARALLEL DO PRIVATE(actaux,actaux2)
	do i=1,dimse1
		actaux=0.0
		actaux2=0.0
		do j=1,dimse
			actaux=actaux+(excitonvec(j,i)*hopt1(j))
			actaux2=actaux2+(excitonvec(j,i)*hopt2(j)*fdeh(j))
		end do
		activity(i)=actaux*conjg(actaux2)
	end do
	!$OMP END PARALLEL DO
end subroutine dielbseptemp






















!INCLUDE "./subprograms/bse_diel-tool.f90"
!ifort bse_diel-tool.f90 -o spec_bse_diel.x -qopenmp -mkl
subroutine bsedielraw(nthreads,outputfolder,renorm,params,ngrid,nc,nv,ebse0,ebsef,numbse,sme)
	integer :: i,j,erro
	double precision,allocatable,dimension(:) :: ebse
	double precision,allocatable,dimension(:) :: exxf,exyf,exzf 
	double precision,allocatable,dimension(:) :: eyyf,eyzf,ezzf 
	double precision :: flag
	double precision :: vc
	double precision :: elux
	double precision :: rxx,rxy,rxz
	double precision :: ryy,ryz,rzz
	double precision :: ixx,ixy,ixz
	double precision :: iyy,iyz,izz
	character(len=2) :: aread
	double precision :: spinf
	double precision,allocatable,dimension(:,:) :: dielfxx,dielfyy,dielfzz
	double precision,allocatable,dimension(:,:) :: dielfxy,dielfxz,dielfyz	
	integer :: nthreads
	integer,dimension(3) :: ngrid
	integer :: nc,nv
	double precision :: edos0,edosf,numdos
	double precision :: ebse0,ebsef,numbse
	double precision :: sme,exc,rk
	double precision,dimension(3) :: mshift
	double precision :: ktol
	character(len=70) :: params   !parametros TB
	character(len=70) :: orbw     !peso orbitais lcount
	character(len=70) :: kpaths    !kpath
	character(len=70) :: kpathsbse    !kpath
	character(len=70) :: outputfolder    !pasta saida
	character(len=70) :: calcparms
	character(len=70) :: meshtype
	character(len=5) :: coultype
	double precision,dimension(3) :: ediel
	logical :: renorm

 !if(process==0) print *,'MODULE_BSE bsedielraw'
    call MSG('MODULE_BSE: start bsedielraw')
 
	call hamiltonian_input_read(200,params)

	if ( systype .eq. "NP" ) then
	spinf = 2.0
	else
	spinf = 1.0
	end if

	dimbse=ngrid(1)*ngrid(2)*ngrid(3)*nc*nv

 if(Process==0) then
	OPEN(UNIT=100, FILE=trim(outputfolder)//"/bse_opt_diel.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse_opt_diel input file"
	OPEN(UNIT=300, FILE=trim(outputfolder)//"/bse_diel_xx.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel_xx output file"
	OPEN(UNIT=301, FILE=trim(outputfolder)//"/bse_diel_xy.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel_xy output file"
	OPEN(UNIT=302, FILE=trim(outputfolder)//"/bse_diel_xz.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel_xz output file"
	OPEN(UNIT=303, FILE=trim(outputfolder)//"/bse_diel_yy.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel_yy output file"
	OPEN(UNIT=304, FILE=trim(outputfolder)//"/bse_diel_yz.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel_yz output file"
	OPEN(UNIT=305, FILE=trim(outputfolder)//"/bse_diel_zz.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel_zz output file"
 endif

	call vcell3D(rlat,vc)

	allocate(ebse(dimbse))
	allocate(exxf(dimbse),exyf(dimbse),exzf(dimbse),eyyf(dimbse))
	allocate(eyzf(dimbse),ezzf(dimbse))

	if(Process==0) read(100,*) aread
	call P_sendA(aread)

 if(Process==0) then
	do j=1,dimbse
	 read(100,*) ebse(j),exxf(j),eyyf(j),ezzf(j),exyf(j),exzf(j),eyzf(j)
	end do
 endif
    call P_sendR1(ebse,dimbse)
    call P_sendR1(exxf,dimbse)
    call P_sendR1(eyyf,dimbse)
    call P_sendR1(ezzf,dimbse)
    call P_sendR1(exyf,dimbse)
    call P_sendR1(exzf,dimbse)
    call P_sendR1(eyzf,dimbse)
 
	allocate(dielfxx(int(numbse),3),dielfyy(int(numbse),3),dielfzz(int(numbse),3))
	allocate(dielfxy(int(numbse),3),dielfxz(int(numbse),3),dielfyz(int(numbse),3))

 if(Process==0) then
	write(300,*) "#","  ","energy","  ","real","  ","imag"
	write(301,*) "#","  ","energy","  ","real","  ","imag"
	write(302,*) "#","  ","energy","  ","real","  ","imag"
	write(303,*) "#","  ","energy","  ","real","  ","imag"
	write(304,*) "#","  ","energy","  ","real","  ","imag"
	write(305,*) "#","  ","energy","  ","real","  ","imag"
 endif

	!$omp do 
	!ordered
	do j=1,int(numbse)

		elux= ebse0 + dble(((ebsef-ebse0)*(j-1))/(numbse-1.))

		call optdiel(1,vc,dimbse,ngrid,elux,ebse,exxf,sme,rxx,ixx)
		call optdiel(1,vc,dimbse,ngrid,elux,ebse,eyyf,sme,ryy,iyy)
		call optdiel(1,vc,dimbse,ngrid,elux,ebse,ezzf,sme,rzz,izz)


		call optdiel(0,vc,dimbse,ngrid,elux,ebse,exyf,sme,rxy,ixy)
		call optdiel(0,vc,dimbse,ngrid,elux,ebse,exzf,sme,rxz,ixz)
		call optdiel(0,vc,dimbse,ngrid,elux,ebse,eyzf,sme,ryz,iyz)
		
		dielfxx(j,1) = elux
		dielfxx(j,2) = rxx
		dielfxx(j,3) = ixx
		
		dielfyy(j,1) = elux
		dielfyy(j,2) = ryy
		dielfyy(j,3) = iyy
		
		dielfzz(j,1) = elux
		dielfzz(j,2) = rzz
		dielfzz(j,3) = izz
		
		dielfxy(j,1) = elux
		dielfxy(j,2) = rxy
		dielfxy(j,3) = ixy
		
		dielfxz(j,1) = elux
		dielfxz(j,2) = rxz
		dielfxz(j,3) = ixz
		
		dielfyz(j,1) = elux
		dielfyz(j,2) = ryz
		dielfyz(j,3) = iyz				
		
																			

		! $omp ordered
		!write(300,*) elux,real(spinf*rxx),real(spinf*ixx)
		!write(301,*) elux,real(spinf*rxy),real(spinf*ixy)
		!write(302,*) elux,real(spinf*rxz),real(spinf*ixz)
		!write(303,*) elux,real(spinf*ryy),real(spinf*iyy)
		!write(304,*) elux,real(spinf*ryz),real(spinf*iyz)
		!write(305,*) elux,real(spinf*rzz),real(spinf*izz)
	
		!call flush(300)
		! $omp end ordered


	end do
	!$omp end  do
	
	if (renorm) then
	call imagrenorm(dimbse,int(numbse),ebse,dielfxx)
	call imagrenorm(dimbse,int(numbse),ebse,dielfyy)
	call imagrenorm(dimbse,int(numbse),ebse,dielfzz)
	call imagrenorm(dimbse,int(numbse),ebse,dielfxy)
	call imagrenorm(dimbse,int(numbse),ebse,dielfxz)
	call imagrenorm(dimbse,int(numbse),ebse,dielfyz)
	
	else
		continue
	end if					
	
 if(Process==0) then
	do i=1,int(numbse)
		write(300,*) dielfxx(i,1),real(spinf*dielfxx(i,2)),real(spinf*dielfxx(i,3))
		write(301,*) dielfxy(i,1),real(spinf*dielfxy(i,2)),real(spinf*dielfxy(i,3))
		write(302,*) dielfxz(i,1),real(spinf*dielfxz(i,2)),real(spinf*dielfxz(i,3))
		write(303,*) dielfyy(i,1),real(spinf*dielfyy(i,2)),real(spinf*dielfyy(i,3))
		write(304,*) dielfyz(i,1),real(spinf*dielfyz(i,2)),real(spinf*dielfyz(i,3))
		write(305,*) dielfzz(i,1),real(spinf*dielfzz(i,2)),real(spinf*dielfzz(i,3))
		!write(301,*) elux,real(spinf*rxy),real(spinf*ixy)
		!write(302,*) elux,real(spinf*rxz),real(spinf*ixz)
		!write(303,*) elux,real(spinf*ryy),real(spinf*iyy)
		!write(304,*) elux,real(spinf*ryz),real(spinf*iyz)
		!write(305,*) elux,real(spinf*rzz),real(spinf*izz)
	end do
 endif
 
	deallocate(ebse)
	deallocate(exxf,exyf,exzf,eyyf)
	deallocate(eyzf,ezzf)
	deallocate(rvec,hopmatrices)
	deallocate(ihopmatrices,ffactor)
	deallocate(dielfxx,dielfyy,dielfzz)
	deallocate(dielfxy,dielfxz,dielfyz)	

 if(Process==0) then
	close(100)
	close(200)
	close(300)
	close(301)
	close(302)
	close(303)
	close(304)
	close(305)
 endif
end subroutine bsedielraw




!INCLUDE "./subprograms/bse_diel-tool-pol.f90"
!ifort bse_diel-tool.f90 -o spec_bse_diel.x -qopenmp -mkl
subroutine bsedielrawpol(nthreads,outputfolder,renorm,params,ngrid,nc,nv,ebse0,ebsef,numbse,sme)
	integer :: i,j,erro
	integer :: dimbse     
	double precision,allocatable,dimension(:) :: esp
	double precision,allocatable,dimension(:) :: exxf,eyyf,ezzf 
	double precision,allocatable,dimension(:) :: espf,esmf 
	double precision :: flag
	double precision :: vc
	double precision :: elux
	double precision :: rxx,ryy,rzz
	double precision :: rsp,rsm
	double precision :: ixx,iyy,izz
	double precision :: isp,ism
	double precision,allocatable,dimension(:,:) :: dielfxx,dielfyy,dielfzz
	double precision,allocatable,dimension(:,:) :: dielfsp,dielfsm	
	character(len=2) :: aread
	double precision :: spinf
	integer :: nthreads
	integer,dimension(3) :: ngrid
	integer :: nc,nv
	double precision :: edos0,edosf,numdos
	double precision :: ebse0,ebsef,numbse
	double precision :: sme,exc,rk
	double precision,dimension(3) :: mshift
	double precision :: ktol
	character(len=70) :: params   !parametros TB
	character(len=70) :: orbw     !peso orbitais lcount
	character(len=70) :: kpaths    !kpath
	character(len=70) :: kpathsbse    !kpath
	character(len=70) :: outputfolder    !pasta saida
	character(len=70) :: calcparms
	character(len=70) :: meshtype
	character(len=5) :: coultype
	double precision,dimension(3) :: ediel
	logical :: renorm

   if(Process==0) print *,'MODULE_BSE   bsedielrawpol'

	call hamiltonian_input_read(200,params)
	if ( systype .eq. "NP" ) then
	spinf = 2.0
	else
	spinf = 1.0
	end if

	dimbse=ngrid(1)*ngrid(2)*ngrid(3)*nc*nv

 if(Process==0) then
	OPEN(UNIT=100, FILE=trim(outputfolder)//"/bse_opt_diel-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse_opt_diel-pol input file"
!	call OMP_SET_NUM_THREADS(nthreads)
	OPEN(UNIT=300, FILE=trim(outputfolder)//"/bse_diel-pol_x.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel-pol_x output file"
	OPEN(UNIT=301, FILE=trim(outputfolder)//"/bse_diel-pol_y.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel-pol_y output file"
	OPEN(UNIT=302, FILE=trim(outputfolder)//"/bse_diel-pol_z.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel-pol_z output file"
	OPEN(UNIT=303, FILE=trim(outputfolder)//"/bse_diel-pol_sp.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel-pol_sp output file"
	OPEN(UNIT=304, FILE=trim(outputfolder)//"/bse_diel-pol_sm.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel-pol_sm output file"
 endif
 
	call vcell3D(rlat,vc)

	allocate(esp(dimbse))
	allocate(exxf(dimbse),eyyf(dimbse),ezzf(dimbse))
	allocate(espf(dimbse),esmf(dimbse))
	if(Process==0) read(100,*) aread
	call P_sendA(aread)

 if(Process==0) then
	do j=1,dimbse
	read(100,*) esp(j),exxf(j),eyyf(j),ezzf(j),espf(j),esmf(j)
	end do
 endif
    call P_sendR1(esp,dimbse)
    call P_sendR1(exxf,dimbse)
    call P_sendR1(eyyf,dimbse)
    call P_sendR1(ezzf,dimbse)
    call P_sendR1(espf,dimbse)
    call P_sendR1(esmf,dimbse)

	allocate(dielfxx(int(numbse),3),dielfyy(int(numbse),3),dielfzz(int(numbse),3))
	allocate(dielfsp(int(numbse),3),dielfsm(int(numbse),3))

 if(Process==0) then
	write(300,*) "#","  ","energy","  ","real","  ","imag"
	write(301,*) "#","  ","energy","  ","real","  ","imag"
	write(302,*) "#","  ","energy","  ","real","  ","imag"
	write(303,*) "#","  ","energy","  ","real","  ","imag"
	write(304,*) "#","  ","energy","  ","real","  ","imag"
 endif

	!$omp do 
	!ordered
	do j=1,int(numbse)
		elux= ebse0 + dble(((ebsef-ebse0)*(j-1))/(numbse-1.))
		call optdiel(1,vc,dimbse,ngrid,elux,esp,exxf,sme,rxx,ixx)
		call optdiel(1,vc,dimbse,ngrid,elux,esp,eyyf,sme,ryy,iyy)
		call optdiel(1,vc,dimbse,ngrid,elux,esp,ezzf,sme,rzz,izz)
		call optdiel(1,vc,dimbse,ngrid,elux,esp,espf,sme,rsp,isp)
		call optdiel(1,vc,dimbse,ngrid,elux,esp,esmf,sme,rsm,ism)
		dielfxx(j,1) = elux
		dielfxx(j,2) = rxx
		dielfxx(j,3) = ixx
		dielfyy(j,1) = elux
		dielfyy(j,2) = ryy
		dielfyy(j,3) = iyy
		dielfzz(j,1) = elux
		dielfzz(j,2) = rzz
		dielfzz(j,3) = izz
		dielfsp(j,1) = elux
		dielfsp(j,2) = rsp
		dielfsp(j,3) = isp
		dielfsm(j,1) = elux
		dielfsm(j,2) = rsm
		dielfsm(j,3) = ism
	end do
	!$omp end  do
	
	if (renorm) then
	call imagrenorm(dimbse,int(numbse),esp,dielfxx)
	call imagrenorm(dimbse,int(numbse),esp,dielfyy)
	call imagrenorm(dimbse,int(numbse),esp,dielfzz)
	call imagrenorm(dimbse,int(numbse),esp,dielfsp)
	call imagrenorm(dimbse,int(numbse),esp,dielfsm)
	end if		
	
 if(Process==0) then
	do i=1,int(numbse)
		write(300,*) dielfxx(i,1),real(spinf*dielfxx(i,2)),real(spinf*dielfxx(i,3))
		write(301,*) dielfyy(i,1),real(spinf*dielfyy(i,2)),real(spinf*dielfyy(i,3))
		write(302,*) dielfzz(i,1),real(spinf*dielfzz(i,2)),real(spinf*dielfzz(i,3))
		write(303,*) dielfsp(i,1),real(spinf*dielfsp(i,2)),real(spinf*dielfsp(i,3))
		write(304,*) dielfsm(i,1),real(spinf*dielfsm(i,2)),real(spinf*dielfsm(i,3))
	end do
 endif

	deallocate(esp)
	deallocate(exxf,eyyf,ezzf)
	deallocate(espf,esmf)
	deallocate(rvec,hopmatrices)
	deallocate(ihopmatrices,ffactor)
	deallocate(dielfxx,dielfyy,dielfzz)
	deallocate(dielfsp,dielfsm)

 if(Process==0) then
	close(100)
	close(200)
	close(300)
	close(301)
	close(302)
	close(303)
	close(304)
 endif
end subroutine bsedielrawpol




!INCLUDE "./subprograms/bse_kpath-tool.f90"
!gfortran -mcmodel=large bse_kpath-tool.f90 -o bse_kpath-tool.x -llapack95 -lopenblas -fopenmp 
subroutine bsebnds(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,ez,w1,r0,lc,rk,meshtype,bsewf,excwf0,excwff)
!	use omp_lib
!	use module_hamiltonian   !_input_variables
!	implicit none
	double precision,parameter:: pi=acos(-1.)
	integer :: dimbse
	double precision,allocatable,dimension(:,:) :: exk !(dimbse,nqpts*npath)
	integer :: ncaux,nvaux
	integer,allocatable,dimension(:,:) :: stt !(ngrid*ngrid*nc*nv,4)
	double precision,allocatable,dimension(:,:) :: kpt,qpt !pontos k do grid (ngrid*ngrid,2)
	double precision,dimension(4) :: q
	double precision,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
	double complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores
	double precision,allocatable,dimension(:,:) :: energy !variavel que guarda os autovalores para cada ponto k
	double complex,allocatable,dimension(:,:,:) :: vector !variavel que guarda os autovetores para cada ponto k
	double precision,allocatable,dimension(:,:) :: energyq !variavel que guarda os autovalores para cada ponto k
	double complex,allocatable,dimension(:,:,:) :: vectorq !variavel que guarda os autovetores para cada ponto k
	double precision,allocatable,dimension(:,:) :: qauxv !(nqpts*npath,4) 
	integer:: counter,c,v,i,j,f,l,h,erro,i2,k
	double precision :: a,r0,ed
	integer :: ngkpt
	double precision:: t0,tf
	integer,dimension(8) :: values,values2
	integer :: nks
	integer :: nkpts
	double precision,allocatable,dimension(:,:) :: ks
	integer,allocatable,dimension(:) :: nocpk,nocpq
	INTEGER   ::       ifail
	double precision,parameter :: ABSTOL=1.0e-6
	INTEGER          INFO
	double precision,allocatable,dimension(:) :: W,RWORK
        INTEGER ::          LWMAX
   	INTEGER ::         LWORK
	INTEGER ::         LIWORK, LRWORK
	INTEGER,allocatable,dimension(:) :: IWORK
        double complex,allocatable,dimension (:) :: WORK
	integer :: nthreads
	integer,dimension(3) :: ngrid
	integer :: nc,nv
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
	double precision :: ez,w1,lc
	logical :: bsewf
	integer :: excwf0,excwff
 if(Process==0) print *,'******* bsebnds'
 
	OPEN(UNIT=500, FILE= kpathsbse,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse-kpath input file"
	!OPEN(UNIT=201, FILE= diein,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada ambiente dieletrico"

	!OUTPUT : criando arquivos de saida
	OPEN(UNIT=300, FILE=trim(outputfolder)//"/log_bse_kpath.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening log_bse_kpath output file"
	OPEN(UNIT=400, FILE=trim(outputfolder)//"/bands_bse.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bands_bse output file"

	call cpu_time(t0)
	call date_and_time(VALUES=values)
!	call OMP_SET_NUM_THREADS(nthreads)

	call hamiltonian_input_read(200,params)
	!ediel(2) = edielh

 if(Process==0) then
	read(500,*) nks
	read(500,*) nkpts
 endif
    call P_sendI(nks)
    call P_sendI(nkpts)

	allocate(ks(nks,3))

 if(Process==0) then
	do i=1,nks
		read(500,*) ks(i,1),ks(i,2),ks(i,3)
	end do
 endif
    call P_sendR2(ks,nks,3)

	ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
	dimbse = ngkpt*nc*nv

	!call alat(systype,rlat,a)
	!lc = rlat3(3)

	!if (systype .eq. "2D") then
	!	ed = (ediel(1)+ediel(3))/2.
	!	r0= ((ediel(2)-1.0)*lc)/(ediel(1)+ediel(3))
	!else
	!	r0 = 1.0
	!	ed = ediel(2)
		
	!end if

	!definindo kpath
	!allocate(qauxv(nkpts*(nks-1),4))
	allocate(qauxv((nks/2)*nkpts,4))	

	call kpathbse(outputfolder,rlat(1,:),rlat(2,:),rlat(3,:),nks,ks,nkpts,qauxv)

	!Informaç��es para o arquivo de log do calculo

 if(Process==0) then
	write(300,*) 'threads:', nthreads
	write(300,*)
	write(300,*) 'grid:',ngrid(1),ngrid(2),ngrid(3)
	write(300,*)
	write(300,*) 'ktol-coulomb:', ktol
	write(300,*)
	write(300,*) 'kmesh shift:',mshift(1),mshift(2),mshift(3)
	write(300,*)
	write(300,*) 'conduction bands, nc:',nc,'  ','valence bands, nv:',nv
	write(300,*)
	write(300,*) 'number of kpoints in the path:','   ',(nks/2)*nkpts
	write(300,*)
	write(300,*)
	write(300,*) 'begin','  ','day',values(3),'',values(5),'hours',values(6),'min',values(7),'seg'
	write(300,*) 
	call flush(300)
 endif

	allocate(kpt(ngkpt,3))
	
	!shift= 0.0
	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpt)
	!call gridgenmhp(ngrid,rlat,kpt)

	allocate(eaux(w90basis),vaux(w90basis,w90basis))
	allocate(energy(ngkpt,nc+nv),vector(ngkpt,nc+nv,w90basis))
	allocate(nocpk(ngkpt))
	allocate(nocpq(ngkpt))

	! $omp parallel do
	do i=1,ngkpt
		call eigsys(nthreads,scs,exc,nocpk(i),ffactor,kpt(i,1),kpt(i,2),kpt(i,3),w90basis,nvec,&
			    rlat,rvec,hopmatrices,&
		             ihopmatrices,efermi,eaux,vaux)
			do j=1,nc+nv
				energy(i,j)= eaux(nocpk(i)-nv+j)
			end do
			do l=1,nc+nv
				do h=1,w90basis
					vector(i,l,h)=vaux(nocpk(i)-nv+l,h)
				end do
			end do
	end do
	! $omp end parallel do

	deallocate(eaux,vaux)

	!definindo os numeros quanticos dos estados
    !write(*,*) "quantum numbers for exciton basis set finished"
	if(Process==0) write(300,*) "quantum numbers for exciton basis set finished"
	if(Process==0) call flush(300)	
	counter=counter-1 !numero total de estados para equação bse
	allocate(qpt(ngkpt,3))
	allocate(exk(dimbse,nkpts*(nks-1)))

	do i=1,(nks/2)*nkpts
		!definindo os pontos q

		q(1)= qauxv(i,1)
		q(2)= qauxv(i,2) 
		q(3)= qauxv(i,3)
		q(4)= qauxv(i,4)

		!gerando o grid k+q
		
		call monhkhorst_packq(q,ngrid(1),ngrid(2),ngrid(3),mshift,&
				      rlat(1,:),rlat(2,:),rlat(3,:),qpt)


		allocate(eaux(w90basis),vaux(w90basis,w90basis))

		allocate(energyq(ngkpt,nc+nv),vectorq(ngkpt,nc+nv,w90basis))

	allocate (stt(ngkpt*nc*nv,4))


	! $omp parallel do
		do i2=1,ngkpt

            !write(*,*) i2

			call eigsys(nthreads,scs,exc,nocpq(i2),ffactor,qpt(i2,1),qpt(i2,2),qpt(i2,3),w90basis,nvec,&
			             rlat,rvec,hopmatrices,ihopmatrices,efermi,eaux,vaux)

				do j=1,nc+nv
					energyq(i2,j)= eaux(nocpq(i2)-nv+j) 
				end do

				do l=1,nc+nv
					do h=1,w90basis
						vectorq(i2,l,h)=vaux(nocpq(i2)-nv+l,h)
					end do
				end do
		end do
	! $omp end parallel do

		call quantumnumbers2(w90basis,ngkpt,nc,nv,nocpk,nocpq,stt)

		deallocate(eaux,vaux)


		allocate(hbse(dimbse,dimbse))
		allocate(W(dimbse))

		hbse=0.

	!$omp parallel do 
        !collapse(2)

		do i2=1,dimbse
			do j=i2,dimbse
hbse(i2,j)= matrizelbsekq(coultype,ktol,w90basis,ediel,lc,ez,w1,r0,ngrid,q,rlat,stt(i2,:),energyq(stt(i2,4),stt(i2,3))&
          ,energy(stt(i2,4),stt(i2,2)),vectorq(stt(i2,4)&
          ,stt(i2,3),:) ,vector(stt(i2,4),stt(i2,2),:),kpt(stt(i2,4),:),stt(j,:)&
          ,energyq(stt(j,4),stt(j,3)),energy(stt(j,4),stt(j,2))&
          ,vectorq(stt(j,4),stt(j,3),:),vector(stt(j,4),stt(j,2),:),kpt(stt(j,4),:))
			end do
		end do

	!$omp end parallel do

    		!call LA_HEEVR( hbse, W, JOBZ='N', UPLO='U', ABSTOL=ABSTOL, INFO=INFO ) 



      	!LWORK = 2*dimbse+dimbse**2
      	!LIWORK = 3 + 5*dimbse
      	!LRWORK = 1 + 5*dimbse + 2*dimbse**2

	!allocate(WORK(LWORK))
	!allocate (IWORK(LIWORK))
	!allocate (RWORK(LRWORK))

 

      	!CALL ZHEEVD( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
      	!LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      	!LRWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      	!LIWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

      	!CALL ZHEEVD( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK,&
	!	      RWORK, LRWORK, IWORK, LIWORK,INFO )
      	!IF( INFO.GT. 0 ) THEN
        !WRITE(*,*)'The algorithm failed to compute eigenvalues.'
        ! STOP
      	!END IF   


	allocate (RWORK(3*dimbse-2))
        LWMAX = 2*dimbse-1
	allocate(WORK(LWMAX))
	
	if (bsewf) then
	
	LWORK = -1
      	CALL ZHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, INFO )
      	LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      	CALL ZHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK,INFO )
      	IF( INFO.GT. 0 ) THEN
        WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      	END IF 
	
	
	else

     	LWORK = -1
      	CALL ZHEEV( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, INFO )
      	LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      	CALL ZHEEV( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK,INFO )
      	IF( INFO.GT. 0 ) THEN
        WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      	END IF 
      	
      	end if  

		!write(*,*) W(1)


		do i2=1,dimbse

			exk(i2,i)=W(i2)

		!	write(400,*) qauxv(i,1),W(i2)
		!	call flush(400)

		end do


	if (bsewf) then
	      	do i2=excwf0,excwff
      			call excwfi2(outputfolder,ngkpt,kpt,q,nc,nv,nocpk,stt,W(i2),i2,i,hbse(:,i2))
      		end do
	end if

		deallocate(hbse)
		deallocate(W)
		deallocate(energyq,vectorq)
		deallocate(stt)
		deallocate(RWORK,WORK)
		!deallocate (IWORK)


		qpt=0.0
		nocpq=0

		if(Process==0) write(300,*) 'progress:',i,'/',(nks/2)*nkpts
		if(Process==0) call flush(300)
	end do

 if(Process==0) then
	do i2=1,dimbse
		do i=1,(nks/2)*nkpts
			write(400,*) real(qauxv(i,1)),exk(i2,i)
			call flush(400)
		end do
		write(400,*)
	end do
 endif


	deallocate(energy,vector)
	deallocate(kpt)
	deallocate(qpt)
	deallocate(exk)
	deallocate(nocpq,nocpk)
	deallocate(rvec,hopmatrices)
	deallocate(ihopmatrices,ffactor)

 if(Process==0) then
	call cpu_time(tf)
	call date_and_time(VALUES=values2)
	write(300,*)
	write(300,*) 'end','   ','day',values2(3),'',values2(5),'hours',values2(6),'min',values2(7),'seg'
	write(300,*)
	close(200)
	close(300)
	close(400)
	close(500)
 endif
end subroutine bsebnds



!INCLUDE "./subprograms/bse_kpath-tool-temp.f90"
!gfortran -mcmodel=large bse_kpath-tool.f90 -o bse_kpath-tool.x -llapack95 -lopenblas -fopenmp 
subroutine bsebndstemp(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,ez,w1,r0,lc,rk,meshtype,bsewf,excwf0,excwff,st,phavg,ta,temp)
	double precision,parameter:: pi=acos(-1.)
	integer :: dimbse
	double precision,allocatable,dimension(:,:) :: exk !(dimbse,nqpts*npath)
	integer :: ncaux,nvaux
	integer,allocatable,dimension(:,:) :: stt !(ngrid*ngrid*nc*nv,4)
	double precision,allocatable,dimension(:,:) :: kpt,qpt !pontos k do grid (ngrid*ngrid,2)
	double precision,dimension(4) :: q
	double precision,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
	double complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores
	double precision,allocatable,dimension(:,:) :: energy !variavel que guarda os autovalores para cada ponto k
	double complex,allocatable,dimension(:,:,:) :: vector !variavel que guarda os autovetores para cada ponto k
	double precision,allocatable,dimension(:,:) :: energyq !variavel que guarda os autovalores para cada ponto k
	double complex,allocatable,dimension(:,:,:) :: vectorq !variavel que guarda os autovetores para cada ponto k
	double precision,allocatable,dimension(:,:) :: qauxv !(nqpts*npath,4) 
	integer:: counter,c,v,i,j,f,l,h,erro,i2,k
	double precision :: a,r0,ed
	integer :: ngkpt
	!variaveis relacionadas a marcacao do tempo
	double precision:: t0,tf
	integer,dimension(8) :: values,values2
	!variaveis kpoints
	integer :: nks
	integer :: nkpts
	double precision,allocatable,dimension(:,:) :: ks
	integer,allocatable,dimension(:) :: nocpk,nocpq
	INTEGER   ::       ifail
	double precision,parameter :: ABSTOL=1.0e-6
	INTEGER          INFO
	double precision,allocatable,dimension(:) :: W,RWORK
        INTEGER ::          LWMAX
   	INTEGER ::         LWORK
	INTEGER ::         LIWORK, LRWORK
	INTEGER,allocatable,dimension(:) :: IWORK
        double complex,allocatable,dimension (:) :: WORK
	integer :: nthreads
	integer,dimension(3) :: ngrid
	integer :: nc,nv
	double precision :: edos0,edosf,numdos
	double precision :: ebse0,ebsef,numbse
	double precision :: sme,exc,rk
	double precision,dimension(3) :: mshift
	double precision :: ktol
	character(len=70) :: params   !parametros TB
	character(len=70) :: orbw     !peso orbitais lcount
	character(len=70) :: kpaths    !kpath
	character(len=70) :: kpathsbse    !kpath
	character(len=70) :: outputfolder    !pasta saida
	character(len=70) :: calcparms
	character(len=70) :: meshtype
	character(len=5) :: coultype
	character(len=2) :: ta
	double precision,dimension(3) :: ediel
	double precision :: ez,w1,lc
	logical :: bsewf
	integer :: excwf0,excwff
	double precision :: st,phavg,temp
    double precision :: tcor

 if(Process==0) print *,'****** bsebndstemp'

 if(Process==0) then
	OPEN(UNIT=500, FILE= kpathsbse,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse-kpath input file"
	!OPEN(UNIT=201, FILE= diein,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada ambiente dieletrico"

	!OUTPUT : criando arquivos de saida
	OPEN(UNIT=300, FILE=trim(outputfolder)//"/log_bse_kpath.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening log_bse_kpath output file"
	OPEN(UNIT=400, FILE=trim(outputfolder)//"/bands_bse.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bands_bse output file"
	call cpu_time(t0)
	call date_and_time(VALUES=values)
!	call OMP_SET_NUM_THREADS(nthreads)
 endif

	call hamiltonian_input_read(200,params)
	!ediel(2) = edielh

 if(Process==0) then
	read(500,*) nks
	read(500,*) nkpts
 endif
    call P_sendI(nks)
    call P_sendI(nkpts)

	allocate(ks(nks,3))

 if(Process==0) then
	do i=1,nks
		read(500,*) ks(i,1),ks(i,2),ks(i,3)
	end do
 endif
    call P_sendR2(ks,nks,3)


	!termino leitura parametros


	!parametros do calculo


	!termino parametros calculo 




	ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
	dimbse = ngkpt*nc*nv

	!call alat(systype,rlat,a)
	!lc = rlat3(3)

	!if (systype .eq. "2D") then
	!	ed = (ediel(1)+ediel(3))/2.
	!	r0= ((ediel(2)-1.0)*lc)/(ediel(1)+ediel(3))
	!else
	!	r0 = 1.0
	!	ed = ediel(2)
		
	!end if

	!definindo kpath
	!allocate(qauxv(nkpts*(nks-1),4))
	allocate(qauxv((nks/2)*nkpts,4))	

	call kpathbse(outputfolder,rlat(1,:),rlat(2,:),rlat(3,:),nks,ks,nkpts,qauxv)

	!Informações para o arquivo de log do calculo

 if(Process==0) then
	write(300,*) 'threads:', nthreads
	write(300,*)
	write(300,*) 'temperature', temp
	write(300,*)
	write(300,*) 'grid:',ngrid(1),ngrid(2),ngrid(3)
	write(300,*)
	write(300,*) 'ktol-coulomb:', ktol
	write(300,*)
	write(300,*) 'kmesh shift:',mshift(1),mshift(2),mshift(3)
	write(300,*)
	write(300,*) 'conduction bands, nc:',nc,'  ','valence bands, nv:',nv
	write(300,*)
	write(300,*) 'number of kpoints in the path:','   ',(nks/2)*nkpts
	write(300,*)
	write(300,*)
	write(300,*) 'begin','  ','day',values(3),'',values(5),'hours',values(6),'min',values(7),'seg'
	write(300,*) 
	call flush(300)
 endif

	allocate(kpt(ngkpt,3))
	
	!shift= 0.0
	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpt)
	!call gridgenmhp(ngrid,rlat,kpt)

	allocate(eaux(w90basis),vaux(w90basis,w90basis))
	allocate(energy(ngkpt,nc+nv),vector(ngkpt,nc+nv,w90basis))

	allocate(nocpk(ngkpt))
	allocate(nocpq(ngkpt))

	select case (ta)
	
	case("FA")
	
		tcor = 0.0
	
	case("VE")
	
		tcor = gapcortemp(st,phavg,temp)
	
	case("BE")
	
		tcor = gapcortemp2(st,phavg,temp)
	
	case default
	
		tcor= 0.00
	
	end select

	! $omp parallel do
	do i=1,ngkpt
		call eigsys(nthreads,scs+tcor,exc,nocpk(i),&
			    ffactor,kpt(i,1),kpt(i,2),kpt(i,3),w90basis,nvec,&
			    rlat,rvec,hopmatrices,&
		            ihopmatrices,efermi,eaux,vaux)
			do j=1,nc+nv
				energy(i,j)= eaux(nocpk(i)-nv+j)
			end do
			do l=1,nc+nv
				do h=1,w90basis
					vector(i,l,h)=vaux(nocpk(i)-nv+l,h)
				end do
			end do
	end do
	! $omp end parallel do

	deallocate(eaux,vaux)


	!definindo os numeros quanticos dos estados




	if(Process==0) write(300,*) "quantum numbers for exciton basis set finished"
	if(Process==0) call flush(300)	

	counter=counter-1 !numero total de estados para equação bse

	allocate(qpt(ngkpt,3))
	allocate(exk(dimbse,nkpts*(nks-1)))



	do i=1,(nks/2)*nkpts

		
		!definindo os pontos q

		q(1)= qauxv(i,1)
		q(2)= qauxv(i,2) 
		q(3)= qauxv(i,3)
		q(4)= qauxv(i,4)

		!gerando o grid k+q
		
		call monhkhorst_packq(q,ngrid(1),ngrid(2),ngrid(3),mshift,&
				      rlat(1,:),rlat(2,:),rlat(3,:),qpt)


		allocate(eaux(w90basis),vaux(w90basis,w90basis))

		allocate(energyq(ngkpt,nc+nv),vectorq(ngkpt,nc+nv,w90basis))

	allocate (stt(ngkpt*nc*nv,4))


	! $omp parallel do
		do i2=1,ngkpt



			call eigsys(nthreads,scs+tcor,exc,nocpq(i2),&
			            ffactor,qpt(i2,1),qpt(i2,2),qpt(i2,3),w90basis,nvec,&
			             rlat,rvec,hopmatrices,&
		                      ihopmatrices,efermi,eaux,vaux)

				do j=1,nc+nv
	
					energyq(i2,j)= eaux(nocpq(i2)-nv+j) 

				end do

	
				do l=1,nc+nv


					do h=1,w90basis

						vectorq(i2,l,h)=vaux(nocpq(i2)-nv+l,h)

	
					end do
			

				end do



	
		end do
	! $omp end parallel do

		call quantumnumbers2(w90basis,ngkpt,nc,nv,nocpk,nocpq,stt)

		deallocate(eaux,vaux)


		allocate(hbse(dimbse,dimbse),W(dimbse))

		hbse=0.

	!$omp parallel do 
        !collapse(2)

		do i2=1,dimbse



			do j=i2,dimbse




hbse(i2,j)= matrizelbsekqtemp(coultype,ktol,w90basis,ediel,lc,ez,w1,r0,ngrid,q,rlat,stt(i2,:),energyq(stt(i2,4),stt(i2,3))&
          ,energy(stt(i2,4),stt(i2,2)),vectorq(stt(i2,4)&
          ,stt(i2,3),:) ,vector(stt(i2,4),stt(i2,2),:),kpt(stt(i2,4),:),stt(j,:)&
          ,energyq(stt(j,4),stt(j,3)),energy(stt(j,4),stt(j,2))&
          ,vectorq(stt(j,4),stt(j,3),:),vector(stt(j,4),stt(j,2),:),kpt(stt(j,4),:),temp)


			end do



		end do

	!$omp end parallel do

    		!call LA_HEEVR( hbse, W, JOBZ='N', UPLO='U', ABSTOL=ABSTOL, INFO=INFO ) 



      	!LWORK = 2*dimbse+dimbse**2
      	!LIWORK = 3 + 5*dimbse
      	!LRWORK = 1 + 5*dimbse + 2*dimbse**2

	!allocate(WORK(LWORK))
	!allocate (IWORK(LIWORK))
	!allocate (RWORK(LRWORK))

 

      	!CALL ZHEEVD( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
      	!LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      	!LRWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      	!LIWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

      	!CALL ZHEEVD( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK,&
	!	      RWORK, LRWORK, IWORK, LIWORK,INFO )
      	!IF( INFO.GT. 0 ) THEN
        !WRITE(*,*)'The algorithm failed to compute eigenvalues.'
        ! STOP
      	!END IF   


	allocate (RWORK(3*dimbse-2))
        LWMAX = 2*dimbse-1
	allocate(WORK(LWMAX))
	
	if (bsewf) then
	
	LWORK = -1
      	CALL ZHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, INFO )
      	LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      	CALL ZHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK,INFO )
      	IF( INFO.GT. 0 ) THEN
        WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      	END IF 
	
	
	else

     	LWORK = -1
      	CALL ZHEEV( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, INFO )
      	LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      	CALL ZHEEV( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK,INFO )
      	IF( INFO.GT. 0 ) THEN
        WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      	END IF 
      	
      	end if  

		!write(*,*) W(1)


		do i2=1,dimbse

			exk(i2,i)=W(i2)

		!	write(400,*) qauxv(i,1),W(i2)
		!	call flush(400)

		end do


	if (bsewf) then
	      	do i2=excwf0,excwff
      			call excwfi2(outputfolder,ngkpt,kpt,q,nc,nv,nocpk,stt,W(i2),i2,i,hbse(:,i2))
      		end do
	end if

		deallocate(hbse,W)
		deallocate(energyq,vectorq)
		deallocate(stt)
		deallocate(RWORK,WORK)
		!deallocate (IWORK)


		qpt=0.0
		nocpq=0

		write(300,*) 'progress:',i,'/',(nks/2)*nkpts
		call flush(300)

	end do



	do i2=1,dimbse

		do i=1,(nks/2)*nkpts

			write(400,*) real(qauxv(i,1)),exk(i2,i)
			call flush(400)

		end do

		write(400,*)

	end do



	deallocate(energy,vector)
	deallocate(kpt)
	deallocate(qpt)
	deallocate(exk)
	deallocate(nocpq,nocpk)
	deallocate(rvec,hopmatrices)
	deallocate(ihopmatrices,ffactor)

 if(Process==0) then
	call cpu_time(tf)
	call date_and_time(VALUES=values2)
	write(300,*)
	write(300,*) 'end','   ','day',values2(3),'',values2(5),'hours',values2(6),'min',values2(7),'seg'
	write(300,*)
	close(200)
	close(300)
	close(400)
	close(500)
 endif
end subroutine bsebndstemp



!#ifdef FAST1
!INCLUDE "./subprograms/bse_solver-tool-diel-f.f90"
!INCLUDE "./subprograms/bse_solver-tool-diel-temp-f.f90"
!#else
!INCLUDE "./subprograms/bse_solver-tool-diel.f90"
subroutine bsesolver(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,ez,w1,r0,lc,rk,meshtype,bsewf,excwf0,excwff,&
		     dtfull,cpol,tmcoef)
	double precision,parameter:: pi=acos(-1.)
	double precision,allocatable,dimension(:,:) :: eigv
	double complex,allocatable,dimension(:,:,:) :: vector
	double precision,allocatable,dimension(:,:) :: kpt
	integer,allocatable,dimension(:,:) :: stt
	double precision,allocatable,dimension(:) :: eaux
	double complex,allocatable,dimension(:,:) :: vaux
	integer:: counter,c,v,i,j,f,l,h,k,kl,erro
	integer :: ncaux,nvaux
	double complex,allocatable,dimension(:) :: hrx,hry,hrz,hrsp,hrsm 
	double precision,allocatable,dimension(:) :: actxx,actyy,actzz,actxy,actxz,actyz,actsp,actsm
	double complex,parameter :: imag=(0.d0,1.d0)
	double precision :: a,r0,ed,ec,ev,egap
	integer :: ngkpt
	real,allocatable,dimension(:,:) :: vecres
	double precision:: t0,tf
	integer,dimension(8) :: values,values2
	integer,allocatable,dimension(:) :: nocpk
	INTEGER   ::       ifail
	double precision,parameter :: ABSTOL=1.0e-6
!	COMPLEX*16,allocatable,dimension(:,:) :: hbse
!	INTEGER ::         LIWORK, LRWORK
!	INTEGER,allocatable,dimension(:) :: IWORK
	integer :: nthreads
	integer,dimension(3) :: ngrid
	integer :: nc,nv
	double precision :: edos0,edosf,numdos
	double precision :: ebse0,ebsef,numbse
	double precision :: sme,exc,rk
	double precision,dimension(3) :: mshift
	double precision :: ktol
	character(len=70) :: params
	character(len=70) :: orbw
	character(len=70) :: kpaths
	character(len=70) :: kpathsbse
	character(len=70) :: outputfolder
	character(len=70) :: calcparms
	character(len=70) :: meshtype
	character(len=5) :: coultype
	double precision,dimension(3) :: ediel
	logical :: bsewf
	integer :: excwf0,excwff
	double precision :: ez,w1,lc
	logical :: cpol,dtfull
	logical :: tmcoef
	integer                  :: dimbse1      ! number of BSE solutions
	integer                  :: ist,ngr      ! for parallel calculation

!  if(Process==0) print *,'  bsesolver'

  call MSG('MODULE_BSE: open files')
 if(Process==0) then
 ! print *,'bsesolver: open files'
	OPEN(UNIT=300, FILE=trim(outputfolder)//"/log_bse-diel.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening log_bse-diel output file"
	OPEN(UNIT=301, FILE=trim(outputfolder)//"/bse_opt_diel.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse_opt_diel output file"
    	OPEN(UNIT=302, FILE=trim(outputfolder)//"/bse_opt_diel-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening output file"	


	OPEN(UNIT=401, FILE=trim(outputfolder)//"/sp_optics.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse_opt_diel-pol output file"
    	
    	
    	if (tmcoef) then
    	OPEN(UNIT=402, FILE=trim(outputfolder)//"/tm_coef.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening tm_coef output file"
    	OPEN(UNIT=404, FILE=trim(outputfolder)//"/tm_coef-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening tm_coef-pol output file"    	
    	else
    	 continue
    	end if    
    	
    	OPEN(UNIT=403, FILE=trim(outputfolder)//"/sp_optics-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening sp_optics-pol output file"	

	call cpu_time(t0)
	call date_and_time(VALUES=values)
!!	call OMP_SET_NUM_THREADS(nthreads)
!#ifdef gnu
!	call OPENBLAS_SET_NUM_THREADS(nthreads)
!#endif
 endif


	!inicio leitura parametros
	
	!if(Process==0) print *,' params=',params
  call MSG('MODULE_BSE: hamiltonian_input_read')
	call hamiltonian_input_read(200,params)
	!if(Process==0) print *,' end hamiltonian_input_read '

	!ediel(2) = edielh


	!termino leitura parametros

	!parametros do calculo


	!dtfull = .false.
	!cpol = .true.


	 if(Process==0) write(301,*) "#","  ", "exciton energy","  ","xx","  ","yy","  ","zz"," ","xy","  ","xz","  ","yz"
	 if(Process==0) write(302,*) "#","  ", "exciton energy","  ","xx","  ","yy","  ","zz"," ","sp","  ","sm"	


	!termino parametros calculo 

	ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
	dimbse = ngkpt*nc*nv


	if(Process==0) print *,'   dimbse=',dimbse

	!call alat(systype,rlat,a)

	!if (systype .eq. "2D") then
	!	ed = (ediel(1)+ediel(3))/2.
	!	r0= ((ediel(2)-1.0)*lc)/(ediel(1)+ediel(3))
	!else
	!	r0 = 1.0
	!	ed = ediel(2)
		
	!end if

 if(Process==0) then
	!Informações para o arquivo de log do calculo
	write(300,*)
	write(300,*)
	write(300,*)
	write(300,*) 'threads:', nthreads
	write(300,*)
	write(300,*) 'grid:',ngrid(1),ngrid(2),ngrid(3)
	write(300,*)
	write(300,"(A14,1E15.4)") 'ktol-coulomb:', ktol
	write(300,*)
	write(300,"(A13,3F15.4)") 'kmesh shift:',mshift(1),mshift(2),mshift(3)
	!write(300,*) 'angulo theta:',beta,'rad'

	write(300,*)
	write(300,*) 'conduction bands, nc:',nc,'  ','valence bands, nv:',nv
	write(300,*)
	!write(300,*) 'constante dielétrica:',ed,'  ','r0:',r0

	!write(300,*)
	!write(300,*) 'intensidade termo exchange:',exc
	!write(300,*)
	!write(300,*) 'intensidade termo exchange l:',excl
	!write(300,*)
	!write(300,*) 'intensidade campo eletrico:',ez
	!write(300,*)
	!write(300,*) 'direção da magnetização do termo de exchange:',mag(1),'x',mag(2),'y',mag(3),'z'
	!write(300,*)

	!write(300,*) 'arquivo de saida 1:',logbse
	!write(300,*) 'arquivo de saida 2:',bseoptx
	!write(300,*) 'arquivo de saida 3:',bseopty
	!write(300,*) 'arquivo de saida 4:',bseoptsp
	!write(300,*) 'arquivo de saida 5:',bseoptsm

	
	write(300,*)
	write(300,*) 'begin','  ','day',values(3),'',values(5),'hours',values(6),'min',values(7),'seg'
	write(300,*) 

	call flush(300)
 endif
 
	allocate(kpt(ngkpt,3))

	!shift = 0.0
	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpt)


	allocate(eaux(w90basis),vaux(w90basis,w90basis))
	allocate(eigv(ngkpt,nc+nv),vector(ngkpt,nc+nv,w90basis))
	allocate(nocpk(ngkpt))

	allocate (stt(ngkpt*nc*nv,4))
	allocate(hrx(dimbse),hry(dimbse),hrz(dimbse),hrsp(dimbse),hrsm(dimbse))

	allocate(hbse(dimbse,dimbse))
	allocate(Ws(dimbse))


	allocate(actxx(dimbse),actyy(dimbse),actzz(dimbse),actxy(dimbse))
	allocate(actxz(dimbse),actyz(dimbse))
	allocate(actsp(dimbse),actsm(dimbse))	

	!allocate(orbweight(w90basis))
	
	!do i=1,w90basis
	!	read(203,*) orbweight(i)
	!end do

	!allocate(lcount(ngkpt,w90basis))

	egap = 50.0

	! $omp parallel do
	do i=1,ngkpt
		call eigsys(nthreads,scs,exc,nocpk(i),ffactor,kpt(i,1),kpt(i,2),kpt(i,3),w90basis,nvec,&
			    rlat,rvec,hopmatrices,&
		             ihopmatrices,efermi,eaux,vaux)
			do j=1,nc+nv
				eigv(i,j)= eaux(nocpk(i)-nv+j)
			end do
			if ((eigv(i,nv+1)-eigv(i,nv)) .le. egap) then
				egap = eigv(i,nv+1)-eigv(i,nv)
			end if
			do l=1,nc+nv
				do h=1,w90basis
					vector(i,l,h)=vaux(nocpk(i)-nv+l,h)
				end do
			end do
	end do
	! $omp end parallel do


    !print *,'calc eigsys'
     if(Process==0) call print_calc_t('calc eigsys',t0,tf)


 if(Process==0) then
    print *
    print *,'******** BSESOLVER:'
	print *,'Efermi=',efermi
	do i=1,ngkpt
	 print *,'k=',i,'  EFermi=',nocpk(i)
	enddo
	print *,'conduction bands, nc:',nc,'  ','valence bands, nv:',nv
    print *,'eigv:'
	do i=1,ngkpt
     print 11,i,(eigv(i,j),j=1,nc+nv)
    enddo
11  format(I3,4F10.5)    
 endif


	deallocate(eaux,vaux)
	 if(Process==0) write(300,*) 'direct gap:', egap
	 if(Process==0) write(300,*) 'eigenvalues and eigenvectors calculated'
	 if(Process==0) call flush(300)


	!definindo os numeros quanticos dos estados


	!allocate (stt(ngkpt*nc*nv,4))

	call quantumnumbers2(w90basis,ngkpt,nc,nv,nocpk,nocpk,stt)

    !if(Process==0) print *,'calc quantumnumbers2'
    !if(Process==0) call print_calc_t('calc quantumnumbers2',t0,tf)
    call MSGT('MODULE_BSE: quantumnumbers2')

	!allocate(hrx(dimbse),hry(dimbse),hrz(dimbse))

	 if(Process==0) write(300,*) 'quantum numbers for exciton basis set finished'
	 if(Process==0) call flush(300)


	
	allocate(vecres(dimbse,15))	
	
 if(Process==0) then
	write(401,*) "#","  ", "energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(403,*) "#","  ", "energy","  ","xx","  ","yy","  ","zz","  ","sp","  ","sm"
	
	if (tmcoef) then	
	write(402,*) "number of kpoints:",ngkpt
	write(402,*) "number of conduction states",nc
	write(402,*) "number of valence states",nv
	write(402,*) "#"," ", "kx", " ", "ky"," ", "kz"," ","nocp"," ", "nc"," ", "nv","  ", "energy","  ","xx","  ",&
		      "yy","  ","zz","  ","xy","  ","xz","  ","yz"
		      
	write(404,*) "number of kpoints:",ngkpt
	write(404,*) "number of conduction states",nc
	write(404,*) "number of valence states",nv
	write(404,*) "#"," ", "kx", " ", "ky"," ", "kz"," ","nocp"," ", "nc"," ", "nv","  ", "energy","  ","xx","  ",&
		      "yy","  ","zz","  ","sp","  ","sm"		   
	end if	
 endif

	! $omp parallel do &
	! $OMP DEFAULT (SHARED)	

	do i=1,dimbse

		!ec = eigv(stt(i,4),stt(i,3))
		!ev = eigv(stt(i,4),stt(i,2))

		call optsp(eigv(stt(i,4),stt(i,2)),vector(stt(i,4),stt(i,2),:),&
		     eigv(stt(i,4),stt(i,3)),vector(stt(i,4),stt(i,3),:),&
		     kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),ffactor,sme,&
		     w90basis,nvec,rlat,rvec,hopmatrices,&
		     ihopmatrices,hrx(i),hry(i),hrz(i))
		     
		     hrsp(i) = (hrx(i)+cmplx(0.,1.)*hry(i))*(1.0/sqrt(2.))
		     hrsm(i) = (hrx(i)-cmplx(0.,1.)*hry(i))*(1.0/sqrt(2.))
		     
		vecres(i,1) = real(kpt(stt(i,4),1))
		vecres(i,2) = real(kpt(stt(i,4),2))
		vecres(i,3) = real(kpt(stt(i,4),3))
		vecres(i,4) = real(nocpk(stt(i,4)))
		vecres(i,5) = real(nocpk(stt(i,4))-nv+stt(i,3))
		vecres(i,6) = real(nocpk(stt(i,4))-nv+stt(i,2))		
		vecres(i,7) = real(eigv(stt(i,4),stt(i,3))-eigv(stt(i,4),stt(i,2)))		
		vecres(i,8) = real(hrx(i)*conjg(hrx(i)))		
		vecres(i,9) = real(hry(i)*conjg(hry(i)))		
		vecres(i,10) = real(hrz(i)*conjg(hrz(i)))
		vecres(i,11) = real(hrx(i)*conjg(hry(i)))
		vecres(i,12) = real(hrx(i)*conjg(hrz(i)))
		vecres(i,13) = real(hry(i)*conjg(hrz(i)))		 
		vecres(i,14) = real(hrsp(i)*conjg(hrsp(i)))
		vecres(i,15) = real(hrsm(i)*conjg(hrsm(i)))				     


	end do
	
    !if(Process==0) print *,'calc optsp'
!    if(Process==0) call print_calc_t('calc optsp',t0,tf)
  call MSGT('MODULE_BSE: optsp')


	! $omp end parallel do
	
	call Bubblem(7,15,vecres, dimbse)

    !if(Process==0) print *,'MODULE_BSE Bubblem'
!    if(Process==0) call print_calc_t('MODULE_BSE Bubblem',t0,tf)
  call MSGT('MODULE_BSE: Bubblem')

 if(Process==0) then
	do i=1,dimbse
	
		write(401,"(7F15.4)") vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),vecres(i,11),vecres(i,12),vecres(i,13)
		write(403,"(6F15.4)") vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),vecres(i,14),vecres(i,15)
		
		if (tmcoef) then
	        write(402,"(3F10.4,3I10.0,7F10.4)") vecres(i,1),vecres(i,2),vecres(i,3),int(vecres(i,4)),int(vecres(i,5)),&
						  int(vecres(i,6)),vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),&
						  vecres(i,11),vecres(i,12),vecres(i,13)
	        
	        write(404,"(3F10.4,3I10.0,6F10.4)") vecres(i,1),vecres(i,2),vecres(i,3),int(vecres(i,4)),int(vecres(i,5)),&
						  int(vecres(i,6)),vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),&
						  vecres(i,14),vecres(i,15)						  
						  
		end if
	
	end do	

	write(300,*) 'single particle optics finished'
	call flush(300)
 endif
 
	deallocate(vecres)
	
	!go to 789

	!allocate(hbse(dimbse,dimbse),W(dimbse))

	hbse=0.0

	!$omp parallel do 
!collapse(2)

!   call P_calc_group(dimbse,ist,ngr)                   ! calculate groups to calculate in MPI parallel
!    do i=ist+1,ist+ngr                                ! not even distribution !!!!!!!!!
    do i=1,dimbse
        do j=i,dimbse
  hbse(i,j)=matrizelbse(coultype,ktol,w90basis,ediel,lc,ez,w1,r0,ngrid,rlat,stt(i,:),eigv(stt(i,4)&
            ,stt(i,3)),eigv(stt(i,4),stt(i,2)),vector(stt(i,4)&
            ,stt(i,3),:) ,vector(stt(i,4),stt(i,2),:),kpt(stt(i,4),:),stt(j,:),eigv(stt(j,4),stt(j,3))&
            ,eigv(stt(j,4),stt(j,2)),vector(stt(j,4),stt(j,3),:),vector(stt(j,4),stt(j,2),:),kpt(stt(j,4),:))
        end do
    end do
!    do j=1,dimbse
!     call P_receiveZ(hbse(1,j),dimbse,ist,ngr)             ! collect calculated parts hbse(ist+1:ist+ngr,j) from all processes
!    enddo

    do i=1,dimbse    !SLEPc
     do j=1,i-1
      hbse(i,j)=conjg(hbse(j,i))
     end do
    end do
    
  ltest = .true.
!  if(Process==0) then
!   if(ltest) then
!    print *,'bsesolver: hbse 1st row'
!    do i=1,1
!     print 1,hbse(1:dimbse,i)
!    enddo
!   endif
!  endif

    !if(Process==0) print *,'calc hbse'
!    if(Process==0) call print_calc_t('MODULE_BSE calculate hbse matrix',t0,tf)
  call MSGT('MODULE_BSE: calculate hbse matrix')


   dimbse1 = 100
   allocate(HE(dimbse,dimbse1))
   HE = (0.d0,0.d0)
!   if(Process==0) print *,'call slepc_solver'

   call slepc_solver(hbse,Ws,HE,dimbse,dimbse1)          !(dimbse,dimbse1) !
!   call slepc_check_results(hbse,Ws,HE,dimbse,dimbse1)

!  if(Process==0) then
!    print *,'dimbse1=',dimbse1
!   if(ltest) then
!    print *,'bsesolver: HE first solution'
!    do i=1,1
!     print 1,HE(1:dimbse,i)
!    enddo
!1   format(100(2E12.3,3x))
!   endif
!  endif

   hbse = (0.d0,0.d0)
   do i=1,dimbse1
    hbse(:,i) = HE(:,i)                  ! write the dimbse1 solutions
   enddo
   deallocate(HE)

!   call mkl_solver

    !if(Process==0) print *,'calc SLEPc'
!    if(Process==0) call print_calc_t('MODULE_BSE SLEPc calculation',t0,tf)
  call MSGT('MODULE_BSE: SLEPc total calculation time')

 if(Process==0) then
	if (bsewf) then
	      	do i=excwf0,excwff
      			call excwf2(outputfolder,ngkpt,kpt,nc,nv,nocpk,stt,Ws(i),i,hbse(:,i))
      		end do
	end if
 endif

     !if(Process==0) print *,'calc excwf2'
!     if(Process==0) call print_calc_t('MODULE_BSE excwf2',t0,tf)
  call MSGT('MODULE_BSE: excwf2')

	 if(Process==0) write(300,*) 'exciton Hamiltonian diagonalized'
	 if(Process==0) call flush(300)
	 if(Process==0) write(300,*) "exciton ground state",Ws(1)

	deallocate(eigv,vector)
	deallocate(rvec,hopmatrices)
	deallocate(ihopmatrices,ffactor)

 ltest = .true.
	call dielbsev(nthreads,dimbse,dimbse1,hbse,hrx,actxx)                     ! change dimbse to dimbse1
	if(Process==0) write(300,*) 'xx tensor component'
	if(Process==0) call flush(300)
		
	call dielbsev(nthreads,dimbse,dimbse1,hbse,hry,actyy)
	if(Process==0) write(300,*) 'yy tensor component'
	if(Process==0) call flush(300)
		
	call dielbsev(nthreads,dimbse,dimbse1,hbse,hrz,actzz)
	if(Process==0) write(300,*) 'zz tensor component'
	if(Process==0) call flush(300)
	
	if (dtfull) then
	
	call dielbsep(nthreads,dimbse,dimbse1,hbse,hrx,hry,actxy)
	if(Process==0) write(300,*) 'xy tensor component'
	if(Process==0) call flush(300)
	
	call dielbsep(nthreads,dimbse,dimbse1,hbse,hrx,hrz,actxz)
	if(Process==0) write(300,*) 'xz tensor component'
	if(Process==0) call flush(300)
		
	call dielbsep(nthreads,dimbse,dimbse1,hbse,hry,hrz,actyz)
	if(Process==0) write(300,*) 'yz tensor component'
	if(Process==0) call flush(300)
	
	else
	
	actxy = 0.0000
	if(Process==0) write(300,*) 'xy tensor component set to 0'
	if(Process==0) call flush(300)
	
	actxz = 0.0000
	if(Process==0) write(300,*) 'xz tensor component set to 0'
	if(Process==0) call flush(300)
		
	actyz = 0.0000
	if(Process==0) write(300,*) 'yz tensor component set to 0'
	if(Process==0) call flush(300)
	
	end if
	
	if (cpol) then
	
	call dielbsev(nthreads,dimbse,dimbse1,hbse,hrsp,actsp)
	if(Process==0) write(300,*) 'sp polarization'
	if(Process==0) call flush(300)
		
	call dielbsev(nthreads,dimbse,dimbse1,hbse,hrsm,actsm)
	if(Process==0) write(300,*) 'sm polarization'
	if(Process==0) call flush(300)
	
	else 
	
	actsp = 0.0000
	if(Process==0) write(300,*) 'sp polarization set to 0'
	if(Process==0) call flush(300)
		
	actsm = 0.0000
	if(Process==0) write(300,*) 'sm polarization set to 0'
	if(Process==0) call flush(300)
	
	end if	

 if(Process==0) then
!    call print_calc_t('MODULE_BSE dielbse',t0,tf)
  call MSGT('MODULE_BSE: dielbse')

	write(300,*) 'optics finished'
	call flush(300)

	!$omp do ordered
	do i=1,dimbse                 ! 1                                 ! change dimbse to dimbse1
		!$omp ordered
		write(301,"(7F15.4)") Ws(i),actxx(i),actyy(i),actzz(i),actxy(i),actxz(i),actyz(i)
		call flush(301)
		
		if (cpol) then
		write(302,"(6F15.4)") Ws(i),actxx(i),actyy(i),actzz(i),actsp(i),actsm(i)
		call flush(302)		
		end if

		!$omp end ordered
	end do
	!$omp end do
 endif
 
	deallocate(hrx,hry,hrz,hrsp,hrsm)
	deallocate(actxx,actxy,actxz,actyy,actyz,actzz)
	deallocate(actsp,actsm)
	deallocate(hbse,Ws,stt,nocpk)
	deallocate(kpt)
	
!	deallocate(HE)

 if(Process==0) then
	call cpu_time(tf)
	call date_and_time(VALUES=values2)
	write(300,*)
	write(300,*) 'end','   ','day',values2(3),'',values2(5),'hours',values2(6),'min',values2(7),'seg'
	write(300,*)
	close(200)
	close(300)
	close(301)
	close(302)
	close(401)
	close(402)
	close(403)
	close(404)			
 endif
end subroutine bsesolver



!subroutine bsesolver_test
!  integer     :: i,j
!  dimbse = 7744
!  allocate(hbse(dimbse,dimbse))
!  allocate(Ws(dimbse))
!  hbse = (0.d0,0.d0)
!  do i=1,dimbse
!   hbse(i,i) = (1.d0,0.d0)
!  enddo 
!  call slepc_solver
!end subroutine bsesolver_test





   subroutine mkl_solver
    real(8),allocatable         :: RWORK(:)
    complex(8),allocatable      :: WORK(:)
    INTEGER                     :: LWMAX
    INTEGER                     :: INFO
    INTEGER                     :: LWORK
    allocate(RWORK(3*dimbse-2))
    LWMAX = 2*dimbse-1
    allocate(WORK(LWMAX))

!MKL
!   !$omp end parallel do
!    print *,'calc exciton Hamiltonian matrix'
!    call print_calc_t(t0,tf)
!	write(300,*) 'exciton Hamiltonian matrix finished'
!	call flush(300)	
     	LWORK = -1
      	CALL ZHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, Ws, WORK, LWORK, RWORK, INFO )
      	LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      	CALL ZHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, Ws, WORK, LWORK, RWORK,INFO )
      	IF( INFO.GT. 0 ) THEN
        WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      	END IF   
!MKL
!    print *,'calc ZHEEV'
    deallocate(WORK)
    deallocate(RWORK)
   end subroutine mkl_solver





!INCLUDE "./subprograms/bse_solver-tool-diel-temp.f90"
subroutine bsesolvertemp(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,ez,w1,r0,lc,rk,meshtype,bsewf,excwf0,excwff,&
		     dtfull,cpol,tmcoef,st,phavg,ta,temp)
	double precision,parameter:: pi=acos(-1.)
	integer :: dimbse !=ngrid*ngrid*nc*nv ! dimensão da matriz bse
	double precision,allocatable,dimension(:,:) :: eigv
	double complex,allocatable,dimension(:,:,:) :: vector
	double precision,allocatable,dimension(:,:) :: kpt !pontos k do grid
	integer,allocatable,dimension(:,:) :: stt
	double precision,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
	double complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores
	integer:: counter,c,v,i,j,f,l,h,k,kl,erro
	integer :: ncaux,nvaux
	double complex,allocatable,dimension(:) :: hrx,hry,hrz,hrsp,hrsm 
	double precision,allocatable,dimension(:) :: actxx,actyy,actzz,actxy,actxz,actyz,actsp,actsm
	double complex,parameter :: imag=cmplx(0.0,1.0)
	double precision :: a,r0,ed,ec,ev,egap
	integer :: ngkpt
	real,allocatable,dimension(:,:) :: vecres
	double precision:: t0,tf
	integer,dimension(8) :: values,values2
	integer,allocatable,dimension(:) :: nocpk
	INTEGER   ::       ifail
	double precision,parameter :: ABSTOL=1.0e-6
	INTEGER          INFO
	double precision,allocatable,dimension(:) :: W,RWORK
	COMPLEX*16,allocatable,dimension(:,:) :: hbse
        INTEGER ::          LWMAX
   	INTEGER ::         LWORK
	INTEGER ::         LIWORK, LRWORK
	INTEGER,allocatable,dimension(:) :: IWORK
        double complex,allocatable,dimension (:) :: WORK
	integer :: nthreads
	integer,dimension(3) :: ngrid
	integer :: nc,nv
	double precision :: edos0,edosf,numdos
	double precision :: ebse0,ebsef,numbse
	double precision :: sme,exc,rk
	double precision,dimension(3) :: mshift
	double precision :: ktol
	character(len=70) :: params   !parametros TB
	character(len=70) :: orbw     !peso orbitais lcount
	character(len=70) :: kpaths    !kpath
	character(len=70) :: kpathsbse    !kpath
	character(len=70) :: outputfolder    !pasta saida
	character(len=70) :: calcparms
	character(len=70) :: meshtype
	character(len=5) :: coultype
	character(len=2) :: ta
	double precision,dimension(3) :: ediel
	logical :: bsewf
	integer :: excwf0,excwff
	double precision :: ez,w1,lc
	logical :: cpol,dtfull
	logical :: tmcoef	
	double precision :: st,phavg,temp
	double precision :: tcor
	double precision,allocatable,dimension(:) :: fdeh	

	!call input_read

 if(Process==0) then
	! INPUT : lendo os parametros do modelo de tight-binding
	!OPEN(UNIT=203, FILE= orbw,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada orb weight"

	!OUTPUT
	OPEN(UNIT=300, FILE=trim(outputfolder)//"/log_bse-diel.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening log_bse-diel output file"
	OPEN(UNIT=301, FILE=trim(outputfolder)//"/bse_opt_diel.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse_opt_diel output file"
    	OPEN(UNIT=302, FILE=trim(outputfolder)//"/bse_opt_diel-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening output file"	


	OPEN(UNIT=401, FILE=trim(outputfolder)//"/sp_optics.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse_opt_diel-pol output file"
    	
    	
    	if (tmcoef) then
    	OPEN(UNIT=402, FILE=trim(outputfolder)//"/tm_coef.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening tm_coef output file"
    	OPEN(UNIT=404, FILE=trim(outputfolder)//"/tm_coef-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening tm_coef-pol output file"    	
    	else
    	 continue
    	end if    
    	
    	OPEN(UNIT=403, FILE=trim(outputfolder)//"/sp_optics-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening sp_optics-pol output file"	
	call cpu_time(t0)
	call date_and_time(VALUES=values)
!	call OMP_SET_NUM_THREADS(nthreads)
	!call OPENBLAS_SET_NUM_THREADS(nthreads)
 endif


	!inicio leitura parametros

	call hamiltonian_input_read(200,params)

	!ediel(2) = edielh


	!termino leitura parametros

	!parametros do calculo


	!dtfull = .false.
	!cpol = .true.

 if(Process==0) then
	write(301,*) "#","  ", "exciton energy","  ","xx","  ","yy","  ","zz"," ","xy","  ","xz","  ","yz"
	write(302,*) "#","  ", "exciton energy","  ","xx","  ","yy","  ","zz"," ","sp","  ","sm"	
 endif

	!termino parametros calculo 

	ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
	dimbse = ngkpt*nc*nv

	!call alat(systype,rlat,a)

	!if (systype .eq. "2D") then
	!	ed = (ediel(1)+ediel(3))/2.
	!	r0= ((ediel(2)-1.0)*lc)/(ediel(1)+ediel(3))
	!else
	!	r0 = 1.0
	!	ed = ediel(2)
		
	!end if

 if(Process==0) then
	!Informações para o arquivo de log do calculo
	write(300,*)
	write(300,*)
	write(300,*)
	write(300,*) 'threads:', nthreads
	write(300,*)
	write(300,*) 'temperature', temp
	write(300,*)
	write(300,*) 'grid:',ngrid(1),ngrid(2),ngrid(3)
	write(300,*)
	write(300,"(A14,1E15.4)") 'ktol-coulomb:', ktol
	write(300,*)
	write(300,"(A13,3F15.4)") 'kmesh shift:',mshift(1),mshift(2),mshift(3)
	!write(300,*) 'angulo theta:',beta,'rad'

	write(300,*)
	write(300,*) 'conduction bands, nc:',nc,'  ','valence bands, nv:',nv
	write(300,*)
	!write(300,*) 'constante dielétrica:',ed,'  ','r0:',r0

	!write(300,*)
	!write(300,*) 'intensidade termo exchange:',exc
	!write(300,*)
	!write(300,*) 'intensidade termo exchange l:',excl
	!write(300,*)
	!write(300,*) 'intensidade campo eletrico:',ez
	!write(300,*)
	!write(300,*) 'direção da magnetização do termo de exchange:',mag(1),'x',mag(2),'y',mag(3),'z'
	!write(300,*)
	!write(300,*) 'arquivo de saida 1:',logbse
	!write(300,*) 'arquivo de saida 2:',bseoptx
	!write(300,*) 'arquivo de saida 3:',bseopty
	!write(300,*) 'arquivo de saida 4:',bseoptsp
	!write(300,*) 'arquivo de saida 5:',bseoptsm
	write(300,*)
	write(300,*) 'begin','  ','day',values(3),'',values(5),'hours',values(6),'min',values(7),'seg'
	write(300,*) 
	call flush(300)
 endif
 
	allocate(kpt(ngkpt,3))

	!shift = 0.0
	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpt)


	allocate(eaux(w90basis),vaux(w90basis,w90basis))
	allocate(eigv(ngkpt,nc+nv),vector(ngkpt,nc+nv,w90basis))
	allocate(nocpk(ngkpt))

	allocate (stt(ngkpt*nc*nv,4))
	allocate(hrx(dimbse),hry(dimbse),hrz(dimbse),hrsp(dimbse),hrsm(dimbse))
	allocate(fdeh(dimbse))

	allocate(hbse(dimbse,dimbse),W(dimbse))

	allocate (RWORK(3*dimbse-2))
        LWMAX = 2*dimbse-1
	allocate(WORK(LWMAX))

	allocate(actxx(dimbse),actyy(dimbse),actzz(dimbse),actxy(dimbse))
	allocate(actxz(dimbse),actyz(dimbse))
	allocate(actsp(dimbse),actsm(dimbse))	

	!allocate(orbweight(w90basis))
	
	!do i=1,w90basis
	!	read(203,*) orbweight(i)
	!end do

	!allocate(lcount(ngkpt,w90basis))

	egap = 50.0
	
	select case (ta)
	
	case("FA")
	
		tcor = 0.0
	
	case("VE")
	
		tcor = gapcortemp(st,phavg,temp)
	
	case("BE")
	
		tcor = gapcortemp2(st,phavg,temp)
	
	case default
	
		tcor= 0.00
	
	end select



	! $omp parallel do
	do i=1,ngkpt


		call eigsys(nthreads,scs+tcor,exc,nocpk(i),ffactor,kpt(i,1),kpt(i,2),kpt(i,3),w90basis,nvec,&
			    rlat,rvec,hopmatrices,&
		             ihopmatrices,efermi,eaux,vaux)


			do j=1,nc+nv
	
				eigv(i,j)= eaux(nocpk(i)-nv+j)

			end do
			
			!if (ntype .eq. 2) then

			!do kl = 1,w90basis

			!	call layercont(w90basis,vaux(kl,:),orbweight,lcount(i,kl))


			!end do

			!else 

			!	continue

			!end if

			if ((eigv(i,nv+1)-eigv(i,nv)) .le. egap) then

				egap = eigv(i,nv+1)-eigv(i,nv)

			else

				continue

			end if

			


			do l=1,nc+nv


				do h=1,w90basis

					vector(i,l,h)=vaux(nocpk(i)-nv+l,h)


				end do
			

			end do

	
	end do
	! $omp end parallel do

	deallocate(eaux,vaux)
	if(Process==0) write(300,*) 'direct gap:', egap
	if(Process==0) write(300,*) 'eigenvalues and eigenvectors calculated'
	if(Process==0) call flush(300)


	!definindo os numeros quanticos dos estados


	!allocate (stt(ngkpt*nc*nv,4))

	call quantumnumbers2(w90basis,ngkpt,nc,nv,nocpk,nocpk,stt)



	!allocate(hrx(dimbse),hry(dimbse),hrz(dimbse))

	if(Process==0) write(300,*) 'quantum numbers for exciton basis set finished'
	if(Process==0) call flush(300)


	
	allocate(vecres(dimbse,15))	

 if(Process==0) then
	write(401,*) "#","  ", "energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(403,*) "#","  ", "energy","  ","xx","  ","yy","  ","zz","  ","sp","  ","sm"
	if (tmcoef) then	
	 write(402,*) "number of kpoints:",ngkpt
	 write(402,*) "number of conduction states",nc
	 write(402,*) "number of valence states",nv
	 write(402,*) "#"," ", "kx", " ", "ky"," ", "kz"," ","nocp"," ", "nc"," ", "nv","  ", "energy","  ","xx","  ",&
		      "yy","  ","zz","  ","xy","  ","xz","  ","yz"
	 write(404,*) "number of kpoints:",ngkpt
	 write(404,*) "number of conduction states",nc
	 write(404,*) "number of valence states",nv
	 write(404,*) "#"," ", "kx", " ", "ky"," ", "kz"," ","nocp"," ", "nc"," ", "nv","  ", "energy","  ","xx","  ",&
		      "yy","  ","zz","  ","sp","  ","sm"		   
	end if	
 endif

	! $omp parallel do &
	! $OMP DEFAULT (SHARED)	

	do i=1,dimbse

		ec = eigv(stt(i,4),stt(i,3))
		ev = eigv(stt(i,4),stt(i,2))
		fdeh(i) = fermidisteh(ec,ev,temp)

		call optsp(eigv(stt(i,4),stt(i,2)),vector(stt(i,4),stt(i,2),:),&
		     eigv(stt(i,4),stt(i,3)),vector(stt(i,4),stt(i,3),:),&
		     kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),ffactor,sme,&
		     w90basis,nvec,rlat,rvec,hopmatrices,&
		     ihopmatrices,hrx(i),hry(i),hrz(i))
		     
		     hrsp(i) = (hrx(i)+cmplx(0.,1.)*hry(i))*(1.0/sqrt(2.))
		     hrsm(i) = (hrx(i)-cmplx(0.,1.)*hry(i))*(1.0/sqrt(2.))
		     
		vecres(i,1) = real(kpt(stt(i,4),1))
		vecres(i,2) = real(kpt(stt(i,4),2))
		vecres(i,3) = real(kpt(stt(i,4),3))
		vecres(i,4) = real(nocpk(stt(i,4)))
		vecres(i,5) = real(nocpk(stt(i,4))-nv+stt(i,3))
		vecres(i,6) = real(nocpk(stt(i,4))-nv+stt(i,2))		
		vecres(i,7) = real(eigv(stt(i,4),stt(i,3))-eigv(stt(i,4),stt(i,2)))		
		vecres(i,8) = real(hrx(i)*conjg(hrx(i)))*fdeh(i)		
		vecres(i,9) = real(hry(i)*conjg(hry(i)))*fdeh(i)		
		vecres(i,10) = real(hrz(i)*conjg(hrz(i)))*fdeh(i)
		vecres(i,11) = real(hrx(i)*conjg(hry(i)))*fdeh(i)
		vecres(i,12) = real(hrx(i)*conjg(hrz(i)))*fdeh(i)
		vecres(i,13) = real(hry(i)*conjg(hrz(i)))*fdeh(i)		 
		vecres(i,14) = real(hrsp(i)*conjg(hrsp(i)))*fdeh(i)
		vecres(i,15) = real(hrsm(i)*conjg(hrsm(i)))*fdeh(i)			     


	end do
	


	! $omp end parallel do
	
	call Bubblem(7,15,vecres, dimbse)
	
	do i=1,dimbse
	
		write(401,"(7F15.4)") vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),vecres(i,11),vecres(i,12),vecres(i,13)
		write(403,"(6F15.4)") vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),vecres(i,14),vecres(i,15)
		
		if (tmcoef) then
	        write(402,"(3F10.4,3I10.0,7F10.4)") vecres(i,1),vecres(i,2),vecres(i,3),int(vecres(i,4)),int(vecres(i,5)),&
						  int(vecres(i,6)),vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),&
						  vecres(i,11),vecres(i,12),vecres(i,13)
	        
	        write(404,"(3F10.4,3I10.0,6F10.4)") vecres(i,1),vecres(i,2),vecres(i,3),int(vecres(i,4)),int(vecres(i,5)),&
						  int(vecres(i,6)),vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),&
						  vecres(i,14),vecres(i,15)						  
						  
		end if
	
	end do	

	if(Process==0) write(300,*) 'single particle optics finished'
	if(Process==0) call flush(300)

	deallocate(vecres)
	
	!go to 789

	!allocate(hbse(dimbse,dimbse),W(dimbse))

	hbse=0.0

	!$omp parallel do 
!collapse(2)

	do i=1,dimbse



		do j=i,dimbse



  hbse(i,j)= matrizelbsetemp(coultype,ktol,w90basis,ediel,lc,ez,w1,r0,ngrid,rlat,stt(i,:),eigv(stt(i,4)&
  	    ,stt(i,3)),eigv(stt(i,4),stt(i,2)),vector(stt(i,4)&
            ,stt(i,3),:) ,vector(stt(i,4),stt(i,2),:),kpt(stt(i,4),:),stt(j,:),eigv(stt(j,4),stt(j,3))&
  	    ,eigv(stt(j,4),stt(j,2)) &
            ,vector(stt(j,4),stt(j,3),:),vector(stt(j,4),stt(j,2),:),kpt(stt(j,4),:),temp)



		end do



	end do

	!$omp end parallel do

	if(Process==0) write(300,*) 'exciton Hamiltonian matrix finished'
	if(Process==0) call flush(300)	



      	!LWORK = 2*dimbse+dimbse**2
      	!LIWORK = 3 + 5*dimbse
      	!LRWORK = 1 + 5*dimbse + 2*dimbse**2


!	allocate(WORK(LWORK))
!	allocate (IWORK(LIWORK))
!	allocate (RWORK(LRWORK))

 !     	CALL ZHEEVD( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK,& 
!	      RWORK, LRWORK, IWORK, LIWORK,INFO )
 !     	IF( INFO.GT. 0 ) THEN
  !      WRITE(*,*)'The algorithm failed to compute eigenvalues.'
  !       STOP
  !    	END IF   

	!allocate (RWORK(3*dimbse-2))
        !LWMAX = 2*dimbse-1
	!allocate(WORK(LWMAX))

     	LWORK = -1
      	CALL ZHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, INFO )
      	LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      	CALL ZHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK,INFO )
      	IF( INFO.GT. 0 ) THEN
        WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      	END IF   

	!allocate(pinter(dimbse),pintra(dimbse))

	!if (ntype .eq. 2) then

	!do i=1,dimbse

	!	call excitonil(i,w90basis,nc,nv,ngrid*ngrid,stt,dimbse,hbse(:,i),lcount,pinter(i),pintra(i))

	!	write(306,*) W(i),pinter(i),pintra(i)
	!	call flush(306)

	!end do 

	!else 

	!	continue

	!end if 
	
	if (bsewf) then
	      	do i=excwf0,excwff
      			call excwf2(outputfolder,ngkpt,kpt,nc,nv,nocpk,stt,W(i),i,hbse(:,i))
      		end do
	end if

	if(Process==0) write(300,*) 'exciton Hamiltonian diagonalized'
	if(Process==0) call flush(300)
	if(Process==0) write(300,*) "exciton ground state",W(1)

	deallocate(eigv,vector)
	deallocate(rvec,hopmatrices)
	deallocate(ihopmatrices,ffactor)
	deallocate(RWORK,WORK)
	!deallocate (IWORK)

	!allocate(actxx(dimbse),actyy(dimbse),actzz(dimbse),actxy(dimbse))
	!allocate(actxz(dimbse),actyz(dimbse))

 if(Process==0) then
	call dielbseptemp(nthreads,dimbse,dimbse,hbse,hrx,hrx,fdeh,actxx)
	write(300,*) 'xx tensor component'
	call flush(300)
		
	call dielbseptemp(nthreads,dimbse,dimbse,hbse,hry,hry,fdeh,actyy)
	write(300,*) 'yy tensor component'
	call flush(300)
		
	call dielbseptemp(nthreads,dimbse,dimbse,hbse,hrz,hrz,fdeh,actzz)
	write(300,*) 'zz tensor component'
	call flush(300)
	
	if (dtfull) then
	
	call dielbseptemp(nthreads,dimbse,dimbse,hbse,hrx,hry,fdeh,actxy)
	write(300,*) 'xy tensor component'
	call flush(300)
	
	call dielbseptemp(nthreads,dimbse,dimbse,hbse,hrx,hrz,fdeh,actxz)
	write(300,*) 'xz tensor component'
	call flush(300)
		
	call dielbseptemp(nthreads,dimbse,dimbse,hbse,hry,hrz,fdeh,actyz)
	write(300,*) 'yz tensor component'
	call flush(300)
	
	else
	
	actxy = 0.0000
	write(300,*) 'xy tensor component set to 0'
	call flush(300)
	
	actxz = 0.0000
	write(300,*) 'xz tensor component set to 0'
	call flush(300)
		
	actyz = 0.0000
	write(300,*) 'yz tensor component set to 0'
	call flush(300)
	
	end if
	
	if (cpol) then
	
	call dielbseptemp(nthreads,dimbse,dimbse,hbse,hrsp,hrsp,fdeh,actsp)
	write(300,*) 'sp polarization'
	call flush(300)
		
	call dielbseptemp(nthreads,dimbse,dimbse,hbse,hrsm,hrsm,fdeh,actsm)
	write(300,*) 'sm polarization'
	call flush(300)
	
	else 
	
	actsp = 0.0000
	write(300,*) 'sp polarization set to 0'
	call flush(300)
		
	actsm = 0.0000
	write(300,*) 'sm polarization set to 0'
	call flush(300)
	
	end if	

	write(300,*) 'optics finished'
	call flush(300)

	!$omp do ordered
	do i=1,dimbse
		!$omp ordered
		write(301,"(7F15.4)") W(i),actxx(i),actyy(i),actzz(i),actxy(i),actxz(i),actyz(i)
		call flush(301)
		
		if (cpol) then
		write(302,"(6F15.4)") W(i),actxx(i),actyy(i),actzz(i),actsp(i),actsm(i)
		call flush(302)		
		end if

		!$omp end ordered
	end do
	!$omp end do
 endif
	
	deallocate(hrx,hry,hrz,hrsp,hrsm)
	deallocate(actxx,actxy,actxz,actyy,actyz,actzz)
	deallocate(actsp,actsm)

	deallocate(hbse,W,stt,nocpk)
	deallocate(fdeh)	

	deallocate(kpt)

!789     continue


 if(Process==0) then
	call cpu_time(tf)
	call date_and_time(VALUES=values2)
	write(300,*)
	write(300,*) 'end','   ','day',values2(3),'',values2(5),'hours',values2(6),'min',values2(7),'seg'
	write(300,*)
	close(200)
	close(300)
	close(301)
	close(302)
	close(401)
	close(402)
	close(403)
	close(404)
 endif
end subroutine bsesolvertemp


!#endif


!INCLUDE "./subprograms/exciton_lifetime.f90"

subroutine lifetime(sysdim,numbse,ngrid,rlat,nc,nv,outputfolder)
	integer :: nc,nv
	double precision :: numbse
	integer,dimension(3) :: ngrid
	double precision,dimension(3,3) :: rlat
	character(len=70) :: outputfolder
	character(len=70) :: cflag
	character(len=5) :: sysdim	
	
	integer :: erro
	double precision,allocatable,dimension(:,:) :: refr
	double precision,allocatable,dimension(:,:) :: fosc
	integer :: dimbse,i,j
	
	double precision :: resx,resy,resz
	double precision :: refraux,foscaux
	double precision,dimension(4) :: lft
	
	
	!arquivo de entrada
 if(Process==0) then
	OPEN(UNIT=100, FILE=trim(outputfolder)//"/bse_opt_diel.dat",STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse_opt_diel input file"
	!OPEN(UNIT=544, FILE=trim(outputfolder)//"/bse_indice_refracao.dat",STATUS='old', IOSTAT=erro)
	!if (erro/=0) stop "Erro na abertura do arquivo de entrada indice de refracao"
	
	!arquivo de saida
	OPEN(UNIT=300, FILE=trim(outputfolder)//"/exciton_lifetime.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening exciton_lifetime output file"
 endif
	
	dimbse=ngrid(1)*ngrid(2)*ngrid(3)*nc*nv	
	
	allocate(fosc(dimbse,7),refr(int(numbse),7))
	
	!lendo os inputs
	if(Process==0) read(100,*) cflag
	call P_sendA(cflag)

 if(Process==0) then
	do i=1,dimbse
	 read(100,*) fosc(i,1),fosc(i,2),fosc(i,3),fosc(i,4),fosc(i,5),fosc(i,6),fosc(i,7)
	end do
 endif
    call P_sendR2(fosc,dimbse,7)

	!read(544,*) cflag
	
	
	!do i=1,int(numbse)
	
	 !read(544,*) refr(i,1),refr(i,2),refr(i,3),refr(i,4),refr(i,5),refr(i,6),refr(i,7)

	 
	!end do
	
	if(Process==0) write(300,*) "#energy(eV) sum-lft(s)  x-lft(s) y-lft(s) z-lft(s)"
	
	do i=1,dimbse
		!call interp1D(int(numbse),7,2,refr,fosc(i,1),resx)
		!call interp1D(int(numbse),7,3,refr,fosc(i,1),resy)	
		!call interp1D(int(numbse),7,4,refr,fosc(i,1),resz)	
		foscaux = fosc(i,2)+fosc(i,3)+fosc(i,4)
		call exclft(sysdim,ngrid,rlat,foscaux,fosc(i,1),lft(1))
		call exclft(sysdim,ngrid,rlat,fosc(i,2),fosc(i,1),lft(2))
		call exclft(sysdim,ngrid,rlat,fosc(i,3),fosc(i,1),lft(3))	
		call exclft(sysdim,ngrid,rlat,fosc(i,4),fosc(i,1),lft(4))		
		!write(*,*) refraux
		if(Process==0) write(300,"(1F15.4,4E15.4)") fosc(i,1),lft(1),lft(2),lft(3),lft(4)	
	end do
	
	deallocate(fosc,refr)

 if(Process==0) then
	close(100)
	!close(544)
	close(300)
 endif
end subroutine lifetime



subroutine lifetimepol(sysdim,numbse,ngrid,rlat,nc,nv,outputfolder)
	integer :: nc,nv
	double precision :: numbse
	integer,dimension(3) :: ngrid
	double precision,dimension(3,3) :: rlat
	character(len=70) :: outputfolder
	character(len=70) :: cflag
	character(len=5) :: sysdim	
	integer :: erro
	double precision,allocatable,dimension(:,:) :: refr
	double precision,allocatable,dimension(:,:) :: fosc
	integer :: dimbse,i,j
	double precision :: resx,resy,resz
	double precision :: refraux,foscaux
	double precision,dimension(6) :: lft
	
	
	!arquivo de entrada

 if(Process==0) then
	OPEN(UNIT=100, FILE=trim(outputfolder)//"/bse_opt_diel-pol.dat",STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse_opt_diel-pol input file"
	!OPEN(UNIT=544, FILE=trim(outputfolder)//"/bse_indice_refracao.dat",STATUS='old', IOSTAT=erro)
	!if (erro/=0) stop "Erro na abertura do arquivo de entrada indice de refracao"
	
	!arquivo de saida
	OPEN(UNIT=300, FILE=trim(outputfolder)//"/exciton_lifetime-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening exciton_lifetime-pol output file"
 endif
	
	dimbse=ngrid(1)*ngrid(2)*ngrid(3)*nc*nv	
	
	allocate(fosc(dimbse,6),refr(int(numbse),6))
	
	!lendo os inputs
	if(Process==0) read(100,*) cflag
	call P_sendA(cflag)

 if(Process==0) then
	do i=1,dimbse
	 read(100,*) fosc(i,1),fosc(i,2),fosc(i,3),fosc(i,4),fosc(i,5),fosc(i,6)
	end do
 endif
    call P_sendR2(fosc,dimbse,6)

	!read(544,*) cflag
	
	
	!do i=1,int(numbse)
	
	 !read(544,*) refr(i,1),refr(i,2),refr(i,3),refr(i,4),refr(i,5),refr(i,6),refr(i,7)

	 
	!end do
	
	if(Process==0) write(300,*) "#energy(eV) sum-lft(s)  x-lft(s) y-lft(s) z-lft(s) sp-lft(s) sm-lft(s)"
	
	do i=1,dimbse
		!call interp1D(int(numbse),7,2,refr,fosc(i,1),resx)
		!call interp1D(int(numbse),7,3,refr,fosc(i,1),resy)	
		!call interp1D(int(numbse),7,4,refr,fosc(i,1),resz)	
		foscaux = fosc(i,2)+fosc(i,3)+fosc(i,4)
		call exclft(sysdim,ngrid,rlat,foscaux,fosc(i,1),lft(1))
		call exclft(sysdim,ngrid,rlat,fosc(i,2),fosc(i,1),lft(2))
		call exclft(sysdim,ngrid,rlat,fosc(i,3),fosc(i,1),lft(3))	
		call exclft(sysdim,ngrid,rlat,fosc(i,4),fosc(i,1),lft(4))
		call exclft(sysdim,ngrid,rlat,fosc(i,5),fosc(i,1),lft(5))
		call exclft(sysdim,ngrid,rlat,fosc(i,6),fosc(i,1),lft(6))					
		!write(*,*) refraux
		if(Process==0) write(300,"(1F15.4,6E15.4)") fosc(i,1),lft(1),lft(2),lft(3),lft(4),lft(5),lft(6)	
	end do

	deallocate(fosc,refr)

 if(Process==0) then
	close(100)
	!close(544)
	close(300)
 endif
end subroutine lifetimepol



!INCLUDE "./subprograms/diel-pp.f90"
!ifort diel-pp.f90 -o diel-pp.x -qopenmp -mkl
subroutine spoptprop(numbse,outputfolder)
	integer :: i,j,erro,dimpp
	double precision,allocatable,dimension(:) :: rxx,rxy,rxz
	double precision,allocatable,dimension(:) :: ryy,ryz,rzz
	double precision,allocatable,dimension(:) :: ixx,ixy,ixz
	double precision,allocatable,dimension(:) :: iyy,iyz,izz
	double precision,allocatable,dimension(:) :: energy
	double precision :: flag
	character(len=2) :: aread
	double precision :: refr_xx,refr_xy,refr_xz,refr_yy,refr_yz,refr_zz !indice refracao
	double precision :: ext_xx,ext_xy,ext_xz,ext_yy,ext_yz,ext_zz !indice extincao
	double precision :: refl_xx,refl_xy,refl_xz,refl_yy,refl_yz,refl_zz !reflectibilidade
	double precision :: abs_xx,abs_xy,abs_xz,abs_yy,abs_yz,abs_zz !coeficiente de absorcao
	double precision :: els_xx,els_xy,els_xz,els_yy,els_yz,els_zz !energy loss function
	character(len=70) :: outputfolder
	double precision ::  numbse

	!call input_read

	dimpp=int(numbse)

	!arquivo de entrada

 if(Process==0) then
	OPEN(UNIT=100, FILE=trim(outputfolder)//"/sp_diel_xx.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel_xx input file"
	OPEN(UNIT=101, FILE=trim(outputfolder)//"/sp_diel_xy.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel_xy input file"
	OPEN(UNIT=102, FILE=trim(outputfolder)//"/sp_diel_xz.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel_xz input file"
	OPEN(UNIT=103, FILE=trim(outputfolder)//"/sp_diel_yy.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel_yy input file"
	OPEN(UNIT=104, FILE=trim(outputfolder)//"/sp_diel_yz.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel_yz input file"
	OPEN(UNIT=105, FILE=trim(outputfolder)//"/sp_diel_zz.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel_zz input file"

	!arquivos de saida
	OPEN(UNIT=200, FILE=trim(outputfolder)//"/sp_refractive_index.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_refractive_index output file"
	OPEN(UNIT=300, FILE=trim(outputfolder)//"/sp_extinction_coef.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_extinction_coef output file"
	OPEN(UNIT=400, FILE=trim(outputfolder)//"/sp_reflectibility.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_reflectibility output file"
	OPEN(UNIT=500, FILE=trim(outputfolder)//"/sp_absorption_coef.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_absorption_coef output file"
	OPEN(UNIT=600, FILE=trim(outputfolder)//"/sp_en_loss_func.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_en_loss_func output file"	

	read(100,*) aread
	read(101,*) aread
	read(102,*) aread
	read(103,*) aread
	read(104,*) aread
	read(105,*) aread
 endif
    call P_sendA(aread)
   
	allocate(energy(dimpp))
	allocate(rxx(dimpp),ryy(dimpp),rzz(dimpp))
	allocate(rxy(dimpp),rxz(dimpp),ryz(dimpp))
	allocate(ixx(dimpp),iyy(dimpp),izz(dimpp))
	allocate(ixy(dimpp),ixz(dimpp),iyz(dimpp))

 if(Process==0) then
	do i=1,dimpp
		read(100,*) energy(i),rxx(i),ixx(i)
		read(101,*) flag,rxy(i),ixy(i)
		read(102,*) flag,rxz(i),ixz(i)
		read(103,*) flag,ryy(i),iyy(i)
		read(104,*) flag,ryz(i),iyz(i)
		read(105,*) flag,rzz(i),izz(i)
	end do
 endif
    call P_sendR1(energy,dimpp)
    call P_sendR1(rxx,dimpp)
    call P_sendR1(ixx,dimpp)
    call P_sendR1(rxy,dimpp)
    call P_sendR1(ixy,dimpp)
    call P_sendR1(rxz,dimpp)
    call P_sendR1(ixz,dimpp)
    call P_sendR1(ryy,dimpp)
    call P_sendR1(iyy,dimpp)
    call P_sendR1(ryz,dimpp)
    call P_sendR1(iyz,dimpp)
    call P_sendR1(rzz,dimpp)
    call P_sendR1(izz,dimpp)
    call P_sendR(flag)


 if(Process==0) then
	write(200,*) "#","  ","energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(300,*) "#","  ","energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(400,*) "#","  ","energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(500,*) "#","  ","energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(600,*) "#","  ","energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"	
 endif
 
 
	do i=1,dimpp
		call refracao(rxx(i),ixx(i),refr_xx)
		call refracao(ryy(i),iyy(i),refr_yy)
		call refracao(rzz(i),izz(i),refr_zz)
		call refracao(rxy(i),ixy(i),refr_xy)
		call refracao(rxz(i),ixz(i),refr_xz)
		call refracao(ryz(i),iyz(i),refr_yz)
		if(Process==0) write(200,"(7F15.4)") energy(i),refr_xx,refr_yy,refr_zz,refr_xy,refr_xz,refr_yz
		call extincao(rxx(i),ixx(i),ext_xx)
		call extincao(ryy(i),iyy(i),ext_yy)
		call extincao(rzz(i),izz(i),ext_zz)
		call extincao(rxy(i),ixy(i),ext_xy)
		call extincao(rxz(i),ixz(i),ext_xz)
		call extincao(ryz(i),iyz(i),ext_yz)

		if(Process==0) write(300,"(7F15.4)") energy(i),ext_xx,ext_yy,ext_zz,ext_xy,ext_xz,ext_yz

		call reflectibilidade(rxx(i),ixx(i),refl_xx)
		call reflectibilidade(ryy(i),iyy(i),refl_yy)
		call reflectibilidade(rzz(i),izz(i),refl_zz)
		call reflectibilidade(rxy(i),ixy(i),refl_xy)
		call reflectibilidade(rxz(i),ixz(i),refl_xz)
		call reflectibilidade(ryz(i),iyz(i),refl_yz)

		if(Process==0) write(400,"(7F15.4)") energy(i),refl_xx,refl_yy,refl_zz,refl_xy,refl_xz,refl_yz

		call abscoef(rxx(i),ixx(i),energy(i),abs_xx)
		call abscoef(ryy(i),iyy(i),energy(i),abs_yy)
		call abscoef(rzz(i),izz(i),energy(i),abs_zz)
		call abscoef(rxy(i),ixy(i),energy(i),abs_xy)
		call abscoef(rxz(i),ixz(i),energy(i),abs_xz)
		call abscoef(ryz(i),iyz(i),energy(i),abs_yz)

		if(Process==0) write(500,"(7F15.4)") energy(i),abs_xx,abs_yy,abs_zz,abs_xy,abs_xz,abs_yz

		call enloss(rxx(i),ixx(i),els_xx)
		call enloss(ryy(i),iyy(i),els_yy)
		call enloss(rzz(i),izz(i),els_zz)
		call enloss(rxy(i),ixy(i),els_xy)
		call enloss(rxz(i),ixz(i),els_xz)
		call enloss(ryz(i),iyz(i),els_yz)

		if(Process==0) write(600,"(7F15.4)") energy(i),els_xx,els_yy,els_zz,els_xy,els_xz,els_yz		
	end do


	deallocate(energy)
	deallocate(rxx,ryy,rzz)
	deallocate(rxy,rxz,ryz)
	deallocate(ixx,iyy,izz)
	deallocate(ixy,ixz,iyz)

 if(Process==0) then
	close(100)
	close(101)
	close(102)
	close(103)
	close(104)
	close(105)
	close(200)
	close(300)
	close(400)
	close(500)
	close(600)	
 endif
end subroutine spoptprop


!INCLUDE "./subprograms/diel-pp-pol.f90"
!ifort diel-pp.f90 -o diel-pp.x -qopenmp -mkl
!bseoptproppol
subroutine spoptproppol(numbse,outputfolder)
	integer :: i,j,erro,dimpp
	double precision,allocatable,dimension(:) :: rxx,ryy,rzz
	double precision,allocatable,dimension(:) :: rsp,rsm
	double precision,allocatable,dimension(:) :: ixx,iyy,izz
	double precision,allocatable,dimension(:) :: isp,ism
	double precision,allocatable,dimension(:) :: energy
	double precision :: flag
	character(len=2) :: aread
	double precision :: refr_xx,refr_sp,refr_sm,refr_yy,refr_zz !indice refracao
	double precision :: ext_xx,ext_sp,ext_sm,ext_yy,ext_zz !indice extincao
	double precision :: refl_xx,refl_sp,refl_sm,refl_yy,refl_zz !reflectibilidade
	double precision :: abs_xx,abs_sp,abs_sm,abs_yy,abs_zz !coeficiente de absorcao
	double precision :: els_xx,els_sp,els_sm,els_yy,els_zz !energy loss function	
	character(len=70) :: outputfolder
	double precision ::  numbse

	!call input_read

	dimpp=int(numbse)

	!arquivo de entrada

 if(Process==0) then
	OPEN(UNIT=100, FILE=trim(outputfolder)//"/sp_diel-pol_x.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel-pol_x input file"
	OPEN(UNIT=101, FILE=trim(outputfolder)//"/sp_diel-pol_y.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel-pol_y input file"
	OPEN(UNIT=102, FILE=trim(outputfolder)//"/sp_diel-pol_z.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel-pol_z input file"
	OPEN(UNIT=103, FILE=trim(outputfolder)//"/sp_diel-pol_sp.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel-pol_sp input file"
	OPEN(UNIT=104, FILE=trim(outputfolder)//"/sp_diel-pol_sm.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel-pol_sm input file"
	!arquivos de saida
	OPEN(UNIT=200, FILE=trim(outputfolder)//"/sp_refractive_index-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_refractive_index-pol output file"
	OPEN(UNIT=300, FILE=trim(outputfolder)//"/sp_extinction_coef-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_extinction_coef-pol output file"
	OPEN(UNIT=400, FILE=trim(outputfolder)//"/sp_reflectibility-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_reflectibility-pol output file"
	OPEN(UNIT=500, FILE=trim(outputfolder)//"/sp_absorption_coef-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_absorption_coef-pol output file"
	OPEN(UNIT=600, FILE=trim(outputfolder)//"/sp_en_loss_func-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_en_loss_func-pol output file"		
	read(100,*) aread
	read(101,*) aread
	read(102,*) aread
	read(103,*) aread
	read(104,*) aread
 endif
    call P_sendA(aread)

	allocate(energy(dimpp))
	allocate(rxx(dimpp),ryy(dimpp),rzz(dimpp))
	allocate(rsp(dimpp),rsm(dimpp))
	allocate(ixx(dimpp),iyy(dimpp),izz(dimpp))
	allocate(isp(dimpp),ism(dimpp))

 if(Process==0) then
	do i=1,dimpp
		read(100,*) energy(i),rxx(i),ixx(i)
		read(101,*) flag,ryy(i),iyy(i)
		read(102,*) flag,rzz(i),izz(i)
		read(103,*) flag,rsp(i),isp(i)
		read(104,*) flag,rsm(i),ism(i)
	end do
 endif
    call P_sendR1(energy,dimpp)
    call P_sendR1(rxx,dimpp)
    call P_sendR1(ixx,dimpp)
    call P_sendR1(ryy,dimpp)
    call P_sendR1(iyy,dimpp)
    call P_sendR1(rzz,dimpp)
    call P_sendR1(izz,dimpp)
    call P_sendR1(rsp,dimpp)
    call P_sendR1(isp,dimpp)
    call P_sendR1(rsm,dimpp)
    call P_sendR1(ism,dimpp)
    call P_sendR(flag)
 

 if(Process==0) then
	write(200,*) "#","  ","energy","  ","x","  ","y","  ","z","  ","sp","  ","sm"
	write(300,*) "#","  ","energy","  ","x","  ","y","  ","z","  ","sp","  ","sm"
	write(400,*) "#","  ","energy","  ","x","  ","y","  ","z","  ","sp","  ","sm"
	write(500,*) "#","  ","energy","  ","x","  ","y","  ","z","  ","sp","  ","sm"
	write(600,*) "#","  ","energy","  ","x","  ","y","  ","z","  ","sp","  ","sm"	
 endif

	do i=1,dimpp
		call refracao(rxx(i),ixx(i),refr_xx)
		call refracao(ryy(i),iyy(i),refr_yy)
		call refracao(rzz(i),izz(i),refr_zz)
		call refracao(rsp(i),isp(i),refr_sp)
		call refracao(rsm(i),ism(i),refr_sm)

		if(Process==0) write(200,"(6F15.4)") energy(i),refr_xx,refr_yy,refr_zz,refr_sp,refr_sm

		call extincao(rxx(i),ixx(i),ext_xx)
		call extincao(ryy(i),iyy(i),ext_yy)
		call extincao(rzz(i),izz(i),ext_zz)
		call extincao(rsp(i),isp(i),ext_sp)
		call extincao(rsm(i),ism(i),ext_sm)

		if(Process==0) write(300,"(6F15.4)") energy(i),ext_xx,ext_yy,ext_zz,ext_sp,ext_sm

		call reflectibilidade(rxx(i),ixx(i),refl_xx)
		call reflectibilidade(ryy(i),iyy(i),refl_yy)
		call reflectibilidade(rzz(i),izz(i),refl_zz)
		call reflectibilidade(rsp(i),isp(i),refl_sp)
		call reflectibilidade(rsm(i),ism(i),refl_sm)


		if(Process==0) write(400,"(6F15.4)") energy(i),refl_xx,refl_yy,refl_zz,refl_sp,refl_sm

		call abscoef(rxx(i),ixx(i),energy(i),abs_xx)
		call abscoef(ryy(i),iyy(i),energy(i),abs_yy)
		call abscoef(rzz(i),izz(i),energy(i),abs_zz)
		call abscoef(rsp(i),isp(i),energy(i),abs_sp)
		call abscoef(rsm(i),ism(i),energy(i),abs_sm)

		if(Process==0) write(500,"(6F15.4)") energy(i),abs_xx,abs_yy,abs_zz,abs_sp,abs_sm
		
		call enloss(rxx(i),ixx(i),els_xx)
		call enloss(ryy(i),iyy(i),els_yy)
		call enloss(rzz(i),izz(i),els_zz)
		call enloss(rsp(i),isp(i),els_sp)
		call enloss(rsm(i),ism(i),els_sm)

		if(Process==0) write(600,"(6F15.4)") energy(i),els_xx,els_yy,els_zz,els_sp,els_sm		
	end do


	deallocate(energy)
	deallocate(rxx,ryy,rzz)
	deallocate(rsp,rsm)
	deallocate(ixx,iyy,izz)
	deallocate(isp,ism)

 if(Process==0) then
	close(100)
	close(101)
	close(102)
	close(103)
	close(104)
	close(200)
	close(300)
	close(400)
	close(500)
	close(600)
 endif
end subroutine spoptproppol



!INCLUDE "./subprograms/diel-pp-bse.f90"
!ifort diel-pp-bse.f90 -o diel-pp-bse.x -qopenmp -mkl
subroutine bseoptprop(numbse,outputfolder)
	integer :: i,j,erro,dimpp
	double precision,allocatable,dimension(:) :: rxx,rxy,rxz
	double precision,allocatable,dimension(:) :: ryy,ryz,rzz
	double precision,allocatable,dimension(:) :: ixx,ixy,ixz
	double precision,allocatable,dimension(:) :: iyy,iyz,izz
	double precision,allocatable,dimension(:) :: energy
	double precision :: flag
	character(len=2) :: aread
	double precision :: refr_xx,refr_xy,refr_xz,refr_yy,refr_yz,refr_zz !indice refracao
	double precision :: ext_xx,ext_xy,ext_xz,ext_yy,ext_yz,ext_zz !indice extincao
	double precision :: refl_xx,refl_xy,refl_xz,refl_yy,refl_yz,refl_zz !reflectibilidade
	double precision :: abs_xx,abs_xy,abs_xz,abs_yy,abs_yz,abs_zz !coeficiente de absorcao
	double precision :: els_xx,els_xy,els_xz,els_yy,els_yz,els_zz !energy loss function	
	character(len=70) :: outputfolder
	double precision ::  numbse

	!call input_read

	dimpp=int(numbse)

	!arquivo de entrada

 if(Process==0) then
	OPEN(UNIT=100, FILE=trim(outputfolder)//"/bse_diel_xx.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel_xx input file"
	OPEN(UNIT=101, FILE=trim(outputfolder)//"/bse_diel_xy.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel_xy input file"
	OPEN(UNIT=102, FILE=trim(outputfolder)//"/bse_diel_xz.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel_xz input file"
	OPEN(UNIT=103, FILE=trim(outputfolder)//"/bse_diel_yy.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel_yy input file"
	OPEN(UNIT=104, FILE=trim(outputfolder)//"/bse_diel_yz.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel_yz input file"
	OPEN(UNIT=105, FILE=trim(outputfolder)//"/bse_diel_zz.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel_zz input file"

	!arquivos de saida
	OPEN(UNIT=200, FILE=trim(outputfolder)//"/bse_refractive_index.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_refractive_index output file"
	OPEN(UNIT=300, FILE=trim(outputfolder)//"/bse_extinction_coef.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_extinction_coef output file"
	OPEN(UNIT=400, FILE=trim(outputfolder)//"/bse_reflectibility.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_reflectibility output file"
	OPEN(UNIT=500, FILE=trim(outputfolder)//"/bse_absorption_coef.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_absorption_coef output file"
	OPEN(UNIT=600, FILE=trim(outputfolder)//"/bse_en_loss_func.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_en_loss_func output file"	
	read(100,*) aread
	read(101,*) aread
	read(102,*) aread
	read(103,*) aread
	read(104,*) aread
	read(105,*) aread
 endif
    call P_sendA(aread)

	allocate(energy(dimpp))
	allocate(rxx(dimpp),ryy(dimpp),rzz(dimpp))
	allocate(rxy(dimpp),rxz(dimpp),ryz(dimpp))
	allocate(ixx(dimpp),iyy(dimpp),izz(dimpp))
	allocate(ixy(dimpp),ixz(dimpp),iyz(dimpp))

 if(Process==0) then
	do i=1,dimpp
		read(100,*) energy(i),rxx(i),ixx(i)
		read(101,*) flag,rxy(i),ixy(i)
		read(102,*) flag,rxz(i),ixz(i)
		read(103,*) flag,ryy(i),iyy(i)
		read(104,*) flag,ryz(i),iyz(i)
		read(105,*) flag,rzz(i),izz(i)
	end do
 endif
    call P_sendR1(energy,dimpp)
    call P_sendR1(rxx,dimpp)
    call P_sendR1(ixx,dimpp)
    call P_sendR1(rxy,dimpp)
    call P_sendR1(ixy,dimpp)
    call P_sendR1(rxz,dimpp)
    call P_sendR1(ixz,dimpp)
    call P_sendR1(ryy,dimpp)
    call P_sendR1(iyy,dimpp)
    call P_sendR1(ryz,dimpp)
    call P_sendR1(iyz,dimpp)
    call P_sendR1(rzz,dimpp)
    call P_sendR1(izz,dimpp)
    call P_sendR(flag)


 if(Process==0) then
	write(200,*) "#","  ","energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(300,*) "#","  ","energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(400,*) "#","  ","energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(500,*) "#","  ","energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(600,*) "#","  ","energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"	
 endif

	do i=1,dimpp

		call refracao(rxx(i),ixx(i),refr_xx)
		call refracao(ryy(i),iyy(i),refr_yy)
		call refracao(rzz(i),izz(i),refr_zz)
		call refracao(rxy(i),ixy(i),refr_xy)
		call refracao(rxz(i),ixz(i),refr_xz)
		call refracao(ryz(i),iyz(i),refr_yz)

		if(Process==0) write(200,"(7F15.4)") energy(i),refr_xx,refr_yy,refr_zz,refr_xy,refr_xz,refr_yz

		call extincao(rxx(i),ixx(i),ext_xx)
		call extincao(ryy(i),iyy(i),ext_yy)
		call extincao(rzz(i),izz(i),ext_zz)
		call extincao(rxy(i),ixy(i),ext_xy)
		call extincao(rxz(i),ixz(i),ext_xz)
		call extincao(ryz(i),iyz(i),ext_yz)

		if(Process==0) write(300,"(7F15.4)") energy(i),ext_xx,ext_yy,ext_zz,ext_xy,ext_xz,ext_yz

		call reflectibilidade(rxx(i),ixx(i),refl_xx)
		call reflectibilidade(ryy(i),iyy(i),refl_yy)
		call reflectibilidade(rzz(i),izz(i),refl_zz)
		call reflectibilidade(rxy(i),ixy(i),refl_xy)
		call reflectibilidade(rxz(i),ixz(i),refl_xz)
		call reflectibilidade(ryz(i),iyz(i),refl_yz)

		if(Process==0) write(400,"(7F15.4)") energy(i),refl_xx,refl_yy,refl_zz,refl_xy,refl_xz,refl_yz

		call abscoef(rxx(i),ixx(i),energy(i),abs_xx)
		call abscoef(ryy(i),iyy(i),energy(i),abs_yy)
		call abscoef(rzz(i),izz(i),energy(i),abs_zz)
		call abscoef(rxy(i),ixy(i),energy(i),abs_xy)
		call abscoef(rxz(i),ixz(i),energy(i),abs_xz)
		call abscoef(ryz(i),iyz(i),energy(i),abs_yz)

		if(Process==0) write(500,"(7F15.4)") energy(i),abs_xx,abs_yy,abs_zz,abs_xy,abs_xz,abs_yz
		
		call enloss(rxx(i),ixx(i),els_xx)
		call enloss(ryy(i),iyy(i),els_yy)
		call enloss(rzz(i),izz(i),els_zz)
		call enloss(rxy(i),ixy(i),els_xy)
		call enloss(rxz(i),ixz(i),els_xz)
		call enloss(ryz(i),iyz(i),els_yz)

		if(Process==0) write(600,"(7F15.4)") energy(i),els_xx,els_yy,els_zz,els_xy,els_xz,els_yz		
	end do


	deallocate(energy)
	deallocate(rxx,ryy,rzz)
	deallocate(rxy,rxz,ryz)
	deallocate(ixx,iyy,izz)
	deallocate(ixy,ixz,iyz)

 if(Process==0) then
	close(100)
	close(101)
	close(102)
	close(103)
	close(104)
	close(105)
	close(200)
	close(300)
	close(400)
	close(500)
	close(600)
 endif
end subroutine bseoptprop


!INCLUDE "./subprograms/diel-pp-bse-pol.f90"
!ifort diel-pp.f90 -o diel-pp.x -qopenmp -mkl
!bseoptproppol
subroutine bseoptproppol(numbse,outputfolder)
	integer :: i,j,erro,dimpp
	double precision,allocatable,dimension(:) :: rxx,ryy,rzz
	double precision,allocatable,dimension(:) :: rsp,rsm
	double precision,allocatable,dimension(:) :: ixx,iyy,izz
	double precision,allocatable,dimension(:) :: isp,ism
	double precision,allocatable,dimension(:) :: energy
	double precision :: flag
	character(len=2) :: aread
	double precision :: refr_xx,refr_sp,refr_sm,refr_yy,refr_zz !indice refracao
	double precision :: ext_xx,ext_sp,ext_sm,ext_yy,ext_zz !indice extincao
	double precision :: refl_xx,refl_sp,refl_sm,refl_yy,refl_zz !reflectibilidade
	double precision :: abs_xx,abs_sp,abs_sm,abs_yy,abs_zz !coeficiente de absorcao
	double precision :: els_xx,els_sp,els_sm,els_yy,els_zz !energy loss function	
	character(len=70) :: outputfolder
	double precision ::  numbse

	!call input_read
	dimpp=int(numbse)
	!arquivo de entrada

 if(Process==0) then
	OPEN(UNIT=100, FILE=trim(outputfolder)//"/bse_diel-pol_x.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel-pol_x input file"
	OPEN(UNIT=101, FILE=trim(outputfolder)//"/bse_diel-pol_y.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel-pol_y input file"
	OPEN(UNIT=102, FILE=trim(outputfolder)//"/bse_diel-pol_z.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel-pol_z input file"
	OPEN(UNIT=103, FILE=trim(outputfolder)//"/bse_diel-pol_sp.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel-pol_sp input file"
	OPEN(UNIT=104, FILE=trim(outputfolder)//"/bse_diel-pol_sm.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel-pol_sm input file"
	!arquivos de saida
	OPEN(UNIT=200, FILE=trim(outputfolder)//"/bse_refractive_index-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_refractive_index-pol output file"
	OPEN(UNIT=300, FILE=trim(outputfolder)//"/bse_extinction_coef-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_extinction_coef-pol output file"
	OPEN(UNIT=400, FILE=trim(outputfolder)//"/bse_reflectibility-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_reflectibility-pol output file"
	OPEN(UNIT=500, FILE=trim(outputfolder)//"/bse_absorption_coef-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_absorption_coef-pol output file"
	OPEN(UNIT=600, FILE=trim(outputfolder)//"/bse_en_loss_func-pol.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_en_loss_func-pol output file"	
	read(100,*) aread
	read(101,*) aread
	read(102,*) aread
	read(103,*) aread
	read(104,*) aread
 endif
    call P_sendA(aread)

	allocate(energy(dimpp))
	allocate(rxx(dimpp),ryy(dimpp),rzz(dimpp))
	allocate(rsp(dimpp),rsm(dimpp))
	allocate(ixx(dimpp),iyy(dimpp),izz(dimpp))
	allocate(isp(dimpp),ism(dimpp))

 if(Process==0) then
	do i=1,dimpp
		read(100,*) energy(i),rxx(i),ixx(i)
		read(101,*) flag,ryy(i),iyy(i)
		read(102,*) flag,rzz(i),izz(i)
		read(103,*) flag,rsp(i),isp(i)
		read(104,*) flag,rsm(i),ism(i)
	end do
 endif
    call P_sendR1(energy,dimpp)
    call P_sendR1(rxx,dimpp)
    call P_sendR1(ixx,dimpp)
    call P_sendR1(ryy,dimpp)
    call P_sendR1(iyy,dimpp)
    call P_sendR1(rzz,dimpp)
    call P_sendR1(izz,dimpp)
    call P_sendR1(rsp,dimpp)
    call P_sendR1(isp,dimpp)
    call P_sendR1(rsm,dimpp)
    call P_sendR1(ism,dimpp)
    call P_sendR(flag)
 
 if(Process==0) then
	write(200,*) "#","  ","energy","  ","x","  ","y","  ","z","  ","sp","  ","sm"
	write(300,*) "#","  ","energy","  ","x","  ","y","  ","z","  ","sp","  ","sm"
	write(400,*) "#","  ","energy","  ","x","  ","y","  ","z","  ","sp","  ","sm"
	write(500,*) "#","  ","energy","  ","x","  ","y","  ","z","  ","sp","  ","sm"
	write(600,*) "#","  ","energy","  ","x","  ","y","  ","z","  ","sp","  ","sm"	
 endif

	do i=1,dimpp
		call refracao(rxx(i),ixx(i),refr_xx)
		call refracao(ryy(i),iyy(i),refr_yy)
		call refracao(rzz(i),izz(i),refr_zz)
		call refracao(rsp(i),isp(i),refr_sp)
		call refracao(rsm(i),ism(i),refr_sm)

		if(Process==0) write(200,"(6F15.4)") energy(i),refr_xx,refr_yy,refr_zz,refr_sp,refr_sm

		call extincao(rxx(i),ixx(i),ext_xx)
		call extincao(ryy(i),iyy(i),ext_yy)
		call extincao(rzz(i),izz(i),ext_zz)
		call extincao(rsp(i),isp(i),ext_sp)
		call extincao(rsm(i),ism(i),ext_sm)

		if(Process==0) write(300,"(6F15.4)") energy(i),ext_xx,ext_yy,ext_zz,ext_sp,ext_sm

		call reflectibilidade(rxx(i),ixx(i),refl_xx)
		call reflectibilidade(ryy(i),iyy(i),refl_yy)
		call reflectibilidade(rzz(i),izz(i),refl_zz)
		call reflectibilidade(rsp(i),isp(i),refl_sp)
		call reflectibilidade(rsm(i),ism(i),refl_sm)

		if(Process==0) write(400,"(6F15.4)") energy(i),refl_xx,refl_yy,refl_zz,refl_sp,refl_sm

		call abscoef(rxx(i),ixx(i),energy(i),abs_xx)
		call abscoef(ryy(i),iyy(i),energy(i),abs_yy)
		call abscoef(rzz(i),izz(i),energy(i),abs_zz)
		call abscoef(rsp(i),isp(i),energy(i),abs_sp)
		call abscoef(rsm(i),ism(i),energy(i),abs_sm)

		if(Process==0) write(500,"(6F15.4)") energy(i),abs_xx,abs_yy,abs_zz,abs_sp,abs_sm
		
		call enloss(rxx(i),ixx(i),els_xx)
		call enloss(ryy(i),iyy(i),els_yy)
		call enloss(rzz(i),izz(i),els_zz)
		call enloss(rsp(i),isp(i),els_sp)
		call enloss(rsm(i),ism(i),els_sm)

		if(Process==0) write(600,"(6F15.4)") energy(i),els_xx,els_yy,els_zz,els_sp,els_sm		
	end do

	deallocate(energy)
	deallocate(rxx,ryy,rzz)
	deallocate(rsp,rsm)
	deallocate(ixx,iyy,izz)
	deallocate(isp,ism)

 if(Process==0) then
	close(100)
	close(101)
	close(102)
	close(103)
	close(104)
	close(200)
	close(300)
	close(400)
	close(500)
	close(600)
 endif
end subroutine bseoptproppol



!INCLUDE "./subprograms/sp_diel-tool.f90"
!ifort rpa_diel-tool.f90 -o spec_rpa_diel.x -qopenmp -mkl
subroutine spdielraw(nthreads,outputfolder,renorm,params,ngrid,nc,nv,ebse0,ebsef,numbse,sme)
	integer :: i,j,erro
	integer :: dimbse     
	double precision,allocatable,dimension(:) :: esp
	double precision,allocatable,dimension(:) :: exxf,exyf,exzf 
	double precision,allocatable,dimension(:) :: eyyf,eyzf,ezzf 
	double precision :: flag
	double precision :: vc
	double precision :: elux
	double precision :: rxx,rxy,rxz
	double precision :: ryy,ryz,rzz
	double precision :: ixx,ixy,ixz
	double precision :: iyy,iyz,izz
	character(len=2) :: aread
	double precision :: spinf
	double precision,allocatable,dimension(:,:) :: dielfxx,dielfyy,dielfzz
	double precision,allocatable,dimension(:,:) :: dielfxy,dielfxz,dielfyz	
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
	logical :: renorm	


	!call input_read

	dimbse=ngrid(1)*ngrid(2)*ngrid(3)*nc*nv

	!arquivo de entrada

 if(Process==0) then
	OPEN(UNIT=100, FILE=trim(outputfolder)//"/sp_optics.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening sp_optics input file"
 endif
!	call OMP_SET_NUM_THREADS(nthreads)

	call hamiltonian_input_read(200,params)
	!ediel(2) = edielh

	if ( systype .eq. "NP" ) then

	spinf = 2.0

	else

	spinf = 1.0

	end if

	!arquivos de saida
 if(Process==0) then
	OPEN(UNIT=300, FILE=trim(outputfolder)//"/sp_diel_xx.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel_xx output file"
	OPEN(UNIT=301, FILE=trim(outputfolder)//"/sp_diel_xy.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel_xy output file"
	OPEN(UNIT=302, FILE=trim(outputfolder)//"/sp_diel_xz.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel_xz output file"
	OPEN(UNIT=303, FILE=trim(outputfolder)//"/sp_diel_yy.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel_yy output file"
	OPEN(UNIT=304, FILE=trim(outputfolder)//"/sp_diel_yz.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel_yz output file"
	OPEN(UNIT=305, FILE=trim(outputfolder)//"/sp_diel_zz.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel_zz output file"
 endif

	!parametros do calculo


	!termino parametros calculo


	!call vcell(systype,rlat,vc)
	call vcell3D(rlat,vc)

	allocate(esp(dimbse))
	allocate(exxf(dimbse),exyf(dimbse),exzf(dimbse),eyyf(dimbse))
	allocate(eyzf(dimbse),ezzf(dimbse))

	if(Process==0) read(100,*) aread
	call P_sendA(aread)

 if(Process==0) then
	do j=1,dimbse
	 read(100,*) esp(j),exxf(j),eyyf(j),ezzf(j),exyf(j),exzf(j),eyzf(j)
	end do
 endif
    call P_sendR1(esp,dimbse)
    call P_sendR1(exxf,dimbse)
    call P_sendR1(eyyf,dimbse)
    call P_sendR1(ezzf,dimbse)
    call P_sendR1(exyf,dimbse)
    call P_sendR1(exzf,dimbse)
    call P_sendR1(eyzf,dimbse)

	allocate(dielfxx(int(numbse),3),dielfyy(int(numbse),3),dielfzz(int(numbse),3))
	allocate(dielfxy(int(numbse),3),dielfxz(int(numbse),3),dielfyz(int(numbse),3))

 if(Process==0) then
	write(300,*) "#","  ","energy","  ","real","  ","imag"
	write(301,*) "#","  ","energy","  ","real","  ","imag"
	write(302,*) "#","  ","energy","  ","real","  ","imag"
	write(303,*) "#","  ","energy","  ","real","  ","imag"
	write(304,*) "#","  ","energy","  ","real","  ","imag"
	write(305,*) "#","  ","energy","  ","real","  ","imag"
 endif
 
	!$omp do 
	!ordered
	do j=1,int(numbse)

		elux= ebse0 + dble(((ebsef-ebse0)*(j-1))/(numbse-1.))

		call optdiel(1,vc,dimbse,ngrid,elux,esp,exxf,sme,rxx,ixx)
		call optdiel(1,vc,dimbse,ngrid,elux,esp,eyyf,sme,ryy,iyy)
		call optdiel(1,vc,dimbse,ngrid,elux,esp,ezzf,sme,rzz,izz)

		call optdiel(0,vc,dimbse,ngrid,elux,esp,exyf,sme,rxy,ixy)
		call optdiel(0,vc,dimbse,ngrid,elux,esp,exzf,sme,rxz,ixz)
		call optdiel(0,vc,dimbse,ngrid,elux,esp,eyzf,sme,ryz,iyz)
		
		dielfxx(j,1) = elux
		dielfxx(j,2) = rxx
		dielfxx(j,3) = ixx
		
		dielfyy(j,1) = elux
		dielfyy(j,2) = ryy
		dielfyy(j,3) = iyy
		
		dielfzz(j,1) = elux
		dielfzz(j,2) = rzz
		dielfzz(j,3) = izz
		
		dielfxy(j,1) = elux
		dielfxy(j,2) = rxy
		dielfxy(j,3) = ixy
		
		dielfxz(j,1) = elux
		dielfxz(j,2) = rxz
		dielfxz(j,3) = ixz
		
		dielfyz(j,1) = elux
		dielfyz(j,2) = ryz
		dielfyz(j,3) = iyz			

		!!$omp ordered
		!write(300,*) real(elux),real(spinf*rxx),real(spinf*ixx)
		!write(301,*) real(elux),real(spinf*rxy),real(spinf*ixy)
		!write(302,*) real(elux),real(spinf*rxz),real(spinf*ixz)
		!write(303,*) real(elux),real(spinf*ryy),real(spinf*iyy)
		!write(304,*) real(elux),real(spinf*ryz),real(spinf*iyz)
		!write(305,*) real(elux),real(spinf*rzz),real(spinf*izz)
	
		!call flush(300)
		!!$omp end ordered
	end do
	!$omp end  do
	
	if (renorm) then
	 call imagrenorm(dimbse,int(numbse),esp,dielfxx)
	 call imagrenorm(dimbse,int(numbse),esp,dielfyy)
	 call imagrenorm(dimbse,int(numbse),esp,dielfzz)
	 call imagrenorm(dimbse,int(numbse),esp,dielfxy)
	 call imagrenorm(dimbse,int(numbse),esp,dielfxz)
	 call imagrenorm(dimbse,int(numbse),esp,dielfyz)
	end if		

 if(Process==0) then
	do i=1,int(numbse)
		write(300,*) dielfxx(i,1),real(spinf*dielfxx(i,2)),real(spinf*dielfxx(i,3))
		write(301,*) dielfxy(i,1),real(spinf*dielfxy(i,2)),real(spinf*dielfxy(i,3))
		write(302,*) dielfxz(i,1),real(spinf*dielfxz(i,2)),real(spinf*dielfxz(i,3))
		write(303,*) dielfyy(i,1),real(spinf*dielfyy(i,2)),real(spinf*dielfyy(i,3))
		write(304,*) dielfyz(i,1),real(spinf*dielfyz(i,2)),real(spinf*dielfyz(i,3))
		write(305,*) dielfzz(i,1),real(spinf*dielfzz(i,2)),real(spinf*dielfzz(i,3))
														
		!write(301,*) elux,real(spinf*rxy),real(spinf*ixy)
		!write(302,*) elux,real(spinf*rxz),real(spinf*ixz)
		!write(303,*) elux,real(spinf*ryy),real(spinf*iyy)
		!write(304,*) elux,real(spinf*ryz),real(spinf*iyz)
		!write(305,*) elux,real(spinf*rzz),real(spinf*izz)
	end do	
 endif

	deallocate(esp)
	deallocate(exxf,exyf,exzf,eyyf)
	deallocate(eyzf,ezzf)
	deallocate(rvec,hopmatrices)
	deallocate(ihopmatrices,ffactor)
	deallocate(dielfxx,dielfyy,dielfzz)
	deallocate(dielfxy,dielfxz,dielfyz)	

 if(Process==0) then
	close(100)
	close(200)
	close(300)
	close(301)
	close(302)
	close(303)
	close(304)
	close(305)	
 endif
end subroutine spdielraw







!INCLUDE "./subprograms/sp_diel-tool-pol.f90"
!ifort rpa_diel-tool.f90 -o spec_rpa_diel.x -qopenmp -mkl
subroutine spdielrawpol(nthreads,outputfolder,renorm,params,ngrid,nc,nv,ebse0,ebsef,numbse,sme)
	integer :: i,j,erro
	integer :: dimbse     
	double precision,allocatable,dimension(:) :: esp
	double precision,allocatable,dimension(:) :: exxf,eyyf,ezzf 
	double precision,allocatable,dimension(:) :: espf,esmf 
	double precision :: flag
	double precision :: vc
	double precision :: elux
	double precision :: rxx,ryy,rzz
	double precision :: rsp,rsm
	double precision :: ixx,iyy,izz
	double precision :: isp,ism
	double precision,allocatable,dimension(:,:) :: dielfxx,dielfyy,dielfzz
	double precision,allocatable,dimension(:,:) :: dielfsp,dielfsm	
	character(len=2) :: aread
	double precision :: spinf
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
	logical :: renorm

	!call input_read

	dimbse=ngrid(1)*ngrid(2)*ngrid(3)*nc*nv

	!arquivo de entrada

 if(Process==0) then
	OPEN(UNIT=100, FILE=trim(outputfolder)//"/sp_optics-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening sp_optics-pol input file"
 endif
 
!	call OMP_SET_NUM_THREADS(nthreads)

	call hamiltonian_input_read(200,params)
	!ediel(2) = edielh

	if ( systype .eq. "NP" ) then

	spinf = 2.0

	else

	spinf = 1.0

	end if

	!arquivos de saida

 if(Process==0) then
	OPEN(UNIT=300, FILE=trim(outputfolder)//"/sp_diel-pol_x.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel-pol_x output file"
	OPEN(UNIT=301, FILE=trim(outputfolder)//"/sp_diel-pol_y.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel-pol_y output file"
	OPEN(UNIT=302, FILE=trim(outputfolder)//"/sp_diel-pol_z.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel-pol_z output file"
	OPEN(UNIT=303, FILE=trim(outputfolder)//"/sp_diel-pol_sp.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel-pol_sp output file"
	OPEN(UNIT=304, FILE=trim(outputfolder)//"/sp_diel-pol_sm.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel-pol_sm output file"
 endif


	!parametros do calculo
	!termino parametros calculo
	!call vcell(systype,rlat,vc)
	call vcell3D(rlat,vc)

	allocate(esp(dimbse))
	allocate(exxf(dimbse),eyyf(dimbse),ezzf(dimbse))
	allocate(espf(dimbse),esmf(dimbse))

	if(Process==0) read(100,*) aread
	call P_sendA(aread)
		
 if(Process==0) then
	do j=1,dimbse
	 read(100,*) esp(j),exxf(j),eyyf(j),ezzf(j),espf(j),esmf(j)
	end do
 endif
    call P_sendR1(esp,dimbse)
    call P_sendR1(exxf,dimbse)
    call P_sendR1(eyyf,dimbse)
    call P_sendR1(ezzf,dimbse)
    call P_sendR1(espf,dimbse)
    call P_sendR1(esmf,dimbse)
    
    
	allocate(dielfxx(int(numbse),3),dielfyy(int(numbse),3),dielfzz(int(numbse),3))
	allocate(dielfsp(int(numbse),3),dielfsm(int(numbse),3))

 if(Process==0) then
	write(300,*) "#","  ","energy","  ","real","  ","imag"
	write(301,*) "#","  ","energy","  ","real","  ","imag"
	write(302,*) "#","  ","energy","  ","real","  ","imag"
	write(303,*) "#","  ","energy","  ","real","  ","imag"
	write(304,*) "#","  ","energy","  ","real","  ","imag"
 endif

	!$omp do 
	!ordered
	do j=1,int(numbse)

		elux= ebse0 + dble(((ebsef-ebse0)*(j-1))/(numbse-1.))

		call optdiel(1,vc,dimbse,ngrid,elux,esp,exxf,sme,rxx,ixx)
		call optdiel(1,vc,dimbse,ngrid,elux,esp,eyyf,sme,ryy,iyy)
		call optdiel(1,vc,dimbse,ngrid,elux,esp,ezzf,sme,rzz,izz)
		
		call optdiel(1,vc,dimbse,ngrid,elux,esp,espf,sme,rsp,isp)
		call optdiel(1,vc,dimbse,ngrid,elux,esp,esmf,sme,rsm,ism)
		
		dielfxx(j,1) = elux
		dielfxx(j,2) = rxx
		dielfxx(j,3) = ixx
		
		dielfyy(j,1) = elux
		dielfyy(j,2) = ryy
		dielfyy(j,3) = iyy
		
		dielfzz(j,1) = elux
		dielfzz(j,2) = rzz
		dielfzz(j,3) = izz
		
		dielfsp(j,1) = elux
		dielfsp(j,2) = rsp
		dielfsp(j,3) = isp
		
		dielfsm(j,1) = elux
		dielfsm(j,2) = rsm
		dielfsm(j,3) = ism		


		!!$omp ordered
		!write(300,*) real(elux),real(spinf*rxx),real(spinf*ixx)
		!write(301,*) real(elux),real(spinf*ryy),real(spinf*iyy)
		!write(302,*) real(elux),real(spinf*rzz),real(spinf*izz)
		
		!write(303,*) real(elux),real(spinf*rsp),real(spinf*isp)
		!write(304,*) real(elux),real(spinf*rsm),real(spinf*ism)

	
		!call flush(300)
		!!$omp end ordered
	end do
	!$omp end  do
	
	if (renorm) then
	 call imagrenorm(dimbse,int(numbse),esp,dielfxx)
	 call imagrenorm(dimbse,int(numbse),esp,dielfyy)
	 call imagrenorm(dimbse,int(numbse),esp,dielfzz)
	 call imagrenorm(dimbse,int(numbse),esp,dielfsp)
	 call imagrenorm(dimbse,int(numbse),esp,dielfsm)
	end if	
	
 if(Process==0) then
	do i=1,int(numbse)
		write(300,*) dielfxx(i,1),real(spinf*dielfxx(i,2)),real(spinf*dielfxx(i,3))
		write(301,*) dielfyy(i,1),real(spinf*dielfyy(i,2)),real(spinf*dielfyy(i,3))
		write(302,*) dielfzz(i,1),real(spinf*dielfzz(i,2)),real(spinf*dielfzz(i,3))
		write(303,*) dielfsp(i,1),real(spinf*dielfsp(i,2)),real(spinf*dielfsp(i,3))
		write(304,*) dielfsm(i,1),real(spinf*dielfsm(i,2)),real(spinf*dielfsm(i,3))
		!write(301,*) elux,real(spinf*rxy),real(spinf*ixy)
		!write(302,*) elux,real(spinf*rxz),real(spinf*ixz)
		!write(303,*) elux,real(spinf*ryy),real(spinf*iyy)
		!write(304,*) elux,real(spinf*ryz),real(spinf*iyz)
		!write(305,*) elux,real(spinf*rzz),real(spinf*izz)
	end do			
 endif

	deallocate(esp)
	deallocate(exxf,eyyf,ezzf)
	deallocate(espf,esmf)
	deallocate(rvec,hopmatrices)
	deallocate(ihopmatrices,ffactor)
	deallocate(dielfxx,dielfyy,dielfzz)
	deallocate(dielfsp,dielfsm)	

 if(Process==0) then
	close(100)
	close(200)
	close(300)
	close(301)
	close(302)
	close(303)
	close(304)
 endif
end subroutine spdielrawpol




!INCLUDE "./subprograms/sp_opt_bz-tool.f90"
!funciona apenas para semicondutores
!gfortran -mcmodel=large sp_opt_bz-tool.f90 -o sp_opt_bz-tool.x -llapack95 -lopenblas -fopenmp
subroutine spoptpolbz(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,rk,meshtype)
	double precision,parameter:: pi=acos(-1.)
	integer :: dimrpa !=ngrid*ngrid*nc*nv ! dimensão da matriz bse
	double precision,allocatable,dimension(:,:) :: eigv
	double complex,allocatable,dimension(:,:,:) :: vector
	double precision,allocatable,dimension(:,:) :: kpt !pontos k do grid
	integer,allocatable,dimension(:,:,:) :: stto
	double precision,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
	double complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores
	integer:: counter,c,v,i,j,f,l,h,k,kl,erro
	integer :: ncaux,nvaux
	double complex :: hxsp,hysp,hzsp
	double precision :: auxx,auxy,auxz,auxsp,auxsm
	double complex,parameter :: imag=cmplx(0.0,1.0)
	integer :: nk
	double precision,allocatable,dimension(:,:) :: output
	double precision,parameter :: norm=(1.0/sqrt(2.))
	integer :: ngkpt
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

	! INPUT : lendo os parametros do modelo de tight-binding
	!OPEN(UNIT=201, FILE= diein,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada ambiente dieletrico"

 if(Process==0) then
	!OUTPUT
	OPEN(UNIT=301, FILE=trim(outputfolder)//"/bz_act_x.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bz_act_x output file"
	OPEN(UNIT=302, FILE= trim(outputfolder)//"/bz_act_y.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bz_act_y output file"
	OPEN(UNIT=303, FILE=trim(outputfolder)//"/bz_act_z.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bz_act_z output file"
	OPEN(UNIT=304, FILE=trim(outputfolder)//"/bz_act_sp.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bz_act_sp output file"
	OPEN(UNIT=305, FILE=trim(outputfolder)//"/bz_act_sm.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bz_act_sm output file"
	OPEN(UNIT=306, FILE=trim(outputfolder)//"/dichroism.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening dichroism output file"
 endif



!	call OMP_SET_NUM_THREADS(nthreads)



	!inicio leitura parametros

	call hamiltonian_input_read(200,params)
	!ediel(2) = edielh


	!termino leitura parametros

	!parametros do calculo


	!termino parametros calculo

	ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
	dimrpa = nc*nv

	allocate(kpt(ngkpt,3))

	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpt)


	allocate(eaux(w90basis),vaux(w90basis,w90basis))
	allocate(eigv(ngkpt,w90basis),vector(ngkpt,w90basis,w90basis))
	allocate(nocpk(ngkpt))


	! $omp parallel do
	do i=1,ngkpt


		call eigsys(nthreads,scs,exc,nocpk(i),ffactor,kpt(i,1),kpt(i,2),kpt(i,3),w90basis,nvec,&
			    rlat,rvec,hopmatrices,&
		             ihopmatrices,efermi,eaux,vaux)


			do j=1,w90basis
	
				eigv(i,j)= eaux(j)

			end do
			


			do l=1,w90basis


				do h=1,w90basis

					vector(i,l,h)=vaux(l,h)


				end do
			

			end do

	
	end do
	! $omp end parallel do

	deallocate(eaux,vaux)



	!definindo os numeros quanticos dos estados




	allocate (stto(ngkpt,nc*nv,3))

do i=1,ngkpt

	counter=1
	
	ncaux=(nocpk(i)+1)+nc-1

	nvaux=nocpk(i)-nv+1

	do c=(nocpk(i)+1),ncaux


		do v=nvaux,nocpk(1)


			stto(i,counter,1) = counter !numero do estado

			stto(i,counter,2) = v   !numero da banda da valencia

			stto(i,counter,3) = c   !numero da banda da condução

			counter=counter+1
 			

		end do



	end do

end do

	counter=counter-1 !numero total de estados para equação bse

	!allocate(hoptxf(nk,dimrpa),hoptyf(nk,dimrpa),hoptspf(nk,dimrpa),hoptsmf(nk,dimrpa))

	!allocate(output(nk,7))

	!$omp do ordered
	do j=1,ngkpt

	auxx = 0.0
	auxy = 0.0
	auxz = 0.0
	auxsp = 0.0
	auxsm = 0.0

	! $omp critical
	do i=1,dimrpa

		call optspbz(vector(j,stto(j,i,2),:),vector(j,stto(j,i,3),:),&
		     kpt(j,1),kpt(j,2),kpt(j,3),&
		     ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		     ihopmatrices,hxsp,hysp,hzsp)

		hxsp= hxsp/( cmplx(eigv(j,stto(j,i,3))-eigv(j,stto(j,i,2)),sme,8) )
		hysp= hysp/( cmplx(eigv(j,stto(j,i,3))-eigv(j,stto(j,i,2)),sme,8) )
		hzsp= hzsp/( cmplx(eigv(j,stto(j,i,3))-eigv(j,stto(j,i,2)),sme,8) )
		
		auxx = auxx+real(hxsp*conjg(hxsp),8)	
		auxy = auxy+real(hysp*conjg(hysp),8)
		auxz = auxz+real(hzsp*conjg(hzsp),8)
		auxsp = auxsp&
		+real(norm*(hxsp+cmplx(0.,1.)*hysp)*conjg(norm*(hxsp+cmplx(0.,1.)*hysp)),8)
	auxsm = auxsm&
	+real(norm*(hxsp-cmplx(0.,1.)*hysp)*conjg(norm*(hxsp-cmplx(0.,1.)*hysp)),8)		



	end do
	! $omp end critical
		!output(j,1) = kpt(j,1)
		!output(j,2) = kpt(j,2)
		!output(j,3) = auxx
		!output(j,4) = auxy
		!output(j,5) = auxsp
		!output(j,6) = auxsm
		!output(j,7) = (auxsp-auxsm)/(auxsp+auxsm)
		!$omp ordered
   if(Process==0) then
		write(301,*) real(kpt(j,1)),real(kpt(j,2)),real(kpt(j,3)),auxx
		write(302,*) real(kpt(j,1)),real(kpt(j,2)),real(kpt(j,3)),auxy
		write(303,*) real(kpt(j,1)),real(kpt(j,2)),real(kpt(j,3)),auxz
		write(304,*) real(kpt(j,1)),real(kpt(j,2)),real(kpt(j,3)),auxsp
		write(305,*) real(kpt(j,1)),real(kpt(j,2)),real(kpt(j,3)),auxsm
		write(306,*) real(kpt(j,1)),real(kpt(j,2)),real(kpt(j,3)),(auxsp-auxsm)/(auxsp+auxsm)
   endif
		!$omp end ordered
	end do
	!$omp end do


	!do j=1,nk
	!	write(301,*) output(j,1),output(j,2),output(j,3)
	!	write(302,*) output(j,1),output(j,2),output(j,4)
	!	write(303,*) output(j,1),output(j,2),output(j,5)
	!	write(304,*) output(j,1),output(j,2),output(j,6)
	!	write(305,*) output(j,1),output(j,2),output(j,7)

	!end do



	deallocate(eigv,vector)
	deallocate(rvec,hopmatrices,ihopmatrices,ffactor)
	!deallocate(hoptxf,hoptyf,hoptspf,hoptsmf)
	deallocate(stto)
	deallocate(kpt)
	!deallocate(output)




 if(Process==0) then
	close(200)
	close(300)
	close(301)
	close(302)
	close(303)
	close(304)
	close(305)
	close(306)
 endif
end subroutine spoptpolbz


!INCLUDE "./subprograms/sp_solver-tool-diel.f90"
!gfortran -mcmodel=large rpa_optics-tool-diel.f90 -o rpa_opt-tool-diel.x -llapack95 -lopenblas -fopenmp
subroutine spoptics(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,rk,meshtype,tmcoef)
	double precision,parameter:: pi=acos(-1.)
	integer :: dimsp !=ngrid*ngrid*nc*nv ! dimensão da matriz bse
	double precision,allocatable,dimension(:,:) :: eigv
	double complex,allocatable,dimension(:,:,:) :: vector
	double precision,allocatable,dimension(:,:) :: kpt !pontos k do grid
	integer,allocatable,dimension(:,:) :: stt
	double precision,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
	double complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores
	!double precision,allocatable,dimension(:) :: esp
	integer:: counter,c,v,i,j,f,l,h,k,kl,erro
	integer :: ncaux,nvaux
	double complex:: matrizelbse
	double complex :: hxsp,hysp,hzsp,hrsp,hrsm
	!double complex :: exx,exy,exz,eyy,eyz,ezz
	double precision :: ec,ev
	double complex,parameter :: imag=cmplx(0.0,1.0)
	integer :: ngkpt
	!real,allocatable,dimension(:) :: auxx,auxy,auxz,auyy,auyz,auzz
	real,allocatable,dimension(:,:) :: vecres
	integer,allocatable,dimension(:) :: nocpk
	!definicoes diagonalizacao LA_HEEVR
	INTEGER   ::       ifail
	double precision,parameter :: ABSTOL=1.0e-6
	INTEGER          INFO
	double precision,allocatable,dimension(:) :: W,RWORK
	COMPLEX*16,allocatable,dimension(:,:) :: hbse
        INTEGER ::          LWMAX
   	INTEGER ::         LWORK
        double complex,allocatable,dimension (:) :: WORK
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
	logical :: tmcoef

	!fim modificacoes versao 2.1
	!call input_read

	! INPUT : lendo os parametros do modelo de tight-binding
	!OPEN(UNIT=201, FILE= diein,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada ambiente dieletrico"

 if(Process==0) then
	!OUTPUT
	OPEN(UNIT=301, FILE=trim(outputfolder)//"/sp_optics.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening sp_optics output file"
    	
    	
    	if (tmcoef) then
    	OPEN(UNIT=302, FILE=trim(outputfolder)//"/tm_coef.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening tm_coef output file"
    	OPEN(UNIT=304, FILE=trim(outputfolder)//"/tm_coef-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening tm_coef-pol output file"    	
    	else
    	 continue
    	end if    
    	
    	OPEN(UNIT=303, FILE=trim(outputfolder)//"/sp_optics-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening sp_optics-pol output file"	
    	
 endif  	



!	call OMP_SET_NUM_THREADS(nthreads)



	!inicio leitura parametros

	call hamiltonian_input_read(200,params)
	!ediel(2) = edielh


	!termino leitura parametros

	!parametros do calculo


	!termino parametros calculo

	ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
	dimsp = ngkpt*nc*nv


	allocate(kpt(ngkpt,3))
	allocate(nocpk(ngkpt))

	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpt)

	allocate(eaux(w90basis),vaux(w90basis,w90basis))
	allocate(eigv(ngkpt,nc+nv),vector(ngkpt,nc+nv,w90basis))


	! $omp parallel do
	do i=1,ngkpt


		call eigsys(nthreads,scs,exc,nocpk(i),ffactor,kpt(i,1),kpt(i,2),kpt(i,3),w90basis,nvec,&
			    rlat,rvec,hopmatrices,&
		             ihopmatrices,efermi,eaux,vaux)


			do j=1,nc+nv
	
				eigv(i,j)= eaux(nocpk(i)-nv+j)

			end do
			


			do l=1,nc+nv


				do h=1,w90basis

					vector(i,l,h)=vaux(nocpk(i)-nv+l,h)

				end do
			

			end do

	
	end do
	! $omp end parallel do
	deallocate(eaux,vaux)



	!definindo os numeros quanticos dos estados



	allocate (stt(ngkpt*nc*nv,4))

	call quantumnumbers2(w90basis,ngkpt,nc,nv,nocpk,nocpk,stt)

	!allocate(esp(dimsp))
	

	!allocate(auxx(dimsp),auxy(dimsp),auxz(dimsp))
	!allocate(auyy(dimsp),auyz(dimsp),auzz(dimsp))
	
	allocate(vecres(dimsp,15))

 if(Process==0) then
	write(301,*) "#","  ", "energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(303,*) "#","  ", "energy","  ","xx","  ","yy","  ","zz","  ","sp","  ","sm"
	if (tmcoef) then	
	write(302,*) "number of kpoints:",ngkpt
	write(302,*) "number of conduction states",nc
	write(302,*) "number of valence states",nv
	write(302,*) "#"," ", "kx", " ", "ky"," ", "kz"," ","nocp"," ", "nc"," ", "nv","  ", "energy","  ","xx","  ",&
		      "yy","  ","zz","  ","xy","  ","xz","  ","yz"
		      
	write(304,*) "number of kpoints:",ngkpt
	write(304,*) "number of conduction states",nc
	write(304,*) "number of valence states",nv
	write(304,*) "#"," ", "kx", " ", "ky"," ", "kz"," ","nocp"," ", "nc"," ", "nv","  ", "energy","  ","xx","  ",&
		      "yy","  ","zz","  ","sp","  ","sm"		   
	end if
 endif
	! $omp do ordered
	!$omp do 

	do i=1,dimsp
		ec = eigv(stt(i,4),stt(i,3))
		ev = eigv(stt(i,4),stt(i,2))
	       call optsp(eigv(stt(i,4),stt(i,2)),vector(stt(i,4),stt(i,2),:),&
		     eigv(stt(i,4),stt(i,3)),vector(stt(i,4),stt(i,3),:),&
		     kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),ffactor,sme,&
		     w90basis,nvec,rlat,rvec,hopmatrices,&
		     ihopmatrices,hxsp,hysp,hzsp)  
		
		hrsp = (hxsp+cmplx(0.,1.)*hysp)*(1.0/sqrt(2.))
		hrsm = (hxsp-cmplx(0.,1.)*hysp)*(1.0/sqrt(2.))
		
		vecres(i,1) = real(kpt(stt(i,4),1))
		vecres(i,2) = real(kpt(stt(i,4),2))
		vecres(i,3) = real(kpt(stt(i,4),3))
		vecres(i,4) = real(nocpk(stt(i,4)))
		vecres(i,5) = real(nocpk(stt(i,4))-nv+stt(i,3))
		vecres(i,6) = real(nocpk(stt(i,4))-nv+stt(i,2))		
		vecres(i,7) = real(eigv(stt(i,4),stt(i,3))-eigv(stt(i,4),stt(i,2)))		
		vecres(i,8) = real(hxsp*conjg(hxsp))		
		vecres(i,9) = real(hysp*conjg(hysp))		
		vecres(i,10) = real(hzsp*conjg(hzsp))
		vecres(i,11) = real(hxsp*conjg(hysp))
		vecres(i,12) = real(hxsp*conjg(hzsp))
		vecres(i,13) = real(hysp*conjg(hzsp))		 
		vecres(i,14) = real(hrsp*conjg(hrsp))
		vecres(i,15) = real(hrsm*conjg(hrsm))		

		! $omp ordered
		!write(301,"(7F15.4)") esp(i),auxx(i),auyy(i),auzz(i),auxy(i),auxz(i),auyz(i)
		
		!if (tmcoef) then
		!write(302,"(3F10.4,3I10.0,7F10.4)") kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),nocpk(stt(i,4)),stt(i,3),&
		!				  stt(i,2),esp(i),auxx(i),auyy(i),auzz(i),auxy(i),auxz(i),auyz(i)
						  
		!else
		 !continue
		!end if 						  
		! $omp end ordered

		!write(2077,*) i,"/",dimrpa

	end do

	!$omp end do

	call Bubblem(7,15,vecres, dimsp)
	
 if(Process==0) then
	do i=1,dimsp
		write(301,"(7F15.4)") vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),vecres(i,11),vecres(i,12),vecres(i,13)
		write(303,"(6F15.4)") vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),vecres(i,14),vecres(i,15)
		if (tmcoef) then
	        write(302,"(3F10.4,3I10.0,7F10.4)") vecres(i,1),vecres(i,2),vecres(i,3),int(vecres(i,4)),int(vecres(i,5)),&
						  int(vecres(i,6)),vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),&
						  vecres(i,11),vecres(i,12),vecres(i,13)
	        
	        write(304,"(3F10.4,3I10.0,6F10.4)") vecres(i,1),vecres(i,2),vecres(i,3),int(vecres(i,4)),int(vecres(i,5)),&
						  int(vecres(i,6)),vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),&
						  vecres(i,14),vecres(i,15)						  
		end if
	end do
 endif

	deallocate(eigv,vector)
	deallocate(rvec,hopmatrices,ihopmatrices,ffactor)
	!deallocate(auxx,auxy,auxz,auyy,auyz,auzz)
	deallocate(nocpk)
	deallocate(stt)
	!deallocate(esp)
	deallocate(kpt)
	deallocate(vecres)

 if(Process==0) then
	close(200)
	close(300)
	close(301)
	close(303)	
	if (tmcoef) then
	 close(302)
	 close(304)
	end if
 endif
end subroutine spoptics


!INCLUDE "./subroutines/diel-pp-subs.f90"
subroutine refracao(e1,e2,refrac)
	double precision :: e1,e2
	double precision :: aux,refrac
	aux = sqrt((e1*e1)+(e2*e2))
	refrac = sqrt((aux+e1)/2.0)
end subroutine refracao



subroutine extincao(e1,e2,extinc)
	double precision :: e1,e2
	double precision :: aux,extinc
	aux = sqrt((e1*e1)+(e2*e2))
	extinc = sqrt((aux-e1)/2.0)
end subroutine extincao



subroutine reflectibilidade(e1,e2,reflec)
	double precision :: e1,e2
	double precision :: reflec,refrac,extinc
	double precision :: aux1,aux2,aux3
	call extincao(e1,e2,extinc)
	call refracao(e1,e2,refrac)
	aux1 = (refrac-1.0)*(refrac-1.0)
	aux2 = (refrac+1.0)*(refrac+1.0)
	aux3 = extinc*extinc
	reflec = (aux1+aux3)/(aux2+aux3)
end subroutine reflectibilidade



subroutine abscoef(e1,e2,efoton,absc)
	double precision :: e1,e2
	double precision :: absc
	double precision :: aux,efoton
	double precision,parameter :: hc=19.746E-06
	aux = sqrt((e1*e1)+(e2*e2))
	absc=(sqrt(2.0)*efoton)*(1.0/hc)*sqrt(aux-e1)
end subroutine abscoef



subroutine enloss(e1,e2,eloss)
	double precision :: e1,e2
	double precision :: eloss
	double precision :: aux
	aux = (e1*e1)+(e2*e2)
	eloss= e2/aux
end subroutine enloss








  end module module_bse





