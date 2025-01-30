

module module_input
    use MPI_P
	use module_hamiltonian   
	implicit none 
	integer                        :: nthreads
	integer                        :: ngrid(3)
	integer                        :: nc,nv
	integer                        :: excwf0,excwff
	real(8)                        :: rk
	real(8)                        :: numdos
	real(8)                        :: ebse0,ebsef,numbse
	real(8)                        :: sme,exc,cshift
	real(8)                        :: mshift(3), ediel(3)
	real(8)                        :: ktol
	real(8)                        :: ez,w,lc,r0
	real(8)                        :: st,phavg,temp
	real(8)                        :: ctemp,tmax,eg,egd
	real(8)                        :: egs,ebgs
	character(len=70)              :: params   !parametros TB
	character(len=70)              :: orbw     !peso orbitais lcount
	character(len=70)              :: kpaths    !kpath
	character(len=70)              :: kpathsbse    !kpath
	character(len=70)              :: outputfolder    !pasta saida
	character(len=70)              :: calcparms
	character(len=70)              :: meshtype
	character(len=5)               :: coultype
	character(len=5)               :: sysdim
	character(len=1)               :: dft
	character(len=2)               :: ta
	character(len=6)               :: ses
	logical                        :: bandscalc,doscalc
	logical                        :: bse,bsepol,bsekpath
	logical                        :: spdiel,spdielpol,sppolbz
	logical                        :: spec
	logical                        :: berryk,berrybz
	logical                        :: pponly
	logical                        :: bsewf
	logical                        :: tmcoef
	logical                        :: dtfull,cpol
	logical                        :: bset,bsetbnd
	logical                        :: pce 
	logical                        :: renorm
	logical                        :: KBE                               ! KBE

 Contains

subroutine input_read
	character(len=70)               :: a,b
	character(len=70)               :: add
	integer                         :: erro

!default values
	nthreads = 1
	ngrid(1)= 1
	ngrid(2)= 1
	ngrid(3)= 1
	rk= 0.0
	nc = 1
	nv = 1
	numdos = 6001
	ebse0 = 0.0
	ebsef = 3.0
	numbse = 6001
	sme = 0.08
	cshift = 0.01
	exc = 0.000
	mshift(1) = 0.0
	mshift(2) = 0.0
	mshift(3) = 0.0
	ktol = 0.001 
	outputfolder = "./"
	calcparms = "./"
	meshtype = "MKH"
	coultype = "V3D"
	orbw = "non declared file"
	kpaths = "non declared file"
	kpathsbse = "non declared file"  
	bandscalc = .false.
	doscalc = .false.
	bse = .false.
    kbe = .false.                            ! KBE
	bsepol = .false.
	bsekpath = .false.
	spdiel = .false.
	spdielpol = .false.
	sppolbz = .false.
	spec = .false.
	berryk = .false.
	berrybz = .false.
	pponly = .false.
	bsewf = .false.
	tmcoef = .false.
	dtfull = .false.
	cpol = .false.
	excwf0 = 1
	excwff = 2
	sysdim = "3D"	
	dft = "V"
	ta = "FA"
	bset = .false.
	bsetbnd = .false.
	renorm = .true.
	ediel(1) = 1.0
	ediel(2) = 1.0
	ediel(3) = 1.0
	ez = 1.0
	w = 0.0
	lc = 1.0
	r0 = 1.0	
	st = 0.0 
	phavg = 0.0 
	temp = 0.0
	pce = .false.
	ses = "AM15G"
	ctemp = 298.15 
	tmax = 5E-06
	eg = 0.00
	egd = 0.00
	egs = 0.00
	ebgs = 0.00
!end default values

 if(Process==0) then
	do 

	read(*,*,iostat=erro) a,b
	if (erro/=0) exit

	select case (a)
	case ("RNMD=")
		read(b,*) renorm 	
	case ("PCE=")
		read(b,*) pce 	
	case ("SES=")
		read(b,*) ses
	case ("CTEMP=")
		read(b,*) ctemp
	case ("THMAX=")
		read(b,*) tmax	
	case ("EG=")
		read(b,*) eg	
	case ("EGD=")
		read(b,*) egd	
	case ("EGS=")
		read(b,*) egs	
	case ("EBGS=")
		read(b,*) ebgs
	case ("TA=")
		read(b,*) ta 	
	case ("BSET=")
		read(b,*) bset 	
	case ("BSET_BND=")
		read(b,*) bsetbnd 	
	case ("TEMP=")
		read(b,*) temp 	
	case ("PHAVG=")
		read(b,*) phavg 	
	case ("ST=")
		read(b,*) st 	
	case ("DFT=")
		read(b,*) dft 	
	case ("CSHIFT=")
		read(b,*) cshift 	
	case ("DTDIAG=")
		read(b,*) dtfull 
	case ("CPOL=")
		read(b,*) cpol 			
	case ("SYSDIM=")
		read(b,*) sysdim 	
	case ("R_0=")
		read(b,*) r0	
	case ("LC=")
		read(b,*) lc	
	case ("EDIEL_Z=")
		read(b,*) ez
	case ("W_COUL=")
		read(b,*) w
	case ("TMCOEF=")
		read(b,*) tmcoef
	case ("EXC_WF_F=")
		read(b,*) excwff
	case ("EXC_WF_I=")
		read(b,*) excwf0
	case ("BSE_WF=")
		read(b,*) bsewf
	case ("PP_ONLY=")
		read(b,*) pponly
	case ("BERRY_BZ=")
		read(b,*) berrybz
	case ("BERRY=")
		read(b,*) berryk
	case ("SPEC=")
		read(b,*) spec
	case ("OPT_BZ=")
		read(b,*) sppolbz
	!case ("DIEL_POL=")
	!	read(b,*) spdielpol
	case ("DIEL=")
		read(b,*) spdiel
	case ("BSE_BND=")
		read(b,*) bsekpath
	!case ("BSE_POL=")
	!	read(b,*) bsepol
	case ("BSE=")
		read(b,*) bse
	case ("KBE=")
		read(b,*) kbe
		print *,'KBE=',kbe
		lkbe = kbe
	case ("DOS=")
		read(b,*) doscalc
	case ("BANDS=")
		read(b,*) bandscalc
	case ("RK=")
		read(b,*) rk
	case ("COULOMB_POT=")
		coultype = b
	case ("MESH_TYPE=")
		meshtype = trim(adjustl(b))
	case ("OUTPUT=")
		outputfolder = b
	case ("CALC_DATA=")
		calcparms = b
	case ("NGX=")
		read(b,*) ngrid(1) 
	case ("NGY=")
		read(b,*) ngrid(2) 
	case ("NGZ=")
		read(b,*) ngrid(3) 
	case ("NBANDSC=")
		read(b,*) nc 
	case ("NBANDSV=")
		read(b,*) nv 
	!case ("ENDOSI=")
	!	read(b,*) edos0 
	!case ("ENDOSF=")
	!	read(b,*) edosf 
	case ('NEDOS=')
		read(b,*) numdos 
	case ("ENSPECI=")
		read(b,*) ebse0 
	case ("ENSPECF=")
		read(b,*) ebsef 
	case ("NESPEC=")
		read(b,*) numbse 
	case ("SIGMA=")
		read(b,*) sme 		
	case ("KTOL=")
		read(b,*) ktol 
	case ("PARAMS_FILE=")
		params = b
	case ("KPATH_FILE=")
		kpaths = b
	case ("KPATH_BSE=")
		kpathsbse = b
	case ("ORB_W=")
		orbw = b
	case ("EDIEL_T=")
		read(b,*) ediel(1) 
	case ("EDIEL=")
		read(b,*) ediel(2) 		
	case ("EDIEL_B=")
		read(b,*) ediel(3) 
	case ("EXC_CORR=")
		read(b,*) exc 
	case ("NTHREADS=")
		read(b,*) nthreads 
	case ("SHIFT_1=")
		read(b,*) mshift(1) 
	case ("SHIFT_2=")
		read(b,*) mshift(2) 
	case ("SHIFT_3=")
		read(b,*) mshift(3) 
	case default
	 continue
	end select
	end do
	print *,'input: ngrid=',ngrid
	print *,'input: rk=',rk
	print *,'input: meshtype=',meshtype
 endif

   call P_sendI(nthreads)                       ! send variables to all Processes
   call P_sendI(nc)
   call P_sendI(nv)
   call P_sendI(excwf0)
   call P_sendI(excwff)
   call P_sendI1(ngrid,3)
   call P_sendR(numdos)
   call P_sendR(ktol)
   call P_sendR(ez)
   call P_sendR(w)
   call P_sendR(lc)
   call P_sendR(r0)
   call P_sendR(st)
   call P_sendR(phavg)
   call P_sendR(temp)
   call P_sendR(ctemp)
   call P_sendR(tmax)
   call P_sendR(eg)
   call P_sendR(egd)
   call P_sendR(egs)
   call P_sendR(ebgs)
   call P_sendR(ebse0) 
   call P_sendR(ebsef)
   call P_sendR(numbse)
   call P_sendR(sme)
   call P_sendR(rk)
   call P_sendR(exc)
   call P_sendR(cshift)
   call P_sendR1(mshift,3) 
   call P_sendR1(ediel,3) 
   call P_sendA(params)
   call P_sendA(orbw)
   call P_sendA(kpaths)
   call P_sendA(kpathsbse)
   call P_sendA(outputfolder)
   call P_sendA(calcparms)
   call P_sendA(meshtype)
   call P_sendA(coultype)
   call P_sendA(sysdim)
   call P_sendA(dft)
   call P_sendA(ta)
   call P_sendA(ses)
   call P_sendL(bandscalc)
   call P_sendL(doscalc)
   call P_sendL(bse)
   call P_sendL(bsepol)
   call P_sendL(bsekpath)
   call P_sendL(spdiel)
   call P_sendL(spdielpol)
   call P_sendL(sppolbz)
   call P_sendL(spec)
   call P_sendL(berryk)
   call P_sendL(berrybz)
   call P_sendL(pponly)
   call P_sendL(bsewf)
   call P_sendL(tmcoef)
   call P_sendL(dtfull)
   call P_sendL(cpol)
   call P_sendL(bset)
   call P_sendL(bsetbnd)
   call P_sendL(pce) 
   call P_sendL(renorm)
   call P_sendL(KBE)
end subroutine input_read



  subroutine print_input_param
   print *,'nthreads=',nthreads
   print *,'nc=',nc
   print *,'nv=',nv
   print *,'excwf0=',excwf0
   print *,'excwff=',excwff
   print *,'ngrid=',ngrid
   print *,'numdos=',numdos
   print *,'ktol=',ktol
   print *,'ez=',ez
   print *,'w=',w
   print *,'lc=',lc
   print *,'r0=',r0
   print *,'st=',st
   print *,'phavg=',phavg
   print *,'temp=',temp
   print *,'ctemp=',ctemp
   print *,'tmax=',tmax
   print *,'eg=',eg
   print *,'egd=',egd
   print *,'egs=',egs
   print *,'ebgs=',ebgs
   print *,'ebse0=',ebse0 
   print *,'ebsef=',ebsef
   print *,'numbse=',numbse
   print *,'sme=',sme
   print *,'rk=',rk
   print *,'exc=',exc
   print *,'cshift=',cshift
   print *,'mshift=',mshift 
   print *,'ediel=',ediel 
   print *,'params=',params
   print *,'orbw=',orbw
   print *,'kpaths=',kpaths
   print *,'kpathsbse=',kpathsbse
   print *,'outputfolder=',outputfolder
   print *,'calcparms=',calcparms
   print *,'meshtype=',meshtype
   print *,'coultype=',coultype
   print *,'sysdim=',sysdim
   print *,'dft=',dft
   print *,'ta=',ta
   print *,'ses=',ses
   print *,'bandscalc=',bandscalc
   print *,'doscalc=',doscalc
   print *,'bse=',bse
   print *,'bsepol=',bsepol
   print *,'bsekpath=',bsekpath
   print *,'spdiel=',spdiel
   print *,'spdielpol=',spdielpol
   print *,'sppolbz=',sppolbz
   print *,'spec=',spec
   print *,'berryk=',berryk
   print *,'berrybz=',berrybz
   print *,'pponly=',pponly
   print *,'bsewf=',bsewf
   print *,'tmcoef=',tmcoef
   print *,'dtfull=',dtfull
   print *,'cpol=',cpol
   print *,'bset=',bset
   print *,'bsetbnd=',bsetbnd
   print *,'pce=',pce
   print *,'renorm=',renorm
   print *,'KBE=',KBE
  end subroutine print_input_param





subroutine param_out(unitout,nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,bandscalc,doscalc,bse,bsepol,bsekpath,spec,&
		     spdiel,spdielpol,sppolbz,berryk,berrybz,pponly,bsewf,excwf0,excwff,&
		     tmcoef,ez,w,lc,r0,sysdim,dtfull,cpol,cshift,dft,bset,bsetbnd,&
		     st,phavg,temp,ta,pce,ses,ctemp,tmax,eg,egd,egs,ebgs,renorm)

	integer :: unitout
	real(8),dimension(3) :: ediel
	integer :: nthreads
	integer,dimension(3) :: ngrid
	integer :: nc,nv,excwf0,excwff
	real(8) :: numdos
	real(8) :: ebse0,ebsef,numbse
	real(8) :: sme,exc,cshift
	real(8),dimension(3) :: mshift
	real(8) :: ktol
	character(len=70) :: params   !parametros TB
	character(len=70) :: orbw     !peso orbitais lcount
	character(len=70) :: kpaths    !kpath
	character(len=70) :: kpathsbse    !kpath
	character(len=70) :: outputfolder    !pasta saida
	character(len=70) :: calcparms
	character(len=70) :: meshtype
	character(len=5) :: coultype
	character(len=5) :: sysdim
	character(len=1) :: dft
	character(len=2) :: ta	
	character(len=6) :: ses		
	logical :: bandscalc,doscalc
	logical :: bse,bsepol,bsekpath
	logical :: spdiel,spdielpol,sppolbz
	logical :: spec,tmcoef
	logical :: berryk,berrybz,pponly,bsewf
	logical :: dtfull,cpol	
	logical :: bset,bsetbnd
	logical :: pce,renorm	
	real(8) :: ez,w,lc,r0
	real(8) :: st,phavg,temp
	real(8) :: ctemp,tmax,eg,egd
	real(8) :: egs,ebgs		

 if(Process==0) then
	write(unitout,"(A10,I0)") "NTHREADS= ",nthreads
	write(unitout,"(A8,A2)") "SYSDIM= ",sysdim
	write(unitout,"(A5,A1)") "DFT= ",dft		
	write(unitout,*)
	write(unitout,*) "INPUT/OUTPUT FILES"
	write(unitout,*)
	write(unitout,"(A8,A70)") "OUTPUT= ", outputfolder
	write(unitout,"(A11,A70)") "CALC_DATA= ", calcparms
	write(unitout,"(A13,A70)") "PARAMS_FILE= ", params
	write(unitout,"(A12,A70)") "KPATH_FILE= ", kpaths
	write(unitout,"(A11,A70)") "KPATH_BSE= ", kpathsbse
	write(unitout,"(A7,A70)") "ORB_W= ", orbw
	write(unitout,*)	
	write(unitout,*) "CONTROL"
	write(unitout,*)
	write(unitout,"(A7,L1)")  "BANDS= ",bandscalc
	write(unitout,"(A5,L1)")  "DOS= ",doscalc	
	write(unitout,"(A5,L1)")  "BSE= ",bse
	write(unitout,"(A6,L1)")  "BSET= ",bset	
	write(unitout,"(A9,L1)")  "BSE_BND= ",bsekpath
	write(unitout,"(A10,L1)")  "BSET_BND= ",bsetbnd	
	write(unitout,"(A6,L1)")  "DIEL= ",spdiel
	write(unitout,"(A8,L1)")  "OPT_BZ= ",sppolbz
	write(unitout,"(A6,L1)")  "SPEC= ",spec
	write(unitout,"(A10,L1)") "BERRY_BZ= ",berrybz
	write(unitout,"(A7,L1)")  "BERRY= ",berryk
	write(unitout,"(A9,L1)")  "PP_ONLY= ",pponly
	write(unitout,"(A8,L1)") "BSE_WF= ",bsewf	
	write(unitout,"(A8,L1)") "TMCOEF= ",tmcoef
	write(unitout,"(A8,L1)") "DTDIAG= ",dtfull
	write(unitout,"(A6,L1)") "CPOL= ",cpol
	write(unitout,"(A5,L1)") "PCE= ",pce			
	write(unitout,*)
	write(unitout,*) "K-MESH"
	write(unitout,*)
	write(unitout,"(A5,I0)") "NGX= ", ngrid(1)
	write(unitout,"(A5,I0)") "NGY= ", ngrid(2)
	write(unitout,"(A5,I0)") "NGZ= ", ngrid(3)
	write(unitout,"(A9,F8.4)") "SHIFT_1= ", mshift(1)
	write(unitout,"(A9,F8.4)") "SHIFT_2= ", mshift(2)
	write(unitout,"(A9,F8.4)") "SHIFT_3= ", mshift(3)
	write(unitout,*)
	write(unitout,*) "DOS"
	write(unitout,*)
	write(unitout,"(A7,I0)") "NEDOS= ", int(numdos)
	write(unitout,"(A7,F8.4)") "SIGMA= ", sme
	write(unitout,*)
	write(unitout,*) "BSE/OPTICAL PROPERTIES"
	write(unitout,*)
	write(unitout,"(A9,I0)") "NBANDSC= ", nc
	write(unitout,"(A9,I0)") "NBANDSV= ", nv
	write(unitout,"(A9,F8.4)") "ENSPECI= ", ebse0
	write(unitout,"(A9,F8.4)") "ENSPECF= ", ebsef
	write(unitout,"(A8,I0)") "NESPEC= ", int(numbse)
	write(unitout,"(A6,F8.4)") "KTOL= ", ktol
	write(unitout,"(A10,I0)") "EXC_WF_I= ",excwf0
	write(unitout,"(A10,I0)") "EXC_WF_F= ",excwff
	write(unitout,"(A13,A5)") "COULOMB_POT= ",coultype
	write(unitout,"(A8,F8.4)") "CSHIFT= ", cshift	
	write(unitout,"(A6,L1)") "RNMD= ",renorm		
	write(unitout,*)
	write(unitout,*) "PARAMETERS FOR COULOMB POTENTIALS"
	write(unitout,*)
	write(unitout,"(A9,F8.4)") "EDIEL_Z= ",ez
	write(unitout,"(A8,F8.4)") "W_COUL= ",w
	write(unitout,"(A4,F8.4)") "LC= ",lc
	write(unitout,"(A5,F8.4)") "R_0= ",r0
	write(unitout,"(A9,F8.4)") "EDIEL_T= ", ediel(1)
	write(unitout,"(A9,F8.4)") "EDIEL_B= ", ediel(3)
	write(unitout,"(A7,F8.4)") "EDIEL= ", ediel(2)	
	write(unitout,*)	
	write(unitout,*) "PARAMETERS FOR BSE with Temperature Effects"
	write(unitout,*)
	write(unitout,"(A4,A2)") "TA= ",ta			
	write(unitout,"(A4,F8.4)") "ST= ",st
	write(unitout,"(A7,F8.4)") "PHAVG= ",phavg
	write(unitout,"(A6,F8.4)") "TEMP= ",temp
	write(unitout,*)	
	write(unitout,*) "PARAMETERS FOR PCE"
	write(unitout,*)
	write(unitout,"(A5,A6)") "SES= ",ses
	write(unitout,"(A7,F8.4)") "CTEMP= ",ctemp
	write(unitout,"(A7,E15.4)") "THMAX= ",tmax				
	write(unitout,"(A4,F8.4)") "EG= ",eg
	write(unitout,"(A5,F8.4)") "EGD= ",egd
	write(unitout,"(A5,F8.4)") "EGS= ",egs
	write(unitout,"(A6,F8.4)") "EBGS= ",ebgs			
	write(unitout,*)
	write(unitout,*)
    call flush(unitout)
 endif
end subroutine param_out


end module module_input





