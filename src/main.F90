


program main
    use MPI_P                                                             ! MPI module 
    use module_input
    use module_hamiltonian
    use module_bse
    use module_pce
    use module_tools
    use module_berry
    implicit none
    integer :: erro
    double precision:: t0,tf
    integer,dimension(8) :: values,values2
    double precision,dimension(3,3) :: rlatv

    call P_start                                                   ! start parallel calculation

    call input_read

    call MSG('MAIN: open file '//trim(adjustl(params)))
    if(Process==0) then
!     print *,'open file ',params
     OPEN(UNIT=2055, FILE= trim(adjustl(params)),STATUS='old', IOSTAT=erro)
        if (erro/=0) stop "Error opening wannier90 hr input file"   

     read(2055,*) systype
     read(2055,*) scs
     read(2055,*) efermi
     read(2055,*) rlatv(1,1),rlatv(1,2),rlatv(1,3)
     read(2055,*) rlatv(2,1),rlatv(2,2),rlatv(2,3)
     read(2055,*) rlatv(3,1),rlatv(3,2),rlatv(3,3)
     close(2055)
    endif
    call P_sendA(systype)                                 ! send array to all other processes
    call P_sendR(scs)                                     ! send variable to all other processes
    call P_sendR(efermi)                                  ! send variable to all other processes
    call P_sendR2(rlatv,3,3)                               ! send array to all other processes

if(Process==1) then
   call print_input_param
   print *,'  systype=',systype
   print *,'  scs=',scs
   print *,'  efermi=',efermi
   print *,'Process=',Process,'  rlatv=',rlatv
endif
    if (meshtype .eq. "RK3D") then

        call rkmesh(rk,rlatv,ngrid)

    else if (meshtype .eq. "RK2D") then

        call rkmesh2D(rk,rlatv,ngrid)
    end if
call P_wait

   call MSG('MAIN: write log.dat file')
   if(Process==0) then
    !print *,'open log.dat file'
    OPEN(UNIT=2077, FILE= trim(calcparms)//"/log.dat",STATUS='unknown', IOSTAT=erro)
        if (erro/=0) stop "Error opening log output file "

    call param_out(2077,nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
             ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
             exc,mshift,coultype,bandscalc,doscalc,bse,bsepol,bsekpath,spec,&
             spdiel,spdielpol,sppolbz,berryk,berrybz,pponly,bsewf,excwf0,excwff,&
             tmcoef,ez,w,lc,r0,sysdim,dtfull,cpol,cshift,dft,bset,bsetbnd,&
             st,phavg,temp,ta,pce,ses,ctemp,tmax,eg,egd,egs,ebgs,renorm)

    call cpu_time(t0)
    call date_and_time(VALUES=values)

    write(2077,*)
    write(2077,*) 'Begin','  ','day',values(3),'',values(5),'hours',values(6),'min',values(7),'seg'
    write(2077,*)
   endif

    if (pponly) go to 131

    call MSG('MAIN: call bandstool')
    if (bandscalc) then

     call bandstool(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
             ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
             exc,mshift,coultype,rk,dft)
     if(Process==0) write(2077,*) "Band Structure finished"
     if(Process==0)call flush(2077)
    end if


    call MSG('MAIN: call dostool')
    if (doscalc) then
     call dostool(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
             ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
             exc,mshift,coultype,rk,meshtype,dft)
     if(Process==0)write(2077,*) "DOS finished"
     if(Process==0) call flush(2077)
    end if

    call MSG('MAIN: call bsesolver')
    if (bse .and. bset) then
     if(Process==0) write(2077,*) "ERROR: Choose BSE=T or BSET=T not both"
     stop
    elseif (bse) then
     call bsesolver(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
             ebse0,ebsef,numbse,cshift,ktol,params,kpaths,kpathsbse,orbw,ediel, &
             exc,mshift,coultype,ez,w,r0,lc,rk,meshtype,bsewf,excwf0,excwff,dtfull,cpol,tmcoef)


     if(Process==0) write(2077,*) "BSE and Single particle dielectric calculation finished"
     if(Process==0) call flush(2077)

    elseif (bset) then
    call MSG('MAIN: call bsesolvertemp')
    call bsesolvertemp(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
             ebse0,ebsef,numbse,cshift,ktol,params,kpaths,kpathsbse,orbw,ediel, &
             exc,mshift,coultype,ez,w,r0,lc,rk,meshtype,bsewf,excwf0,excwff,dtfull,&
             cpol,tmcoef,st,phavg,ta,temp)
     if(Process==0) write(2077,*) "BSE and Single particle dielectric calculation, with temperature, finished"  
     if(Process==0) call flush(2077)            
    end if  
    
    if (bsekpath .and. bsetbnd) then
     if(Process==0) write(2077,*) "ERROR: Choose BSE_BND=T or BSET_BND=T not both"
     stop
    elseif (bsekpath) then
     call MSG('MAIN: call bsebnds')
     call bsebnds(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
             ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
             exc,mshift,coultype,ez,w,r0,lc,rk,meshtype,bsewf,excwf0,excwff)    
     if(Process==0) write(2077,*) "BSE exciton band structure finished"
     if(Process==0) call flush(2077)             
    elseif (bsetbnd) then
     call MSG('MAIN: call bsebndstemp')
     call bsebndstemp(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
             ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
             exc,mshift,coultype,ez,w,r0,lc,rk,meshtype,bsewf,excwf0,excwff,&
             st,phavg,ta,temp)  
     if(Process==0) write(2077,*) "BSE exciton band structure, with temperature, finished"
     if(Process==0) call flush(2077)
    end if      

    if (spdiel) then
     call MSG('MAIN: call spoptics')
     call spoptics(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
             ebse0,ebsef,numbse,cshift,ktol,params,kpaths,kpathsbse,orbw,ediel, &
             exc,mshift,coultype,rk,meshtype,tmcoef)
     if(Process==0) write(2077,*) "Single particle dielectric calculation finished"
     if(Process==0) call flush(2077)
    end if

    if (sppolbz) then
     call MSG('MAIN: call spoptpolbz')
     call spoptpolbz(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
             ebse0,ebsef,numbse,cshift,ktol,params,kpaths,kpathsbse,orbw,ediel, &
             exc,mshift,coultype,rk,meshtype)
     if(Process==0) write(2077,*) "Single particle optical activity in BZ finished"
     if(Process==0) call flush(2077)
    end if

    !calculo curvatura de berry

    if (berryk) then
     call MSG('MAIN: call berrycurv')
     call berrycurv(nthreads,outputfolder,params,kpaths,sme)
     if(Process==0) write(2077,*) "Berry curvature in kpath finished"
     if(Process==0) call flush(2077)
    end if

    if (berrybz) then
     call MSG('MAIN: call berrycurvbz')
     call berrycurvbz(nthreads,outputfolder,params,sme,ngrid,mshift)
     if(Process==0) write(2077,*) "Berry curvature in BZ finished"
     if(Process==0) call flush(2077)
    end if
    
    
131 continue

    !calculo espectro
    if ((spec) .and. (bse .or. bset)) then
     call MSG('MAIN: call bsedielraw')
     call bsedielraw(nthreads,outputfolder,renorm,params,ngrid,nc,nv,ebse0,ebsef,numbse,cshift) !calculate dielectric constants with BSE
     call MSG('MAIN: call bseoptprop')
     call bseoptprop(numbse,outputfolder) !calculate abs coeficient and other properties with BSE
    if (.not.bset) then
     call MSG('MAIN: call lifetime')
     call lifetime(sysdim,numbse,ngrid,rlat,nc,nv,outputfolder)
    end if
     if(Process==0) write(2077,*) "BSE dielectric properties calculated"
     if(Process==0) call flush(2077)

    if (cpol) then
     call MSG('MAIN: call bsedielrawpol')
     call bsedielrawpol(nthreads,outputfolder,renorm,params,ngrid,nc,nv,ebse0,ebsef,numbse,cshift) !calculate dielectric constants with BSE for                                                     different light polarization
     call MSG('MAIN: call bseoptproppol')
     call bseoptproppol(numbse,outputfolder) !calculate abs coeficient and other properties with BSE
     if (.not.bset) then
      call MSG('MAIN: call lifetimepol')
      call lifetimepol(sysdim,numbse,ngrid,rlat,nc,nv,outputfolder)
     end if
     if(Process==0) write(2077,*) "BSE absorption spectrum with light polarization calculated"
     if(Process==0) call flush(2077)
     end if
    end if

    if ((spec) .and. ((spdiel) .or. (bse .or. bset) )) then
     call MSG('MAIN: call spdielraw')
     call spdielraw(nthreads,outputfolder,renorm,params,ngrid,nc,nv,ebse0,ebsef,numbse,cshift) !calculate dielectric constants
     call MSG('MAIN: call spoptprop')
     call spoptprop(numbse,outputfolder) !calculate abs coeficient and other properties
     if(Process==0) write(2077,*) "Single particle dielectric properties calculated"
     if(Process==0) call flush(2077)
    end if

    if ((spec) .and. ((cpol) .or. (bse .or. bset) )) then
     call MSG('MAIN: call spdielrawpol')
     call spdielrawpol(nthreads,outputfolder,renorm,params,ngrid,nc,nv,ebse0,ebsef,numbse,cshift) !calculate dielectric constants different                                                        light polarizations
     call MSG('MAIN: call spoptproppol')
     call spoptproppol(numbse,outputfolder) !calculate abs coeficient and other properties
     if(Process==0) write(2077,*) "Single particle absorption spectrum with light polarization calculated"
     if(Process==0) call flush(2077)
    end if

    if (pce) then
     call MSG('MAIN: call pcecalc IPA')
     call pcecalc(outputfolder,numbse,"IPA",ctemp,ses,tmax,eg,egd)
     if(Process==0) write(2077,*) "PCE with single particle absortion spectrum"
     if(Process==0) call flush(2077)
     call MSG('MAIN: call pcecalc BSE')
     call pcecalc(outputfolder,numbse,"BSE",ctemp,ses,tmax,egs,ebgs)    
     if(Process==0) write(2077,*) "PCE with excitonic effects (BSE)"
     if(Process==0) call flush(2077)
    end if

 if(Process==0) then
    call cpu_time(tf)
    call date_and_time(VALUES=values2)
    write(2077,*)
    write(2077,*) 'End','   ','day',values2(3),'',values2(5),'hours',values2(6),'min',values2(7),'seg'
    write(2077,*)
    !write(2077,*) "Total Time:",(tf-t0)/nthreads,"s"   
    close(2077)
 endif

   call MSG('MAIN:           NORMAL TERMINATION')
   call P_stop                                          ! stop parallel calculation

end program main







