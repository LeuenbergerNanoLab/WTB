


  module module_coulomb
   use module_math
   implicit none
  contains
  
!INCLUDE "./subroutines/coulomb_pot.f90"
!potencial 2D Keldysh (Ridolfi) DOI:https://doi.org/10.1103/PhysRevB.97.205409
real(8) function v2dk(kpt1,kpt2,ediel,rlat,ngrid,lc,tolr)

!	implicit none

	double precision,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	double precision,parameter:: pi=acos(-1.)
	double precision,parameter :: alpha1 = 1.76
	double precision,parameter :: alpha2 = 1.0
	double precision,parameter :: alpha3 = 0d0

	integer,dimension(3) :: ngrid
	double precision,dimension(3) :: kpt1,kpt2
	double precision,dimension(3,3) :: rlat
	double precision :: tolr
	double precision,dimension(3) :: ediel
	double precision :: lc
	double precision :: gridaux1,ed


	double precision :: a0,modk,vc
	double precision :: r0,vbz,auxi
!	double precision :: v2dk

	call alat2D(rlat,a0)
	call modvec(kpt1,kpt2,modk)
	call vcell2D(rlat,vc)

	r0= ((ediel(2)-1.0)*lc)/(ediel(1)+ediel(3))

	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))

	gridaux1 = dble(ngrid(1)*ngrid(2))

	auxi = (2.*pi*r0)/(a0*sqrt(gridaux1))

	ed = (ediel(1)+ediel(3))/(2.0)


	if (modk .lt. tolr) then

		v2dk = vbz*(cic/ed)*(a0*sqrt(gridaux1)/(2.*pi))*(alpha1+auxi*alpha2+alpha3*auxi**2)
	else 

		v2dk = vbz*(cic/ed)*(1./(modk*(1+(r0*modk))))
	end if

	


end function v2dk

!potencial 3D tradicional

real(8) function vcoul(kpt1,kpt2,rlat,ngrid,tolr)

!	implicit none

	double precision,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	double precision,parameter:: pi=acos(-1.)

	integer,dimension(3) :: ngrid
	double precision,dimension(3) :: kpt1,kpt2
	double precision,dimension(3,3) :: rlat
	double precision,dimension(3) :: ediel

	double precision :: modk,ed,tolr,vbz,vc
!	double precision :: vcoul

	call modvec(kpt1,kpt2,modk)
	call vcell3D(rlat,vc)

	ed = 1.0

	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))

	if (modk .lt. tolr) then

		vcoul = 0.0
	else 

		vcoul = vbz*(cic/ed)*(1.0/(modk*modk))
	end if

end function vcoul

!potencial 3D considerando ambiente dielétrico

real(8) function v3diel(kpt1,kpt2,ediel,rlat,ngrid,tolr)

!	implicit none

	double precision,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	double precision,parameter:: pi=acos(-1.)

	integer,dimension(3) :: ngrid
	double precision,dimension(3) :: kpt1,kpt2
	double precision,dimension(3,3) :: rlat
	double precision :: tolr
	double precision,dimension(3) :: ediel

	double precision :: modk,ed,vbz,vc
!	double precision :: v3diel

	call modvec(kpt1,kpt2,modk)
	call vcell3D(rlat,vc)

	ed = ediel(2)

	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))

	if (modk .lt. tolr) then

		v3diel = 0.0
	else 

		v3diel = vbz*(cic/ed)*(1.0/(modk*modk))
	end if



end function v3diel 

!potencial 2D truncado (DOI: 10.1103/PhysRevB.73.205119)

real(8) function v2dt(kpt1,kpt2,ngrid,rlat,tolr)

	implicit none

	double precision,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	double precision,parameter:: pi=acos(-1.)

	double precision,dimension(3) :: kpt1,kpt2,vkpt
	double precision,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid

	double precision :: vc,vbz,modk,factor
	double precision :: gpar,gz,rc
!	double precision :: v2dt
	double precision :: tolr

	double precision :: aux1,aux2,aux3,aux4,aux5

	call vcell3D(rlat,vc)
	call modvec(kpt1,kpt2,modk)

	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))

	vkpt=kpt1-kpt2

	gz= abs(vkpt(3))

	gpar= sqrt((vkpt(1)*vkpt(1))+(vkpt(2)*vkpt(2)))

	rc = 0.5*rlat(3,3)

	!factor= 4.0*pi
	factor=1.0


	if ( (gpar .lt. tolr) .and. (gz .lt. tolr) ) then

		!v2dt = (vbz*cic)*(-2.0*pi*rc*rc)
		v2dt = (vbz*cic)*(-0.5*rc*rc)		

	else if ((gpar .lt. tolr) .and. (gz .ge. tolr)) then

		v2dt = (vbz*cic)*((factor)/(modk*modk))*(1.0-cos(gz*rc)-(gz*rc*sin(gz*rc)))

	else 
		aux1= gz/gpar
		aux2= gpar*rc
		aux3= gz*rc
		
		aux4= aux1*sin(aux3)
		aux5= cos(aux3)

		v2dt = (vbz*cic)*((factor)/(modk*modk))*(1.0+(dexp(-aux2)*(aux4-aux5)))

	end if



end function v2dt

!potencial 0D truncado (DOI: 10.1103/PhysRevB.73.205119)

real(8) function v0dt(kpt1,kpt2,ngrid,rlat,tolr)

	implicit none

	double precision,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	double precision,parameter:: pi=acos(-1.)

	double precision,dimension(3) :: kpt1,kpt2
	double precision,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid

	double precision :: vc,modk,cr,vbz,factor
	double precision,dimension(3) :: vsize
	double precision :: tolr

	call vcell3D(rlat,vc)
	call modvec(kpt1,kpt2,modk)

	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))

	call vecsize(rlat(1,:),vsize(1))
	call vecsize(rlat(2,:),vsize(2))
	call vecsize(rlat(3,:),vsize(3))

	cr = MIN(vsize(1),vsize(2),vsize(3))
	cr = 0.5*cr

	!factor= 4.0*pi
	factor=1.0


	if (modk .lt. tolr) then

		!v0dt = (vbz*cic)*(2.0*pi)*cr*cr
		v0dt = (vbz*cic)*(0.5)*cr*cr		
	else 

		v0dt = (vbz*cic)*((factor)/(modk*modk))*(1.0-cos(cr*modk))
	end if


end function v0dt

!potencial 2D truncado v2 (DOI: 10.1103/PhysRevB.73.233103)


real(8) function v2dt2(kpt1,kpt2,ngrid,rlat,lc,tolr)

	implicit none

	double precision,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	double precision,parameter:: pi=acos(-1.)

	double precision,dimension(3) :: kpt1,kpt2,vkpt
	double precision,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid

	double precision :: vc,modk,vbz,qxy,lc
	double precision,dimension(3) :: vsize
	double precision :: factor
	double precision :: tolr

	!factor= 4.0*pi
	factor=1.0

	call vcell3D(rlat,vc)
	call modvec(kpt1,kpt2,modk)

	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))

	vkpt=kpt1-kpt2

	qxy= sqrt((vkpt(1)*vkpt(1))+(vkpt(2)*vkpt(2)))


	if (modk .lt. tolr) then 

		v2dt2 = 0.0
	else 

		v2dt2 = (vbz*cic)*((factor)/(modk*modk))*(1.0-dexp(-0.5*qxy*lc)*cos(0.5*lc*vkpt(3)))
	end if
	


end function v2dt2

!potencial Rytova-Keldysh (DOI: 10.1103/PhysRevB.98.125308)

real(8) function v2drk(kpt1,kpt2,ngrid,rlat,ediel,lc,ez,w,r0,tolr)

!	implicit none

	double precision,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	double precision,parameter:: pi=acos(-1.)

	double precision,dimension(3) :: kpt1,kpt2,vkpt
	double precision,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid
	double precision :: modk
	double precision,dimension(3) :: ediel
	double precision :: vc,lc,vbz,tolr

	
	!parameters
	double precision :: w,ew
	double precision :: ez,epar
	double precision :: eta,kappa
	double precision :: et,eb,pb,pt
	double precision :: r0
	double precision :: aux1,aux2,aux3
	

	call vcell2D(rlat,vc)
	call modvec(kpt1,kpt2,modk)
	
	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))
	!vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3)))
	
	!r0= ((ediel(2)-1.0)*rlat(3,3))/(ediel(1)+ediel(3))
	!r0 = 40.0
	!lc = 6.5


	epar = ediel(2)
	et = ediel(1)
	eb = ediel(3)
	
	eta = sqrt(epar/ez)
	kappa = sqrt(epar*ez)
	
	pb = (eb-kappa)/(eb+kappa)
	pt = (et-kappa)/(et+kappa)
	

	
	if (modk .lt. tolr) then

		v2drk = 0.0
	else 
	
		aux1 = (1.0-(pb*pt*dexp(-2.0*modk*eta*lc)))*kappa
		aux2 = (1.0-(pt*dexp(-eta*modk*lc)))*(1.0-(pb*dexp(-eta*modk*lc)))
		aux3 = r0*modk*dexp(-modk*w)
	
		ew = (aux1/aux2)+aux3
		
		v2drk = (vbz*cic)*dexp(-modk*w)*(1.0/ew)*(1.0/modk)
	
	end if
	
!novos inputs : w,ez,r0

end function v2drk

!potential ohno (DOI: 10.1103/PhysRevB.89.235410)

real(8) function v2dohono(kpt1,kpt2,ngrid,rlat,ediel,w,ez,tolr)

!	implicit none

	double precision,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	double precision,parameter:: pi=acos(-1.)

	double precision,dimension(3) :: kpt1,kpt2,vkpt
	double precision,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid
	double precision :: modk
	double precision,dimension(3) :: ediel
	double precision :: vc,lc,vbz,w,ez,tolr
	
	
	call vcell2D(rlat,vc)
	call modvec(kpt1,kpt2,modk)
	
	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))	
	!vbz= 1.0/vc !/((ngrid(1)*ngrid(2)*ngrid(3)))	
	
	if (modk .lt. tolr) then

		v2dohono = 0.0
	else 
		
		!v2dohono = (vbz*cic)*(1.0/(4.0*pi*pi))*dexp(-w*modk)*(1.0/(ez*modk))
		v2dohono = (vbz*cic)*dexp(-w*modk)*(1.0/(ez*modk))
	
	end if
	
	!novos inputs : w

end function v2dohono

!1D truncated potential  (DOI: 10.1103/PhysRevB.73.205119)

real(8) function v1dt(kpt1,kpt2,ngrid,rlat,tolr,lc)

!	implicit none
	double precision,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	double precision,parameter:: pi=acos(-1.)
	
	double precision,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid
	double precision,dimension(3) :: kpt1,kpt2,vkpt
	double precision :: tolr,lc
	double precision :: modk,vc,vbz
	
	double precision :: c1,c2,qyz,res
	double precision :: j0,j1,k0,k1

	call vcell3D(rlat,vc)
	call modvec(kpt1,kpt2,modk)
	
	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))
	vkpt=kpt1-kpt2
	
	qyz=sqrt( (vkpt(2)*vkpt(2)) +(vkpt(3)*vkpt(3)) )
	
	
	if ( (abs(vkpt(1)) .lt. tolr) .and. (qyz .ge. tolr) ) then
	
		call integralJ0pot(lc,qyz,res)
	
		v1dt = (vbz*cic)*res
	
	else if  ( (abs(vkpt(1)) .lt. tolr) .and. (qyz .lt. tolr) ) then
	

		v1dt = (vbz*cic)*((lc*lc)/16.0)*(2.0*log(lc)-1.0) 
	
	else

		call caljy0( qyz*lc, j0, 0 )
		call caljy1( qyz*lc, j1, 0 )
		call calck0( abs(vkpt(1))*lc, k0, 1 )
		call calck1( abs(vkpt(1))*lc, k1, 1 )
		
		c1 = qyz*lc*j1*k0
		
		c2 = abs(vkpt(1))*lc*j0*k1
		
		v1dt =	(vbz*cic)*(1.0+c1+c2)
	
	end if


end function v1dt



  
  end module module_coulomb
  


