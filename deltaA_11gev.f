	real function deltaA_11gev()

	include '?'

	real eg
	real egamma

	real ebeam,pbeam,lambda,elaser,gamma
	real egamma_max,a,rho
	real first,second,term1,term2,term3,asy,delta
	real Me,re,hbarc,e,c,pi,alpha ! some constants
	real betabar,s
	real betacm,gammacm,elasercm,ebeamcm,betaecm,costhcm


	parameter(hbarc=1.9732858E-11) ! MeV-cm
	parameter(Me=0.51099006)
	parameter(re=2.818e-13)
	parameter(c=2.9979e10)
	parameter(e=1.6022e-19)
	parameter(pi=3.141592654)
	parameter(alpha=1.0/137.036)


	lambda=532e-7  ! cm
	ebeam = 11000.0 ! MeV


	elaser = hbarc*2*pi/lambda
	gamma = ebeam/Me
	a = 1.0/(1.0+4.0*gamma*elaser/Me)
	egamma_max = 4.0*gamma**2*a*elaser

	rho = Egam_dep/egamma_max
	egamma=rho*egamma_max

	pbeam = sqrt(ebeam**2-Me**2)


	first = rho**2*(1.0-a)**2/(1.0-rho*(1.-a))
	second = (1-rho*(1+a))/(1-rho*(1-a))

	term1 = 1.0/(first+1.0+second**2)
	term2 = 1.0-rho*(1+a)
	term3 = 1.0-1.0/(1.0-rho*(1.-a))**2

* start new stuff to calcualate deltaA
* From Denner and Dittmaier Nuc.Phys.B 540 (1999) 58-86
	betabar=sqrt(1.0-Me**2/ebeam**2)
	betacm=(betabar*ebeam-elaser)/(ebeam+elaser) !eqn 2.10
	s=Me**2+2.0*elaser*ebeam*(1.0+betabar) ! eqn. 2.11
c	gammacm=1.0/sqrt(1.0-betacm**2)
	gammacm=(ebeam+elaser)/sqrt(s)
	elasercm=elaser*gammacm*(1.0+betacm) ! eqn. 2.12
	ebeamcm=sqrt(elasercm**2+Me**2)      ! eqn. 2.2
	betaecm=elasercm/ebeamcm             ! eqn. 2.2
	costhcm=(gammacm*elasercm-egamma)/
     >       (elasercm*gammacm*betacm) ! eqn. 2.13
	if(rho.le.1.0) then
	   if(rho.lt.elaser/egamma_max) then
	      delta=0.0
	      asy=0.0
	   else
	      delta=alpha/pi *(3.0*costhcm-1.)/(4.*(betaecm+costhcm)) ! eqn 3.9
	      asy = term1*term2*term3
	   endif
	else
	   delta=0.0
	   asy=0.0
	endif

C asyBorn = asy
C asyRad = asy*(1+delta)
C delta_A = asy*(1+delta)-asy = asy*delta

	deltaA_11gev = asy*delta

	return 
	end
