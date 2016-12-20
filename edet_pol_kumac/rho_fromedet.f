      function rho_fromedet(x) 
      DOUBLE PRECISION FITPAD(24),FITFUN
      COMMON/HCFITD/FITPAD,FITFUN
       
      real en,me,pi,hbarc,lambda,alpha,ro,k,gamma,a
      real cr1,cr2,cr3,cr,rho,a1,a2,a3,as,distance
      real maxdist,cedge,p2,numberofparam,rhotest
      real rhc1,rhc2,rhc3,rhc4
 
      common/KCWORK/VECTOR(100)
      
      dimension x(1)

      maxdist = VECTOR(1)
      cedge = VECTOR(2)
      p2 = VECTOR(3)
      numberofparam = VECTOR(4)      
      en = VECTOR(5)
      me = VECTOR(6)
      pi = VECTOR(7)
      hbarc = VECTOR(8)
      lambda = VECTOR(9)
      rhc1 = VECTOR(10)
      rhc2 = VECTOR(11)
      rhc3 = VECTOR(12)
      rhc4 = VECTOR(13)
      nstrips = VECTOR(14)
      strwid = VECTOR (15)
      
      alpha=1.0/137.035999
c      ro=alpha*hbarc/me   
      ro=2.8179402894
  
      k=2*pi*hbarc/(lambda*1e6)
      gamma=en/me
      a=1.0/(1.0+4*k*gamma/me)
       
      distance=maxdist-(strwid*1000.)*(cedge-X(1)-0.5)

      rho=(rhc1)+(rhc2)*distance
      rho=rho+(rhc3)*distance**2
      rho=rho+(rhc4)*distance**3
 
* correction
*      rho=rho+(-0.56950E-2-0.55302E-4*X(1))
*      rho=rho+(-0.31289E-2-0.67064E-4*X(1))

      rho_fromedet=rho

       
      END                
