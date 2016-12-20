      function asum_fit_rho(x) 
      DOUBLE PRECISION FITPAD(24),FITFUN
      COMMON/HCFITD/FITPAD,FITFUN
       
      real en,me,pi,hbarc,lambda,alpha,ro,k,gamma,a
      real cr1,cr2,cr3,cr,rho,a1,a2,a3,as,distance
      real maxdist,cedge,p2,numberofparam,rhotest
      real rhc1,rhc2,rhc3,rhc4
 
      common/KCWORK/VECTOR(100)
      
      dimension x(1)

      numberofparam = VECTOR(1)      
      en = VECTOR(2)
      me = VECTOR(3)
      pi = VECTOR(4)
      hbarc = VECTOR(5)
      lambda = VECTOR(6)

      
      alpha=1.0/137.035999
c      ro=alpha*hbarc/me  
      ro=2.8179402894
  
      k=2*pi*hbarc/(lambda*1e6)
      gamma=en/me
      a=1.0/(1.0+4*k*gamma/me)
       
      rho=x(1)
 
*old   rho=-0.11426E-04+0.58985E-01*distance
*      rho=rho-0.13994E-03*distance**2
*      rho=rho+0.30729E-06*distance**3
  
      cr1=(rho*rho*(1.0-a)*(1.0-a))/(1.0-rho*(1.0-a)) 
      cr2=((1.0-rho*(1.0+a))/(1.0-rho*(1.0-a)))**2    
      cr3=2.0*pi*ro*ro*a			      
      cr=cr3*(cr1+1.0+cr2)			      

      a1=2.0*pi*ro*ro*a/cr			      
      a2=(1.0-rho*(1.0+a))			      
      a3=1.0/(1-rho*(1-a))**2			      
      as=a1*a2*(1-a3)				      
        
c      write(*,*)x(1),distance,rho,rhot,rhc1,rhc2,rhc3,rhc4
     
        asum_fit_rho=fitpad(1)*as       
c        asum_fit=0.85*as       

      END                
