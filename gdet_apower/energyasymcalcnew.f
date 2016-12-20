      real function energyasymcalcnew(nentries)

      include '?'
C X = asymmetry, W = weight = cross section weighted energy


      integer count,ncontribute
      real W,X,WX,X2,W2,W2X,W2X2
      real sumW,sumX,sumWX,sumX2,sumW2,sumW2X2,sumW2X

      real asym,easym,easym2

      data count/0/
      data ncontribute/0/
      data sumW/0/
      data sumWX/0/
      data sumX/0/
      data sumX2/0/
      data sumWX/0/
      data sumW2X/0/
      data sumW2X2/0/


      count = count+1

c      if(Egam_dep.gt.1.0E-6) then
         ncontribute=ncontribute+1
         X=apower
c         X=apower
         W=Egam_dep*1000.0*wXsect*wLumin/20000.0
*         W=Eparticle(1)*1000.0*wXsect*wLumin/20000.0
         WX=W*X
         X2=X*X
c    stuff for calculating errors
         W2=W*W
         W2X=W*W*X
         W2X2=W*W*X*X
         
         sumX=sumX+X
         sumX2=sumX2+X2
         sumW=sumW+W
         sumWX=sumWX+WX
         sumW2=sumW2+W2
         sumW2X=sumW2X+W2X
         sumW2X2=sumW2X2+W2X2

      
         asym=sumWX/sumW
c      endif
      if(count.eq.nentries) then
         easym2 = sumW2X2/sumW**2 + sumWX**2*sumW2/sumW**4 
     >   - 2.*sumWX*sumW2X/sumW**3 

         if(easym2.gt.0.) then
            easym=sqrt(easym2)
         else
            easym=0
         endif
         write(6,*) 'number of contributing events', ncontribute
         write(6,*) 'Analyzing power is', asym,'+/-',easym
         write(6,*) 'fractional unc. is', 100.0*easym/asym,'%'
      endif

      return
      end

