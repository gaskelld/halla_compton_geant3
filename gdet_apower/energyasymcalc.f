      real function energyasymcalc(nentries)

      include '?'

      integer count,iasy,ifill
      real egplustot,egminustot
      real asym(2),asym_burst
      real egtot,asy_egtot
      real eg,egp,egm,aeg
      real eg_burst,egp_burst,egm_burst


      data egplustot/0/
      data egminustot/0/
      data egtot/0/
      data egp_burst/0/
      data egm_burst/0/
      data eg_burst/0/
      data count/0/
      data iasy/0/

      ifill=1000

      count = count+1
      iasy = iasy+1

      eg=Egam_dep*1000.0*wXsect*wLumin/50000.0
      aeg=apower*eg
      egp=Egam_dep*1000.0*(1.0+(apower/0.85))*wXsect*wLumin/50000.0
      egm=Egam_dep*1000.0*(1.0-(apower/0.85))*wXsect*wLumin/50000.0

      egplustot=egplustot+egp
      egminustot=egminustot+egm
      egtot = egtot+eg

      asym(1)=(egplustot-egminustot)/(egplustot+egminustot)
      asy_egtot = asy_egtot+aeg
      asym(2)=asy_egtot/egtot

      eg_burst=eg_burst+eg
      egp_burst=egp_burst+egp
      egm_burst=egm_burst+egm

      if(iasy.eq.ifill) then
         asy_burst=(egp_burst-egm_burst)/(egp_burst+egm_burst)
         call hfill(9000,asy_burst,0.0,1.0)
         iasy=0
         egp_burst=0.0
         egm_burst=0.0
         eg_burst=0.0
      endif


      if(count.eq.1e5) then
c         write(6,*) egplustot,egminustot,egtot,egtot2,asym(1),asym(2)
      endif

      if(count.eq.nentries) then
         write(6,*) egplustot,egminustot,egtot,asym(1),asym(2)
      endif



      return
      end

