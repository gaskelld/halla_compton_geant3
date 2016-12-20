	subroutine gufld (location, field)

c This routine returns the magnetic field as a 3-vector evaluated at a
c given point
c Inputs:
c	location(3): coordinates of point in MRS (cm)
c Outputs:
c	field(3): components of field (kG)
c
c PW 20-Oct-93
c  - original code
c
c Richard Jones June 7, 2003
c  - replaced the use of fortran "record" statement with standard fortran-77,
c    and changed column alignment to conform to the fortran source form.
c  - replaced trig intrinsics based on degrees (sind,cosd,tand..) with
c    local versions because g77 does not support them
c
c Patricia Solvignon April 8, 2010
c  - add option to use field map instead of a single value for the dipole field


	implicit none

	include 'user.inc'


*	structure /point/
*	real x, y, z
*	end structure

*	record /point/ location, field

        integer k
	real location(3),field(3)
        real posdip(3)
        real ymap_limit,zmap_limit
        real xphy_limit,yphy_limit
        real xdip(4),ydip(4),zdip(4)
	real Iscale
        logical fieldmap       ! true  --> use field mad values
                               ! false --> use one value for the field everywhere in the dipole
c
c-----------------
c
c
	if (dipole_map.gt.0.1) then
	   fieldmap=.true.
c	   Iscale = (dipole_current*52.335+65.660)/
c     >              (110.0*52.335+65.660)
c	   Iscale=14.734568/14.97348*1.00477
	   Iscale=14.734568/14.97348*dipole_current
c	   write(6,*) 'scaling field map by',Iscale
	else
	   fieldmap=.false.
	endif
        
	field(1) = 0				! initially no field
	field(2) = 0
	field(3) = 0

	if(dipole_map.gt.0.1 .and. dipole_map.lt.1.5) then
	   zmap_limit = 70.0
	   ymap_limit=5.25
	elseif(dipole_map.gt.1.5) then
	   zmap_limit = 105.0
	   ymap_limit = 2.0
	endif

        xphy_limit = 22.0/2.0
        yphy_limit = dipole_height/2.0

	if(fieldmap) then
	   if(yphy_limit.gt.ymap_limit) then
	      write(6,*) 'Danger Will Robinson: Y-range for field map'
	      write(6,*) 'larger than physical dipole size. Field map'
	      write(6,*) 'routine may not work properly - I am outta here'
	      stop
	   endif
	endif

        xdip(1)  = 0.0
        xdip(2)  = 0.0
        xdip(3)  = 0.0
        xdip(4)  = 0.0

        ydip(1)  = chicane_drop-dipole_h/2.0
        ydip(2)  = 0.0+dipole_h/2.0
        ydip(3)  = 0.0+dipole_h/2.0
        ydip(4)  = chicane_drop-dipole_h/2.0

c        zdip(1)  = -chicane_length/2.0 + dipole_length/2.0
        zdip(1)  = - chicane_mid/2. - dipole_length 
     &             - chicane_drift - dipole_length/2.
        zdip(2)  = -chicane_mid/2.0 - dipole_length/2.0
        zdip(3)  = chicane_mid/2.0 + dipole_length/2.0
c        zdip(4)  = chicane_length/2.0 - dipole_length/2.0
        zdip(4)  = chicane_mid/2. + dipole_length 
     &             + chicane_drift + dipole_length/2.


        if (fieldmap) then
           posdip(1) = location(1)
           posdip(2) = location(2)

           k=0
 100	   k=k+1
           if (location(3).ge.(zdip(k)-zmap_limit).and.
     &         location(3).le.(zdip(k)+zmap_limit)) then
               if (posdip(1).lt.(xdip(k)-xphy_limit).or.
     &             posdip(1).gt.(xdip(k)+xphy_limit).or.
     &             posdip(2).lt.(ydip(k)-yphy_limit).or.
     &             posdip(2).gt.(ydip(k)+yphy_limit)) then

                  field(1) = 0.0
               else      
                  posdip(3) = location(3) - zdip(k)
		  posdip(2) = location(2) - ydip(k)
c                 if (posdip(3).gt.0.0) posdip(3) = -posdip(3)
		  if(dipole_map.gt.1.5) then
		     if(k.eq.2.or.k.eq.3) then
			posdip(3)=-posdip(3)
		     endif
		     call get_field_shim(posdip,field)
		  else
		     call get_field(posdip,field)
		  endif
                 field(1) = Iscale*field(1)
                 field(2) = Iscale*field(2)
                 field(3) = Iscale*field(3)

c                 if (k.eq.2.or.k.eq.3) field(1) = -1.0*field(1)
                 if (k.eq.2.or.k.eq.3) field(1) = -1.0*field(1)

               endif
           else
              if (k.lt.4) then
                 goto 100
              else 
                 field(1) = 0.0
                 goto 101
              endif
           endif

 101	   continue

           goto 11
        else
           goto 10
        endif


 10	continue
        if ((location(3).le.-chicane_length/2.0 + dipole_length+10.0)
     1                                 .or.
     2     (location(3).ge.chicane_length/2.0 - dipole_length-10.0))
     3  then
	   field(1) = dipole_field
	else
	   field(1) = -dipole_field
        endif

 11	continue

c	write(6,*) 'Field:', location(1), location(2), location(3), 
c     >      field(1)

	end

c----------
        subroutine get_field(posdip,field)
C Note that "x" in the field map refers to the dispersive direction, i.e., "y" in the lab system.
        implicit none
        integer i,j,iy,jz,ifound
        real posdip(3),field(3)
        real f(16,101)
        real A,B,C,D
        real field_tmp(101)
	logical first
        character line*357
        real ydipole(16),zdipole(101)
        data ydipole/-5.25,-4.55,-3.85, -3.15, -2.45, -1.75, -1.05,
     &               -0.35, 0.35, 1.05,  1.75,  2.45,  3.15,  3.85,
     &                4.55, 5.25/

        data zdipole/-70.0,-69.0,-68.0,-67.0,-66.0,-65.0,-64.0,-63.0,
     &               -62.0,-61.0,-60.0,-59.0,-58.0,-57.0,-56.0,-55.0,
     &               -54.0,-53.0,
     &               -52.0,-51.0,-50.0,-49.0,-48.0,-47.0,-46.0,-45.0,
     &               -44.0,-43.0,-42.0,-41.0,-40.0,-38.0,-36.0,-34.0,
     &               -32.0,-30.0,-28.0,-26.0,-24.0,-22.0,-20.0,-18.0,
     &               -16.0,-14.0,-12.0,-10.0, -8.0, -6.0, -4.0, -2.0,
     &                 0.0,  2.0,  4.0,  6.0,  8.0, 10.0, 12.0, 14.0,
     &                16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0,
     &                32.0, 34.0, 36.0, 38.0, 40.0, 41.0, 42.0, 43.0,
     &                44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0,
     &                52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0,
     &                60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0,
     &                68.0, 69.0, 70.0/ 

	data first/.true./
C No amnesia between calls!!!

        save

c       Find corresponding lines in dipole_fieldmap.dat

        do i=1,16
           if (posdip(2).le.ydipole(i).and.
     &         posdip(2).gt.ydipole(i-1)) then
              iy = i
              goto 901              
           endif
        enddo

 901	continue

	if (first) then ! read in field map first time through
	   write(6,*) 'reading in field map'
	   open(unit=62,file='dipole_fieldmap.dat',status='old')
	   do j=1,101
	      read(62,*)f(1,j), f(2,j), f(3,j), f(4,j), f(5,j), f(6,j),
     &               f(7,j), f(8,j), f(9,j),f(10,j),f(11,j),f(12,j),
     &               f(13,j),f(14,j),f(15,j),f(16,j)
	   enddo
	   close(62)
	   first=.false.
	endif

c	write(6,*) 'posdip',posdip(1),posdip(2),posdip(3)
	do j=1,101
	   if (posdip(3).le.zdipole(j).and.
     &         posdip(3).gt.zdipole(j-1)) then 
              jz = j
              goto 902
           endif
        enddo

 902	continue

        ! interpolate wrt y  (simple linear interpolation)
        A = (ydipole(iy) - posdip(2))  / (ydipole(iy)-ydipole(iy-1))
        B = (posdip(2) - ydipole(iy-1))/ (ydipole(iy)-ydipole(iy-1))
        field_tmp(jz-1) = A*f(iy-1,jz-1) + B*f(iy,jz-1)
        field_tmp(jz)   = A*f(iy-1,jz)   + B*f(iy,jz)

        ! interpolate wrt z  (simple linear interpolation)
        C = (zdipole(jz) - posdip(3))  / (zdipole(jz)-zdipole(jz-1))
        D = (posdip(3) - zdipole(jz-1))/ (zdipole(jz)-zdipole(jz-1))
        field(1) = C*field_tmp(jz-1) + D*field_tmp(jz)

c        field(1) = field(1)*1E-3
        field(1) = field(1)*1E-3
        field(2) = 0.0
        field(3) = 0.0

        goto 1000 ! jump the screen printing of the following control values 
        write(6,*)'************************************************'
        write(6,*)'interp in x:', 
     &              ydipole(iy-1),'<',posdip(2),'<',ydipole(iy)
        write(6,*)'field at (jz-1):', 
     &              f(iy-1,jz-1),field_tmp(jz-1),f(iy,jz-1)
        write(6,*)'field at (jz):',   
     &              f(iy-1,jz), field_tmp(jz), f(iy,jz)
        write(6,*)'interp coeff A & B:',A,B
        write(6,*)'-------'
        write(6,*)'interp in z:', 
     &              zdipole(jz-1),'<',posdip(3),'<',zdipole(jz)
        write(6,*)'field at:',  
     &              field_tmp(jz-1),field(1)*1E3,field_tmp(jz)
        write(6,*)'interp coeff C & D:',C,D

 1000	continue
        return
        end

c----------
        subroutine get_field_shim(posdip,field)
C Note that "x" in the field map refers to the dispersive direction, i.e., "y" in the lab system.
        implicit none
        integer i,j,k,iy,jz,ifound
        real posdip(3),field(3)
        real f(41,1051)
        real A,B,C,D
        real field_tmp(1051)
	real ydum,xdum,zdum,bx,bz
	logical first
        character line*357
        real ydipole(41),zdipole(1051)

	data first/.true./
C No amnesia between calls!!!

        save

c initialize grid
	if(first) then
	   do i=1,41
	      ydipole(i)=-2.0+(i-1)*0.1
	   enddo
	   
	   do i=1,1051
	      zdipole(i)=-105.0+(i-1)*0.2
	   enddo
	endif


c       Find corresponding lines in dipole_fieldmap.dat

        do i=1,41
           if (posdip(2).le.ydipole(i).and.
     &         posdip(2).gt.ydipole(i-1)) then
              iy = i
              goto 901              
           endif
        enddo

 901	continue

	if (first) then ! read in field map first time through
	   write(6,*) 'reading in shim field map'
	   open(unit=62,file='dipole_shim_short.txt',status='old')
	   do j=1,41
	      do k=1,1051
		 if(k.lt.126) then
		    f(j,k)=0.0
		 else
		    read(62,*)ydum,xdum,zdum,bx,f(j,k),bz
		    f(j,k)=-1.0*f(j,k)
		 endif
c		 write(6,*) j,k,ydum,zdum,f(j,k)
	      enddo
	   enddo

	   close(62)
	   first=.false.
	endif


c	write(6,*) 'posdip',posdip(1),posdip(2),posdip(3)
	do j=2,1051
	   if (posdip(3).le.zdipole(j).and.
     &         posdip(3).gt.zdipole(j-1)) then 
              jz = j
              goto 902
           endif
        enddo

 902	continue

        ! interpolate wrt y  (simple linear interpolation)
        A = (ydipole(iy) - posdip(2))  / (ydipole(iy)-ydipole(iy-1))
        B = (posdip(2) - ydipole(iy-1))/ (ydipole(iy)-ydipole(iy-1))
        field_tmp(jz-1) = A*f(iy-1,jz-1) + B*f(iy,jz-1)
        field_tmp(jz)   = A*f(iy-1,jz)   + B*f(iy,jz)

        ! interpolate wrt z  (simple linear interpolation)
        C = (zdipole(jz) - posdip(3))  / (zdipole(jz)-zdipole(jz-1))
        D = (posdip(3) - zdipole(jz-1))/ (zdipole(jz)-zdipole(jz-1))
        field(1) = C*field_tmp(jz-1) + D*field_tmp(jz)

c        field(1) = field(1)*1E-3
        field(1) = field(1)*1E-3
        field(2) = 0.0
        field(3) = 0.0

c	write(6,*) 'posdip',posdip(1),posdip(2),posdip(3),field(1)

        goto 1000 ! jump the screen printing of the following control values 
        write(6,*)'************************************************'
	write(6,*)'y and z indicies:  ',iy,jz
	write(6,*) f(iy-1,jz-1),f(iy-1,jz),f(iy,jz-1),f(iy,jz)
        write(6,*)'interp in x:', 
     &              ydipole(iy-1),'<',posdip(2),'<',ydipole(iy)
        write(6,*)'field at (jz-1):', 
     &              f(iy-1,jz-1),field_tmp(jz-1),f(iy,jz-1)
        write(6,*)'field at (jz):',   
     &              f(iy-1,jz), field_tmp(jz), f(iy,jz)
        write(6,*)'interp coeff A & B:',A,B
        write(6,*)'-------'
        write(6,*)'interp in z:', 
     &              zdipole(jz-1),'<',posdip(3),'<',zdipole(jz)
        write(6,*)'field at:',  
     &              field_tmp(jz-1),field(1)*1E3,field_tmp(jz)
        write(6,*)'interp coeff C & D:',C,D

 1000	continue
        return
        end
