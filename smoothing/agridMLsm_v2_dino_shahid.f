c--- program agridMLsm.v2.f
c --- f95 or gfortran. Needs iosubs
c--- calculates gridded 10**a-values from catalog
c--- this version uses unequal sampling periods
c--- written by Art Frankel  6/94
c -- v2 mods: Steve Harmsen 1/2007.
c-- to run: agridMLsm.v2 parfile
c-- try: agridMLsm.v2 WUSagrid.in
	parameter (nsrcmx=10000,ngrmx=800000,nmagmx=12)
c nsrcmx = maximum # of source records in input catalog
c ngrmx = maximum #  of source grid points in output matrix
c integer*4  to insure that iosubs gets correct information.
      integer*4 readn,nsrc
      real magmin
      logical isok
      character*72 name,name2
      real, dimension (nmagmx) :: ymag
      integer, dimension (nmagmx) :: iyr
      real, dimension (ngrmx) :: asum,asumnew
c      write(6,*) "enter name of input file"
	if(iargc().lt.1)stop"Usage: agridMLsm input.file"
900   format(a)
      call getarg(1,name)
101      inquire(file=name,exist=isok)
      if(isok)then
      open(unit=1,file=name,status='old')
      else
      write(6,50)'Input file not found, please try again: '
50	format(a,$)
      read *,name
      goto 101
      endif
c      open(unit=2,file=name2,status='new')
      write(6,50) "for sites: min lat, max lat of box out: "
      read(1,*) ymin, ymax
	write(6,*)ymin,ymax
      write(6,50) "for sites: min lon, max lon of box out: "
      read(1,*) xmin, xmax
	write(6,*)xmin,xmax
      write(6,50) "increments dlat,dlon for source cell: "
      read(1,*) dsy, dsx
	write(6,*)dsy, dsx
c v2 change: make the nsx, nsy calculations more robust with nint
      nsx= nint((xmax-xmin)/dsx) +1
      nsy= nint((ymax-ymin)/dsy) +1
      nsrc= nsx*nsy
      write(6,*) "Source box: ",nsx,nsy,nsrc
      if(nsrc.gt.ngrmx)stop 'please increase ngrmx to accomodate'
      write(6,*) "enter number of magnitude categories"
      read(1,*) nmag
      if(nmag.gt.nmagmx)stop'number of mags exceeds nmagmx'
      write(6,*) "enter min mag and min year of completeness"
      do 11 i=1,nmag
 11   read(1,*) ymag(i),iyr(i)
c
c
      write(6,*) "enter end year of catalog"
      read(1,*)iyrend 
      write(6,50) "enter name of catalog file: "
      read(1,900) name
      write(6,900)name
      inquire(file=name,exist=isok)
      if(.not.isok)write(6,*)'file not found: ',name
      open(unit=2, file=name,status='old')
      write(6,*) "enter minimum mag,dmag, bval"
      read(1,*) magmin,dmag,bval
      dmag= dmag/2.
c-- above for conversion from cumulative to incremental a value
c-- Herrmann's formula is +-dmag
c
      sum= 0.
	asum = 0.	!use vector arithmetic
      atotal=0.
      do 1 i=1,nsrcmx
      read(2,*,end=999,err=999) iyear,xlat,xlon,depth,xmag
      if(depth.gt.50) goto 1
c      if((depth.lt.30).or.(depth.gt.50)) goto 1
      do 12 ii=1,nmag
      if((xmag.ge.ymag(ii)).and.(iyear.ge.iyr(ii))) then
c 2nd place where nint may improve the robustness of the calculation. v2 mod.
         iy= nint((ymax-xlat)/dsy)
         ix= nint((xlon-xmin)/dsx)
      if((iy.lt.0).or.(iy.ge.nsy).or.(ix.lt.0).or.(ix.ge.nsx)) go to 1
c      write(6,*) ix,iy
         n= iy*nsx + ix +1
         asum(n) = asum(n)+1.
         atotal= atotal+1
c      write(6,*) iy,ix,n,xlat,xlon
         go to 1
         endif
 12   continue
 1    continue
 999  continue
      fac= alog10(10.**(bval*dmag)-10.**(-bval*dmag))
      write(6,*) 'fac ',fac
      beta= bval*2.303
      a1= 0.
      a2= 0.
      do 14 i=1,nmag
      a1= a1+ exp(-beta*(ymag(i)+0.5))
 14   a2= a2 + (iyrend-iyr(i))*exp(-beta*(ymag(i)+0.5))
      fac2= a1/a2
      write(6,*) 'fac2=a1/a2 ',fac2
c----above is Herrmann's formula to convert from cumul. to interval a-value
      atotal2=0.
      do 3 i=1,nsrc
      if(asum(i).ne.0.) then
c         write(6,*) asum(i)
         asum(i)= asum(i)*fac2
         atotal2= atotal2+asum(i)
         asum(i)= alog10(asum(i)) + bval*magmin
         asum(i)= asum(i)+ fac
         asum(i)= 10.**asum(i)
         sum= sum+asum(i)
c         write(6,*) asum(i)
c         read(5,*) idum
         endif
   3  continue
      write(6,*) 'atotal, atotal2,sum ',atotal, atotal2, sum
c      read(5,*) idum
      write(6,50) "enter correlation distance in km: "
      read(1,*) a
      write(6,*)a
      a2= a*a
      coef= 3.14159/180.
      avelat= (ymin+ymax)/2.
c-------smoothing part
c--- go out to 3 times correlation distance
c 3rd place where nint may improve the robustness of the calculation. v2 mod.
      lx= nint(3.*a/(111.1*cos(avelat*coef)*dsx))
      ly= nint(3.*a/(111.1*dsy))
      write(6,*) 'correlation distances grid units: ',lx,ly
      do 10 iy=1,nsy
c      write(6,*) iy
      do 10 ix=1,nsx
      indx1= ix-lx
      indy1= iy-ly
      n= (iy-1)*nsx+ix
      rx= xmin+(ix-1)*dsx
      ry= ymax-(iy-1)*dsy
      sum1= 0.
      sum2= 0.
      do 20 iyy= indy1,indy1+2*ly
      do 20 ixx= indx1,indx1+2*lx
       if((ixx.lt.1).or.(ixx.gt.nsx)) go to 20
       if((iyy.lt.1).or.(iyy.gt.nsy)) go to 20
       n2= (iyy-1)*nsx+ixx
       sx= xmin+(ixx-1)*dsx
       sy= ymax-(iyy-1)*dsy
       call delaz(sy,sx,ry,rx,dist,az,baz)
       dist2= dist*dist
c----- Gaussian smoothing
       sum1= asum(n2)*exp(-dist2/a2) + sum1
       sum2= exp(-dist2/a2)+sum2
 20    continue
c---- normalize by area under Gaussian
      asumnew(n)= sum1/sum2
c      write(6,*) asumnew(n)
 10   continue
      sum=0.
      do 35 i=1,nsrc
      sum= sum+ asumnew(i)
 35   continue
      write(6,*) 'summed mean rate of M0: ',sum
      read(1,900) name2
c---output smoothed file with 10**a values
      open(99,file='grid_rate_0_50km_50_0.1_grid.xyz',status='unknown')
      do 100 iy=1,nsy
      do 100 ix=1,nsx
      n= (iy-1)*nsx+ix
      rx= xmin+(ix-1)*dsx
      ry= ymax-(iy-1)*dsy
      write(99,313) rx,ry,asumnew(n)
 
 100  continue
c 313  format(F7.2,3X,F7.2,3X,E10.4)       
 313  format(F7.2,',',F7.2,',',E10.4)       
c      call openw(name2)
c      call putbuf2(asumnew,nsrc,readn)
c      write(6,801) nsrc,readn
 801  format('Write: ndata=',i8,' readn=',i8)
      end

      subroutine delaz(sorlat,sorlon,stnlat,stnlon,delta,az,baz)
      coef= 3.14159/180.
      xlat= sorlat*coef
      xlon= sorlon*coef
      st0= cos(xlat)
      ct0= sin(xlat) 
      phi0= xlon
      xlat= stnlat*coef
      xlon= stnlon*coef
      ct1= sin(xlat)
      st1= cos(xlat)
      sdlon= sin(xlon-phi0)
      cdlon= cos(xlon-phi0)
      cdelt= st0*st1*cdlon+ct0*ct1
      x= st0*ct1-st1*ct0*cdlon
      y= st1*sdlon
      sdelt= sqrt(x*x+y*y)
      delta= atan2(sdelt,cdelt)
      delta= delta/coef
      az= atan2(y,x)
      az= az/coef
      x= st1*ct0-st0*ct1*cdlon
      y= -sdlon*st0
      baz= atan2(y,x)
      baz= baz/coef
      delta= delta*111.11
      return
      end subroutine delaz        

