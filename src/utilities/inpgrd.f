        SUBROUTINE inpgrd(ixr,iyr,ng)

!***********************************************************************
!    Copyright (C) 1987-2018 by Nicholas Kouwen  
        
!    This file is part of WATFLOOD (R)      
        
!    WATFLOOD(R) is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    any later version.

!    WATFLOOD is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.

!    You should have received a copy of the GNU Lesser General Public License
!    along with WATFLOOD.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

! - THIS SUBROUTINE INPUTS RAIN GAUGE GRID DATA.
!
! Only the first section is needed flr ragmet. Just the extend of the domain is needed
! domain is needed to make the .met and .tem files.
! The domain can be anything larger than the .shd domain. 
!
! FIX: in the future, check that the .met domain .ge. than .shd domain 
!
! Modified by Tricia Stadnyk - September 2000
! Converted common blocks to modules and added dynamically allocated
! run-time arrays as part of Fortran 90 conversion.
!
! - list of arguments:
!
!   o - sta( )   cmplx    gauge locations (grid square numbers)
!   o - ewg( )   real*4   east-west utm coordinates of gauges
!   o - sng( )   real*4   north-south utm coordinates of gauges
!   o - ixr      int      number of east-west grid squares
!   o - iyr      int      number of north-south grid squares
!   o - ng       int      number of rain gauges
!   o - radw     int      utm coordinate of west side of grid
!   o - rads     int      utm coordinate of south side of grid
!       rgrd     int      spl grid size in km
!   o - gname( ) char*12  rain gauge names
!   i - fln( )   char*12  file names
!
!***********************************************************************

c      USE area2
c      USE area12
c      USE area16
c	use areawfo

      use area_watflood

      INTEGER       :: iAllocateStatus
      character(10) :: junk1
      character(20) :: junk

!      include 'debug.for'

! OPEN INPUT FILE

      open(unit=33,file=fln(3),recl=flen,status='old',iostat=ios)
      if(ios.eq.9)then
        print*,' file ',fln(3)(1:40),' not found   ios=',ios
	  print*,'It looks like there was no bsnm.pdl file'
	  print*,'If you have already run bsn.exe, you need to rename'
	  print*,'new.pdl to bsnm.pdl (where bsnm = your basin nime)'
        stop 'program aborted in inpgrd.for @40'
      endif
      if(ios.ne.0)then
        print*,' problems opening pdl file '
        print*,fln(3)(1:72)
        print*,' ios=',ios
	  print*,' filename =',fln(3)(1:40)
        stop 'program aborted in inpgrd.for @ 44'
      endif

!             if(ng.eq.1)then

!       NG IS USED AS FLAG HERE:
!		ng = 	1 READ RAIN GAUGES
!		ng = 	0 READ SNOW COURSE LOCATIONS
!		ng = -1 READ CLIMATE STATION LOCATIONS

!       INPUT RAIN GAUGE GRID DATA

         if(iopt.ge.1)write(51,*)'Opened ',fln(3)(1:50)
d        write(*,*)'Opened ',fln(3)(1:50)
!5000   format(' reading the rain gage grid file: ',a30/)

      nnprint=0
      newformat=0
!     rev. 9.1.58  Jul.  12/04  - NK: New header for the .shd file
      read(33,3001,iostat=ios)junk       !    #
      print*
       if(iopt.ge.1)write(51,*)'First line pdl file =  ',junk
      write(*,*)'First line pdl file =  ',junk

	if(junk(1:1).eq.'#')then
!       new format
!       new format
!       new format
        print*,'New format ',fln(3)(1:30),' pdl file found'
        print*,'Please look in ???_info.txt for error hints'
	  print*
	  newformat=1
         if(iopt.ge.1)write(51,*)' '
         if(iopt.ge.1)write(51,*)'Echo data from    ',fln(3)(1:40) 
         if(iopt.ge.1)write(51,*)' '
!        read(33,3002,iostat=ios)junk,junk1
        read(33,*,iostat=ios)junk,junk1
	   if(iopt.ge.1)write(51,*)junk,junk1
        if(ios.ne.0)then
          print*,'Problems reading line 2 in file: ',fln(3)
          print*
	    stop 'Program aborted in inpgrd @ 119'
	  endif
c        read(33,3004,iostat=ios)junk,coordsys1
        read(33,*,iostat=ios)junk,coordsys1
	   if(iopt.ge.1)write(51,3004)junk,coordsys1
        if(ios.ne.0)then
          print*,'Problems reading line 5 in file: ',fln(3)
          print*
	    stop 'Program aborted in inpgrd @ 119'
	  endif
c        read(33,3004,iostat=ios)junk,datum1
        read(33,*,iostat=ios)junk,datum1
	   if(iopt.ge.1)write(51,3004)junk,datum1
        if(ios.ne.0)then
          print*,'Problems reading line 6 in file: ',fln(3)
          print*
	    stop 'Program aborted in inpgrd @ 135'
	  endif
c        read(33,3004,iostat=ios)junk,zone1
!       fixed for reading zone  June 8/15  NK
        if(coordsys1(1:3).eq.'UTM')then
          read(33,*,iostat=ios)junk,zone1
	  else 
          read(33,*,iostat=ios)junk
          zone1=0
        endif
	   if(iopt.ge.1)write(51,3004)junk,zone1
        if(ios.ne.0)then
          print*,'Problems reading line 7 in file: ',fln(3)
          print*
	    stop 'Program aborted in inpgrd @ 149'
	  endif
        read(33,3002,iostat=ios)junk       !    #
	   if(iopt.ge.1)write(51,3002)junk
        if(ios.ne.0)then
          print*,'Problems reading line 8 in file: ',fln(3)
          print*
	    stop 'Program aborted in inpgrd @ 156'
	  endif
        read(33,*,iostat=ios)junk,xorigin
c        read(33,3003,iostat=ios)junk,xorigin
	   if(iopt.ge.1)write(51,3003)junk,xorigin
        if(ios.ne.0)then
          print*,'Problems reading line 9 in file: ',fln(3)
          print*
	    stop 'Program aborted in inpgrd @ 164'
	  endif
c        read(33,*,iostat=ios)junk,yorigin
        read(33,*,iostat=ios)junk,yorigin
	   if(iopt.ge.1)write(51,3003)junk,yorigin
        if(ios.ne.0)then
          print*,'Problems reading line 9 in file: ',fln(3)
          print*
	    stop 'Program aborted in inpgrd @ 172'
	  endif
c        read(33,3003,iostat=ios)junk       !    #
        read(33,*,iostat=ios)junk       !    #
	   if(iopt.ge.1)write(51,3003)junk
        if(ios.ne.0)then
          print*,'Problems reading line 10 in file: ',fln(3)
          print*
	    stop 'Program aborted in inpgrd @ 180'
	  endif
c        read(33,3001,iostat=ios)junk,xcount
        read(33,*,iostat=ios)junk,xcount
	   if(iopt.ge.1)write(51,3001)junk,xcount
        if(ios.ne.0)then
          print*,'Problems reading line 11 in file: ',fln(3)
          print*
	    stop 'Program aborted in inpgrd @ 188'
	  endif
c        read(33,3001,iostat=ios)junk,ycount
        read(33,*,iostat=ios)junk,ycount
	   if(iopt.ge.1)write(51,3001)junk,ycount
        if(ios.ne.0)then
          print*,'Problems reading line 12 in file: ',fln(3)
          print*
	    stop 'Program aborted in inpgrd @ 196'
	  endif
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!       eventually, imax & jnmax should be replaced everywhere in the code
c        imax=ycount
c	  jmax=xcount
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
c        read(33,3003,iostat=ios)junk,xdelta
        read(33,*,iostat=ios)junk,xdelta
	   if(iopt.ge.1)write(51,3003)junk,xdelta
        if(ios.ne.0)then
          print*,'Problems reading line 13 in file: ',fln(3)
          print*
	    stop 'Program aborted in inpgrd @ 209'
	  endif
c        read(33,3003,iostat=ios)junk,ydelta
        read(33,*,iostat=ios)junk,ydelta
	   if(iopt.ge.1)write(51,3003)junk,ydelta
        if(ios.ne.0)then
          print*,'Problems reading line 14 in file: ',fln(3)
          print*
	    stop 'Program aborted in inpgrd @ 217'
	  endif

        if(iopt.ge.1)write(51,*)'Finished reading header ',fln(3)(1:40)

!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!       eventually, imax & jnmax should be replaced everywhere in the code
c        imax=ycount
c	  jmax=xcount

!       eventually, grdn & grde should be replaced everywhere 
!       in the code by xdelta & ydelta   <<< fix
        grde=xdelta/1000.
        grdn=ydelta/1000.  

!       added ll separation Jul. 27/04  nk
        if(coordsys1.eq.'LATLONG   ')then
          iymin=int(yorigin*60.0)
          iymax=int((yorigin+ycount*ydelta)*60.0) 
          jxmin=int(xorigin*60.0)  
          jxmax=int((xorigin+xcount*xdelta)*60.0)  
        else
          jxmin=int(xorigin/1000.)
          jxmax=jxmin+grde*(xcount-1)
          iymin=int(yorigin/1000.)
	    iymax=iymin+grdn*(ycount-1)
        endif

        read(33,3001,iostat=ios)junk       !    #
	   if(iopt.ge.1)write(51,3001)junk
        if(ios.ne.0)then
          print*,'Problems reading line 26 in file: ',fln(3)
          print*
	    stop 'Program aborted in inpgrd @ 250'
	  endif

        ixr=xcount
        iyr=ycount
         if(iopt.ge.1)write(51,1203)iyr,ixr

      else     !   if(junk.eq.........
!       old format
!       old format
!       old format
        rewind 33
        print*,'New format ',fln(3),' pdl file NOT found'
        print*,'Please look in ???_info.txt for error hints'
	  print*
         if(iopt.ge.1)write(51,*),'Old format found'
	   if(iopt.ge.1)write(51,*)
        read(33,1000,iostat=ios)rgrd,rads,radn,radw,rade,
     *    latdegmin,latminmin,latdegmax,latminmax,
     *    londegmin,lonminmin,londegmax,lonminmax,grdn,grde
         if(iopt.ge.1)write(51,1000)rgrd,rads,radn,radw,rade,
     *    latdegmin,latminmin,latdegmax,latminmax,
     *    londegmin,lonminmin,londegmax,lonminmax,grdn,grde
!        write(*,1000)rgrd,rads,radn,radw,rade,
!     *    latdegmin,latminmin,latdegmax,latminmax,
!     *    londegmin,lonminmin,londegmax,lonminmax,grdn,grde

        print*,'fln(3)= ',fln(3)
        print*,'iostat= ',ios
        print*
        if(ios.ne.0)then
          print*,' Problem reading ',fln(3)
          print*,' The above line = values found'
          stop 'Program aborted in inpgrd @ 235'
        endif

!       FOR LAT-LONG COORDINATES RGRD,RADS,RADN,RADW,RADE ARE SET EQUAL
!       TO -1 AND THE ALTERNATE LAT LONG BORDERS ARE USED INSTEAD

        if(rgrd.le.0)then
          llflg='y'
        else
          llflg='n'
        endif

        if(llflg.eq.'y')then
           if(iopt.ge.1)write(51,1004)
     *    latdegmin,latminmin,latdegmax,latminmax,
     *    londegmin,lonminmin,londegmax,lonminmax,grdn,grde

!         CALCULATE THE GRID CORNERS IN MINUTES 
          rads=latdegmin*60+latminmin
          radn=latdegmax*60+latminmax
          radw=londegmin*60+lonminmin
          rade=londegmax*60+lonminmax
          rgrdn=int(grdn)
          rgrde=int(grde)
        else
*          write(*,1000)rgrd,rads,radn,radw,rade
           if(iopt.ge.1)write(51,1000)rgrd,rads,radn,radw,rade

!         FOR UTM GRID SPACING IS EQUAL IN NS & EW DIRECTIONS
          grdn=float(rgrd)
          grde=float(rgrd)
          rgrdn=float(rgrd)
          rgrde=float(rgrd)
        endif

        ixr=int(float(rade-radw)/grde)+1
        iyr=int(float(radn-rads)/grdn)+1
         if(iopt.ge.1)write(51,1203)iyr,ixr

      endif   !   if(junk.eq.........

       if(iopt.ge.1)write(51,*)'finished reading the pdl file header'
      print*,'finished reading the pdl file header'
      print*

!     ragmet.f and tmp.f only needs the grid information. 
!     The station locations & remainder
!     of the data is in the raing\yyyymmdd.rag file


	print*,'ng=',ng

      if(ng.eq.2.or.ng.eq.-1)then        ! FOR RAGMET ONLY
!         added ng.eq.-1  Nov 18/07 nk
        close(unit=33, status='keep')
         if(iopt.ge.1)write(51,*)'closed u 33,file name: ',fln(3)(1:30)
	   if(iopt.ge.1)write(51,*)
	   if(iopt.ge.1)write(51,*)'-----------------------------------------'
	   if(iopt.ge.1)write(51,*)
d       print*,'closed unit 33'
d	  print*,'file name:  ',fln(3)(1:50)
        print*
        return
      endif

!     *****************************************************************

!     for EVENTS.FOR

      if(ng.eq.1)then
        read(33,1001,iostat=ios)ng
         if(iopt.ge.1)write(51,*)'no of temperature gages = ',ng

!     TS - ALLOCATION OF AREA16A ARRAYS ADDED
      if(iall.eq.0)then
      allocate(smcrag(ng),xsta(ng),ysta(ng),sng(ng),sngmin(ng),ewg(ng),
     *         ewgmin(ng),gname(ng),stat=iAllocateStatus)
      iall=iall+1
      if(iAllocateStatus.ne.0) STOP
     *   '**Allocation failed for AREA16A arrays in inpgrda l118**'
      endif

        if(llflg.ne.'y')then
!         FOR UTM 
          do n=1,ng
            read(33,1100,iostat=ios)sng(n),ewg(n),gname(n)
!            write(*,1101)sng(n),ewg(n),gname(n)
             if(iopt.ge.1)write(51,1101)sng(n),ewg(n),gname(n)
          end do
        else
!         FOR LAT-LONG
          do n=1,ng
            read(33,1200,iostat=ios)
     *        sng(n),sngmin(n),ewg(n),ewgmin(n),gname(n)
!            write(*,1201)sng(n),sngmin(n),ewg(n),ewgmin(n),gname(n)
             if(iopt.ge.1)write(51,1201)sng(n),sngmin(n),
     *                            ewg(n),ewgmin(n),gname(n)
            sng(n)=sng(n)*60+sngmin(n)
            ewg(n)=ewg(n)*60+ewgmin(n)
          end do
        endif

        if(ios.ne.0)then
          print*,' Problem reading ',fln(3)
          print*,' The above line = values found'
          stop 'Program aborted in inpgrd @ 331'
        endif

        do n=1,ng
!          ysta(n)=float(sng(n)-(rads+rgrdn/2))/grdn
!          xsta(n)=float(ewg(n)-(radw+rgrde/2))/grde
          ysta(n)=(float(sng(n))-(float(rads)+rgrdn/2.0))/grdn
          xsta(n)=(float(ewg(n))-(float(radw)+rgrde/2.0))/grde
!          write(*,1202)n,sng(n),ewg(n),gname(n),ysta(n),xsta(n)
           if(iopt.ge.1)write(51,1202)n,sng(n),ewg(n),
     *                          gname(n),ysta(n),xsta(n)
        end do

      elseif(ng.eq.0)then

         if(iopt.ge.1)write(51,5020)fln(3)
 5020   format(' '/' reading the snow course location file: ',a30/)
         
!       READ OVER RAIN GAUGE GRID DATA
         
!        read(33,1000,iostat=ios)
        read(33,1001,iostat=ios)ng
         if(iopt.ge.1)write(51,1001)ng

!     TS - ALLOCATION OF AREA16A ARRAYS ADDED
        if(iall.eq.0)then
          allocate(smcrag(ng),xsta(ng),ysta(ng),sng(ng),sngmin(ng),
     *         ewg(ng),ewgmin(ng),gname(ng),stat=iAllocateStatus)
          iall=iall+1
          if(iAllocateStatus.ne.0) STOP
     *     '**Allocation failed for AREA16A arrays in inpgrda l118**'
        endif

        do n=1,ng
          read(33,*,iostat=ios)
        end do
         
        read(33,1003,iostat=ios)ng,nclass
         if(iopt.ge.1)write(51,1003)ng,nclass

        if(llflg.ne.'y')then
!         FOR UTM 
          do n=1,ng
            read(33,1100)sng(n),ewg(n),gname(n)
*            write(*,1101)sng(n),ewg(n),gname(n)
             if(iopt.ge.1)write(51,1101)sng(n),ewg(n),gname(n)
          end do
        else
!         FOR LAT-LONG
          do n=1,ng
            read(33,1200,iostat=ios)
     *        sng(n),sngmin(n),ewg(n),ewgmin(n),gname(n)
!            write(*,1201)sng(n),sngmin(n),ewg(n),ewgmin(n),gname(n)
             if(iopt.ge.1)write(51,1201)sng(n),sngmin(n),
     *                               ewg(n),ewgmin(n),gname(n)
            sng(n)=sng(n)*60+sngmin(n)
            ewg(n)=ewg(n)*60+ewgmin(n)
          end do
        endif

        if(ios.ne.0)then
          print*,' Problem reading ',fln(3)
          print*,' The above line = values found'
          stop 'Program aborted in inpgrd @ 392'
        endif

        do n=1,ng
!          ysta(n)=float(sng(n)-(rads+rgrdn/2))/grdn
!          xsta(n)=float(ewg(n)-(radw+rgrde/2))/grde
          ysta(n)=(float(sng(n))-(float(rads)+rgrdn/2.0))/grdn
          xsta(n)=(float(ewg(n))-(float(radw)+rgrde/2.0))/grde
           if(iopt.ge.1)write(51,1202)n,sng(n),ewg(n),
     *                                    gname(n),ysta(n),xsta(n)
        end do

      elseif(ng.eq.-1)then
!	  SKIP OVER THE RAIN AND SNOW LOCATIONS AND READ 
!	  THE CLIMATE STATION LOCATIONS:

!         read(33,1000)

!	PASS OVER THE RAIN GAUGES

        read(33,1001,iostat=ios)ng
!        write(*,1001)ng
         if(iopt.ge.1)write(51,1001)ng

        do n=1,ng
          read(33,*,iostat=ios)
        end do

!	PASS OVER THE SNOW COURSES 

        read(33,1003,iostat=ios)ng
!        write(*,1003)ng
         if(iopt.ge.1)write(51,1003)ng
        do n=1,ng
          read(33,*)
        end do

!	AND FINALLY READ THE TEMPERATURE STUFF

        read(33,1001,iostat=ios)ng
!        write(*,1001)ng
         if(iopt.ge.1)write(51,1001)ng

!     TS - ALLOCATION OF AREA16A ARRAYS ADDED
        if(iall.eq.0)then
          allocate(smcrag(ng),xsta(ng),ysta(ng),sng(ng),sngmin(ng),
     *         ewg(ng),ewgmin(ng),gname(ng),stat=iAllocateStatus)
        iall=iall+1
        if(iAllocateStatus.ne.0) STOP
     *   '**Allocation failed for AREA16A arrays in inpgrda l118**'
        endif

        if(llflg.ne.'y')then
!         FOR UTM 
          do n=1,ng
            read(33,1100,iostat=ios)sng(n),ewg(n),gname(n)
!            write(*,1101)sng(n),ewg(n),gname(n)
             if(iopt.ge.1)write(51,1101)sng(n),ewg(n),gname(n)
          end do
        else
!         FOR LAT-LONG
          do n=1,ng
            read(33,1200,iostat=ios)
     *        sng(n),sngmin(n),ewg(n),ewgmin(n),gname(n)
!            write(*,1201)sng(n),sngmin(n),ewg(n),ewgmin(n),gname(n)
             if(iopt.ge.1)write(51,1201)sng(n),sngmin(n),
     *                            ewg(n),ewgmin(n),gname(n)
            sng(n)=sng(n)*60+sngmin(n)
            ewg(n)=ewg(n)*60+ewgmin(n)
          end do
        endif

        do n=1,ng
!          ysta(n)=float(sng(n)-(rads+rgrdn/2))/grdn
!          xsta(n)=float(ewg(n)-(radw+rgrde/2))/grde
          ysta(n)=(float(sng(n))-(float(rads)+rgrdn/2.0))/grdn
          xsta(n)=(float(ewg(n))-(float(radw)+rgrde/2.0))/grde
!          write(*,1202)n,sng(n),ewg(n),gname(n),ysta(n),xsta(n)
           if(iopt.ge.1)write(51,1202)n,sng(n),ewg(n),gname(n),
     *      ysta(n),xsta(n)
        end do
        
      endif

!        if(ios.ne.0)then
!          print*,' Problem reading ',fln(3)
!          print*,' The above line = values found'
!          stop 'Program aborted in inpgrd @ 477'
!        endif

      close(unit=33, status='keep')
	write(51,*)'closed unit 33, file name:  ',fln(3)(1:30)
	write(51,*)
	write(51,*)'-----------------------------------------------'
	write(51,*)
d      print*,'closed unit 33'
d      print*,'file name:  ',fln(3)(1:50)
d      print*

      return

! FORMATS

1000  format(13i5,2f5.1)
1001  format(i5)
1002  format(' ',3i5,2x,a12,2f5.0/)
1003  format(2i5)
1004  format(8i5,2f8.1)
1100  format(2i5,1x,a12)
1101  format(5x,2i5,1x,a12)
1200  format(4i5,1x,a12)
1201  format(5x,4i5,1x,a12)
1202  format(3i10,1x,a12,2f10.2)
1203  format('iyr=',i5,'  ixr=',i5)

 3001 format(a20,i12)
 3002 format(2a20)
 3003 format(a20,f12.3)
 3004 format(a20,a10)


      END SUBROUTINE inpgrd

