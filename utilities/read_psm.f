      SUBROUTINE read_psm(nhg,conv,ng,fileformat)

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
!***********************************************************************
! - THIS SUBROUTINE INPUTS soil moisture DATA.

!   converted to 64 bit Jul 27/11 NK
!   Replaced ntype by classcount-1

! - list of arguments:

!   o - rain( , )real*4    rain gauge precipitation data
!   i - ng       int       number of rain gauges
!   i - nhg      int       number of hours of data
!   i - conv     real*4    conversion factor
!   i - fln( )   char*12   file names
!   o - sta( )   cmplx    gauge locations (grid square numbers)
!   o - ewg( )   real*4   east-west utm coordinates of gauges
!   o - sng( )   real*4   north-south utm coordinates of gauges
!   o - ng       int      number of rain gauges
!   o - radw     int      utm coordinate of west side of grid
!   o - rads     int      utm coordinate of south side of grid
!   o - rads     int      utm coordinate of south side of grid
!   o - gname( ) char*12  rain gauge names
!   i - fln( )   char*30  file names

!***********************************************************************

c      USE area2  
c      use area3
c      USE area12
c      USE area16
c      use areawfo

      use area_watflood

!      include 'debug.for'

      logical       :: exists
      INTEGER       :: iAllocateStatus,ng,nhg,classcount_local
	REAL          :: conv
      CHARACTER(10):: junk1
	character(20) :: junk
      character(30) :: newfilename
      character(10) ::  fileformat
      character(4)  :: yyyy
!character(2)  :: yy,mm,dd,hh

! OPEN INPUT FILE




      INQUIRE(FILE=fln(39),EXIST=exists)
      IF(exists)THEN
        open(unit=269,file=fln(39),status='unknown',iostat=ios)
        if(ios.ne.0)then
          print*,'Problems opening unit=269 file name = ',fln(39)
          print*
          stop 'Program aborted in rdsm.for @ 46'
        endif
        print*,'Opened fln(39)=',fln(39)(1:40)
      else
        print*,'***',fln(39),'*** not found'
        print*
        stop 'Program aborted in rdsm.for @ 49'
      endif

!#######################################################################
!     programmd for pt2 file type 14/07/06 nk

      call find_filetype(39)

      print*
      print*,'filetype=',filetype

      if(filetype.eq.'pt2')then
        do while(junk1.ne.':Name     ')
          read(269,*)junk1
          write(51,*)junk1
        end do
        read(269,*)junk1
        write(51,*)junk1
        read(269,*,iostat=ios)junk,coordsys2
          write(*,*)junk,coordsys2
	  if(coordsys2.eq.'UTM       ')then
          read(269,*,iostat=ios)junk,zone2
          write(*,*)junk,zone2
	  endif
        read(269,*,iostat=ios)junk,datum2
          write(*,*)junk,datum2
	  read(269,*)junk1
	  write(51,*)junk1
	  read(269,*)junk1
	  write(51,*)junk1
	  read(269,*)junk1
	  write(51,*)junk1
	  read(269,*)junk,conv
	  write(51,*)junk,conv
	  read(269,*)junk
	  write(51,*)junk
        classcount_local=0

        do while(junk1(1:10).ne.':EndHeader')
          read(269,*)junk1
          if(junk1(1:10).eq.':Attribute')then
	      classcount_local=classcount_local+1
            write(51,*)junk1,classcount_local
          endif
        end do

        ng=256
	  classcount_local=classcount_local/2-1
	  write(51,*)'# classes found =',classcount_local

!       added Nov. 17, 2010  nk
!       check for allocation adequacy
	  if(classcount_local.gt.classcount)then
	    print*,'Class count in the psm file',classcount_local
		print*,'is greater than the number of'
	    print*,'classes in the shed file',classcount
          stop 'Program aborted in read_psm @ 114'
	  endif

        ng=256
! TS - ALLOCATION OF AREA16A ARRAYS.
        if(iall.eq.0)then
          allocate(smcrag(ng),xsta(ng),ysta(ng),
     *       sng(ng),sngmin(ng),ewg(ng),smc_class(ng,classcount),
     *       ewgmin(ng),gname(ng),stat=iAllocateStatus)
          iall=iall+1
          if(iAllocateStatus.ne.0) STOP
     *     '***Allocation of AREA16A arrays failed in INPRAGA***'
        endif

! TS - ALLOCATION OF ARRAY RAIN().
        allocate(rrain(ng,nhg),stat=iAllocateStatus)
        if(iAllocateStatus.ne.0) STOP
     *     '***Allocation of AREA16A arrays failed in INPRAGA***'

        is=1
        do while(ios.eq.0)
          read(269,*,iostat=ios)xsta(is),ysta(is),gname(is),
     *                    (smc_class(is,ii),ii=1,classcount_local)
	    write(51,*)is
	    write(51,*)xsta(is),ysta(is),gname(is),
     *                    (smc_class(is,ii),ii=1,classcount_local)
          is=is+1
        end do

        ng=is-2

!       convert to grid coordinates:
        do n=1,ng
          ysta(n)=(ysta(n)-yorigin)/ydelta        
          xsta(n)=(xsta(n)-xorigin)/xdelta
        end do

!      write(51,*)'ios=',ios,'       ',ng

!     CLOSE INPUT FILE

      close(unit=39)
	write(51,*)'Closed unit 269. Filename =  ',fln(39)
	write(51,*)
	write(*,*)'Closed unit 269. Filename =  '
	write(*,*)fln(39)(1:72)

      return

      endif


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$





! INPUT RAIN GAUGE GRID DATA

!     GRID ON THE YYMMDD.RAG FILE IS NOW IGNORED
!     THE GRID SIZE IS TAKEN FROM THE basin/bsnm.pdl FILE
!     SO THE GRID SIZE CAN BE EASILY CHANGED


      write(51,*)'Opened  ',fln(39)(1:40),'  in rdsm.for'
	write(51,*)

      read(269,*)junk
      write(51,*)junk

        if(iopt.eq.2)print*,' In ragmet @ checkpoint a1'
        read(269,*)junk,fileformat        ! :FileType
        write(51,*)junk,fileformat        ! :FileType
!       want a different variable name for coordsys as in shed
        read(269,*)junk,coordsys1
        write(51,*)junk,coordsys1
        newformat=1

        read(269,*)junk,datum1
        write(51,*)junk,datum1
        read(269,*)junk,zone1
        write(51,*)junk,zone1
        read(269,*,iostat=ios)junk        !#
        write(51,*)junk
        read(269,*,iostat=ios)junk,startdate
        write(51,*)junk,startdate
        read(269,*,iostat=ios)junk,starttime
        write(51,*)junk,starttime
        read(269,*,iostat=ios)junk        !#
        read(269,*,iostat=ios)junk,ng     
        write(51,*)junk,ng
        read(269,*,iostat=ios)junk,ii  ! not implemented yet
        classcount=ii+1
        write(51,*)junk,classcount-1
        read(269,*,iostat=ios)junk,conv   ! conversion to cms
        write(51,*)junk,conv
        read(269,*,iostat=ios)junk        !#
        write(51,*)junk

! TS - ALLOCATION OF AREA16A ARRAYS.
        if(iall.eq.0)then
          allocate(smcrag(ng),xsta(ng),ysta(ng),
     *       sng(ng),sngmin(ng),ewg(ng),smc_class(ng,classcount),
     *       ewgmin(ng),gname(ng),stat=iAllocateStatus)
          iall=iall+1
          if(iAllocateStatus.ne.0) STOP
     *     '***Allocation of AREA16A arrays failed in INPRAGA***'
        endif

! TS - ALLOCATION OF ARRAY RAIN().
        allocate(rrain(ng,classcount),stat=iAllocateStatus)
        if(iAllocateStatus.ne.0) STOP
     *     '***Allocation of AREA16A arrays failed in INPRAGA***'


!       READ THE GAUGE LOCATIONS NAMES AND FUNCTION COEFFICIENTS
        read(269,*,iostat=ios)(gname(n),n=1,ng)
        write(51,*)(gname(n),n=1,ng)
        if(ios.ne.0)then
          print*,' error reading precip station name n=',n
          print*,  (sng(n),n=1,ng)
          STOP ' program aborted in ragmet.for @ 166'
        endif
        read(269,*,iostat=ios)(xsta(n),n=1,ng)
        write(51,*)(xsta(n),n=1,ng)
        if(ios.ne.0)then
          print*,' error reading x coordinate n=',n
          print*,(xsta(n),n=1,ng)
          STOP ' program aborted in ragmet.for @ 178'
        endif
        read(269,*,iostat=ios)(ysta(n),n=1,ng)
        write(51,*)(ysta(n),n=1,ng)
        if(ios.ne.0)then
          print*,' error reading y coordinate n=',n
          print*,(ysta(n),n=1,ng)
          STOP ' program aborted in ragmet.for @ 172'
        endif
        read(269,*,iostat=ios)junk         !#
        write(51,5005)junk                   !#
        read(269,*,iostat=ios)junk         !: endHeader
        write(51,5005)junk                        !: endHeader
        if(iopt.eq.2)print*,' In rdsm @ checkpoint a2'
        write(51,*)
!       convert to grid coordinates:
        do n=1,ng
          ysta(n)=(ysta(n)-yorigin)/ydelta        
          xsta(n)=(xsta(n)-xorigin)/xdelta
        end do



      write(51,*)
      write(51,98001)'gname   ',(gname(n),n=1,ng)
      write(51,98002)'ysta    ',(ysta(n),n=1,ng)
      write(51,98002)'xsta    ',(xsta(n),n=1,ng)
      write(51,*)



!       INPUT THE INITIAL SOIL MOISTURE
!       read(269,5000)(smcrag(i),i=1,ng)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! fix fix

!       this is just to read in the new format
!       still needs to be programmed for multiple classes
	  write(51,*)'Initial soil moisture at precip station:'
        do ii=1,classcount-1
          read(269,*,iostat=ios)(smc_class(i,ii),i=1,ng)
          write(51,*)(smc_class(i,ii),i=1,ng)              ! writes once only
        end do
        if(ios.ne.0)then
	    print*,'ios=',ios,' end of file found with only',ii
	    print*,'lines of soil moisture data'
	    print*,classcount-1,' lines of data were expected (# classes)'
	    print*
          do  jj=ii,classcount-1
	      do i=1,ng
              smc_class(i,jj)=smc_class(i,ii)
	      end do
	    end do
        endif   
        write(51,*)
        write(51,*)'Echo soil moisture data:'

6000    format('+', i5,'/',i5)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        deltat=1
        startdate='unknown   '
        starttime='unknown   '


! CLOSE INPUT FILE

      close(unit=39)
	write(51,*)'Closed unit 269. Filename =  ',fln(39)(1:40)
	write(51,*)
	write(*,*)'Closed unit 269. Filename =  ',fln(39)(1:40)

!  FORMATS

1000  format(10f8.2)
1100  format(2i5,1x,a12)
1101  format(5x,2i5,1x,a12)
1200  format(4i5,1x,a12)
1201  format(5x,4i5,1x,a12)
1202  format(3i10,1x,a12,2f10.2)
2000  format(2i5,f5.0)
2001  format(2i5,f5.2)
2002  format(a25, '     ignored - in inprag.for')
2003  format(a25)
2100  format(2i5,1x,a12)
!	note this format is to comply with the vbasic i/o menus
5000  format(1200f8.0)
 5005 format(a20,i5)
98001 format(a10,999a10)
98002 format(a10,999f10.2)


      return

      END SUBROUTINE read_psm



