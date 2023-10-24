        PROGRAM moist

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
! RAGMETA - This program generates a precipitation data file from
!           rain gauge data.
!
! Modified by Tricia Stadnyk - January 2000
!  Converted common blocks to modules and added dynamically allocated
!  run-time arrays as part of Fortran 90 conversion.
!
!        subroutines called:
!                rdevt
!                inpraga
!                inpgrda
!                weighta
!                distra
!                outmeta
!
!  VERSION 2      Dec. 2000  - TS: Added dynamic memory allocation
!
!        variable list:
!     NG = NUMBER OF RAIN GAGE STATIONS
!     nhg = NUMBER OF HOURS OF RECORDED DATA
!
!   I - R( , )  REAL*4     radar precipitation data
!   I - nhg     INT        number of hours of data
!   I - xcount     INT        number of east-west grid squares
!   I - ycount     INT        number of north-south grid squares
!   I - SMCRAG  REAL*4     initial soil moisture content
!   I - MO      INT        month of year
!   I - CONV    REAL*4     conversion factor
!   I - DATE    CHAR*12    event identifier
!   I - fln     CHAR*30    file names

!***********************************************************************


c        USE area1
c        USE area2
c        USE area3
c        USE area4
c        USE area5
c        USE area6
c        USE area10
c        USE area11
c        USE area12
c        USE area16
c        USE areaet
c        USE areamelt
c	  USE areawet
c        use areawfo

      use area_watflood

      use ef_module

        INTEGER         :: iAllocateStatus,iDeallocateStatus
        INTEGER       :: ng,nhg
      	REAL          :: conv
        CHARACTER(14)   :: date
        DIMENSION       :: smc5(16)
        INTEGER(kind=2) :: nrad,I
        INTEGER(2) status1
        CHARACTER(2) buf
        CHARACTER(1) stopflg
        character(10) ::  fileformat,time

        iall=0

! GET THE DISTANCE WEIGHTING PARAMETER NW FROM THE PROGRAM ARGUMENT
! USED IN WEIGHT

!     USE DFLIB


      CALL GETARG(1, buf, status1)
      print*,status1
      print*, buf
      if(buf.eq.'-2')then
        stopflg='y'
      else
        stopflg='n'
      endif
      if(status1.ne.1)buf=' '
      print*, buf
      if(buf.ne.' ')then
        open(unit=99,file='junk',status='unknown')
        write(99,9901)buf
9901    format(a1)
        rewind 99
        read(99,9902)nw
9902    format(i2)
      else
        nw=2
      endif


       
        write(6,6001)
        write(6,6002)
        write(6,6003)
        write(6,6004) 
        write(6,6005) 
        write(6,6006) 
        write(6,6007) 
        write(6,6008) 
        write(6,6009) 
        write(6,6002)
        write(6,6001)

 6001 format(1x,'*******************************************************
     ***********************')
 6002 format(1x,'*                                                      
     *                     *')
 6003 format(1x,'*  ww  ww         ww  a    tttttttt ffffff ll       ooo
     *      ooo    ddddd   *')
 6004 format(1x,'*   ww  ww       ww  aaa   tttttttt ffffff ll      oooo
     *o    ooooo   dddddd  *')
 6005 format(1x,'*    ww  ww     ww  aa aa     tt    ff     ll     oo   
     *oo  oo   oo  dd   dd *')
 6006 format(1x,'*     ww  ww   ww  aaa aaa    tt    ffff   ll     oo   
     *oo  oo   oo  dd   dd *')
 6007 format(1x,'*      ww  ww ww  aaaaaaaaa   tt    ffff   ll     oo   
     *oo  oo   oo  dd   dd *')
 6008 format(1x,'*       ww  www  aa       aa  tt    ff     llllll  oooo
     *o    ooooo   dddddd  *')
 6009 format(1x,'*        ww  w  aa         aa tt    ff     llllll   ooo
     *      ooo    ddddd   *')


	print*,'********************************************************'
	print*,'*                                                      *'
	print*,'*    MM       MM   OOOOOO   II   SSSSSSS   TTTTTTTT    *'
	print*,'*    MMM     MMM  OOOOOOOO  II  SSSSSSSSS  TTTTTTTT    *'
	print*,'*    MMMM   MMMM  OO    OO  II  SS            TT       *'
	print*,'*    MM MM MM MM  OO    OO  II  SSSSSSSS      TT       *'
	print*,'*    MM  MMM  MM  OO    OO  II   SSSSSSSS     TT       *'
	print*,'*    MM       MM  OO    OO  II         SS     TT       *'
	print*,'*    MM       MM  OOOOOOOO  II  SSSSSSSSS     TT       *'
	print*,'*    MM       MM   OOOOOO   II   SSSSSSS      TT       *'
	print*,'*                                                      *'
	print*,'*                  WATFLOOD (TM)                       *'
	print*,'*         Version 9.9     July 22, 2022                *'
	print*,'*           (c) N. Kouwen, 1972-2015                   *'
	print*,'*                                                      *'
	print*,'********************************************************'
      print*

 

!      print*,'stopflg=    ',stopflg
!      print*,'nw=',nw
!      print*
!      pause


!         include 'fsublib.fi'          
!         character(128) :: arg
!         integer        :: argc,arglen
!         argc=iargc()
!         arglen=igetarg(0,arg)
!         arglen=igetarg(1,arg)
!         i=1
!         if(arglen.le.0)arglen=2

!      data ioflag/2/itogo/0/
      calling_program_flg='moist     '


! PRINT EVENT TITLE

      write(6,1000)

      open(unit=51,file='moist_info.txt',recl=flen,status='unknown',
     *       iostat=ios)
	if(ios.ne.0)then
	  print*,'Problem opening moist_info.txt in working directory'
	  print*
	  stop 'Program aborted in moist @ 110'
	endif       

! INPUT EVENT PARTICULARS

! TS - ALLOCATION STATEMENT TO DYNAMICALLY ALLOCATE ARRAY'S FOR AREA12

      allocate(fln(601),outfln(100),stat=iAllocateStatus)
      if (iAllocateStatus .ne. 0) STOP 
     *    '**Allocation failed for area12**' 

      iopt=1

! INPUT EVENT PARTICULARS

      fln(99)='event/event.evt'

c      call rdevt(date,conv1,scale,smc5,nhg,nhf)
      call read_evt(date,conv1,scale,smc5,nhg,nhf)

      mo=mo1   ! mo1 = in area 2 - used to be in arg list


! INPUT RAIN GAUGE FILE DATA

cc      ng=ioflag

! IOFLAG IS HARDWIRED IN RAGMET DATA STATEMENT
! FOR IOFLAG=2 GRID IS READ FROM basin\nsnm.rag FILE
! STATION LOCATIONS ARE READ FROM THE raing\yymmdd.rag FILE

! FOR IOFLAG=1 GRID IS READ FROM basin\bsnm.rag FILE
! STATION LOCATIONS ARE ALSO READ FROM THE basin\bsnm.rag FILE

! TS - ALLOCATIONS FOR SMCRAG,XSTA,YSTA,SNG,SNGMIN,EWG,EWGMIN,GNAME
!      OCCUR WITHIN SUBROUTINE INPGRDA.

!      call rdgrid(xcount,ycount,ng,npost)

c!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c       call rdshed()
c!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



      call find_filetype(1)

	if(IsFileTypeR2C(fln(1))) then
        print*
	  print*,'Gone to read ',fln(1)(1:40),' with read_shd_ef'
	  call read_shed_ef(31,1)	
	  print*,'Finished reading',fln(1)(1:40)
	  print*
	else
!     rev. 10.5.00 Nov.  07/22  = NK Modificatons for HYPE i/o 
!			call read_shed_hype(31,1)
!			call rdshed(31,1)
			call read_shed_ef(31,1)
	endif

      if(.not.allocated(s))then
       allocate(s(ycount,xcount),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *     'Error with allocation of s(,) in moist @ 213'
	endif

! TS - ALLOCATION FOR RAIN() OCCURS WITHIN SUBROUTINE INPRAGA.

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call read_psm(classcount,conv,ng,fileformat) ! read the point soil moisture
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! TS - ADDED ALLOCATION STATEMENT TO DYNAMICALLY ALLOCATE ARRAYS IN 
!      AREA16
!        Variables used in array allocations are defined in the 
!        following file(s):
!                    xcount,ycount  - inpgrda
!                    ng       - inpraga
!                    nhg      - inpevta, inpraga

      allocate(temp(ng),p(ycount,xcount),psum(ycount,xcount),
     *         pmin(ycount,xcount),w(ycount,xcount,ng),
     *         ntogo(nhg),stat=iAllocateStatus)
      if (iAllocateStatus .ne. 0) STOP  
     *    '**Allocation failed for AREA16A arrays in ragmet @ 153**'  

!     outarray is in areawfo
!     this has to be allocated before calling write_r2c
      allocate(outarray(ycount,xcount),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *    'Error with allocation of outarray in moist'      


! TS - MOVED DO LOOP FROM BEFORE LINE53 TO HERE BECAUSE PRECEEDED
!      ALLOCATION.

      author='moist                               '
      coordsys_temp=coordsys1
      zone_temp=zone1
	datum_temp=datum1
	xorigin_temp=xorigin
	yorigin_temp=yorigin
	xcount_temp=xcount
	ycount_temp=ycount
	xdelta_temp=xdelta
	ydelta_temp=ydelta
	name='Initial Soil Moisture                   '
	attribute_units='mm                                      ' 
      attribute_type='Course                                  '
      unit_conversion=conv
      source_file_name=fln(39)                                    


      open(unit=99,file='junk1',status='unknown')
      write(99,99005)year1,mo1,day1,hour1
99005 format(i4,'/',i2,'/',i2,2x,i2,':00:00')
      rewind 99
	read(99,99006)startdate,starttime
99006 format(2a10)
      close(unit=99,status='delete')

!     write the header
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      call write_r2c(267,37,1,classcount,0,ii,1)   
      call write_r2c(267,37,0,classcount,0,0,1)   
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if(conv.ne.1.0)then
        do ii=1,classcount
          do is=1,ng
            smc_class(is,ii)=smc_class(is,ii)*conv
          end do
        end do
        conv=1.0
      endif

      do ii=1,classcount
c      do ii=1,ntype
!     DISTRIBUTE THE SOIL MOISTURE
         write(51,*)
         write(*,*)
         write(51,1100)ii,classcount
         write(*,1100)ii,classcount
c         write(51,1100)ii,ntype
c         write(*,1100)ii,ntype

        do i=1,ycount
          do j=1,xcount
            p(i,j)=0.0
          end do
        end do

	
        do is=1,ng
          rrain(is,ii)=smc_class(is,ii)
        end do

        iflag=0
        lflag=0
        kflag=0

!       CHECK IF THERE ARE SOIL MOISTURE DATA AT EACH STA:

        do is=1,ng
          if(rrain(is,ii).ge.0.0)lflag=1
          if(rrain(is,ii).gt.0.0)kflag=1
        end do
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call weight(ii,ng,xcount,ycount,0.0)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call distr(ii,xcount,ycount,ng)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       OUTPUT SMC DATA
        rdr='moist'

!       convert local variables to write module variables
        do i=1,ycount
          do j=1,xcount
            outarray(i,j)=p(i,j)
          end do
        end do          

!      SUBROUTINE write_r2c(conv,un,fn,
!     *            no_frames,no_classes,frame_no,class_no,
!     *            no_signf)

!       write the data
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c        call write_r2c(267,37,1,ntype,1,ii,2)
c        call write_r2c(267,37,1,classcount,1,ii,2)
        call write_r2c(267,37,0,classcount,0,ii,2)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	end do

      close(unit=98,status='keep')

! TS - ADDED DEALLOCATION STATEMENT FOR AREA12 AND AREA16 ARRAYS


! FORMATS

1001  format(a30)
1002  format(10x,i5,1x,f5.2,1x,i6,1x,f5.1)
1003  format(' HOUR=',2i5)
1000  format(' Init Soil Moistgure DISTRIBUTION APPLICATION')
1100  format('+',' DISTRIBUTING moisture for class ',i5,' of ',i5,
     +' land cover classes')
1110  format(' ')
1300  format(1200f5.1)
1301  format(1200f5.2)
c6000  format(1200f5.0)
c6001  format(' Hour= -999                    Summary of rainfall:')
      
      print*
      print*,'DISCLAIMER'
	print*,'The WATFLOOD software and other material supplied' 
	print*,'in connection herewith is furnished by N. Kouwen and the' 
	print*,'University of Waterloo and is accepted by the' 
	print*,'user upon the express understanding that N. Kouwen' 
	print*,'or the University of Waterloo make no warranties, either' 
	print*,'express or implied, concerning the accuracy, completeness,' 
	print*,'reliability, usability, performance, or fitness for any' 
	print*,'particular purpose.' 
	print*
	print*,'The material is provided "as is". The entire risk as to' 
	print*,'its quality and performance is with the user.'
      print*
      
      stop 'Normal ending'
      

      if(stopflg.eq.'y')then
        print*,' normal ending'
        print*
        stop 
      else
        stop ' normal ending'
      endif
      
      END PROGRAM moist





