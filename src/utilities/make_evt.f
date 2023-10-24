        PROGRAM make_evt

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
!	include 'debug.for'
 
!*******************************************************************
! EVENTSA - THIS PROGRAM CREATES AN EVENT FILE AND SETS UP NEW DATA
!           FILES : xxxx.met  xxxx.rag  xxxx.str xxxx.rel  xxxx.val
!           OR SETS AN EXISTING EVENT DATA SET AS THE ACTIVE FILES
!
! Modified by Tricia Stadnyk - September 2000
! Converted common blocks to modules and added dynamically allocated
! run-time arrays as part of Fortran 90 conversion.

! Modified Feb. 2006 to make a set of event files. NK

!
!	subroutines called:
!               inpgrda
!               inpsgaa
!               outevta
!
!  VERSION 2      Dec. 2000  - TS: Added dynamic memory allocation
!
!
! - list of arguments:
!
!   iall   = flag to prevent reallocation of inpgrda arrays
!   no     = # of stream gages
!   ng     = # of rain gages  
!   nr     = # of reservoirs
!   ns     = # of snow courses
!   smc    = soil moisture calculated at each rainfall gauge
!   smc5   = soil moistued for each sub-basin
!   weq    = water equivalent in mm
!   ns     = no. of snow courses
!   nclass = no of posts (stations) per snow course
!
!   o - sta( )   cmplx    gauge locations (grid square numbers)
!   o - ewg( )   real*4   east-west utm coordinates of gauges
!   o - sng( )   real*4   north-south utm coordinates of gauges
!   o - ixr      int      number of east-west grid squares
!   o - iyr      int      number of north-south grid squares
!   o - ng       int      number of rain gauges
!   o - radw     int      utm coordinate of west side of grid
!   o - rads     int      utm coordinate of south side of grid
!   o - gname( ) char*12  rain gauge names
!   i - fln( )   char*256 file names
!
!******************************************************************

c      USE area2   
c      USE area10
c      USE area12
c      USE area16

      use area_watflood

      CHARACTER(256) :: fln1,line
      CHARACTER(12) :: date
	character(13) :: start
	character(11) :: ext(50)
      CHARACTER(8)  :: bsnme,evtname,evtType
      CHARACTER(6)  :: dirn
      CHARACTER(5)  :: keepflg5,keepflg6,keepflg7,keepflg8,keepflg14
      CHARACTER(4)  :: yy1
      CHARACTER(2)  :: mm1,dd1,hh1 
      CHARACTER(1)  :: bs(50)              !,qs
      real*4        :: smc5(5),temp1(50,8784),rel(50,8784),weq(50,8784)
      INTEGER       :: iallocatestatus,ideallocatestatus
      real*4        :: smcavg
      integer       :: evtlength,no_months,mohours(24),noevents
	logical       :: firstpass

      DATA mohours/744,672,744,720,744,720,744,744,720,744,720,744,
     *             744,672,744,720,744,720,744,744,720,744,720,744/
	data firstpass/.true./


      iall=0

      keepflg5='y'
      keepflg6='y'
      keepflg7='y'
      keepflg8='y'
      keepflg14='y'
      evtClass='    '

! TS - ALLOCATION STATEMENT TO DYNAMICALLY ALLOCATE ARRAY
      allocate(fln(601),stat=iAllocateStatus)
      if(iAllocateStatus .ne. 0)STOP 
     *    '**Allocation failed for area12 in eventsa l86**' 

	print*,'********************************************************'
	print*,'*                                                      *'
	print*,'*                 WATFLOOD (R)                         *'
	print*,'*                                                      *'
	print*,'*                Program mk_cfg                        *'
	print*,'*                Version 10_2_15                       *'
	print*,'*                                                      *'
	print*,'*           (c) N. Kouwen, 1972-2020                   *'
	print*,'*                                                      *'
	print*,'********************************************************'
      print*
      print*,'Please see file evt_info.txt for information re: this run'


      open(unit=98,file='evt_info.txt',status='unknown')



	print*,'event selection program'
	print*,'warning: no damage yet, but if you enter the name'
	print*,'of an existing event, all old files by that name'
	print*,'and the series of events following'
	print*,'will be over written. enter ^c or ^break to stop'
      
      print*,'NEW 2022 - optional'
      print*,'enter special event class designation '
      print*,'e.g.   CaPA, RDRS etc.  8chr max'
      print*,'for legacy files enter:   na'
      read(*,*)evtType

	print*
	print*,'Enter the no of events to create:'
	read(*,*)noevents
	print*
	print*,'No. of months per event file (1 or 12)'
	print*,'If you are going to run long sequences - e.g. 40 years'
	print*,'use 12 months per event - i.e. yearly events'
	print*,'If you do not, you will have many many files.'
	print*,'Use 12 month events unless you actually'
	print*,'need month-long events for some reason'
	read(*,*)no_months

      ievtflg=0

	print*,'type in start of event - eg. yyyy mm dd hh'
	print*,'please stick with this convention so radar files work'

      read(*,7000)start
      yy1=start(1:4)
	mm1=start(6:7)
	dd1=start(9:10)
	hh1=start(11:13)

      open(unit=90,file='event\events_to_follow.txt',
     *                    	status='unknown',iostat=ios)

!     Got rid of temp files  Feb. 17/18 NK
c      open(unit=99,file='temp0.fls',status='unknown',iostat=ios)
!      write(*,*)(chars(i),i=1,filename_length-4),'_',filetype,'.',ext
      write(line,7012)yy1,mm1,dd1,hh1
c      rewind 99
      read(line,9905)evtname   ! importan to have the EXACT format here
99000 format(a8)
c      rewind 99
	read(line,7013)year1,month1,day1,hour1    ! area2
 7013 format(i4,3i2)
c      close(unit=99,status='keep')

!	write(*,7010)yy1,mm1,dd1,hh1
! 7010 format(a4,3(1x,a2))

!    ***************************************************************
!    ***************************************************************
!    ***************************************************************
!    ***************************************************************
!    ***************************************************************
      do id=1,noevents

      if(id.gt.1)then
        if(no_months.eq.1)then   ! no months per event
          month1=month1+1	 
          month1=mod(month1,12)  ! 12th month converts to 0
	    if(month1.eq.0)month1=12
     	    if(month1.eq.1)year1=year1+1
          open(unit=99,file='temp5.fls',status='unknown',iostat=ios)
          if(month1.le.9)then
            write(99,99013)year1,month1,day1,hour1
99013       format(i4,'0',i1,'01')
          else
            write(99,99014)year1,month1,day1,hour1
99014       format(i4,i2,'01')
          endif
          rewind 99
          read(99,9905)evtname   
	    rewind 99
	    read(99,99015)yy1,mm1,dd1
99015     format(a4,2a2)
        elseif(no_months.eq.12)then
          year1=year1+1
          open(unit=99,file='temp5.fls',status='unknown',iostat=ios)
          if(month1.le.9)then
            write(99,99013)year1,month1,day1,hour1
          else
            write(99,99014)year1,month1,day1,hour1
          endif
!          write(*,99013)year1,month1,day1,hour1
          rewind 99
          read(99,9905)evtname 
!          write(*,9905)evtname 
	    rewind 99
	    read(99,99015)yy1,mm1,dd1
        endif
      endif

!	write(*,9905)evtname
 9905 format(a8)
!	write(*,7014)year1,month1,day1,hour1
 7014 format(i4,3(1x,i2))

	if(no_months.eq.1)then
        evtlength=mohours(month1)
	  if(mod(year1,4).eq.0.and.month1.eq.2)evtlength=evtlength+24
	elseif(no_months.eq.12)then
        evtlength=mohours(month1)
	  do j=month1+1,month1+11
	    evtlength=evtlength+mohours(j)
	  end do
	  if(mod(year1,4).eq.0)evtlength=evtlength+24
	else
	  print*,'No of months per event not eq to 1 or 12'
	  print*
	  stop 'Program aborted in evt @ 130'
	endif

!     convert input & get starting month
      open(unit=99,file='temp1.fls',status='unknown')
      write(99,7021)yy1,mm1,dd1,hh1
 7021 format(a4,1x,3(a2,1x))
      REWIND 99
      read(99,7000)date
 7000 format(a13)
      REWIND 99
      read(99,7022)mo
 7022 format(3x,i2)
      REWIND 99
      close(unit=99,status='keep')


      if(id.eq.1)then
!     create the event name  yyyymmdd
      open(unit=99,file='temp2.fls',status='unknown')
      write(99,7012)yy1,mm1,dd1
 7012 format(a4,3a2)
      REWIND 99
      read( 99,7015)evtname
 7015 format(a8)
      close (unit=99,status='keep')
      endif


!      write(*,7016)evtname
 7016 format('                   event name = ',a8)


      if(ievtflg.eq.0)then
	  print*
	  print*,'What reporting time step would you like in files?'
	  print*,'This should not be shorted than the str file Dt'
	  read(*,*)deltat_report
	  print*
	  print*
        print*,' will you be running the snow melt  routines? y/n'
	  print*,' Note: temperature data needed for this option'
        read(*,5083)snwflg
 5083   format(a1)
        if(snwflg.eq.'y'.or.snwflg.eq.'y')then
            write(6,5084)
            write(6,5085)
            read(*,5086)snconv
        else
          snconv=1.0
        endif
 5084   format(' enter the snow conversion factor')
 5085   format(' e.g.  1.0 is snow wat. eq. in mm, 25. if in inches'/)
 5086   format(f5.0)
        print*
        print*,' will you be disaggregating precipitation? y/n'
	  read(*,5083)smrflg
	  if(smrflg.eq.'y')then
	    print*,'what is your disaggregation rate mm/hr ?'
		read(*,*)smearfactor
	  else
	    smearfactor=0.0
	  endif
        print*,' will you be running the evaporation routines? y/n'
	  print*,' Note: temperature data needed for this option'
        read(*,5083)vapflg
	  print*
        print*,' will you be using reservoir/lake inflow files? y/n'
	  read(*,5083)resinflg
        print*,' will you be using wetland coupling? y/n'
	  read(*,5083)wetflg
	  print*
        print*,' will you be using diversions? y/n'
	  print*,' Note: diversion data needed for this option'
        read(*,5083)divertflg
	  print*
        print*,' will you be running lake evaporation? y/n'
	  print*,' Note: diversion data needed for this option'
        read(*,5083)lakeflg
	  print*
        print*,' will you be using ice factors? y/n'
	  print*,' Note: diversion data needed for this option'
        read(*,5083)iceflg
	  print*
      endif

      if(ievtflg.eq.0)then
        write(6,5011)
        read(*,7009)bsnme
      endif
 5011 format(/' name of shd & par files: eg. gr10k, saug 8 char max'/)
 7009 format(a8)


!	create the file names for the event file
      open(unit=99,file='temp3.fls',status='unknown')

! COUNT THE NUMBER OF CHARACTERS IN THE NAME:

!     write(99,7009)bsnme
!     REWIND 6
!     read(99,6051)(bs(i),i=1,8)
!6051 format(a1)

      bs(1)=bsnme(1:1)
      bs(2)=bsnme(2:2)
      bs(3)=bsnme(3:3)
      bs(4)=bsnme(4:4)
      bs(5)=bsnme(5:5)
      bs(6)=bsnme(6:6)
      bs(7)=bsnme(7:7)
      bs(8)=bsnme(8:8)

!      do 11 i=1,8
!   11 write(6,6011)bs(i)
! 6011 format(10xa1)

      j=8
      do i=1,8
        if(bs(i).eq.' ')then
          j=i-1
          GO TO 10
        endif
      end do
 10   CONTINUE

! J IS THE NUMBER OF CHARACTERS IN THE BASIN NAME
! WRITE THE NEW NAMES TO A FILE WITHOUT PADDED BLANKS:
      ext(1)='_shd.r2c'
      ext(2)='_par.csv'
	ext(3)='_ch_par.r2c'
      ext(4)='.pdl'
      ext(5)='.sdc'
      ext(6)='.wqd'
      ext(7)='_diff.r2c'

      dirn = 'basin\'

! NOW PICK THE FORMAT THAT CORRESPONDS WITH THE CORRECT  
! NUMBER OF CHARACTERS AND WRITE THE .shd .par .rag & .str LINES
      do 9 i=1,7
         GO TO (1,2,3,4,5,6,7,8)j
    1    write(99,6001)dirn,bsnme,ext(i)
         GO TO 9
    2    write(99,6002)dirn,bsnme,ext(i)
         GO TO 9
    3    write(99,6003)dirn,bsnme,ext(i)
         GO TO 9
    4    write(99,6004)dirn,bsnme,ext(i)
         GO TO 9
    5    write(99,6005)dirn,bsnme,ext(i)
         GO TO 9
    6    write(99,6006)dirn,bsnme,ext(i)
         GO TO 9
    7    write(99,6007)dirn,bsnme,ext(i)
         GO TO 9
    8    write(99,6008)dirn,bsnme,ext(i)
    9 CONTINUE

 6001 format(a6,a1,a11)
 6002 format(a6,a2,a11)
 6003 format(a6,a3,a11)
 6004 format(a6,a4,a11)
 6005 format(a6,a5,a11)
 6006 format(a6,a6,a11)
 6007 format(a6,a7,a11)
 6008 format(a6,a8,a11)

      write(99,6065)evtname
      write(99,6066)evtname
      write(99,6067)evtname
      write(99,6068)evtname
      write(99,6069)evtname
      write(99,6070)evtname
      write(99,6071)evtname
      write(99,6072)evtname

 6065 format('raing\',a8,'.rag')
 6066 format('strfw\',a8,'.str')
 6067 format('resrl\',a8,'.rel')
 6068 format('snow1\',a8,'_crs.pt2')
 6069 format('raduc\',a8,'.rad')
 6070 format('radcl\',a8,'_met.r2c')
 6071 format('radar\',a8,'.scn')
 6072 format('raduc\',a8,'.clt')


! WRITE THE SDC FILE (THE 13TH FILE):
      GO TO (11,12,13,14,15,16,17,18)j
   11   write(99,6001)dirn,bsnme,ext(i)
        GO TO 19
   12   write(99,6002)dirn,bsnme,ext(i)
        GO TO 19
   13   write(99,6003)dirn,bsnme,ext(i)
        GO TO 19
   14   write(99,6004)dirn,bsnme,ext(i)
        GO TO 19
   15   write(99,6005)dirn,bsnme,ext(i)
        GO TO 19
   16   write(99,6006)dirn,bsnme,ext(i)
        GO TO 19
   17   write(99,6007)dirn,bsnme,ext(i)
        GO TO 19
   18   write(99,6008)dirn,bsnme,ext(i)
   19 CONTINUE

      write(99,6074)evtname
      write(99,6075)evtname
      write(99,6076)evtname
      write(99,6077)evtname
      write(99,6078)evtname
      write(99,6079)
      write(99,6080)
!      write(99,6081)evtname
      if(evtType(1:2).eq.'na')then
        write(99,6081)evtname
      else
        strLength = LEN_TRIM(evtType)
        write(99,6083)evtname,evtType
      endif  
 6074 format('tempr\',a8,'.tag')
 6075 format('tempr\',a8,'_tem.r2c')
 6076 format('strfw\',a8,'.gwt')
 6077 format('tempr\',a8,'.tmx')
 6078 format('snowg\',a8,'.dsn')
 6079 format('                ')
 6080 format('event\event.evt' )
 6081 format('event\',a8,'.evt')
 6083 format('event\',a8,'_',a<strLength>'.evt')

      REWIND 99
      

      read(99,6082)(fln(i),i=1,23)  ! raise this when adding files
      read(99,6082)fln1             ! reads the name of the file
      close (unit=99,status='keep')
 6082 format(a30)

! CHECK TO SEE IF THIS EVENT EXISTS AND PRINT A DIRE WARNING THAT
! EXISTING FILES WILL BE OVERWRITTEN!!!!!!!



! THIS FILE WILL BE READ BY WATLFD TO SAVE THE EVENT FILE FOR 
! FUTURE USE:

c      open(unit=99,file='temp4.fls',status='unknown')
c      write(99,6071)fln1
c      close (unit=99,status='keep')

c      if(ievtflg.eq.0)then
cc        write(6,5002)
cc         read(*,*)conv
c         write(6,5003)
c         read(*,*)smc5(1)
c      endif
c 5003 format(//' enter the initial soil moisture (0.0-0.33):'//
c     *' enter -1 if you have antecedent precip. data at precip. gauges'/
c     *' or enter average watershed value between .0 and .33  '/)

c      if(smc5(1).gt.0.33)then
c         write(6,5004)
c      end if
c 5004 format(' ','value entered is too large - 0.33 assumed')

      smc5(1)=0.25
      do ii=2,5
         smc5(ii)=smc5(1)
      end do

 110  CONTINUE

! WRITE A NEW EVENT.evt FILE:

	if(firstpass)then
	  print*,'                     if names correspond'
	  print*,'WARNING:'
	  print*,'Existing evt files may be overwritten'
	  print*,'                     if names correspond'
	  print*,'Hit Ctrl C to abort'
	  pause 'Hit enter to continue'
	  Print*
	  print*,'recent change(s):'
	  print*,'8/12/2011 - `ensimflg` has been replaced by `kenueflg`'
	  print*,'in the new evt file'
	  print*,'Pre-existing programs will not accept `kenueflg` &'
	  print*,'will produce an error'
	  print*,'Please update all executables'
	  print*,'New exec`s will recognize ensimflg & kenueflg'
	  print*,'so there is no need to edit your data files'
	  print*
	  pause 'Hit enter to continue'
	  firstpass=.false.
	endif


      scale=1.0 
      fln(99)=fln1
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call write_event(date,conv,scale,smc5,mo,nhr,nhf
      
     *       ,yy1,mm1,dd1,hh1,evtname,evtlength,evtType)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      print*,fln(99)(1:60),'... created'


	write(90,*)fln1(1:72)

	ievtflg=1

	end do   ! id=1,noevents
!    ***************************************************************
!    ***************************************************************
!    ***************************************************************
!    ***************************************************************
!    ***************************************************************

	close(unit=90,status='keep')
	Print*,'A file called events\events_to_follow.txt'
	print*,'has been written'
	print*,'This list can be used in the first event file'
	print*,'of a connected set of events'
	print*



	print*,'********************************************************'
	print*,'*                                                      *'
	print*,'*                  WATFLOOD (R)                        *'
	print*,'*                                                      *'
	print*,'*           Program mk_cfg    Jan. 17 2020             *'
	print*,'*                                                      *'
	print*,'*                Version 10_2_15                       *'
	print*,'*                                                      *'
	print*,'*           (c) N. Kouwen, 1972-2020                   *'
	print*,'*                                                      *'
	print*,'********************************************************'
      print*
      print*,'Please see file evt_info.txt for information re: this run'



	stop 'Normal ending'


      END PROGRAM make_evt


