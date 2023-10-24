        PROGRAM ragmet

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
!  VERSION 3      Dec. 2023  - NK: Added netCDF .nc output
!
!  precip lapse rate - added 2011 NK header

!  same as in DHSVM model
      
!  revised for multiple ID's  nk Dec. 06/09
!  revised for write_r2c arguments Dec. 08/14  NK
!  timestamp algorithm revised for more than 1 year March 19/16 NK
!  modified for leap years  Jan. 05/17 NK          
     
!  Accumulate the climate mean annual precip amouts for each day of the year
!  Addded jan. 31/2020  NK      


!
!        variable list:
!     NG = NUMBER OF RAIN GAGE STATIONS
!     nhg = NUMBER OF HOURS OF RECORDED DATA
!
!   I - R( , )  REAL*4     radar precipitation data
!   I - nhg     INT        number of hours of data
!   I - ixr     INT        number of east-west grid squares
!   I - iyr     INT        number of north-south grid squares
!   I - SMCRAG  REAL*4     initial soil moisture content
!   I - MO      INT        month of year
!   I - CONV    REAL*4     conversion factor
!   I - DATE    CHAR*12    event identifier
!   I - fln     CHAR*30    file names

!***********************************************************************


      use area_watflood
      use ef_module
!      use area_netcdf
      use area_debug
	implicit none

!        include 'debug.for'
        
      INTEGER         :: iAllocateStatus,iDeallocateStatus
      INTEGER         :: ng,nhg,pmax,iyr,ixr,month_no,days_from_start
      INTEGER         :: no_days_year(100),first_day_this_year(100)
      integer         :: cumm_days_year_end(100),ju_start
      integer         :: month_last,day_last,hour_last
!     added aug 18/11 nk      
      integer         :: ioflag,ios,nh,nhf,nhr,j,l,n,old_nh,is,k
      integer         :: ju,max,iallocate,nhdt,iflag,lflag,kflag
      INTEGER         :: un,fn      
      real*4          :: conv1,scale,smc5,time,value
!      
      REAL            :: conv,max_precip,temp_value
      real*4          :: coef1,coef2,coef3,coef4
      CHARACTER(14)   :: date
      character(6)    :: datatype
      DIMENSION       :: smc5(16)
      INTEGER(kind=2) :: nrad,I
      INTEGER(2)      :: status1
      integer         :: no_days,old_nhg,old_ng,attcount
      integer         :: frame_no1,frame_no2,frame_no3
      CHARACTER(3)    :: buf
      CHARACTER(1)    :: stopflg,dataflg,local_smear_flg,answer
      character(10)   ::  fileformat
      character(30)   :: xyz_filename
      logical         :: lpsflg0,lpsflg1,generic
      integer, dimension(:,:), allocatable :: no_hrs_precip
      real,    dimension(:),   allocatable :: lapse_rate_1d
      real,    dimension(:,:), allocatable :: lapse_rate_2d,
     *                                        unused,ptotal
      
      logical, dimension(:),   allocatable :: dataflag
      logical narrflg,exists,elv_max_flg,weirdflg
      real*4, dimension(:,:), allocatable :: elev_grid
      real*4, dimension(:,:), allocatable :: Pclimate  ! climate P sum


      INTEGER        :: ju_mon(24),ju_mon_ly(24)        ! mohours(24),

!     for netCDF addition
      integer         :: edays,esecs,eyears,julian_day,frame_nc
      real,    dimension(:,:,:), allocatable :: outarray6        ! 6 hr time step for FEWS
      
      
c      DATA mohours/744,672,744,720,744,720,744,744,720,744,720,744,
c     *             744,672,744,720,744,720,744,744,720,744,720,744/

      DATA ju_mon/1, 32, 60, 91,121,152,182,213,244,274,305,335,
     *          366,397,425,456,486,517,547,578,609,639,670,700/

!     modified for leap years  Jan. 05/17 NK          
      DATA ju_mon_ly/1, 32, 61, 92,122,153,183,214,245,275,306,336,
     *          367,398,426,457,487,518,548,579,610,640,671,701/

      id=1   ! needed in rdevt

      iall=0
c	debugflg=.true.  ! iopt from the par file will not be changed
c      iopt=1	  
	  
! GET THE DISTANCE WEIGHTING PARAMETER NW FROM THE PROGRAM ARGUMENT
! USED IN WEIGHT

!     USE DFLIB

      buf='   '
      CALL GETARG(1, buf, status1)
      print*,status1
      print*, buf
      
!     buf = nul : write frames with zero precip - for MESH      

!     buf = DDP - distribute daily precip
!     buf = dhp - distribute hourly precip
! unit=201  fln(201)- Point Hourly Precip           yyyymmdd_pcp.tb0
! unit=202  fln(202)- Point Daily Precip            yyyymmdd_pcp.tb0
! unit=203  fln(203)- Gridded Hourly Precip         yyyymmdd_pcp.r2c
! unit=204  fln(204)- Gridded Daily precip          yyyymmdd_pcp.r2c

      rdr='Gauge'
!      allocate(hdrcomment(100),hdrcomment_value(100),stat=iallocate)



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

c        data ioflag/2/itogo/0/

      ioflag=2

! PRINT EVENT TITLE

      write(6,1000)

      open(unit=50,file='debug\ragmet_cumm.txt',recl=flen,
     *       status='unknown',iostat=ios)
      if(ios.ne.0)then
        print*,'Problem opening debug\ragmet_cumm.txt'
        print*,'in \debug directory'
        print*,'Please create bsnm\debug and try again'
        print*
        stop 'Program aborted in ragmet @ 110'
      endif       
      open(unit=51,file='debug\ragmet_info.txt',recl=flen,
     *       status='unknown',iostat=ios)
      if(ios.ne.0)then
        print*,'Problem opening debug\ragmet_info.txt'
        print*,'in \debug directory'
        print*
        stop 'Program aborted in ragmet @ 165'
      endif       
      open(unit=52,file='debug\ragmet_recl.txt',recl=flen,
     *       status='unknown',iostat=ios)
      if(ios.ne.0)then
        print*,'Problem opening debug\ragmet_recl.txt'
        print*,'in \debug directory'
        print*
        stop 'Program aborted in ragmet @ 165'
      endif       
      open(unit=53,file='debug\ragmet_flag.txt',recl=flen,
     *       status='unknown',iostat=ios)
      if(ios.ne.0)then
        print*,'Problem opening debug\ragmet_flag.txt'
        print*,'in \debug directory'
        print*
        stop 'Program aborted in ragmet @ 165'
      endif       
      open(unit=54,file='debug\ragmet_yrly.txt',recl=flen,
     *       status='unknown',iostat=ios)
      if(ios.ne.0)then
        print*,'Problem opening debug\ragmet_yrly.txt'
        print*,'in \debug directory'
        print*
        stop 'Program aborted in ragmet @ 165'
      endif       



! INPUT EVENT PARTICULARS

! TS - ALLOCATION STATEMENT TO DYNAMICALLY ALLOCATE ARRAY'S FOR AREA12

      allocate(fln(999),outfln(999),stat=iAllocateStatus)
      if (iAllocateStatus .ne. 0) STOP 
     *    '**Allocation failed for area12**' 

      
	ni=1

!     Read the first event file and get names & number of others
      fln(99)='event/event.evt'
c     before read_evt      
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call read_evt(date,conv1,scale,smc5,nhg,nhf)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     after read_evt
!     modefied Aug. 1/11 nk
      
	INQUIRE(FILE='basin\ragmet.par',EXIST=exists)
	if(exists)then
	  print*,'found the file `basin\ragmet.par`'
	  print*,'use of this file has been replaced by reading the'
	  print*,'radius of influence and the smoothing distance'
	  print*,'from the bsnm_par.csv file'
	  Print*,'Please delete the `basin\ragmet.par` file'
	  print*,'and ensure the proper values are in the file:'
	  print*,fln(2)(1:60)
	  print*,'Sorry for any inconvenience'
	  print*
	  stop 'ragmet aborted @159'
	endif


      print*,'********************************************************'
      print*,'*                                                      *'
      print*,'* RRRRRR         A       GGGG  MM     MM EEEEEE TTTTTT *'
      print*,'* RRRRRRR       AAA     GGGGGG MMM   MMM EEEEEE TTTTTT *'
      print*,'* RR   RR      AA AA    GG     MMMMMMMMM EE       TT   *'
      print*,'* RRRRRRR     AA   AA   GG     MM MMM MM EEEEE    TT   *'
      print*,'* RRRRRR     AAAAAAAAA  GG GGG MM  M  MM EEP      TT   *'
      print*,'* RR   RR   AA       AA GG  GG MM     MM EEEEEE   TT   *'
      print*,'* RR    RR AA         AA GGGG  MM     MM EEEEEE   TT   *'
      print*,'*                                                      *'
      print*,'*                  WATFLOOD (TM)                       *'
      print*,'*             Version 11   Nov. 20.22                  *'
      print*,'*                                                      *'
      print*,'*  >>>>>>>> read tb0, narr  or rag files <<<<<<<<      *'
      print*,'*                                                      *'
      print*,'*             (c) N. Kouwen, 1972-2022                 *'
      print*,'*                                                      *'
      print*,'********************************************************'
      print*      

c      write(600,61001)'nh','#days','ju','year',
c     *             'month','day','hour'
c61001 format(20a8)

      do id=1,ni     
        print*
        print*,'***********************RAGMET**************************'
        print*,'event #',id,' of',ni
        weirdflg=.false.
 
        if(id.gt.1)then
!         read the next event file
          fln(99)=event_fln(id)
c         before read_evt (2)          
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call read_evt(date,conv,scale,smc5,nhr,nhf)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c         after read_evt (2)
        endif
        mo=mo1   ! mo1 = in area 2 - used to be in arg list

!     rev. 10.4.29 Nov.  23/20  = NK netCDF - nc output for r2c & tb0 inputs
!     Reinstituted Nov. 19/22
        if(netCDFoutflg.eq.'y')then
!         Check to see that the output file name = 'nc'
          If(fltype(10)(1:2).eq.'nc')then
! calculate epoch time at start      https://www.epochconverter.com/   
            Eyears=year1-1970  ! To jan 1        see timer in CHARM
            Edays=Eyears*365+Eyears/4
            jul_day_now=julian_day(year1,month1,day1)
            epoch_min_start=float((Edays+jul_day_now-1)*24*60)
            Esecs=epoch_min_start*60.0
            print*,'NEW:'
            print*,'Epoch time in min since 1970/1/1 00:00 for start:'
            print*,'       year         month        day        hour'
            print*,year1,month1,day1,hour1-1
            print*,'no years ',Eyears
            print*,'Epoch time in minutes ',epoch_min_start
            print*,'Epoch time in seconds ',Esecs
          write(63,*)''
          write(63,*)'Epoch time in min since 1970/1/1 00:00 for start:'
          write(63,*)'       year         month        day        hour'
          write(63,*)year1,month1,day1,hour1-1
          write(63,*)'Julian day',jul_day_now
          write(63,*)'Eyears',Eyears
          write(63,*)'Edays ',Edays
          write(63,*)'Epoch time in minutes ',epoch_min_start
          write(63,*)'Epoch time in seconds ',Esecs
        else
          print*,'ERROR:'
          print*,'Expecting a netCDF .nc format for the output'
          print*,'as in the event\event.evt file'
          print*,'and al additional  event\*.evt files'
          print*,':griddedrainfile              radcl\19800101_met.nc'
          stop 'Ragmet aborted @ 337'
          endif
        endif

        frame_no1=0
        frame_no2=0
        frame_no3=0
        frame_nc=0
        if(id.eq.1)then
!         read the shd file
!         read the shd file
          call find_filetype(1)
          if(IsFileTypeR2C(fln(1))) then
            print*
            print*,'Gone to read ',fln(1)(1:40)
            print*,' with read_shed_ef'
            call read_shed_ef(31,1)   
            print*,'Finished reading',fln(1)(1:40)
            print*
          else
            print*
            print*,'Gone to read ',fln(1)(1:40)
            print*,' with read_shed_hype'
!            call read_shed_hype(31,1)
!            call rdshed(31,1)
            print*,'Finished reading',fln(1)(1:40)
            print*
          endif
          allocate(nsdc(classcount),snocap(classcount),
     *          idump(classcount),stat=iAllocate)
          if(iAllocate.ne.0) STOP
     *     'Error with allocation of areamelta arrays in ragmet @ 275'
          
          n=xcount*ycount
          allocate(Pclimate(n,366),stat=iAllocate)
          if(iAllocate.ne.0) STOP
     *     'Error with allocation of Pclimate array in ragmet @ 280'
          do k=1,366
              do n=1,xcount*ycount
                      Pclimate(n,k)=0.0
              end do
          end do
          
          rlapse=0.0
	    radinfl=1000000.0  ! default - all stations will be used
          smoothdist=0.0         ! default - no smoothing
          call find_filetype(2)
          if(filetype.eq.'csv')then
            print*,'reading  ',fln(2)(1:50)
!           **********************************************************************
	      call read_par_parser(32,2)
!           **********************************************************************
	    endif

	    print*,'ni=',ni
          if(dds_flag.ne.0)iopt=0

!		convert rad of influence from km to # cells
          if(iopt.ge.1)then
            print*,'Echo input parameters:'
            print*,'precip lapse rate /m =',rlapse
	      print*,'radinfl=',radinfl,' km'
	      print*,'smoothdist =',smoothdist ,' km'
	      print*
          endif
	    radinfl=radinfl*1000./al
	    smoothdist=smoothdist*1000./al
          if(iopt.ge.1)then
	      print*,'Converted to cell numbers used in program:'
	      print*,'radinfl=',radinfl,' cells'
	      print*,'smoothdist =',smoothdist ,' cells'
          endif

          allocate(elev_grid(ycount,xcount),stat=iAllocateStatus)
          if (iAllocateStatus .ne. 0) STOP 
     *        '**Allocation failed for elev_grid in tmp @ 163**'
c          allocate(array3D(ycount,xcount,1),stat=iAllocateStatus)
c          if (iAllocateStatus .ne. 0) STOP 
c     *        '**Allocation failed for area3D in ragmet @ 319**'

!         convert to gridded elevations    
c         print*,xcount,ycount
          do i=1,ycount
            do j=1,xcount
              elev_grid(i,j)=0.0
            end do
          end do

          do i=1,ycount
            do j=1,xcount
              n=s(i,j)
              if(n.gt.0)then
                elev_grid(i,j)=elev(s(yyy(n),xxx(n)))
c               write(*,*)n,elev(s(yyy(n),xxx(n)))
              endif
            end do
          end do

          do i=1,ycount
              write(51,87654)(elev_grid(i,j),j=1,xcount)
          end do
87654     format(999f6.0)

!         For the met file, use the domain specified by the pdl file
!         Reason is that we may want a larger domain for grid shifting practice.
!         Note that snw & mosit use the shed file as these fields should 
!         never need to be shifted.
!         change to read pdl once only - Feb. 11/11
!        
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
c         before call_inpgrd
          call inpgrd(ixr,iyr,2)  !2 is a flag to read the grid specs only
c         after call_inpgrd          
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

!         check that the grid in the pdl file matches the grid in the shd file
          if(ixr.ne.xcount.or.iyr.ne.ycount)then
            print*
            print*,'WARNING:'
            print*,'Grid size in the pdl file does not match the'
            print*,'grid in the shd file'
            print*,'xcount shd =',xcount,' in the pdl file =',ixr
            print*,'ycount shd =',ycount,' in the pdl file =',iyr
            print*,'OK if you want the precip grid larger'
            print*,'than the watershed grid'
            print*
          endif
             
      
!         look for mean elevation file & use if present      
          elv_max_flg=.false.
	    INQUIRE(FILE='basin\elv_max.r2c',EXIST=exists)      
          if(exists)then
            print*,'WARNING:'
            print*,'found the file elv_max.r2c '
            print*,'Elev in the shd file will be replaced by those'
            print*,'in the elv_max.r2c file'
            print*
!           SUBROUTINE read_r2c(unitNum,flnNum,hdrflg)
            fln(499)='basin\elv_max.r2c'
c           before read_r2c to read elv_max.r2c
c            call read_static_r2c(499,499,'1')
c           after read_r2c            
            print*,'finished reading the header'
c           before read_r2c to read elv_mean.r2c
c            call read_static_r2c(499,499,'1')
            call read_static_r2c(499,499,attcount)
c           after read_r2c            
!           replace the elevation in the shd file with the mean elevations            
            do i=1,ycount
              do j=1,xcount
                elev_grid(i,j)=inarray3(i,j,1)  !only 1 attribute here
              end do
            end do
            elv_max_flg=.true.
          endif
        endif   ! ID=1
d       print*,'Passed L389'

	  if(smearfactor.eq.0.0.and.smrflg.eq.'y')then
	    smearfactor=a12
	  endif
	  
!       This chech added Feb. 09/14  nk	  
        if(smearfactor.eq.0.0.and.smrflg.eq.'y')then
          print*,'ERROR:'
          print*,'disaggregation flag smrflg = `y`'
          print*,'but a12 in the par file = 0.0'
          print*,'Please set the a12 value or set smrflg = `n`'
          print*
          print*,'RAGMET aborted at 408'
        endif
        
        if(mod(year1,4).eq.0)no_days=366
        
        call find_filetype(5)
d       print*,'Passed L400'
        
        print*
        print*,'filetype=',filetype
        if(filetype.eq.'rag')then
          print*
          print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
          print*,'The input file found is ',fln(5)
          print*,'Check if'
          print*,'The first data line in the .rag file is the '
          print*,'initial soil moisture at selected locations.'
          print*
          print*,'This is no longer allowed as the initial soil'
          print*,'moisture is now entered with the .psm and .gsm '
          print*,'files.'
          print*
          print*,'Please note also that the initial soil moisture'
          print*,'is now required for each land cover class.'
          print*,'Please create the .psm file and run MOIST.exe'
          print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
          print*
!        stop 'Program aborted in RAGMET @ 282'
        endif
d       print*,'Passed L474'

! INPUT RAIN GAUGE FILE DATA

        ng=ioflag

! IOFLAG IS HARDWIRED IN RAGMET DATA STATEMENT
! FOR IOFLAG=2 GRID IS READ FROM basin\nsnm.rag FILE
! STATION LOCATIONS ARE READ FROM THE raing\yymmdd.rag FILE

! FOR IOFLAG=1 GRID IS READ FROM basin\bsnm.rag FILE
! STATION LOCATIONS ARE ALSO READ FROM THE basin\bsnm.rag FILE

! TS - ALLOCATIONS FOR SMCRAG,XSTA,YSTA,SNG,SNGMIN,EWG,EWGMIN,GNAME
!      OCCUR WITHIN SUBROUTINE INPGRDA.

    
! TS -   ALLOCATION FOR RAIN() OCCURS WITHIN SUBROUTINE INPRAGA.
        
!     rev. 10.2.11 Dec.  18/17  - NK: 4 files added for BLEND.exe    
        if(buf.eq.'dhp')then         ! distribute hourly precip
            un=293
            fn=63
            smrflg='n'
        elseif(buf.eq.'ddp')then     ! distrubute daily precip
            un=294
            fn=64
            smrflg='n'
        else                         ! default - daily or hourly
          un=35
          fn=5
        endif

        print*,'Using input file   ',fln(un)(1:40)
        print*,'Creating           ',fln(fn)(1:40)

        if(filetype.eq.'tb0')then
        
c        pause 'to read_rag'
          
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
          call read_rag_ef(un,fn,nhg,conv,ng)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
c         print*,Back from read_rag_ef'

!         reduce precip to sea level              
          lpsflg0=.false.                ! default rlapse .le. 0.0 Mar. 18/15 NK
          lpsflg1=.false.                ! default rlapse .le. 0.0 Mar. 18/15 NK
!         conditional added Feb. 29/12 nk
!         check to see if lapse rate given:
          if(rlapse.gt.0.0)then
c            lpsflg0=.false.
c            lpsflg1=.false.
            do n=1,ng
              if(sta_elv(n).le.0.0)lpsflg0=.true.  !missing elevation
              if(sta_elv(n).gt.0.0)lpsflg1=.true.  !elevation present
            end do
            if(lpsflg0)then
              print*,'WARNING:'   
              print*,'One or more station elevations <= 0.0 were found'
              print*,'while the rlapse .ne. 0.0'
c              print*,'Set temp. lapse rates to 0.0 and continue?    y/n'
              print*,'OK if no data for the station(s)'
              print*,'--------------------------------'
            endif
          endif
        
c        endif
        
        elseif(fln(5)(1:18).eq.'raing\narr\narr380')then
        
          print*,'gone to read_narr'
          datatype='precip'
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
c          call read_narr(35,5,nhg,conv,ng,datatype)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
          narrflg=.true.
c         print*,'narrflg=',narrflg
        
        elseif(filetype.eq.'rag')then
!         still read old formats (non ensim)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
!          call rdrag(nhg,conv,ng,fileformat)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
        else
	    print*
	    print*,'wrong file format found'
	    print*,'file extension =',filetype
	    print*
          stop 'Ragmet aborted @ 429'
        endif
        
! TS -   ADDED ALLOCATION STATEMENT TO DYNAMICALLY ALLOCATE ARRAYS IN 
!        AREA16
!          Variables used in array allocations are defined in the 
!          following file(s):
!                      ixr,iyr  - inpgrda
!                      ng       - inpraga
!                      nhg      - inpevta, inpraga
        
!     revised for multiple ID's  nk Dec. 06/09
	  if(allocated(temp))then
	    if(ng.gt.old_ng.or.nhg.gt.old_nhg)then
	      deallocate(temp,w,ntogo,distance,dataflag)
            allocate(temp(ng),w(iyr,ixr,ng),ntogo(nhg),
     *           distance(iyr,ixr,ng),dataflag(ng),
     *          stat=iAllocateStatus)
	      old_nhg=nhg
	      old_nh=ng
	      print*,'reallocations for'
	      print*,ng,' stations'
	      print*,nhg,' no of hours of data'
	    endif
	  else
          allocate(temp(ng),p(iyr,ixr),psum(iyr,ixr),
     *           ptotal(iyr,ixr),                    
     *           pmin(iyr,ixr),w(iyr,ixr,ng),
     *           ntogo(nhg),distance(iyr,ixr,ng),
     *           no_hrs_precip(iyr,ixr),dataflag(ng),
     *           lapse_rate_1d(ng),lapse_rate_2d(iyr,ixr),
     *           unused(iyr,ixr),
     *          stat=iAllocateStatus)
          if(iAllocateStatus .ne. 0)then 
	      print*,'Allocation failed for temp,p,psum etc.' 
	      print*,'arrays in ragmet '
		  STOP  'Program aborted in ragmet @ 443'
          endif
          p=0.0000000  ! needed to get rid of noise in the netCDF file
	    print*,'allocations for'
	    print*,ng,' stations'
	    print*,nhg,' no of hours of data'
	    print*,iyr,' rows'
	    print*,ixr,' columns'
	    old_ng=ng
	    old_nhg=nhg
!         ` is in areawfo
!         this has to be allocated before calling write_r2c
          allocate(outarray(iyr,ixr),outarray6(iyr,ixr,4),
     *                       stat=iAllocate)
            if(iAllocate.ne.0) STOP
     *      'Error with allocation of outarray in ragmet'      

 	  endif
d       print*,'Finished allocations'
       
!       remove spurious values  -  added Dec. 19/13  nk
        if(deltat.le.1)then
          do n=1,ng
            do j=deltat,nhg,deltat
              if(rrain(n,j).gt.100.0)then
                rrain(n,j)=-1.0
                weirdflg=.true.
              endif
            end do
          end do
        endif
            
!       cALCULATE SEA LEVEL PRECIP
c        if(sta_elv(1).gt.0.0)then
c          if(rlapse.gt.0.0)then
          if(lpsflg1)then
!           only idf elevations are present in the tb0 file          
            do n=1,ng
              do j=deltat,nhg,deltat
                if(rrain(n,j).gt.-1.0)then
c                 cannot use equation same as for temp because it has to
c                 be a multiplier as the total precip changes drastically
c                 over the extent of a watershed. Not so for temp.
                  rrain(n,j)=rrain(n,j)/(1+sta_elv(n)*rlapse)
                endif
              end do
            end do
          endif
c        endif

! TS -   MOVED DO LOOP FROM BEFORE LINE53 TO HERE BECAUSE PRECEEDED
!        ALLOCATION.
               
!       Find out which stations have data & set flag
        do l=1,ng
          dataflag(l)=.false.
          do j=1,nhg
            if(rrain(l,j).ge.0.0)dataflag(l)=.true.
          end do
        end do
                
        if(narrflg)then
          open(unit=99,file='raing\stationsprecip_sta.xyz',
     *                status='unknown',iostat=ios)
          if(ios.eq.0)then
            do l=1,ng
              if(dataflag(l))then
                write(99,*)xsta(l),ysta(l),l,gname(l)
              endif
            end do
            close(unit=99,status='keep')
          endif
        else
          open(unit=99,file='precip_sta.xyz',status='unknown')
          do l=1,ng
            if(dataflag(l))then
              write(99,*)xsta(l),ysta(l),l,gname(l)
            endif
          end do
          close(unit=99,status='keep')
        endif
        
!       moved from read_rag_ef.f so we can print the above Feb. 25/09  NK
        
!       This does the lenght calculations in terms of no of grids.
!       That way it works for utm & latlong - although for lat long the 
!       n-s and e-w distances are not properly accounted for if the grids
!       are not set up as nearly square. But that's the user's problem.
        do n=1,ng
          ysta(n)=(ysta(n)-yorigin)/ydelta        
          xsta(n)=(xsta(n)-xorigin)/xdelta
        end do
        
!        SUBROUTINE write_r2c(un,fn,
!       *            no_frames,no_classes,frame_no,class_no,
!       *            no_signf)
        
        author='ragmet.exe                              '     
        name='Precipitation                           '
        coordsys_temp=coordsys1
        zone_temp=zone1
        datum_temp=datum1
        xorigin_temp=xorigin
        yorigin_temp=yorigin
        xcount_temp=xcount
        ycount_temp=ycount
        xdelta_temp=xdelta
        ydelta_temp=ydelta
        attribute_name='precipitation                           '
        attribute_units='mm                                      ' 
        attribute_type='Gauge                                   '  
        source_file_name=fln(5)
        unit_conversion=conv       

        no_hdrcomments=4
        allocate(hdrcomment(no_hdrcomments),
     *            hdrcomment_value(no_hdrcomments),stat=iAllocateStatus)
        if(iAllocateStatus.ne.0)STOP 
     *      '**Allocation failed for hdrcomments in ragmet @ 777**' 
        hdrcomment(1)='#precip lapse rate   1/m        '
        if(lpsflg1)then
          hdrcomment_value(1)=rlapse
        else
          hdrcomment_value(1)=0.0
        endif
!       added Aug. 1/11 nk:        
        hdrcomment(2)='#radius of influence  km        '
        hdrcomment_value(2)=radinfl*al/1000.
        hdrcomment(3)='#smoothing distance   km        '
        hdrcomment_value(3)=smoothdist*al/1000.
!       added Jan. 20/12 nk        
        if(elv_max_flg)then
          no_hdrcomments=4
          hdrcomment(4)='#elv_max_used                 '
          hdrcomment_value(4)=0.0
        else 
          no_hdrcomments=3
        endif
        
        !     rev. 10.2.11 Dec.  18/17  - NK: 4 files added for BLEND.exe    
        if(buf.eq.'dhp')then         ! distribute hourly precip
            un=295
            fn=65
        elseif(buf.eq.'ddp')then     ! distrubute daily precip
            un=298
            fn=68
        else                         ! default - daily or hourly
          un=40
          fn=10                      ! the gridded precip file for CHARM
        endif

        
!       Write the header (frame_no1 = 0)


        if(FLtype(10).eq.'r2c')then      
d           print*,'writing the header in',fln(un)(1:40)
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call write_r2c(un,fn,nhg,1,0,0,1)   
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        endif                             
                   
        deallocate(hdrcomment,hdrcomment_value)

! Take  n out of data statement in area16 for F90
        do i=1,iyr
          do j=1,ixr
            psum(i,j)=0.0
            pmin(i,j)=999.0
          end do
        end do

        if(id.eq.1)then
          do i=1,iyr
            do j=1,ixr
              ptotal(i,j)=0.0
            end do
          end do
	  endif
        
!       FIND HOW MANY HOURS OF RAINFALL HAVE BEEN RECORDED:
        
        if(conv1.ne.conv)then
          write(6,1111)conv1,conv
 1111     format(' ','WARNING: event and rain file values of conv', 
     *        '     do not match. From event',f6.2,' from raing',f6.2, 
     *        '     The value from the .RAG file is used')
          print*
          pause ' In ragmet @ line ~141. Hit enter to continue'
        endif
        
        write(*,1112)nhg,ng
 1112   format(' file length =',i5,'hours    # gauges =',i5)
        
        
!       MET FILES ARE PRODUCED FOR LENGTH OF RAINFALL RECORD
!       CONVERT RAINGAGE VALUES TO MM
!       ONLY RAIN GAGE FILE CAN BE IN NON-MM AMOUNTS
        
        if(conv.ne.1.0)then
          do nh=1,nhg
            do is=1,ng
              rrain(is,nh)=rrain(is,nh)*conv
            end do
          end do
          conv=1.0
        endif
        
!       DISTRIBUTE RAIN ACCORDING TO WEIGHTS AND GENERATE yyyymmdd_met.r2c FILE 
        time=0.0
        
        write(6,1110)
        year_now=year1
        
!       find # days this year 
        if(mod(year1,4).eq.0)then
          no_days_year(1)=366
        else
          no_days_year(1)=365
        endif 
        
      print*,'days in year 1 =',no_days_year(1)
      
!       find # days in first year    
!         modified for leap years  Jan. 05/17 NK          
        if(mod(year1,4).eq.0)then   !leap year    
          no_days_year(1)=no_days_year(1)-ju_mon_ly(mo1)-day1+2
          ju_start=ju_mon_ly(mo1)+day1-1
        else
          no_days_year(1)=no_days_year(1)-ju_mon(mo1)-day1+2
          ju_start=ju_mon(mo1)+day1-1
        endif
        cumm_days_year_end(1)=no_days_year(1)
c      print*,'no days year 1 =   ',no_days_year(1),cumm_days_year_end(1) 
c      print*,'julian day start =',ju_start
             
!       find # days in each year 
        do i=2,21  ! make this 100
          year_now=year1+i-1
          if(mod(year_now,4).eq.0)then
            no_days_year(i)=366
          else
            no_days_year(i)=365
          endif 
          cumm_days_year_end(i)=cumm_days_year_end(i-1)+no_days_year(i)
          if(iopt.ge.1)print*,'no days ',
     *            year_now,no_days_year(i),cumm_days_year_end(i)
        end do
        
!       Check to see there's precip in the first time increment
!       Need grid in first deltat to ensure met file does not start late.        
                       
        
        
!       time loop
!       time loop
!       time loop
!       time loop
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       timestamp algorithm revised for more than 1 year March 19/16 NK

        month_now=mo1
        year_now=year1
        
        do nh=deltat,nhg,deltat
          days_from_start=(nh-deltat)/24+1
          do i=1,20  
            if(days_from_start.gt.cumm_days_year_end(i))then
              year_now=year1+i
              j=i
            endif
          end do
          if(days_from_start.le.cumm_days_year_end(1))then
            ju=ju_start+days_from_start-1
          else  
            ju=days_from_start-cumm_days_year_end(j)
          endif
          i=12
!         modified for leap years  Jan. 05/17 NK          
          if(mod(year_now,4).eq.0)then
            do while(i.ge.1)
              if(ju.lt.ju_mon_ly(i+1))then 
                mo=i
              endif
              i=i-1  
            end do  
          else
            do while(i.ge.1)
              if(ju.lt.ju_mon(i+1))then 
                mo=i
              endif
              i=i-1  
            end do
          endif
          if(mo.gt.12)mo=mo-12
          month_now=mo        
!         modified for leap years  Jan. 05/17 NK          
          if(mod(year_now,4).eq.0)then
            day_now=ju-ju_mon_ly(mo)+1
          else
            day_now=ju-ju_mon(mo)+1
          endif
          hour_now=mod(nh-deltat,24)+hour1     !+deltat
          if(hour_now.ge.24)hour_now=hour_now-24
          if(hour_now.le.deltat)then
            year_last=year_now
            month_last=month_now
            day_last=day_now
            hour_last=hour_now
          endif
c      write(600,61000)nh,days_from_start,ju,year_now,month_now,
c     *             day_now,hour_now,day1
c61000 format(20i8)

          if(iopt.ge.1)then
            if(mod(nh,nhg/10).eq.0.or.mod(nh,nhg/12).eq.0)
     *        write(*,60000)'distributing hour #',nh
     *        ,ju,year_now,month_now,day_now,hour_now
60000       format(a19,6i8)
          endif
        
          do i=1,iyr
            do j=1,ixr
              p(i,j)=0.0
            end do
          end do
          iflag=0
!         iflag = 1 if weights have to be recalculated
!         this happens when there is missing data at one or more 
!         of the stations and stations have to be ignored 
!         temporarily
          lflag=0
!         lflag = 1 when there is data at at least one gage
!         this can be zero rainfall 
          kflag=-1
!         kflag = -1 when there is non-zero rainfall data
!         when kflag = -1, the rainfall grid is not printed and 
!         nh is set -nh & SPL will assume 0 rainfall for the hour
          
c          if(nh.eq.1)then
          if(nh.eq.deltat)then
          
!           CHECK IF THERE IS ANY RECORDED DATA IN 1st HOUR

            do is=1,ng
              if(rrain(is,nh).ge.0.0) lflag=1
              if(rrain(is,nh).gt.0.0) kflag=1
            end do
        
          else
            do is=1,ng
        
!             RECALCULATE WEIGHTS IF AVAILABILITY OF GAUGES HAS
!             CHANGED

              if(nh.gt.deltat)then
                if(rrain(is,nh).ge.0.0.and.  
     *           rrain(is,nh-deltat).lt.0.0)iflag=1
        
                if(rrain(is,nh).lt.0.0.and.  
     *               rrain(is,nh-deltat).ge.0.0)iflag=1
              endif

!             CHECK IF THERE IS ANY RECORDED DATA THIS HOUR
        
              if(rrain(is,nh).ge.0.0)lflag=1
              if(rrain(is,nh).gt.0.0)kflag=1
             end do
          endif
        
!         find the max recorded precip this delta t
          local_smear_flg=smrflg
          max_precip=0.0
          do is=1,ng
             max_precip=amax1(rrain(is,nh),max_precip)
          end do
          
!         Changed Feb. 09/14 NK
!         No need to disaggregate if time step = 1 hour or there's 
!         less than 1 mm
c          if(max_precip.lt.1.0.and.deltat.le.1)local_smear_flg='n'
          if(smrflg.eq.'y')then
            if(max_precip.le.smearfactor.or.deltat.le.1)then
              local_smear_flg='n'   ! no need to disaggregate
            else
              local_smear_flg='y'
            endif
          endif
          
!         THIS SECTION IS BYPASSED IF THERE IS NO RECORDED DATA
          if(lflag.eq.1)then
!           THIS SECTION IS BYPASSED IF THE STATIONS ARE THE SAME THIS HOUR
c            if(iflag.eq.1.or.nh.eq.1)then
            if(iflag.eq.1.or.nh.eq.deltat.or.nh.eq.hour1)then
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
              call weight(nh,ng,xcount,ycount,0.0)   ! changed Mar. 18/15  NK
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
            endif
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
            call distr(nh,xcount,ycount,ng)   ! changed Mar. 18/15  NK
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
          endif
          
          if(lpsflg1)then
c            if(sta_elv(1).gt.0.0)then
!             adjust the sea level precip for grid elevation
!             adjust the sea level precip for grid elevation
!             adjust the sea level precip for grid elevation
              do i=1,ycount
                do j=1,xcount
                    p(i,j)=p(i,j)*(1.0+elev_grid(i,j)*rlapse)
                  p(i,j)=amax1(0.0,p(i,j))  !prevent -ve precip
                end do
              end do          
          endif

!     Accumulate the climate mean annual precip amouts for each day of the year
!     Addded jan. 31/2020  NK      
          k=nh/deltat
          n=0
          do i=1,iyr
            do j=1,ixr
                n=n+1
!             for daily dt there will be 365(6) values / year              
              Pclimate(n,k)=Pclimate(n,k)+p(i,j)
!              if(i.eq.1.and.j.eq.1)write(400,*)'a'
            end do
          end do

!         DISAGGREGATION
!         Disaggregation  added Apr. 22/08  nk
!         only if the deltat > 1 and flag is 'y' and we are writing an r2c file
!         For netCDF - FEWS  
!         Dissaggredation needs to be for 6 hour intervals alwasy
!         So usual disaggredation is done for the r2c file: 
!         hourly amounts until unused = 0
!         Then the hourly disaggregated amounts are summed at 6 hour intervals
!         which keeps the overall disaggredation the same when using FEWS
!         CHARM will disaggregate the 6 hour amounts
          
          if(local_smear_flg.eq.'y')then
!           find no of hours of precip & rate in each grid
            do i=1,iyr
              do j=1,ixr
                no_hrs_precip(i,j)=min(int(p(i,j)+1),deltat)  ! assumes 1 mm/hr
                unused(i,j)=p(i,j)
              end do
            end do
            
            hh=0
            outarray6=0.0
            do nhdt=1,deltat
              if(mod(nhdt-1,6).eq.0)hh=hh+1
!              print*,nhdt,mod(nhdt,6),hh
              dataflg='n'
              hour_now=mod(nh-deltat,24)+nhdt-1
            outarray=0.0
              
              if(hour_now.gt.24-deltat)then
                year_now=year_last
                month_now=month_last
                day_now=day_last
              endif
              do i=1,ycount
                do j=1,xcount
                  if(p(i,j).le.0.0)then
                    outarray(i,j)=0.0
                    outarray6(i,j,hh)=outarray6(i,j,hh)+outarray(i,j)
                    if(i.eq.1.and.j.eq.1.and.debug_output)
     *              write(63,*)nhdt,hh,outarray6(1,1,hh),p(1,1)
                  elseif(p(i,j).le.smearfactor)then
                    if(unused(i,j).ge.0.0)then
                     outarray(i,j)=unused(i,j)
                     outarray6(i,j,hh)=outarray6(i,j,hh)+outarray(i,j)
                     unused(i,j)=0.0
                     if(outarray(i,j).gt.0.0)dataflg='y'
                    endif
                    if(i.eq.1.and.j.eq.1.and.debug_output)
     *              write(63,*)nhdt,hh,outarray6(1,1,hh),p(1,1)
                  elseif(p(i,j).le.deltat)then
                    if(unused(i,j).gt.0.000001)then
                      outarray(i,j)=amin1(smearfactor,unused(i,j))
                      outarray6(i,j,hh)=outarray6(i,j,hh)+outarray(i,j)
                      unused(i,j)=unused(i,j)-outarray(i,j)
                      dataflg='y'
                    if(i.eq.1.and.j.eq.1.and.debug_output)
     *              write(63,*)nhdt,hh,outarray6(1,1,hh),p(1,1)
                    else
                      outarray(i,j)=0.0
                      outarray6(i,j,hh)=outarray6(i,j,hh)+outarray(i,j)
                    if(i.eq.1.and.j.eq.1.and.debug_output)
     *              write(63,*)nhdt,hh,outarray6(1,1,hh),p(1,1)
                    endif
                  elseif(p(i,j).gt.deltat)then
                    outarray(i,j)=p(i,j)/deltat
                    outarray6(i,j,hh)=outarray6(i,j,hh)+outarray(i,j)
                    dataflg='y'
                    if(i.eq.1.and.j.eq.1.and.debug_output)
     *              write(63,*)nhdt,hh,outarray6(1,1,hh),p(1,1)
                  else
                    print*,'nhdt,i,j,p(i,j)/',nhdt,i,j,p(i,j)
                    print*,'This should never happen!!!'
                    stop
                  endif
                end do
              end do  
              
              do i=1,iyr
                do j=1,ixr
                  if(dataflg.eq.'n'.and.outarray(i,j).gt.0.0)dataflg='y'
c                 outarray(i,j)=amin1(outarray(i,j),999.9)  ! prevents ****.* in file
!                 get the precip each hour & accumulate for reporting in netCDF file 
                end do
              end do     

!             Make sure
              if(deltat.eq.nh)dataflg='y'
           
              if(dataflg.eq.'y')then
                frame_no1=frame_no1+1
                 if(FLtype(10).eq.'r2c')then       
!                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    call write_r2c(un,fn,nhg,1,frame_no1,0,9)
!                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                 endif
c                endif
              endif
            end do  ! nhdt=...
          else
!           OUTPUT PRECIPITATION DATA  -  no disaggredation
!           convert local variables to write module variables
            dataflg='n'
            do i=1,iyr
                do j=1,ixr
                  outarray(i,j)=amin1(p(i,j),999.9)  ! prevents ****.* in file
                  if(dataflg.eq.'n'.and.p(i,j).gt.0.0)dataflg='y'
                end do
            end do  
             
              if(deltat.eq.nh)dataflg='y'

            if(dataflg.eq.'y'.or.buf.eq.'nul')then
              frame_no1=frame_no1+1
              if(FLtype(10).eq.'r2c')then       
!                ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                 call write_r2c(un,fn,nhg,1,frame_no1,0,9)
!                ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              endif
            endif
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         
         endif
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        Calling write_2d_nc
!        nh/deltat - must be the frame #
!        nhg = # hours of data
!        ## file name number
!        'xxx' = 3 character unit designation - eg. 'mm '!         Write to the netCDF format yyyymmdd.nc file
!                Write to the netCDF format yyyymmdd.nc file
!         S/R Arguments
!         (k,nrecs,fn,TEMP_NAME,TEMP_UNITS)
!         k = sequence # 1, 2, 3, 4,.....,nrecs
!         nrecs = 0 - keep writing frames until done
!         .false. = do not close the file
!         fn = file name number
!         TEMP_NAME = Tname = variabel name e.g. pcp, tmp  - 3 chrs
!         TEMP_UNITS = Tunits units e.g. cms, mm , dC
         
!         for FEWS: all deltat = 6 hours.

          do hh=1,4
              frame_nc=frame_nc+1
              do i=1,iyr
                do j=1,ixr
                  outarray(i,j)=outarray6(i,j,hh)
                end do
              end do
              if(debug_output)
     *           write(63,*)nh,frame_nc,hh,outarray6(1,1,hh),p(1,1)
              if(netCDFoutflg.eq.'y')then                  
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!              call write_2D(frame_nc,0,.false.,10,'pcp','mm ',6) 
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            endif    
          
      
      
          end do
          if(debug_output)write(63,*)'********************************'
      
      
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             
          do i=1,iyr
            do j=1,ixr
              psum(i,j) = psum(i,j)+p(i,j)
	        ptotal(i,j)=ptotal(i,j)+p(i,j)
              pmin(i,j) = amax1(pmin(i,j),p(i,j))
!             for daily dt there will be 365(6) values / year              
            end do
          end do
c          endif

        end do   ! end of time loop

!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!       return value to par file setting - i.e. annual value

!       WRITE THE RAINFALL SUMMARY: TOTAL RAINFALL
         
c        write(6,6001)
c        do i=iyr,1,-1
c          write(6,6000)(psum(i,j),j=1,ixr)
c        end do
        
        write(51,6001)
        do i=iyr,1,-1
          write(51,6000)(psum(i,j),j=1,ixr)
        end do
        write(51,6002)
        do i=iyr,1,-1
          write(51,6000)(pmin(i,j),j=1,ixr)
        end do
        
        close(unit=40,status='keep')
        print*
        print*,'Closed Unit 40  file name ',fln(10)(1:40)
        
        close(unit=98,status='keep')
        
        if(weirdflg)then    !added Dec. 19/13  nk
          print*,'WARNING:'
          print*,'Spurious precip values removed for'
          print*,'hourly precip > 100 mm'
          print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          print*
        endif
      
!     close the .nc file              
      if(netCDFoutflg.eq.'y')then                  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!          call write_2D(frame_nc,frame_nc,.true.,10,'pcp','mm ',6) 
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      endif
        
      end do  ! event loop
      
!     Write the header (frame_no = 0)
      no_hdrcomments=0
      fln(599)='debug\ragmet_sum.r2c'
      source_file_name='multiple_events'
      
      do i=1,iyr
        do j=1,ixr
          outarray(i,j)=ptotal(i,j)
        end do
      end do
!      SUBROUTINE write_r2c(un,fn,
!     *            no_frames,no_classes,frame_no,class_no,
!     *            no_signf)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c      call write_r2c(599,599,nhg,1,0,1,1)   
c      call write_r2c(599,599,nhg,1,1,1,10)
            
      call write_r2c(599,599,1,1,0,0,1)   
      call write_r2c(599,599,1,1,1,1,10)
c      SUBROUTINE write_r2c(un,fn,
c     *            no_frames,no_classes,frame_no,class_no,
c     *            no_signf)
      
      
      print*,'Wrote: ',fln(599)(1:60)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
!     Write the climate mean annual precip amouts for each day of the year
!     Addded jan. 31/2020  NK      
!     Write the header (frame_no = 0)
      no_hdrcomments=0
      fln(699)='model\climate_pcp.r2c'
      source_file_name='multiple_events'
d           print*,'writing the header in',fln(un)(1:40)
!     Write the header 
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call write_r2c(699,699,nhg/deltat,1,0,0,1)   
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      frame_no1=0
      do nh=1,nhg/deltat
          
          
          
          year_now=2020
          
          
          
          
          ju=nh+1  ! added 1 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          i=12
!         modified for leap years  Jan. 05/17 NK          
          if(mod(year_now,4).eq.0)then
            do while(i.ge.1)
              if(ju.lt.ju_mon_ly(i+1))then 
                mo=i
              endif
              i=i-1
            end do  
          else
            do while(i.ge.1)
              if(ju.lt.ju_mon(i+1))then 
                mo=i
              endif
              i=i-1  
            end do
          endif
          if(mo.gt.12)mo=mo-12
          month_now=mo        
!         modified for leap years  Jan. 05/17 NK          
          if(mod(year_now,4).eq.0)then
            day_now=ju-ju_mon_ly(mo)+1
          else
            day_now=ju-ju_mon(mo)+1
          endif
! This needs fixing. DA will happen after the time stamp but it should really be before if 
! time stamp is at end of reporting period. Whigh is the case for CaPA  and the forecasts          
          hour_now=0   ! oly one value / day
          n=0
          do i=1,iyr
              do j=1,ixr
                  n=n+1
!                 divide by the # of years of data                  
                  outarray(i,j)=Pclimate(n,nh)/float(ni)
              end do
          end do
          frame_no1=frame_no1+1
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call write_r2c(699,699,nhg/deltat,1,frame_no1,0,9)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end do
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      print*,'Wrote: ',fln(699)(1:60),3
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      
      
      
        
! TS - ADDED DEALLOCATION STATEMENT FOR AREA12 AND AREA16 ARRAYS

! FORMATS

202   format(9999f10.3)
1001  format(a30)
1002  format(10x,i5,1x,f5.2,1x,i6,1x,f5.1)
1003  format(' HOUR=',2i5)
1000  format(' RAIN GAUGE RAINFALL DISTRIBUTION APPLICATION')
1100  format('+',' DISTRIBUTING RAINFALL FOR HOUR ',i5,' OF ',i5,
     +' HOURS OF PRECIPITATION DATA')
1110  format(' ')
1300  format(9999f5.1)
1301  format(9999f5.2)
5000  format(' converted data to sea level:')
5001  format(/' converted data to sea level:')
6000  format(9999f7.0)
6001  format(' Hour= -999                    Summary of rainfall:')
6002  format(' Hour= -999                    Maximum rainfall:')

        print*
        print*,'Please check that the program has operated properly'
        print*,'Compare values in the yyyymmdd_met.r2c file against'
        print*,'data in the yyyymmdd_rag.tb0 file'
        print*,'Event totals should be equal at gauge sites if '
	  print*,'data set for station is complete'
	  print*,'AND damping (smoothdist) is set to 0 !!!!!!!!!!!!'
	    print* 
	  print*,'See ragmet_recl.txt for record length in each event'
	  print*,'See ragmet_flag.txt for data summary'
	  print*,'   p = partial record this gauge'
	  print*,'   n = no data'
	  print*,'   y = full record'
        print*

        print*,'Note change:'
	  print*,'In the ragmet.par file, both radius of influence'
	  print*,'& smoothing distance are now in km!!! Please check your file'
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
      
      print*, 'Normal Ending'

      
      END PROGRAM ragmet




