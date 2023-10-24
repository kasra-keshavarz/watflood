      PROGRAM tmp


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
!*******************************************************************
! TMPDST - 
!
!     subroutines called:
!         inpevta
!         inpgrda
!         inptmpa
!         weighta
!         distra
!         outtema
!
! Modified by Tricia Stadnyk - September 2000
! Converted common blocks to modules and added dynamically allocated
! run-time arrays as part of Fortran 90 conversion.
!
!  VERSION 2      Dec. 2000  - TS: Added dynamic memory allocation
!
!
! - list of arguments:
!
!     ng  = number of rrain gage stations
!     nhg = number of hours of recorded data
!tlapse

!   I - R( , )  REAL*4     radar precipitation data
!   I - nhg      INT        number of hours of data
!   I - IXR     INT        number of east-west grid squares
!   I - IYR     INT        number of northccc-south grid squares
!   I - SMC     REAL*4     initial soil moisture content
!   I - MO      INT        month of year
!   I - CONV    REAL*4     conversion factor
!   I - DATE    CHAR*12    event identifier
!   I - FLN     CHAR*30    file names
!
!*******************************************************************

!     revised for multiple ID's  nk April 07/2010
!     revised for write_r2c arguments Dec. 08/14  NK
!     timestamp algorithm revised for more than 1 year March 19/16 NK
!     modified for leap years  Jan. 05/17 NK          

      use area_watflood
      use ef_module
!      use area_netcdf
      use area_debug

      implicit none

      CHARACTER(12) :: junk
      CHARACTER(14) :: date
      real*4        :: smc5(16)
      INTEGER       :: iAllocateStatus,iDeallocateStatus
      INTEGER       :: iAllocate,attcount
      integer       :: i,j,ioflag,nh,nhg,nhf,ng,ixr,iyr,is,ii,k
      integer       :: ju,iflag,lflag,day_no,no_days,n,l,ios,nhr,m
      integer       :: data_count,ju_last,mo_last,mo_count(12)
      integer       :: frame_no1,frame_no2,frame_no3
      INTEGER         :: no_days_year(100),first_day_this_year(100)
      integer         :: cumm_days_year_end(100),ju_start
      integer         :: month_last,day_last,days_from_start
       
      INTEGER(2)    :: status1
      CHARACTER(1)  :: buf,stopflg,answer
      character(6)    :: datatype
      real*4        :: conv,conv1,scale,cut
      logical narrflg,exists,fix_month,firstpass,elv_mean_flg
      logical       :: lpsflg0,lpsflg1,write292flg

      INTEGER        :: ju_mon(24),ju_mon_ly(24)     !mohours(24),

      real*4, dimension(:,:), allocatable :: elev_grid,mean_temp
      real*4, dimension(:,:), allocatable :: dly_min
      real*4, dimension(:,:), allocatable :: dly_max
      real*4, dimension(:,:,:), allocatable :: dly_min_sum
      real*4, dimension(:,:,:), allocatable :: dly_max_sum
      real*4, dimension(:,:,:), allocatable :: mean_dly_diff_3
      real*4, dimension(:,:), allocatable :: Pclimate  ! climate P sum
      
!     for netCDF addition
      integer         :: edays,esecs,eyears,julian_day,frame_nc
      real,    dimension(:,:), allocatable :: 
     *           outarray00,outarray06,outarray12,outarray18        ! 6 hr time step for FEWS

c      DATA mohours/744,672,744,720,744,720,744,744,720,744,720,744,
c     *             744,672,744,720,744,720,744,744,720,744,720,744/

      DATA ju_mon/1, 32, 60, 91,121,152,182,213,244,274,305,335,
     *          366,397,425,456,486,517,547,578,609,639,670,700/

!     modified for leap years  Jan. 05/17 NK          
      DATA ju_mon_ly/1, 32, 61, 92,122,153,183,214,245,275,306,336,
     *          367,398,426,457,487,518,548,579,610,640,671,701/

      DATA ioflag/1/
      data data_count/0/
      data ju_last/0/
      data mo_last/0/
      data fix_month/.true./
      data firstpass/.true./

!     USE DFLIB

C      CALL GETARG(1, buf, status1)
C      if(status1.ne.1)buf=' '
c      if(buf.eq.'1')then
c        stopflg='y'
c      else
c        stopflg='n'
c      endif


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

c      open(unit=52,file='tmp_recl.txt',recl=flen,status='unknown',
c     *       iostat=ios)
c      open(unit=53,file='tmp_flag.txt',recl=flen,status='unknown',
c     *       iostat=ios)
      open(unit=51,file='debug\tmp_info.txt',recl=16384,
     *      status='unknown', iostat=ios)
      open(unit=52,file='debug\tmp_recl.txt',recl=16384,
     *      status='unknown', iostat=ios)
      open(unit=53,file='debug\tmp_flag.txt',recl=16384,
     *       status='unknown',iostat=ios)

       
      cut=-80.0
      id=1

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
	  stop 'tmp aborted @127'
	endif
	
	rad_influence=1000000.0  ! default - all stations will be used
      min_distance=0.0         ! default - no smoothing
c!	if a value is present in ragmet.par, it will be used as a min
c!	if large, all stations to be used
c	if(exists)then
c	  open(unit=99,file='basin\ragmet.par',status='unknown',
c     *       iostat=ios)
c	  read(99,*)iopt
c	  read(99,*)rad_influence
c	  read(99,*)min_distance
c	  close(unit=99,status='keep')
c	endif


      do j=1,10
      do i=1,100000
        smc5(1)=0.1*(sin(float(i)/100000))**0.25
      end do
      end do

! PRINT EVENT TITLE
      write(*,1000)

      if(firstpass)then
! TS - ALLOCATION STATEMENT TO DYNAMICALLY ALLOCATE ARRAY FOR AREA12A
        allocate(fln(999),stat=iAllocateStatus)
        if (iAllocateStatus .ne. 0) STOP 
     *    '**Allocation failed for area12 in tmpdsta l57**'
        fln(63)='debug\tmp_netCDF_errors.txt'
      endif


      ni=1

! INPUT EVENT PARTICULARS
      fln(99)='event\event.evt'
c      call rdevt(date,conv1,scale,smc5,nhg,nhf)
!          read_evt(date,conv,scale,smc5,nhg,nhf)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call read_evt(date,conv1,scale,smc5,nhg,nhf)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      print*
	print*,'# events=',ni

        print*,'*******************************************************'
        print*,'*                                                     *'
        print*,'*       TTTTTTTTT   MM       MM   PPPPPPP             *'
        print*,'*       TTTTTTTTT   MMM     MMM   PPPPPPPP            *'
        print*,'*          TTT      MMMM   MMMM   PP    PP            *'
        print*,'*          TTT      MM MM MM MM   PPPPPPPP            *'
        print*,'*          TTT      MM  MMM  MM   PPPPPPP             *'
        print*,'*          TTT      MM       MM   PP                  *'
        print*,'*          TTT      MM       MM   PP                  *'
        print*,'*          TTT      MM       MM   PP                  *'
        print*,'*                                                     *'
        print*,'*                  WATFLOOD (R)                      *'
        print*,'*            Version 11   Dec. 05, 2022               *'
        print*,'*                                                     *'
        print*,'*      >>>>>>>> read tb0 or tag files <<<<<<<<        *'
        print*,'*                                                     *'
        print*,'*           (c) N. Kouwen, 1972-2022                  *'
        print*,'*                                                     *'
        print*,'*******************************************************'
      print*      

      do id=1,ni     
        print*
      print*,'************************TMP***************************'
      print*,'event #',id,' of',ni
      print*,'event #',id,' of',ni
      print*,'shed file =',fln(1)(1:50)
      print*,'par file  =',fln(2)(1:50)
      print*
 
        if(id.gt.1)then
!         read the next event file
          fln(99)=event_fln(id)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call read_evt(date,conv,scale,smc5,nhr,nhf)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        endif
     
        print*,'Using input file   ',fln(14)(1:50)
        print*,'Creating           ',fln(15)(1:50)
        print*,'Creating           ',fln(62)(1:50)
        print*
     
        if(mod(year1,4).eq.0)then
          no_days=366
          if(mo1.ge.3)then
            do i=mo1,24
              ju_mon(i)=ju_mon(i)+1
            end do
          endif
        else
          no_days=365
        endif
        
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

c       print*,'no days =',no_days
      
        if(id.eq.1)then
!        read the shd file
!        read the shd file
!        read the shd file
         call find_filetype(1)
      
         if(IsFileTypeR2C(fln(1))) then
           print*,'Gone to read ',fln(1)(1:30),' with read_shd_ef'
           call read_shed_ef(31,1)  
           print*,'Finished reading',fln(1)(1:50)
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
      
         if(firstpass)then
           allocate(elev_grid(ycount,xcount),
     *              mean_dly_diff_3(ycount,xcount,12),
     *              dly_min(ycount,xcount),
     *              dly_max(ycount,xcount),      
     *              dly_min_sum(ycount,xcount,12),
     *              dly_max_sum(ycount,xcount,12),      
     *                            stat=iAllocateStatus)
           if (iAllocateStatus .ne. 0) STOP 
     *      '**Allocation failed for elev_grid in tmp @ 163**'

           allocate(nsdc(classcount),snocap(classcount),
     *          idump(classcount),stat=iAllocate)
           if(iAllocate.ne.0) STOP
     *      'Error with allocation of areamelta arrays in ragmet @ 275'
         endif

          n=xcount*ycount
          allocate(Pclimate(n,2196),stat=iAllocate)  ! for deltat = 4 hours
          if(iAllocate.ne.0) STOP
     *     'Error with allocation of Pclimate array in ragmet @ 280'
          do k=1,2196
              do n=1,xcount*ycount
                      Pclimate(n,k)=0.0
              end do
          end do
  
         call find_filetype(2)

         tlapse=0.0
	   radinfl=1000000.0  ! default - all stations will be used
         smoothdist=0.0         ! default - no smoothing
         if(filetype.eq.'csv')then
               print*,iopt,dds_flag
           print*,'reading  ',fln(2)(1:50)
!          **********************************************************************
c	     call read_par(32,2)
	     call read_par_parser(32,2)
!          **********************************************************************
         endif

         if(dds_flag.ne.0)iopt=0

         if(tlapse.gt.0.0)then
           print*,'temperature lapse rate must be < 0.0 dC/m'
           print*,'Please change the value in the par file'
           print*
           stop 'tmp.exe aborted in tmp @ 288'
         endif

!	   convert rad of influence from km to # cells
         if(iopt.ge.1)then
           print*,'Echo input parameters:'
           print*,'Temperature lapse rate dC/m =',tlapse
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

!convert to gridded elevations    
c        print*,xcount,ycount
         do i=1,ycount
           do j=1,xcount
              elev_grid(i,j)=0.0
           end do
         end do
       
         do i=1,ycount
           do j=1,xcount
             n=s(i,j)
             if(n.gt.0)then
c              elev(s(yyy(n),xxx(n)))=elev(s(yyy(n),xxx(n)))*100.0
               elev_grid(i,j)=elev(s(yyy(n),xxx(n)))
c              write(*,*)n,elev(s(yyy(n),xxx(n)))
             endif
           end do
         end do
         write(51,*)'Grid elevation in m - north = up'
         do i=ycount,1,-1
           write(51,51001)(elev_grid(i,j),j=1,xcount)
51001      format(<xcount>f8.0)                 
         end do
        
        ng=-1     ! flag to say we're doing the temperature
       
!         change to read pdl once only - Feb. 11/11
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call inpgrd(ixr,iyr,ng)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!         check that the grid in the pdl file matches the grid in the shd file
          if(ixr.ne.xcount.or.iyr.ne.xcount)then
            print*,'WARNING:'
            print*,'Grid size in the pdl file does not match the'
            print*,'grid in the shd file'
            print*,'xcount shd =',xcount,' in the pdl file =',ixr
            print*,'ycount shd =',ycount,' in the pdl file =',iyr
            print*,'OK if you want the temperature grid larger'
            print*,'than the watershed grid'
            print*
          endif
       

!         look for mean elevation file & use if present      
          elv_mean_flg=.false.
	    INQUIRE(FILE='basin\elv_means.r2c',EXIST=exists)      
          if(exists)then
            print*,'WARNING:'
            print*,'found the file elv_mean.r2c '
            print*,'Elev in the shd file will be replaced by those'
            print*,'in the elv_means.r2c file'
            print*
!           SUBROUTINE read_r2c(unitNum,flnNum,hdrflg)
            fln(499)='basin\elv_means.r2c'
c            call read_r2c(499,499,'1')
            print*,'finished reading the header'
c            call read_r2c(499,499,'0')
            call read_static_r2c(499,499,attcount)
!           replace the elevation in the shd file with the mean elevations            
            do i=1,ycount
              do j=1,xcount
                elev_grid(i,j)=inarray3(i,j,1)
              end do
            end do
            elv_mean_flg=.true.
          endif
        endif

! INPUT GAUGE FILE DATA

! USED AS A FLAG TO GO TO THE TEMP GRID DATA
        iall=0
       
        if(firstpass)then
! TS - ALLOCATION STATEMENT TO DYNAMICALLY ALLOCATE ARRAY FOR AREA16A
! ALLOCATIONS FOR XSTA,YSTA,EWG,SNG,GNAME OCCUR IN INPGRDA
          allocate(p(iyr,ixr),stat=iAllocateStatus)
          if (iAllocateStatus .ne. 0) STOP 
     *    '**Allocation failed for area16 in tmp l72**'
        endif

!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call find_filetype(14)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        print*,'Found temperature filetype= ',filetype

! INPUT RAIN GAUGE FILE DATA

! TS - ALLOCATION FOR RAIN() OCCURS WITHIN SUBROUTINE INPRAGA.

        if(filetype.eq.'tb0')then
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
          call read_tmp_ef(44,14,nhg,conv,ng)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
          if(tlapse.lt.0.0)then   ! as given in the par file
            lpsflg0=.false.
            lpsflg1=.false.
            do n=1,ng
              if(sta_elv(n).le.0.0)lpsflg0=.true.  !missing elevation
              if(sta_elv(n).gt.0.0)lpsflg1=.true.  !elevation present
            end do
            if(lpsflg0)then
              print*,'WARNING:'   
              print*,'One or more station elevations <= 0.0 were found'
              print*,'while the tlapse .ne. 0.0'
c              print*,'Set temp. lapse rates to 0.0 and continue?    y/n'
              print*,'OK if no data for the station(s)'
              print*,'--------------------------------'
c              read*,answer
c              if(answer.eq.'y')then
c                tlapse=0.0
c              else
c                stop 'TMP aborted @ line 410'
c              endif
            endif
          endif
        
c          if(tlapse.lt.0.0)then
          if(lpsflg1)then
!           convert temperatures to sea level
!           convert temperatures to sea level
!           convert temperatures to sea level
            do n=1,ng
              do j=deltat,nhg,deltat
                rrain(n,j)=rrain(n,j)-sta_elv(n)*tlapse
              end do
            end do
c            if(iopt.gt.0)then
c             write(55,5001)
c             do j=deltat,nhg,deltat
c               write(55,202)(rrain(l,j),l=1,ng)
c               write(*,202)(rrain(l,j),l=1,ng)
c             end do
c            endif
c          else
c!           check to see if lapse rate given:
c            if(tlapse.lt.0.0.and.lpsflg0)then
c              tlapse=0.0
c              print*,'WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
c              print*,'lapse rate given in basin\lapse_rates.txt'
c              print*,'but no reference elevations in the _tag.r2c file'
c              print*,'lapse rate is set to 0.0'
c              print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
c            endif
          endif
        
        elseif(fln(14)(1:18).eq.'tempg\narr\narr267')then

          print*,'gone to read_narr'
          datatype='temp  '
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
c          call read_narr(44,14,nhg,conv,ng,datatype)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
          narrflg=.true.
c         print*,'narrflg=',narrflg
        else
!         still read old formats (non ensim)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
          call rdtmp(ng,nhg,conv,ixr,iyr)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
        endif
        
        if(firstpass)then
c         allocate(w(iyr,ixr,ng),dst(ng),stat=iAllocate)
          allocate(w(iyr,ixr,ng),stat=iAllocate)
           if (iAllocate.ne.0) STOP 
     *      'Warning: error with allocation of dst() in tmp' 
        
!         outarray is in areawfo
!         this has to be allocated before calling write_r2c
          allocate(outarray(ycount,xcount),
     *         outarray00(ycount,xcount),outarray06(ycount,xcount),
     *         outarray12(ycount,xcount),outarray18(ycount,xcount),      ! deltaT = 6
     *               mean_temp(ycount,xcount),stat=iAllocate)
          if(iAllocate.ne.0) STOP
     *      'Error with allocation of outarray in moist'      
          do i=1,ycount
            do j=1,xcount
              mean_temp(i,j)=0.0
              do m=1,12
                dly_max_sum(i,j,m)=0.0
                dly_min_sum(i,j,m)=0.0
              end do
            end do
          end do 
          do m=1,12
            mo_count(m)=0.0
          end do
        endif

!       write the header (frame_no = 0)
       
        author='tmp.exe                               '
        name='Gridded Tempratures                     '
        coordsys_temp=coordsys1
        zone_temp=zone1
        datum_temp=datum1
        xorigin_temp=xorigin
        yorigin_temp=yorigin
        xcount_temp=xcount
        ycount_temp=ycount
        xdelta_temp=xdelta
        ydelta_temp=ydelta
        attribute_name='temperature                             '
        attribute_units='degreeCelcius                           ' 
        attribute_type='Gauge                                   '  
        source_file_name=fln(14)                                   

        no_hdrcomments=4
        allocate(hdrcomment(no_hdrcomments),
     *            hdrcomment_value(no_hdrcomments),stat=iAllocateStatus)
        if(iAllocateStatus.ne.0)STOP 
     *      '**Allocation failed for hdrcomments in tmp @ 458**' 
        hdrcomment(1)='#temperature lapse rate in dC/m '
        if(lpsflg1)then
          hdrcomment_value(1)=tlapse
        else
          hdrcomment_value(1)=0.0
        endif
!       added Aug. 1/11 nk:        
        hdrcomment(2)='#radius of influence km         '
        hdrcomment_value(2)=radinfl*al/1000.
        hdrcomment(3)='#smoothing distance km          '
        hdrcomment_value(3)=smoothdist*al/1000.
!       added Jan. 20/12 nk        
        if(elv_mean_flg)then
          no_hdrcomments=4
          hdrcomment(4)='#elv_means_used                 '
          hdrcomment_value(4)=0.0
        else
          no_hdrcomments=3
        endif

!       Write the header (frame_no = 0,classcount=0) gridded tmp's
!     rev. 10.2.46 Jan.  21/19  - NK: Added netCDF capability
        if(FLtype(15).eq.'r2c')then       !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call write_r2c(45,15,nhg,1,0,0,1)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        endif
       
        deallocate(hdrcomment,hdrcomment_value)

!       Write the header for the daily temperature differences
        author='tmp.exe                               '
        name='Gridded Temprature Differences          '
        attribute_name='dailyTemperatureDifferences             '
        attribute_units='degreeCelcius                           ' 
        attribute_type='Gauge                                   '  
        source_file_name=fln(14)                                   

        no_hdrcomments=0
!       Write the header (frame_no = 0)

        print*,'fln(62)=',fln(62)(1:40)


!       Check if there is a file number in the event file       
        open(unit=99,file=fln(62),status='unknown',iostat=ios)
        print*,'ios=',ios
        close(unit=99)
        if(ios.eq.0)then
!         Write the header (frame_no = 0)  daily diff
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call write_r2c(292,62,nhg,1,0,0,1)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          write292flg=.true.
        else
          write292flg=.false.
        endif
!       convert allows a shifting of the entire temperature field up or down
        if(conv.ne.0.0)then
           do nh=deltat,nhg,deltat
              do is=1,ng
!       conditional added May 13/03 nk.
                 if(rrain(is,nh).gt.-50.0)then
                   rrain(is,nh)=rrain(is,nh)+conv
                 endif
              end do
           end do
           conv=1.0
        endif
      
! DISTRIBUTE TEMP ACCORDING TO WEIGHTS AND GENERATE TEMPR\XXXX.TEM FILE
! LOOP THROUGH EACH TIMESTEP
        rdr='temp'
        
        write(6,*)
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
c      print*
             
!       find # days in each year 
        do i=2,21  ! make this 100
          year_now=year1+i-1
          if(mod(year_now,4).eq.0)then
            no_days_year(i)=366
          else
            no_days_year(i)=365
          endif 
          cumm_days_year_end(i)=cumm_days_year_end(i-1)+no_days_year(i)
c      print*,'no days ',year_now,no_days_year(i),cumm_days_year_end(i)
        end do
                       
        
!       LOOP THROUGH TIME
!       LOOP THROUGH TIME
!       LOOP THROUGH TIME
!       LOOP THROUGH TIME (this event)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       timestamp algorithm revised for more than 1 year March 19/16 NK
        day_now=ju_mon(mo1)+day1-1
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
          endif
c      write(600,61000)nh,days_from_start,ju,year_now,month_now,
c     *             day_now,hour_now
61000     format(20i8)
          if(iopt.ge.1)then
            if(mod(nh,nhg/10).eq.0.or.mod(nh,nhg/12).eq.0)
     *        write(*,60000)'distributing hour #',nh,ju,
     *        year_now,month_now,day_now,hour_now
60000         format(a19,6i8)     
          endif

          do i=1,iyr
            do j=1,ixr
              p(i,j)=0.0
            end do
          end do
          
!         reset these each day:
          if(ju.ne.ju_last)then
            do j=1,xcount
              do i=1,ycount
                dly_max_sum(i,j,mo)=dly_max_sum(i,j,mo)+dly_max(i,j)
                dly_min_sum(i,j,mo)=dly_min_sum(i,j,mo)+dly_min(i,j)
              end do
            end do
            mo_count(mo)=mo_count(mo)+1
c      print*,ju,hour_now,
c     *       dly_max(ycount/2,xcount/2),dly_min(ycount/2,xcount/2),
c     *       dly_max(ycount/2,xcount/2)-dly_min(ycount/2,xcount/2)
            do j=1,xcount
              do i=1,ycount
                dly_max(i,j)=-99.0
                dly_min(i,j)=+99.0
              end do
            end do
c            ju_last=ju
          endif
c      print*,'30',nh,mod(nh,24),nh/24,ju
          
          iflag=0
!         IFLAG = 1 IF WEIGHTS HAVE TO BE RECALCULATED
!         THIS HAPPENS WHEN THERE IS MISSING DATA AT ONE OR MORE OF
!         THE STATIONS AND STATIONS HAVE TO BE IGNORED TEMPORARILY
          lflag=0
!         LFLAG = 1 WHEN THERE IS DATA AT AT LEAST ONE GAGE with data
          if(nh.eq.deltat)then
!           CHECK IF THERE IS ANY RECORDED DATA IN 1ST HOUR
            do is=1,ng
              if(rrain(is,nh).gt.cut)lflag=1
            end do
c      print*,'in first timestep'
          else
            do is=1,ng
!             RECALCULATE WEIGHTS IF AVAILABILITY OF GAUGES HAS CHANGED
              if(rrain(is,nh).gt.cut.and.
     * 			rrain(is,nh-deltat).le.cut) iflag=1
              if(rrain(is,nh).le.cut.and.
     *			rrain(is,nh-deltat).gt.cut) iflag=1
!             CHECK IF THERE IS ANY RECORDED DATA THIS HOUR
              if(rrain(is,nh).gt.cut)lflag=1
            end do
          endif
c      print*,'40',nh,mod(nh,24),nh/24,ju,lflag


      lflag=1


!         THIS SECTION IS BYPASSED IF THERE IS NO RECORDED DATA
          if(lflag.eq.1)then
!           LFLAG = 1 WHEN THERE IS DATA AT AT LEAST ONE GAGE with data
c              if(lflag.eq.1.or.nh.eq.deltat)then
            if(iflag.eq.1.or.nh.eq.deltat)then    ! changed to iflag nk jan11/08
!             THIS SECTION IS BYPASSED IF THE STATIONS ARE THE SAME THIS HOUR
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              call weight(nh,ng,ixr,iyr,cut)
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            endif
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call distr(nh,ixr,iyr,ng)
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
c            if(sta_elv(1).gt.0.0.and.tlapse.lt.0.0)then
            if(lpsflg1)then
!             adjust the sea level temperatures for grid elevation
!             adjust the sea level temperatures for grid elevation
!             adjust the sea level temperatures for grid elevation
              do i=1,ycount
                do j=1,xcount
                  p(i,j)=p(i,j)+elev_grid(i,j)*tlapse
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
        
!           OUTPUT TEMPERATURE DATA
!           convert local variables to write module variables
            do i=1,ycount
              do j=1,xcount
                outarray(i,j)=p(i,j)
                mean_temp(i,j)=mean_temp(i,j)+p(i,j)
                dly_max(i,j)=amax1(p(i,j),dly_max(i,j))
                dly_min(i,j)=amin1(p(i,j),dly_min(i,j))
              end do
            end do          
c            print*,mo,dly_min_sum(ycount/2,xcount/2,mo),
c     *                 dly_max_sum(ycount/2,xcount/2,mo),mo_count(mo)
            data_count=data_count+1
            frame_no1=frame_no1+1  
!           It seems like the first frame has to be 1 to be read
!           so this is what nh-deltat+1 does 
!           Write the gridded temperature for this time step
            if(FLtype(15).eq.'r2c')then       
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               call write_r2c(45,15,nhg,1,frame_no1,0,1)
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          endif

!         Version 11   Dec. 05, 2022
!         netCDF output .nc file added         
          if(netCDFoutflg.eq.'y')then                  
c           write(333,*)frame_nc,nh,p(1,1),mod(nh-4,24),mod(nh-16,24)
!           Write to the netCDF format yyyymmdd.nc file
            if(nh.eq.4)then
              ! midnight 1st day only
              outarray00=p  
              frame_nc=frame_nc+1
              outarray=p
c             write(333,*)frame_nc,nh,outarray(1,1),mod(nh-4,24),mod(nh-16,24)
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!              call write_2D(frame_nc,0,.false.,15,'tmp','tmp',6) 
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            elseif(mod(nh-16,24).eq.0)then
              ! noon
              do i=1,ycount
                 do j=1,xcount
                     outarray12=p                           ! noon
                     outarray06=(outarray00+outarray12)/2   ! 6am
                 end do
              end do
              frame_nc=frame_nc+1
              outarray=outarray06
c             write(333,*)frame_nc,nh-6,outarray(1,1),mod(nh-4,24),mod(nh-16,24)
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!              call write_2D(frame_nc,0,.false.,15,'tmp','tmp',6) 
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              frame_nc=frame_nc+1
              outarray=outarray12
c             write(333,*)frame_nc,nh,outarray(1,1),mod(nh-4,24),mod(nh-16,24)
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!              call write_2D(frame_nc,0,.false.,15,'tmp','tmp',6) 
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            elseif(mod(nh-4,24).eq.0)then
              ! midnight
              do i=1,ycount
                 do j=1,xcount
                     outarray00=p
                     outarray18=(outarray12+outarray00)/2.0
                end do
              end do
              frame_nc=frame_nc+1
              outarray=outarray18
c             write(333,*)frame_nc,nh-6,outarray(1,1),mod(nh-4,24),mod(nh-16,24)
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!              call write_2D(frame_nc,0,.false.,15,'tmp','tmp',6) 
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              frame_nc=frame_nc+1
              outarray=outarray00
c             write(333,*)frame_nc,nh,outarray(1,1),mod(nh-4,24),mod(nh-16,24)
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!              call write_2D(frame_nc,0,.false.,15,'tmp','tmp',6) 
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            endif
          endif
             
c        print*,'60',nh,mod(nh,24),nh/24,ju

!           To keep this program backward compatible, skip writing this
!           file is there is no filename given in the event file
            if(write292flg)then
!             write the dly_diff at the end of each day
              if(mod(nh,24).eq.0)then
                do i=1,ycount
                  do j=1,xcount
                    outarray(i,j)=dly_max(i,j)-dly_min(i,j)
                    if(outarray(i,j).le.0)outarray(i,j)=10.0  
!                   happens when there is only 1 temperature/day
                  end do
                end do        
                frame_no2=frame_no2+1  
c               print*,'70',nh,mod(nh,24),nh/24,ju,frame_no
!               Write the daily difference for this time step
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                call write_r2c(292,62,nhg,1,frame_no2,0,1)
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              endif
            endif
            ju_last=ju
c      print*,'99',nh,mod(nh,24),nh/24,ju
          endif   !lflag.eq.1
          ioflag=0
        end do   ! end of time loop
        firstpass=.false.
        print*,'ycount,xcount',ycount,xcount

        if(netCDFoutflg.eq.'y')then                  
!         close the .nc file      
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!          call write_2D(frame_nc,frame_nc,.true.,15,'tmp','tmp',6) 
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        else
          close(unit=15)
          print*,'Closed ',fln(15)(1:50)
        endif

      end do  ! end of event loop

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

!     write the mean_temp.r2c file
!     Added Dec. 05/12  NK
!     write the header (frame_no = 0)
      author='tmp.exe                               '
      name='Gridded Mean Tempratures                '
      attribute_name='temperature                             '
      attribute_units='degreeCelcius                           ' 
      source_file_name=fln(14)                                   

      no_hdrcomments=4
      allocate(hdrcomment(no_hdrcomments),
     *            hdrcomment_value(no_hdrcomments),stat=iAllocateStatus)
      if(iAllocateStatus.ne.0)STOP 
     *      '**Allocation failed for hdrcomments in tmp @ 458**' 
      hdrcomment(1)='#temperature lapse rate in dC/m '
      hdrcomment_value(1)=tlapse
!     added Aug. 1/11 nk:        
      hdrcomment(2)='#radius of influence km         '
      hdrcomment_value(2)=radinfl*al/1000.0
      hdrcomment(3)='#smoothing distance km          '
      hdrcomment_value(3)=smoothdist*al/1000.0
!     added Jan. 20/12 nk        
      if(elv_mean_flg)then
        no_hdrcomments=4
        hdrcomment(4)='#elv_means_used                 '
        hdrcomment_value(4)=0.0
      else
        no_hdrcomments=3
      endif

!     Write the header (frame_no = 0)
      fln(99)='debug\mean_temperature.r2c'
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call write_r2c(99,99,1,1,0,0,1)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      deallocate(hdrcomment,hdrcomment_value)
      do i=1,ycount
        do j=1,xcount
          outarray(i,j)=mean_temp(i,j)/float(data_count)
        end do
      end do
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call write_r2c(99,99,1,1,1,1,1)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      print*
      print*,'NEW:'
      print*,'The file mean_temp.r2c has been written'
      print*,'in the debug directory'
      print*


!     Create the mean daily differences for each month
!     Create the mean daily differences for each month
!     Create the mean daily differences for each month
      do j=1,xcount
        do i=1,ycount
          do m=1,12
          mean_dly_diff_3(i,j,m)=(dly_max_sum(i,j,m)-dly_min_sum(i,j,m))
     *                                /mo_count(m)
          end do
        end do
      end do
      print*,'NEW:'
      print*,'Sample mean daily diff by month:'
      print*,'        month   north         middle          south'
      do m=1,12
        print*,m,mean_dly_diff_3(ycount,xcount/2,m),
     *           mean_dly_diff_3(ycount/2,xcount/2,m),
     *           mean_dly_diff_3(1,xcount/2,m)
      end do
                          
!     write the basin\diff.r2c file
!     first ask if it's wanted
      author='tmp.exe                               '
      coordsys_temp=coordsys1
      zone_temp=zone1
	datum_temp=datum1
	xorigin_temp=xorigin
	yorigin_temp=yorigin
	xcount_temp=xcount
	ycount_temp=ycount
	xdelta_temp=xdelta
	ydelta_temp=ydelta
	name='climate normals - diff                  '
	attribute_units='dC                                     ' 
      attribute_type='derived                         '
      attribute_count=12
      unit_conversion=conv
      source_file_name='yyyymmdd_tag.tb0 files          '          
      no_hdrcomments=0

      open(unit=99,file='junk1',status='unknown')
      write(99,99005)year1,mo1,day1,hour1
99005 format(i4,'/',i2,'/',i2,2x,i2,':00:00')
      rewind 99
	read(99,99006)startdate,starttime
99006 format(2a10)
      close(unit=99,status='delete')
      fln(60)='basin\new_dly_diff.r2c'
      ii=12     
!     write the header:      
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call write_r2c(290,60,1,0,0,0,1)   
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      do ii=1,12
!       classcount is used as the # of months in the year
!       convert local variables to write module variables
        do i=1,ycount
          do j=1,xcount
            outarray(i,j)=mean_dly_diff_3(i,j,ii)
          end do
        end do          

!      SUBROUTINE write_r2c(un,fn,
!     *            no_frames,no_classes,frame_no,class_no,
!     *            no_signf)

!       write the data (mean daily diff):
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call write_r2c(290,60,1,12,1,ii,2)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	end do
      print*
      print*,'NEW:'
      print*,'The file basin\new_diff.r2c in has been written'
      print*,'Please copy this file to basin\mean_dly_diff.r2c'
      print*,'if you wish to use this file'
      if(.not.write292flg)then
        print*,'WARNING:'
        print*,'New files for daily max/min temperature differences'
        print*,'have NOT been written.'
        print*,'Please make a new set of event files with'
        print*,'make_evt.exe to get these files'
        print*
      else
        print*
        print*,'New files for daily max/min temperature differences'
        print*,'have been written. fln=',fln(62)(1:50) 
      endif
      
      
      
      
!     Write the climate mean annual precip amouts for each day of the year
!     Addded jan. 31/2020  NK      
!     Write the header (frame_no = 0)
      no_hdrcomments=0
      fln(699)='tempr\climate_tmp.r2c'
      source_file_name='multiple_events'
d           print*,'writing the header in',fln(699)(1:40)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call write_r2c(699,699,nhg/deltat,1,0,0,1)   
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      frame_no1=0
        do nh=deltat,nhg,deltat
          year_now=2020
          ju=(nh-4)/24+1
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
          hour_now=mod(nh-4,24)     !+deltat
c          print*,nh,hour_now,day_now,month_now
          n=0
          do i=1,iyr
              do j=1,ixr
                  n=n+1
!                 divide by the # of years of data                  
                  outarray(i,j)=Pclimate(n,nh/deltat)/float(ni)
              end do
          end do
          frame_no1=frame_no1+1
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call write_r2c(699,699,nhg/deltat,1,frame_no1,0,1)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        end do
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      print*,'NEW'
      print*,'Wrote: ',fln(699)(1:60)
      print*
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
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
      

! FORMATS:

  202   format(256f10.3)
 1000 format(' Temperature DISTRIBUTION APPLICATION',/)
 1001 format(A30)
 1002 format(10X,I5,1X,F5.2,1X,I6,1X,F5.1)
 1003 format(' HOUR=',2I5)
1100  format('+',' DISTRIBUTING TEMP FOR HOUR ',i5,' OF ',i5,
     +' HOURS OF TEMPERATURE DATA')
 1113 format(' file length =',i5,'hours - record length =',i5,'hours')
 1300 format(99F5.1)
 1301 format(99f5.2)
 5000   format(' converted data to sea level:')
 5001   format(/' converted data to sea level:')
 6000 format(99F5.0)
    
        print*,'*******************************************************'
        print*,'*                  WATFLOOD (R)                       *'
        print*,'*            Version 11   Dec. 05, 2022               *'
        print*,'*           (c) N. Kouwen, 1972-2022                  *'
        print*,'*******************************************************'


c      if(stopflg.eq.'y')then
c        print*,' normal ending'
c        pause ' Hit enter to exit window'
c        print*
c        stop 
c      else
        stop ' normal ending'
c      endif



      END PROGRAM tmp
