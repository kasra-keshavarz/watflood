      suBROUTINE sub(jan,e1,smc5,conv,scale,icase,smok,optlow,igrdshft)

!***********************************************************************
!    Copyright (C) 1987 by Nicholas Kouwen 
        
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

!  THIS SUBROUTINE ORGANIZES ALL THE CALCULATIONS.

!     REV. 7.51 oct.  08/95   -   revise init channel flow in SUB
!     REV.        aug   07/96 - to read in mica inflow file and
!                       use mica data for optimization when
!                       maxn=2002
!     REV. 7.80 Oct.  29/96   - spl7 added yymmdd.rin for res inflows
!                     - unit = 49   fln = 19
!     REV. 7.9  Dec   18/96   - Todd's Evaporation, changed above
!                       res inflows to unit=99  fln=21
!     REV. 8.1  - Feb.    15/96 - TBC & RSM (to be continued & resume) 
!     REV. 8.61 - Dec.    12/97 - added contflg for statistics cont'n
!     REV. 8.75 - Apr.    27/89 - took da out of the resume file
!     REV. 8.90 - Dec.    04/98 - input to memory for opt runs
!     REV. 8.91 - Dec.    07/98 - read rdevt in sub as well as spl!
!     REV. 8.94 - Feb.    01/99 - crseflg to read resume & snow course
!     REV. 8.82 - July    10/98 - added runoff output option: routeflg
!     REV. 8.82 - July    10/98 - added runoff output option: routeflg
!     REV. 8.83 - Oct.    23/98 - added step to the lst argument list
!     REV. 8.86 - Nov.    02/98 - fixed opt problem found by ted.
!     REV. 8.96.1 May 12/99 - added ireport for reporting interval
!     REV. 8.98   July    15/99 - met grid increased
!     REV. 9.00    Mar.  2000 - TS: CONVERTED TO FORTRAN 90 
!     rev. 9.1.18  Jun.  03/02  - Added sub-watershed modelling capability
!     rev. 9.1.21  Jun.  28/02  - Added wetland storage & outflow to the wfo file
!     rev. 9.1.23  Jul.  23/02  - Added control for nudging in event #1
!     rev. 9.1.28  Sept. 19/02  - Added shedlfg to replace the bsnm.shd file
!     rev  9.1.29  Oct.  24/02  - Added q1, qint & drng to wfo file
!     rev. 9.1.35  Dec.  26/02  - Added wetland & channel heights to the wfo file
!     rev. 9.1.38  Mar.  31/03  - revised str header and routing dt selectable
!     rev. 9.1.42  May     31/03  - Tracer module added - first try
!     rev. 9.1.44  Jun.  11/03  - Added Cumulative precip to the wfo file
!     rev. 9.1.45  Jun.  11/03  - WATROUTE: runoff, recharge and leakage files added
!     rev. 9.1.47  July  11/03  - TS: Tracer s/r call line modified 
!     rev. 9.1.50  Jan.  14/04  - NK: version number added to the wfo_spec.txt file
!     rev. 9.1.51  Jan.  28/04  - NK: added iz.ne.jz conditional to ENSIM output  
!     rev. 9.1.56  Jun.  18/04  - NK: write new rel & rin files to resrl\newfmt folder.
!     rev. 9.1.59  Jul.  15/04  - NK: split rerout into two parts: rdresv & rerout
!     rev. 9.1.67  Oct.  21/04  - NK; added unit 80 for lake_stor & lake_flow
!     rev. 9.1.81  Apr.  04/05  - NK: added sublimation,et and etfs to wfo file
!     rev. 9.2.05  Jul.  15/05  - NK: reversed order of reading resume file 
!     rev. 9.2.07  Jul.  29/05  - NK: soilinit moved from runoff to sub 
!     rev. 9.2.20  Oct.  28/05  - NK: WFO_SPEC - reporting start & finish times 
!     rev. 9.2.25  Dec.  13/05  - NK: ENSIM r2c gridded soil moisture 
!     rev. 9.2.28  Jan.  30/06  - NK: Added low slope a4 for grids with water
!     rev. 9.2.30  Feb.  07/06  - NK: Added class_distribution.txt to output
!     rev. 9.2.35  Mar.  22/06  - NK: Glacier flow bypasses wetlands
!     rev. 9.2.39  May.  09/06  - NK: thr added to route & rerout arg list
!     rev. 9.3.02  Jul.  18/06  - NK: converted runof, rchrg & lkage to r2c
!     rev. 9.3.08  Jan.  15/07  - NK: added lzs_init_new.r2c output to sub.for
!     rev. 9.3.10  Jan.  29/07  - NK: routing pars changed to gridded values
!     rev. 9.4.08  May.  29/07  - NK: changed baseflow argument list
!     rev. 9.4.14  Jul.  09/07  - NK: added lake loss file 
!     rev. 9.5.01  Oct.  15/07  - NK: added wetland continuity check
!     rev. 9.5.03  Dec.  09/07  - NK: added reads for precip isotopes
!     rev. 9.5.12  Feb.  13/08  - NK: added evaporation input file with read_r2c
!     rev. 9.5.18  Mar.  03/08  - NK: added conv to options & sub argument list
!     rev. 9.5.19  Mar.  05/08  - NK: prevented use of tracer * iso models with nudging
!     rev. 9.5.21  Mar.  06/08  - NK: fixed dtmin for first time step each event
!     rev. 9.5.22  Mar.  12/08  - NK: added grdflg to print gridded flow, swe & evap
!     rev. 9.5.30  May.  26/08  - NK: conv back in read_rain & process_rain arg. list
!     rev. 9.5.31  May.  27/08  - NK: moved totsnw(n) computation in sub
!     rev. 9.5.33  Sep.  12/08  - NK: added column labels for grapher in flow_station_location.xyz
!     rev. 9.5.35  Sep.  22/08  - NK: moved flow_sta_location to flowinit
!     rev. 9.5.45  Dec.  16/08  - NK: added various error calculations - user's choice with errflg
!     rev. 9.5.48  Dec.  26/08  - NK: added event_fln() to allow unlimited events
!     rev. 9.5.50  Jan.  05/09  - NK: read evap data for reaches only
!     rev. 9.5.51  Jan.  13/09  - NK: added reading yyyymmdd_ill.pt2 foa all lakes
!     rev. 9.5.56  Mar.  26/09  - NK: Fix bug with month in yearly events
!     rev. 9.5.57  Apr.  13/09  - NK: added ntrlflg for natural lake flows
!     rev. 9.5.62  Sep.  04/09  - NK: new tb0 file for DW routing
!     rev. 9.5.65  Sep.  23/09  - NK: change class frac to whole basin values
!     rev. 9.5.73  Oct.  12/09  - NK: bypass using lake levels when optimizing
!     rev. 9.5.75  Oct.  26/09  - NK: commented "deallocate in sub for watroute reads
!     rev. 9.5.77  Oct.  26/09  - NK: fixed some inits for out of basin gauges
!     rev. 9.5.80  Dec.  20/09  - NK: added swe_locations.txt file for swe input
!     rev. 9.6.05  Apr.  06/10  - NK: added store_error_flag for -ve storage grids
!     rev. 9.7.00  May.  26/10  - NK: dds with pre-emption
!     rev. 9.7.03  Jun.  24/10  - NK: normalized SSE with station Qmean**2
!     rev. 9.7.08  Sep.  21/10  - NK: revised mean squared error weighting for DDS
!     rev. 9.7.17  Jan.  05/11  - NK: Fixed diversions outside sub-basin
!     rev. 9.7.19  Jan.  18/11  - NK: Added sensitivity analysis
!     rev. 9.7.22. Mar.  07/11  - NK: Changed diversion code: give/route take/rerout
!     rev. 9.7.27  May.  26/11  - NK: Add lake_ice_factor
!     rev. 9.8.03  Aug.  08/11  - NK: chech no of mean observed flows in file are ok
!     rev. 9.8.09  Nov.  22/11  - NK: nopt(l)=0 for area_error(l) > 10%
!     rev. 9.8.14  Jan.  27/11  - NK: dds_penalty added for swe not to zero in summer
!     rev. 9.8.16  Mar.  21/11  - NK: reinstate reservoir inflow error for dds
!     rev. 9.8.21  Jun.  18/12  - NK: Added swe observed date & report
!     rev. 9.8.22  Jul.  17/12  - NK: Added resetflg to reset cumm. precip Sept.1
!     rev. 9.8.23  Aug.  03/12  - NK: Added resinid1flg to use resinflg for id=1
!     rev. 9.8.44  Jan.  31/13  - NK: fixed bug in sub.f : uninitialized course_calc(n,j)
!     rev. 9.8.49  Feb.  20/13  - NK: Added n=municipal & irrigation withdrawals
!     rev. 9.8.51  Mar.  11/13  - NK: Link skiplines in s/r stats to value1 in the str file
!     rev. 9.8.91  Oct.  30/13  - NK: Got rid of lzs_init.r2c - data is in flow_init.r2c already
!     rev. 9.8.92  Nov.  06/13  - NK: Changed output file swe.txt to swe.csv
!     rev. 9.8.93  Nov.  12/13  - NK: Added the routing initialization with yyyymmdd_fli.r2c
!     rev. 9.9.06  Jan.  08/14  - NK: Add daily differences to Harfreaves ETHarg.f
!     rev. 9.9.09  Feb.  24/14  - NK: Fixed reading the time stame in r2c frame headers
!     rev. 9.9.11  Mar.  20/14  - NK: Added lake_level_init.pt2 file for a resume
!     rev. 9.9.21  Jul.  27/14  - NK: Added allocation for outarray in sub
!     rev. 9.9.26  Sep.  16/14  - NK: Added precip adjust for forecast & fcstflg
!     rev. 9.9.36  Nov.  03/14  - NK: Revised error message for daily diff choices
!     rev. 9.9.64  Apr.  08/15  - NK: DDS bypass in sub for single runs
!     rev. 10.1.11 Dec.  11/15  - NK: Revised ice factor initialization and calculation   
!     rev. 10.1.58 Dec.  06/16  - NK: corrected tdum > tdum1 for modelflg=i
!
!  id     - storm id number
!  l      - station number
!  ni     - total # of storms
!  nl     - forecast period = total length of the streamflow file
!  no     - total # of stations
!  nr     - the # of hours rainfall is available
!  mhrd   - length of streamflow record = time elapsed since start
!  mhtot  - the simulation length in hours or is the length of rainfall
!         input for simulated forecasts
!  aintvl - precip interval in hours
!  sintvl - precip. interval in seconds
!  tot2   - the total average precipitation over the watershed


!***********************************************************************

      use areacg
      use area_watflood

C///////////////////////// 
C// Added by Dave
      USE EF_module
C// End Dave addition
C/////////////////////////

      implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

        DIMENSION     :: smc5(16)
        CHARACTER(14) :: date
        CHARACTER(3)  :: eofmark
        CHARACTER(1)  :: lineflg,smok,answer
        character(20) :: junk
        REAL(4)       :: optlow,time,time2,tot1,qwert,conv,scale,
     *           smc5,tj1,clock,t,thr,dtmin,dtmax,div,aintvl,sintvl,
     *           tot2,e1,tdum,tdum1,qtemp,diff2,sdlz,dlz,
     *           wfo_spec_version_number,
     *           route_dt,sec_div,hr_div
      integer         :: frame_no1,frame_no2,frame_no3,frame_no4
      integer         :: frame_no5,frame_no6,frame_no7,frame_no8

      real*4          ::  input,runoff,router,output,last,now
      real*4          ::  hhhhh,mmmmm,sssss
      integer         ::  yyyyy,moooo,ddddd
      
  
      INTEGER*4       :: rbin,inta,block,no1,classcount1,na1,
     *             ycount1,xcount1,
     *             iallcnt1,iallcnt2,n1,ii1,jan,m,ios,n,iallocate,
     *             l,ii,juold,julast,jj,kk,lun,nhr,nhf,nfg,
     *             i,j,icase,iz,jz,nchr,mz,ju,mon,
     *             iwest,ieast,isouth,inorth,nocolumns,nno,k,
     *             nu,igrdshft,minsnwflg,oldjan,flgevp22,
     *             noresv1,nj,npick,n_trick,no_frames,
     *             no_dt,iDeallocate,id_last,nd,MC_id
      integer       :: swe_l(99),no_swe_l,no_swe_obs
      integer       :: nrows,ncols
      integer       :: flnNum
      integer       :: n_max,nold,nlines_old,nlines

      CHARACTER(10) :: ctime
      CHARACTER(8)  :: cday
      CHARACTER(12) :: outlet_type
      logical   :: exists,newpafflg,firstpass,errflg_store,msgflg
      logical   :: dataflg,newDataFlag

        CHARACTER(10) :: coordsys
!        INTEGER      :: xcount,ycount
!        REAL     :: xorigin,yorigin,xdelta,ydelta,
      real     :: a66,sum_mean_flow,class_sum
      real :: swe_error,swe_penalty
      real :: ha,fpw,kdn,nratio,xtemp,ytemp,ztemp
      real :: domain_precip,domain_area
      real :: store_live

      data firstpass/.true./
      data id_last/0/
      data no_swe_l/0/
      DATA iallcnt1/0/
      DATA iallcnt2/0/
c      DATA col1/'b','d','f','h','j','l','n','p','r','t','v','x','z'/
c      DATA col0/'c','e','g','i','k','m','o','q','s','u','w','y','a'/
c      DATA col2/' ','a','b','c','d','e','f','g','h','i','j','k','l'
c     *       ,'m','n','o','p','q','r','s','t','v','u','w','x','y'/


!     NOTE: FOR MONTHLY RUNS, DIMENSIONS CAN BE CHANGED FROM 
!         3,8784,500  TO  12,366,3000

!>>>>>>>>>>>>>  AB: STUFF FOR ENSIM
      INTEGER(4) :: wfo_yy,wfo_mm,wfo_dd,wfo_hh,wfo_mi,wfo_ss,
     *             wfo_ms
      INTEGER(4) :: wfo_seq
      real*4     :: dds_penalty

      character(256) line,tmpLine
      CHARACTER(100) sys_command ! Command line input string

!      replaced by areawfo
!      REAL*4, DIMENSION(:,:), ALLOCATABLE :: outwfo

!     DDS error functions
     	real*4    :: log_qhyd,log_qinfl,log_mean_obs           
      real*4    ::  sum_swe
      real*4    ::   score
      real*4    ::   temp_value,sum_sq_error,
     *               temp_value_sum,dds_error

!     WFO IO FUNCTIONS
      INTEGER :: wfo_write_attribute_data
      INTEGER :: wfo_write_timestamp
      
!>>>>>>>>>>>>>

      if(iopt.eq.2)print*,' In sub after definitions'

!     RESET SETS ALL INITIAL VARIABLES

!     CHECK FILES MODE    iopt=99
!     FOR IOPT=99 NL AND MHRD ARE SET TO KT AND THE PROGRAM WILL RUN
!     FOR ONE TIME STEP ONLY - THIS WILL CHECK AND ECHO ALL INPUT FILES
!     >> VERY HANDY FOR CHECKING DATA FILE PRIOR TO LONG RUNS

      id=1  ! just so it's not some value from before

      if(firstpass)then
      
        input=0.0
        runoff=0.0
        router=0.0
        output=0.0
        min_flow_cutoff=0.01    ! error not calculated below this flow

      
!     rev. 10.1.58 Dec.  06/16  - NK: corrected tdum > tdum1 for modelflg=i
        tdum1=1000.*step2/3600.

!       rev. 9.8.21  Jun.  18/12  - NK: Assed swe observed date & report
        if(iopt99)open(unit=951,file=filename(951),iostat=ios)
        if(ios.ne.0)then    ! added Nov. 10/14  nk
          print*
          print*,'Unable to open file',filename(951)(1:40)
          print*,'Possible cause(s):'
          print*,'file in use by another application'
          print*,'or target directory does not exist'
          stop 'Program aborted in sub.f @ 236'
        endif

!     rev. 9.9.00  Dec.  08/13  - NK: Added Lake Evaporation model
        if(lakeflg.eq.'y')then
          allocate(WindSpd(na),WindDir(na),stat=iAll)
	    if(iAll.ne.0) STOP 'Error allocating wind arrays in sub' 
        endif
          
c        if(netCDFflg)then
             allocate(QQsum(na),stat=iAll)
	       if(iAll.ne.0) STOP 'Error allocating QQsum in sub' 
             do n=1,na
                  QQsum(n)=0.0
             end do 
c        endif

!     rev. 9.8.80  Aug   09/13  - NK: Added withdraw.r2c output file in route.f
!       write the header for the withdraw.r2c file 
      if(iopt99)then
        author='watflood                    '
        name='Irrigation withdrawals             '
        coordsys_temp=coordsys1
!       GreenKenue uses LatLong - code below uses LATLONG
        if(coordsys_temp.eq.'LatLong   ')coordsys_temp='LATLONG   '
        if(coordsys_temp.eq.'Cartesian ')coordsys_temp='CARTESIAN '
        zone_temp=zone1
        datum_temp=datum1
        xorigin_temp=xorigin
        yorigin_temp=yorigin
        xcount_temp=xcount
        ycount_temp=ycount
        xdelta_temp=xdelta
        ydelta_temp=ydelta
        attribute_name='withdrawalslow                    '
        attribute_units='cms                           ' 
        attribute_type='Flow                          '
        unit_conversion=1.0   
	    startdate='unknown   '
	    starttime='unknown   '
        source_file_name='various rff,rch,lkg files'     
        no_frames=2      
        frame_no1=0       ! write the header
!       no_frames=2 tricks write_r2c to write frame .....
!       write the header for gridded withdrawal flows
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        fln(23)='results\withdraw.r2c'
        print*,fln(23)(1:40)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call write_r2c(23,23,ni*12,0,0,0,1)   
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      endif
!       added for BC Hydro
        if(iopt99)then
        open(unit=29,file=filename(29),status='unknown',iostat=ios)
        if(ios.ne.0)then    ! added Nov. 10/14  nk
          print*
          print*,'Unable to open file',filename(29)(1:40)
          print*,'Possible cause(s):'
          print*,'file in use by another application'
          print*,'or target directory does not exist'
          stop 'Program aborted in sub.f @ 293'
        else
          write(29,*)'       ID domain_precip'     
        endif
        
        domain_area=0.0 
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          domain_area=domain_area+grid_area(n)
        end do
        domain_precip=0.0
        endif
        
!       rev. 9.7.15  Dec.  14/10  - NK: Create reduced precip & temp files for sub-basins
!       This will only happen is the user has created a folder & file: 'radcl\new_grid\junk'
        open(unit=99,file='radcl\new_grid\junk',
     *             status='unknown',iostat=ios)
        if(ios.eq.0)then
          new_precip_grid_flg=.true.
          print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
          print*,'re-gridded met file will be written in radcl\new_grid'
	    print*,'delete the \new_grid folder to stop this message'
          print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        else
          new_precip_grid_flg=.false.
        endif
        close(unit=99,status='delete')
        open(unit=99,file='tempr\new_grid\junk',
     *             status='unknown',iostat=ios)
c        inquire(FILE='radcl\new_grid\junk',EXIST=exists)
        if(ios.eq.0)then
          new_temp_grid_flg=.true.
          print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
          print*,'re-gridded tem file will be written in tempr\new_grid'
	    print*,'delete the \new_grid folder to stop this message'
          print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        else
          new_temp_grid_flg=.false.
        endif
        close(unit=99,status='delete')
!       end rev. 9.7.15  Dec.  14/10  

!       assume we do a swe report untuil we do not find the swe.tb0 files
        courseflg=.true.

!     rev. 9.8.51  Mar.  11/13  - NK: Link skiplines in s/r stats to value1 in the str file
        skiplines=0
        


!     rev. 9.9.26  Sep.  16/14  - NK: Added precip adjust for forecast & fcstflg
!       this needs to be read along with the event file; include with read_evt()?
!       ? ... or flow_init -- either has correct scope and occurs at start of event
!       these declarations are in area_watflood.f
!       logical :: fcst_exists,fcst_mode ! declare  in header of 
!       integer :: fcst_yr, fcst_hr0, fcst_days2avg
!       real*4 :: fsct_snow_adj, fcst_rain_adj
! --    the rest of this located in main scope; could add to separate sub

        fcst_mode=.false.
        INQUIRE(FILE='fcst_params.txt',EXIST=fcst_exists)
        IF(fcst_exists)THEN
          open(unit=99,file='fcst_params.txt',iostat=ios)    ! 99 is a scratch file #
!           -- unit number selected randomly (fix!);
!           -- can't figure out filename convention
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file  fcst_params.txt'
            print*,'Possible cause(s):'
            print*,'file in use by another application'
            stop 'Program aborted in sub.f @ 358'
          endif
          read(99,*)fcst_yr ! read each line and assign to forecast parameter variables
          read(99,*)fcst_hr0
          read(99,*)fcst_days2avg
          read(99,*)fcst_snow_adj
          read(99,*)fcst_rain_adj
          close(unit=99,status='keep')
!         check if this is the current year
          if(fcst_yr.eq.year_now)then          ! Trish believes 'year1' is the current year; might be 'year' or 'yearnow'?
!                                                 look in timer.f
            fcst_mode=.true.
          end if
        END IF    
!       end rev. 9.9.26  Sep.  16/14  - NK: Added precip adjust for forecast & fcstflg

!     rev. 10.2.01 Oct   08/17  - NK: Moved ruleflg from sub.f to spl.f
c!     rev. 9.9.65  Apr.  `3/15  - NK: Added rule s/r; resrl\rules.txt & ruleflg
c        INQUIRE(FILE='resrl\rules.ts5',EXIST=exists)
c        if(exists)then
c!         this means there are rules for some or all of the lakes & reservoirs
c!         Note:  not necessarily all        
c          ruleflg=.true.
c        else
c!         this means there are no rules at all.        
c          ruleflg=.false.
c        endif      




!     rev. 10.1.08 Dec.  04/15  - NK: Added msg re: replacing "mean_observed_flows.txt"' 
        if(dds_flag.eq.0.and.id.eq.ni)then
          inquire(FILE='mean_observed_flows.txt',EXIST=exists)
          if(exists)then
            print*
            print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
            print*,'found the file mean_observed_flows.txt'
            print*
            print*,'If you wish to replace this file with the current'
            print*,'the old file must be deleted'
            print*
            print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
          endif
        endif
        
        if(iopt99)then
          open(unit=1000,file='debug\benchmark_time.txt',iostat=ios)
          if(ios.ne.0)then
            print*,'Problems opening benchmark_time.txt'
            print*,'needs debug\ directory in the working dir'
            pause 'paused in sub @ 427'
          endif
        endif
        
!     rev. 10.2.46 Jan.  21/19  - NK: Added netCDF capability
!       New here        
        if(.not.allocated(p).and.FLtype(10)(1:2).eq.'nc')then
           allocate(p(ycount,xcount),stat=iAllocate)
           if(iAllocate.ne.0) 
     *       STOP 'Error with allocation of p in sub @ 117'
        endif
        
!       Stuff for FEWS & netCDF files
        deltaT2 = 1  ! Will be reset as soom as we read the first precip file
        deltaT3 = 1  ! Will be reset as soom as we read the first temp file
        
      endif   ! firstpass
 
      if(numa.ne.0)then
c        iopt=0
        trcflg='n'
        ensimflg='n'
        initflg='n'
        modelflg='n'
      endif

!     JAN=2 SUBSEQUENT PASSES
!     JAN=3 LAST PASS  - SET BELOW

c      jan=1

      m=1
      tot1=0.0
      totaltime=0.0       ! used for ensim time series
	mo=mo1              ! added jan 22/11 nk for repeated runs
      wfo_open_flg='n'
!     rev. 9.5.21  Mar.  06/08  - NK: fixed dtmin for first time step each event
      dtmin=a6

!     Write the header in rte.txt
      IF(iopt.ge.1.and.n.eq.nnprint)write(55,5551)
5551  format('   id   time         at     qi1     qi2'
     *          '     qo1     qo2   store1     store2')

d      if(iopt.eq.2)print*,' In sub - gone to rdflow @ 147'

!     get no and nopt for allocations etc.  nk Apr. 8/03 
!     read the str header for allocation and initialization purposes

!      call rdflow('0',date)     
 
      if(IsFileTypeTB0(fln(6)).or.FLtype(6)(1:2).eq.'nc') then
d       print*,'reading flow header'      
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call read_flow_ef('0',date)  !EnSim compatible tb0 file
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      else
         print*,'Old format .str files not accepted'
         print*,'Please create yyyymmdd_str.tb0 files & rerun'
         stop 'Program aborted in sub @ 192'
      endif

!     rev. 9.5.19  Mar.  05/08  - NK: prevented use of tracer * iso models with nudging
      msgflg=.false.
      do l=1,no
        if(nopt(l).eq.2)then
          frcflg='n'
c          trcflg='n'
          if(msgflg)then
            print*,'Flow nudging is turned on in flow station ',l
          else
            print*,'also at station',l
          endif
          msgflg=.true.
        endif
      end do
      if(msgflg)then
        print*
        print*,'WARNING:'
        print*,'Tracer and isotope model data may have issues at'
        print*,'locations downstream from stations with nudging'
        print*,'due to a possible lack of continuity'
c        print*,'trcflg and frcflg set to `n`'
        print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      endif
 
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     Section added to allow for lengthened routing time step for large grids
      if(irdt.gt.kt)then
        a6=float(irdt)*3600.0
        write(51,*)' Warning'
        write(51,*)' Min time step a6 changed to ',a6
        write(51,*)
        write(*,*)
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)' WARNING'
        write(*,*)' deltaT (str file) irdt=',irdt
        write(*,*)' Min time step a6 changed to ',a6
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)
d        pause 'hit enter to continue - in sub @167'
      endif
d     if(iopt.eq.2)print*,' In sub - back from rdflow @ 152'

      if(iallcnt1.eq.0)then
        iallcnt1=1
      endif   ! iallcnt1=0
      if(iallcnt2.eq.0)then
        iallcnt2=1
      endif   ! iallcnt2=0

d     if(iopt.eq.2)print*,' In sub after allocations 1'

      if(snwflg.eq.'y')then
         do n=1,naa
          do ii=1,classcount
!            will print warning if snow never disappears
             snowcmin(n,ii)=1.0e+32
          end do
         end do
      endif
d      if(iopt.eq.2)print*,' In sub after allocations 2'

!     SINGLE RUN USING SOIL MOISTURE GRID IS THE DEFAULT.
!     NO SOIL MOISTURE OPTIMIZATION - THIS CAN BE CHANGED WITH 
!     SETTING ICASE=-1 IN THE PARAMETER FILE & SM GRID WILL BE IGNORED

!     SAVE THE ORIGINAL VALUE:

      flgevp22=flgevp2

d     if(iopt.eq.2)print*,' In sub after allocations 3'

      do n=1,naa
        rechrg(n)=0.0
!       qstream & strloss need to be initialized for watroute
        qstream(n)=0.0
        strloss(n)=0.0
!        rh(n)=.50   ! moved to rdtemp 28/12/04 nk
        
      end do
      juold=0

c!     added mar 28/06  nk
c!     rev. 9.2.35  Mar.  22/06  - NK: Glacier flow bypasses wetlands
c      if(glacier_class_number.ne.0)then
c!       glacier_class_number is assigned in the rdpar file
c        do n=1,naa
c          if(aclass(n,glacier_class_number).gt.0.0)then
c!         there is a glacier in this grid
c          glacier_flag(n)='y'
c          else
c!         there is no glacier in this grid
c          glacier_flag(n)='n'
c          endif
c        end do
c      endif

c      open(873,file='aet.out',status='unknown')

d     if(iopt.eq.2)print*,' In sub before user check'

!     rev. 9.1.46  Jul.  17/03  - WATFLOOD LITE incorporated 
c      if(ichsm.eq.4)then
c!       it's ok to go thru here for ichsm=3 for the first event
c        if(ni.gt.1)then
c          print*,' For WATFLOOD LT, no of events allowed is 1'
c          print*,' Program will abort after first event.'
c          print*
c          ni=1
c!          pause ' Hit enter to continue - in SUB @ 349'
c        endif
c      endif

!     call decread(ha,fpw,kdn,nratio)

 !     NK - ALLOCATION OF AREAwfo ARRAYS
      if(iallcnt1.eq.1)then
        allocate(outwfo(xcount,ycount),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *    'Error with allocation of ensim arrays in sub'    
          iallcnt1=2
      endif
!
!      if(routeflg.eq.'y')then
      if(iallcnt2.eq.1)then
        allocate(outarray(ycount,xcount),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *    'Error with allocation of ensim arrays in sub'    
          iallcnt2=2
      endif
!     endif

!     Initialize all grids for write_r2c
!      if(routeflg.eq.'y')then
        do i=1,ycount
          do j=1,xcount
          outarray(i,j)=0.0
          end do
        end do
!      endif

      if(iopt.eq.2)print*,'In sub before writing header / gridflow.r2c'

      if(iopt99)then
      if(numa.eq.0.and.dds_flag.eq.0)then
        author='watflood                    '
        name='Gridded Channel Flow (SPL)          '
        coordsys_temp=coordsys1
!       GreenKenue uses LatLong - code below uses LATLONG
        if(coordsys_temp.eq.'LatLong   ')coordsys_temp='LATLONG   '
        if(coordsys_temp.eq.'Cartesian ')coordsys_temp='CARTESIAN '
        zone_temp=zone1
        datum_temp=datum1
        xorigin_temp=xorigin
        yorigin_temp=yorigin
        xcount_temp=xcount
        ycount_temp=ycount
        xdelta_temp=xdelta
        ydelta_temp=ydelta
        attribute_name='channel_inflow                    '
        attribute_units='mm                            ' 
        attribute_type='Runoff                        '
        unit_conversion=1.0   
        source_file_name='various rff,rch,lkg files'     
        frame_no2=0
!       write the header for gridded channel flow (gridflow.r2c)
!     rev. 9.9.23  Aug.  10/14  - NK: Uniform arguments for write_r2c
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c        call write_r2c(72,72,no_frames,1,frame_no2,1,1) 
!       use mhtot+1 so the file statys open until the end of the last event.  
        call write_r2c(72,72,0,0,0,0,1)   
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(ensimflg.eq.'n')then
          print*
          print*
          print*,'disabled output to',fln(72)(1:30)
          write(72,*)'to write this file, set ensimflg= y or a'
        endif

        name='Snow Water Equivalent - weighted        '
        attribute_name='Weighted swe                  '
        attribute_type='SWE                           '
        source_file_name=fln(10)   ! met file
        frame_no4=0      
!       write the header for weighted swe
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call write_r2c(61,61,0,0,0,0,1)   
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(ensimflg.eq.'n')then
          print*,'disabled output to',fln(61)(1:30)
          write(*,*)'to write these files, set ensimflg= y or a'
          write(61,*)'to write this file, set ensimflg= y or a'
          print*
          print*
        endif

!     rev. 9.5.09  Feb.  12/08  - NK: added evap.r2c to the output files
        name='Total evaporation - weighted        '
        attribute_name='Weighted evaporation              '
        attribute_type='EVAP                          '
        source_file_name=fln(10)   ! met file
        frame_no5=0
!       write the header for weighted gridded et
!     rev. 9.9.23  Aug.  10/14  - NK: Uniform arguments for write_r2c
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c        call write_r2c(100,100,no_frames,1,frame_no5,1,1)   
        call write_r2c(100,100,0,0,0,0,1)   
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(ensimflg.eq.'n')then
          print*,'disabled output to',fln(100)(1:30)
          write(*,*)'to write these files, set ensimflg= y or a'
          write(100,*)'to write this file, set ensimflg= y or a'
        endif

        name='isostrcon2                  '
        attribute_name='Iso_concentration_stream          '
        attribute_type='Concentration                 '
        source_file_name=fln(10)   ! met file
        fln(200)='results\isostrconc2.r2c'
        frame_no6=0
!       write the header        
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c       call write_r2c(200,200,no_frames,1,frame_no6,1,1)   
c        call write_r2c(200,200,mhtot+8784,1,frame_no6,0,1)   
        call write_r2c(200,200,0,0,0,0,1)   
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      endif
      endif   ! iopt99


!     rev. 9.7.00  May.  26/10  - NK: dds with pre-emption
!     rev. 9.7.00  May.  26/10  - NK: dds with pre-emption
      if(dds_flag.eq.1)then
!       pre-emption 
        inquire(FILE='dds\pre-emption_value.txt',EXIST=exists)
        if(exists)then
!         read the pre-emption value which is the best solution so far
          open(unit=99,file='dds\pre-emption_value.txt',
     *          status='unknown',iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file  dds\pre-emption_value.txt'
            print*,'Possible cause(s):'
            print*,'file in use by another application'
            print*,'or target directory does not exist'
            stop 'Program aborted in sub.f @ 667'
          endif
          read(99,*)pre_emption_value
          print*,'from dds\pre-emption_value.txt'
          print*,'pre_emption_value read=',pre_emption_value
          close(unit=99,status='keep')
        else
!         first evaluation - set the pre-emption value to a large number
          open(unit=99,file='dds\pre-emption_value.txt',
     *          status='unknown',iostat=ios)
          if(ios.ne.0)then
            print*
            print*,'Unable to open file  dds\pre-emption_value.txt'
            print*,'Possible cause(s):'
            print*,'file in use by another application'
            print*,'or target directory does not exist'
c            pause 'hit ctrl C to abort in sub @ 677'
!     rev. 9.9.64  Apr.  08/15  - NK: DDS bypass in sub for single runs
            print*,'It looks like you are making a single run with'
            print*,'a par file generated by DDS or the PS'
            print*,'You can avoid this interruption by changing'
            print*,'ddsflg = 0 in the par file'
            print*
            pause 'For single run, to continue - hit Enter'
            dds_flag=0
          endif
          pre_emption_value=1.0E+35
          write(99,*)pre_emption_value
          close(unit=99,status='keep')
          print*,'first pre_emption_value written=',1.0E+35  ! changed 35 > 32 nk Jan. 22/11
        endif

      endif


!     rev. 9.7.11  Nov.  22/10  - NK: added monthly_climate_deltas.txt file
	climate_delta_flg=.false.     !default
      inquire(FILE='basin\monthly_climate_deltas.txt',EXIST=exists)
	if(exists)then
        open(unit=99,file='basin\monthly_climate_deltas.txt',
     *	  status='unknown',iostat=ios)
        if(ios.ne.0)then    ! added Nov. 10/14  nk
          print*
          print*,'Unable to open file  basin\monthly_climate_deltas.txt'
          print*,'Possible cause(s):'
          print*,'file in use by another application'
          print*,'or target directory does not exist'
          stop 'Program aborted in sub.f @ 700'
        endif
	  read(99,*)(monthly_temperature_delta(j),j=1,12)
	  read(99,*)(monthly_precipitation_delta(j),j=1,12)
	  close(unit=99,status='keep')
	  print*
	  print*,'Montly climate deltas found in file'
	  print*,'basin\monthly_climate_deltas.txt'
	  print*,'All temperatures will be adjusted by dC:'
	  print*,monthly_temperature_delta
	  print*,'All precipition will be adjusted by %:'
	  print*,monthly_precipitation_delta
	  print*
	  do j=1,12
	    monthly_precipitation_delta(j)=
     *               1.0+monthly_precipitation_delta(j)/100.0
	  end do
	  print*,'Do you want to continue with these adjustments?  y/n'
	  read*,answer
	  if(answer.eq.'y')then
	    climate_delta_flg=.true.
	    print*,'All temperatures will be adjusted as per table'
	    print*
	    print*
	  else
	    do j=1,12
	      monthly_temperature_delta(j)=0.0
	      monthly_precipitation_delta(j)=1.0
         end do
	  endif
	endif

!     rev. 9.8.08  Nov.  18/11  - NK: added fratio for interception hight optimization
!     correct the interception capacity by the multiplier fratio
      fratioflg=.false.  ! i.e. it's equal to 1.0
      do i=1,12
        do ii=1,classcount
          h(i,ii)=h(i,ii)*fratio(ii)
          if(fratio(ii).gt.1.00001.or.fratio(ii).lt.0.99999)
     *                                     fratioflg=.true.
        end do
      end do
      
!     rev. 10.2.55 June  09/19  - NK Added read_swe_update.f90 for read in swe adjustment factors
      swe_update=.false.
      inquire(file='snow1\swe_update.xml',EXIST=exists)
      if(exists)then
          call read_swe_update(yy_swe,mm_swe,dd_swe,hh_swe,
     *                                         swe_add,swe_mult)
          if(swe_add.gt.-999.or.swe_mult.gt.-999.0)swe_update=.true.
          if(swe_add.lt.-998.and.swe_mult.lt.-998.0)then
              write(98,*)
     *         'Warning: Both the swe adjustment factors are < -999'
              swe_update=.false.
          endif
          if(swe_update)then
              write(98,*)'INFO: swe update'
              write(98,*)'INFO: swe update for ',yy_swe,mm_swe,dd_swe,hh_swe
              write(98,*)'INFO: swe_add/swe_mult',swe_add,swe_mult
              write(98,*)
          endif
      endif
      
!     rev. 10.2.62 Sep.  09/19  - NK Added read_sm_update.f90 for read in sm adjustment factors
      uzs_update=.false.
      inquire(file='moist\uzs_update.xml',EXIST=exists)
      if(exists)then
          call read_uzs_update(yy_uzs,mm_uzs,dd_uzs,hh_uzs,
     *                                         uzs_add,uzs_mult)
          if(uzs_add.gt.-999.or.uzs_mult.gt.-999.0)uzs_update=.true.
          if(uzs_add.lt.-998.and.uzs_mult.lt.-998.0)then
              write(98,*)
     *         'Warning: Both the uzs adjustment factors are < -999'
              uzs_update=.false.
          endif
          if(uzs_update)then
              write(98,*)'INFO:uzs update for ',
     *                        yy_uzs,mm_uzs,dd_uzs,hh_uzs
              write(98,*)'INFO:uzs_add/uzs_mult',uzs_add,uzs_mult
          endif
      endif

      if(iopt.eq.2)print*,'Before event loop start'

!     * * * * * * * *  EVENT LOOP START  * * * * * * * * * * * * * * *
!     * * * * * * * *  EVENT LOOP START  * * * * * * * * * * * * * * *
!     * * * * * * * *  EVENT LOOP START  * * * * * * * * * * * * * * *
!     * * * * * * * *  EVENT LOOP START  * * * * * * * * * * * * * * *
!     * * * * * * * *  EVENT LOOP START  * * * * * * * * * * * * * * *
!     * * * * * * * *  EVENT LOOP START  * * * * * * * * * * * * * * *
!     * * * * * * * *  EVENT LOOP START  * * * * * * * * * * * * * * *
!     * * * * * * * *  EVENT LOOP START  * * * * * * * * * * * * * * *


      do id=1,ni
     
!         TIMER SETS ALL THE CLOCKS i.e. SUM HOURS, SUM SECONDS, ETC.
!     rev. 10.1.74 Apr.  01/17  - NK: Changed timer to fix 1 day-off problem 
!       This hasd to be moved here so the read_rain, read_temp, read_r2c etc files
!       have a proper initisl clock time        
        time=0.0
        jz=0 
        if(netCDFflg)jz=1         
c        if(netCDFflg)then
c          nnprint=s(ycount/2,xcount/2)         
c          ipr=yyy(nnprint)   
c          jpr=xxx(nnprint)
c          iopt=1
c          iopt99=.true.
c        endif
c        if(.not.iopt99)ensimflg='n'
        
          
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\      
!\\\\\added by Dave Newson
!      rev. 9.9.72  Jul.  21/15  - NK: Dave Newson additions to sub & process_rain
!       DN addition of forecast mode code
        fcst_mode=.false.
        !print *, '*** looking for fcst_params.txt in . ...'
        inquire(FILE="fcst_params.txt",EXIST=exists)
        !inquire(FILE='fcst_params.txt',EXIST=exists)
        
	  !print *,"cwd=",getcwd()
	  
        IF(exists)THEN
          print *, 'FOUND fcst_params.txt! checking forecast year...'
          open(unit=99,file='fcst_params.txt',iostat=ios)    ! 99 is a scratch file #
!           -- unit number selected randomly (fix!);
!           -- can't figure out filename convention
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file  fcst_params.txt'
            print*,'Possible cause(s):'
            print*,'file in use by another application'
            stop 'Program aborted in sub.f @ 358'
          endif
          read(99,*)fcst_yr ! read each line and assign to forecast parameter variables
          read(99,*)fcst_hr0
          read(99,*)fcst_days2avg
          read(99,*)fcst_snow_adj
          read(99,*)fcst_rain_adj
          close(unit=99,status='keep')
!         check if this is the current year
          print *, 'forecast year, current year: ', fcst_yr, year1
          if(fcst_yr.eq.year1)then          ! Trish believes 'year1' is the current year; might be 'year' or 'yearnow'?
!                                                look in timer.f
          fcst_mode=.true.
!         DN 2015-02-07 -- these are temporary debug statements
          !print *, 'in sub at 372, forecast parameters:'
          !print *, fcst_yr,fcst_hr0,fcst_snow_adj,fcst_rain_adj
          print *, 'forecast mode enabled?: ', fcst_mode
          end if
        ELSE 
          if(dds_flag.eq.0)then
            print*
            print*, '!WARNING fcst_params.txt not found!'
            print*
          endif
        END IF
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\      
!\\\\\end Dave Newson addition

d        if(iopt.eq.2)print*,' In sub, passed location  201'
        julast=0
!     rev. 9.9.49  Jan.  06/14  - NK: Added courantflg
        courantflg=.true.    ! All's well Courant criteria not violated
!       
CHECK FOR STOP COMMAND:
        open(unit=99,file='stop.txt',form='formatted',
     *          status='unknown',iostat=ios)
        if(ios.ne.0)print*,'Problems opening stop.txt file ignored'
        read(99,99001,iostat=ios)qwert
        if(ios.ne.0)print*,'Problems reading stop.txt file ignored'
        close(unit=99,status='keep')

d       if(iopt.eq.2)print*,' In sub, passed location  2011'
        if(qwert.gt.0.0)go to 83

!       RESET TO THE ORIGINAL VALUE -  WILL BE CHANGED IF NO DATA
        flgevp2=flgevp22

!     rev. 10.1.79 Apr.  18/17  - NK: Set trcflg=0 for all dds except errflg=10
        if(numa.ne.0.or.dds_flag.ne.0)then
          iopt=0
          trcflg='n'
          ensimflg='n'
          initflg='n'
          modelflg='n'
!         but GW flows are needd to run this dds criteria          
          if(errflg.eq.10)trcflg='y'
        endif

        index=1
!       INDEX = 1 FOR FIRST PASS REROUT ONLY
!       INDEX = 2 FOR SUBSEQUENT PASSES
!       SET   = 1 FOR EACH NEW LINKED EVENT TO READ IN NEXT SET OF 
!              RESERVOIR RELEASES

!     rev. 9.1.59  Jul.  15/04  - NK: split rerout into two parts: rdresv & rerout

d       if(iopt.eq.2)print*,' In sub, passed location  2012'

!       REV. 8.91 - Dec.  07/98 - READ rdevt IN SUB AS WELL AS SPL 
!       READ THE EVENT FILE - HAS TO BE DONE FOR EACH ID
        if(id.le.1)then
!         event file used to be read here as well but no longer needed
!         it is read in spl9 at the beginning
!         except for optimization when it has to be re-read
!         for each new iteration           apr. 30/08 nk
c          if(numa.gt.0)then
c            fln(99)='event\event.evt'
cd           print*,'reading event file  id = 1'    
c            call read_evt(date,conv,scale,smc5,nhr,nhf)
c          endif
c
c          month1=mo1

        else

!     rev. 9.5.48  Dec.  26/08  - NK: added event_fln() to allow unlimited events
          fln(99)=event_fln(id)

!         rdevt WAS CALLED IN SPL9 FOR THE FIRST EVENT
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!          call rdevt(date,conv,scale,smc5,nhr,nhf)
d         print*,'reading event file id > 1'      
          call read_evt(date,conv,scale,smc5,nhr,nhf)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        endif

!       needed of precip fitting: scale is set in options
        if(icase.eq.2)scale=a(1)

        if(numa.eq.0.or.id.eq.1.and.nnn.eq.0)then
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c          if(abs(dds_flag).ne.1)call header()
          call header()
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(abs(dds_flag).eq.1)ensimflg='n'

        endif

!     Changed the order: has to be after opening the next event!!
!     Nov. 24/02 nk

!     make sure we don't write the wfo file when optimizing
!     has to be her because flag is read with each new event
        if(ensimflg.eq.'y'.and.wfo_open_flg.eq.'n')then
!         read the wfo_spec file & do allocations
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d         print*,'reading wfo_spec file'
          call rd_wfo_spec()
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> AB:  ENSIM HEADERS
!         WRITE THE HEADER FOR ENSIM FILES:

!         NK - ALLOCATIONS FOR AREA2 arrays
!         rev  9.1.29  Oct.  24/02    - Added q1, qint & drng to wfo file

!         added nj to the arg lst to allocate attname ...nk  30/01/03
          if(llflg.ne.'y')then
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call write_both_headers('UTM       ',jan)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          else
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call write_both_headers('LATLONG ',jan)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          endif

          wfo_open_flg='y'

!         NOTE: IF YOU CHANGE WHICH DATA YOU OUTPUT YOU WILL HAVE TO EDIT
!           THE ABOVE SUBROUTINE
!           YOU CAN REPLACE 'UTM' WITH 'LATLONG' WHERE NECESSARY
!         ALSO, MAKE THE ARRAY FOR OUTPUTTING
!          xcount     ! number of columns
!          ycount     ! number of rows
           if(iopt.eq.2)print*,' In sub before initializing outwfo()'
          do j=1,xcount
            do i=1,ycount
              outwfo(j,i)=0.0
            end do
          end do
          do n=1,na
            wfo_sum_p(n)=0.0
          end do
          wfo_seq=0
!         LUN IS THE UNIT NUMBER for the watflood.wfo file
          lun=65
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>....>!
        endif ! if(ensimflg.eq......

!       rev. 9.1.28  Sept. 19/02  - Added shdlfg to replace the bsnm_shd.r2c file
        if(shdflg.eq.'y')then      
!         basin/bsnm.shd
          open(unit=31 ,file=fln(1) ,status='old',iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file',fln(1)(1:40)
            print*,'Possible cause(s):'
            print*,'file in use by another application'
            print*,'or target directory does not exist'
            print*,'or wrong # of events listed in event file'
            stop 'Program aborted in sub.f @ 896'
          endif

          if(IsFileTypeR2C(fln(1))) then
d             print*,'reading shd file'      
              call read_shed_ef(31,1) !EnSim compatible r2c file
c	!       **********************************************************************
c	        call read_par_parser(32,2)
c!             **********************************************************************

          else
!             call rdshed()
               print*,'Old format shd files not accepted'
               print*,'Please create EF ????_shd.r2c files & rerun'
               stop 'Program aborted in sub @ 528'
          endif

          write(51,5101)
          write(51,5102)id,fln(1) 
          print*,' New watershed file ',fln(1)
          print*,' read in'
          print*
!          pause 'In sub - read in new shd file'
          close(unit= 31)

        endif     !  if(shdflg.eq......

!     rev. 9.8.93  Nov.  12/13  - NK: Added the routing initialization with yyyymmdd_fli.r2c
        if(fliflg.eq.'y')then
          flnNum=99
          fln(99)=fln(55)
!         There are 2 flowinit possibilities:
!         One with a resume: the flow_init.r2c file will be in the working directory
!         the other when is it done on the fly as an updating tool with
!                    strfw\yyyymmdd_fli.r2c
!             This is on the fly with the file name from the event file
!             ~~~~~~~~~~~~~~~~~~~~~
              call read_flowinit(flnNum)
!             ~~~~~~~~~~~~~~~~~~~~~
        endif
d       if(iopt.eq.2)print*,' In sub, passed location  202'

!     rev. 9.1.52  Mar.  11/04  - NK: continuous water quality modelling
!     wqual/yymmdd.wqd
        if(sedflg.eq.'y')then
          close(unit=256,status='keep',iostat=ios)
          if(ios.ne.0)then
          print*,'Problem closing unit 256 fln=',fln(26)
          print*
          stop ' program aborted in sub @ 645'
          endif
        endif
d       if(iopt.eq.2)print*,' In sub, passed location  2026a'

!       REV 7.9  ADDED GRIDDED RADIATION FILE
        if(ver.ge.7.9)then
          if(flgevp2.eq.3.0)then 
d         if(iopt.eq.2)print*,' In sub, passed location  2044'
          close(unit=49,status='keep')
d         if(iopt.eq.2)print*,' In sub, passed location  2045'
          close(unit=50,status='keep')
          endif
        endif

!       if(resinflg.eq.'y') close(unit=99,status='keep')

d      if(iopt.eq.2)print*,' In sub, passed location  203'

!       INPUT FILES

!       unit    32 = fln(2)  - parameter data file .par
!       unit    36 = fln(6)  - flow file .str
!       unit    37 = fln(7)  - reservoir releases file .rel
!       unit    38 = fln(8)  - snow file .snw
!       unit    40 = fln(10) - precip file .met

!       AND OPEN THE FILES FOR THE NEXT ONE:

d        if(iopt.eq.2)print*,' In sub, passed location  205'

!       rev. 9.1.52  Mar.  11/04  - NK: continuous water quality modelling
!       wqual/yymmdd.wqd
        if(sedflg.eq.'y')then
          open(unit=256,file=fln(40),status='unknown',iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file',fln(40)(1:40)
            print*,'Possible cause(s):'
            print*,'file in use by another application'
            print*,'or target directory does not exist'
            print*
            print*,'This file is optional - used for nutrient loading'
            print*,'needed only if sedflg=y in yymmdd.evt'
            stop 'Program aborted in sub.f @ 990'
          endif
        endif

!       unit 261 = fln(31) - gridded runoff files for watroute .rff
!       unit 262 = fln(32) - gridded recharge for modflow .rch
!       unit 263 = fln(33) - gridded leakage for watroute

!       rev. 9.1.45  Jun.  11/03  - WATROUTE: runoff, recharge and leakage files added 
!       WRITING FILES FOR WATROUTE
        if(modelflg.eq.'l')then
!         route surface  flow & leakage
!         these 2 are added together and then routed
!         open runof\yyyymmdd_rff.r2c
          open(unit=261,file=fln(31),status='unknown',iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file',fln(31)(1:40)
            print*,'Possible cause(s):'
            print*,'file in use by another application'
            print*,'or target directory does not exist'
            stop 'Program aborted in sub.f @ 1011'
          endif
!         open lkage\yyyymmdd_lkg.r2c
          open(unit=263,file=fln(33),status='unknown',iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file',fln(33)(1:40)
            print*,'Possible cause(s):'
            print*,'file in use by another application'
            print*,'or target directory does not exist'
            stop 'Program aborted in sub.f @ 1021'
          endif
        endif

        if(modelflg.eq.'r')then
!         route surface flow to stream and rchrg thru the lzs
!         leakage is computed using the watflood LZ routing module
!         open runof\yyyymmdd_rff.r2c
          open(unit=261,file=fln(31),status='unknown',iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file',fln(31)(1:40)
            print*,'Possible cause(s):'
            print*,'file in use by another application'
            print*,'or target directory does not exist'
            stop 'Program aborted in sub.f @ 1036'
          endif
!         open rchrg\yyyymmdd_rch.r2c
          open(unit=262,file=fln(32),status='unknown',iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file',fln(32)(1:40)
            print*,'Possible cause(s):'
            print*,'file in use by another application'
            print*,'or target directory does not exist'
            stop 'Program aborted in sub.f @ 1046'
          endif
        endif

        if(modelflg.eq.'i')then
!         surface flow only to the stream for routing
!         leakage is not included in the routing
!         open runof\yyyymmdd_rff.r2c
          open(unit=261,file=fln(31),status='unknown',iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file',fln(31)(1:40)
            print*,'Possible cause(s):'
            print*,'file in use by another application'
            print*,'or target directory does not exist'
            stop 'Program aborted in sub.f @ 1061'
          endif
        endif

!       WATROUTE WATROUTE WATROUTE WATROUTE WATROUTE WATROUTE WATROUTE WATROUTE

!       rev. 9.1.45  Jun.  11/03  - WATROUTE: runoff, recharge and leakage files added 
!       Writing FILES FOR WATROUTE
        if(routeflg.eq.'y'.and.numa.eq.0)then
          mhtot=2   ! <<<<<<<<<<<< just to get the right header - reset later
          author='spl.exe                     '    
          name='Gridded Channel Inflow              '
          coordsys_temp=coordsys1
!         GreenKenue uses LatLong - code below uses LATLONG
          if(coordsys_temp.eq.'LatLong   ')coordsys_temp='LATLONG   '
          if(coordsys_temp.eq.'Cartesian ')coordsys_temp='CARTESIAN '
          zone_temp=zone1
          datum_temp=datum1
          xorigin_temp=xorigin
          yorigin_temp=yorigin
          xcount_temp=xcount
          ycount_temp=ycount
          xdelta_temp=xdelta
          ydelta_temp=ydelta
          attribute_name='channel_inflow                  '
          attribute_units='mm                          ' 
          unit_conversion=1.0
          attribute_type='Runoff                      '  
          source_file_name=fln(10)
          frame_no3=0            ! write the headers               
!     rev. 9.9.23  Aug.  10/14  - NK: Uniform arguments for write_r2c
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c         call write_r2c(261,31,mhtot,1,frame_no3,1,1)     
          call write_r2c(261,31,mhtot,0,frame_no3,0,1)     
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          name='Gridded Recharge                '
          attribute_name='recharge                        '
          attribute_units='mm                          ' 
          attribute_type='recharge                        '  
          source_file_name=fln(10)                           
!     rev. 9.9.23  Aug.  10/14  - NK: Uniform arguments for write_r2c
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c         call write_r2c(262,32,mhtot,1,frame_no3,1,1)     
          call write_r2c(262,32,mhtot,0,frame_no3,0,1)     
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          name='Gridded Leakage                 '
          attribute_name='leakage                     '
          attribute_units='mm                          ' 
          attribute_type='leakage                         '  
          source_file_name=fln(10)                           
!     rev. 9.9.23  Aug.  10/14  - NK: Uniform arguments for write_r2c
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c         call write_r2c(263,33,mhtot,1,frame_no3,1,1)     
          call write_r2c(263,33,mhtot,0,frame_no3,0,1)     
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        endif

d       if(iopt.eq.2)print*,' In sub, passed location  206'

        if(ver.ge.7.9)then
          if(flgevp2.eq.3.0)then
!         RADIATION FILE
          open(unit=49,file=fln(19),status='unknown',iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file',fln(19)(1:40)
            print*,'Possible cause(s):'
            print*,'file in use by another application'
            print*,'or target directory does not exist'
            stop 'Program aborted in sub.f @ 1131'
          endif
          open(unit=50,file=fln(20),status='unknown',iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file',fln(20)(1:40)
            print*,'Possible cause(s):'
            print*,'file in use by another application'
            print*,'or target directory does not exist'
            stop 'Program aborted in sub.f @ 1140'
          endif
          endif
        endif

d       if(iopt.eq.2)print*,' In sub, passed location  207'

!       RESET THE CLOCK:   >>>>>>>>>>>> CHECK THIS OUT
        m=1
        tot1=0.0

d       if(iopt.eq.2)print*,' In sub, passed location  209'

! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!       FROM THE STREAMFLOW FILE:

!       REV. 7.32   feb.  07/95 - ADDED NOPT TO SELECT OPT FLOW STA 
!       NOPT IS A VECTOR SIZE = # OF STREAMFLOW STATIONS
!       IF = 0 STATION NOT USED FOR OPTIMIZATION
!       IF = 1 STATION USED FOR OPTIMIZATION
!       if = 2 station is used to "correct" flows for that event only
!       if = 3 in 1st event, station nudges in all events

!       rev. 9.1.23  Jul.  23/02  - Added control for nudging in event #1
        if(id.eq.1)then
!         the first event is the master event for settings    
!         if nopt=3
          do l=1,no
            nopt1(l)=nopt(l)
            if(nopt(l).eq.3)nopt(l)=2   ! set to nudge first event
          end do
        else
!         subsequent events are nudged if nopt1(?)=3
          do l=1,no
            if(nopt1(l).eq.3)then
              nopt(l)=2     ! set to nudge following events
            endif
          end do
        endif

d        if(iopt.eq.2)print*,' In sub, passed location  214'

        if(iopt.eq.99)then
!         THIS OPTION IS TO CHECK ALL INPUT FILES
!          kt=1
!          nl=1
          mhrd=kt
          if(id.le.1)then
          write(*,6300)iopt
          write(98,6300)iopt
          write(*,'(A)',advance='no') 
     *         'In sub: hit any key to continue checking files'
          read(*,*)
          endif
        endif

d        if(iopt.eq.2)print*,' In sub, passed location  722'

!       MHRD CAN BE ENTERED WHEN SAVING THE STR FILE.
!       MHRD DEFAULT VALUE IS THE END OF THE STREAMFLOW DATA.

!*******************
! RAIN (MET FILE)
!*******************
!//////////////////////////////////////////////
!///////////////////////// 
!// Added by Dave
!     rev. 9.5.30  May.  26/08  - NK: conv back in read_rain & process_rain arg. list

c        if(IsFileTypeR2C(fln(10))) then
!         changed argument list  nov. 9/06 nk
!         read the header:
d          if(iopt.eq.2)then
d			print*,'In sub, passed location  904 before read_rain'
d             print*,'If program dies here, check the met file'
d		 endif 

        if(modelflg.eq.'n')then
!       rev. 9.9.35  Oct.  20/14  - NK: Added keyword & file checks
          inquire(FILE=fln(10),EXIST=exists)
          if(.not.exists)then
            print*
            print*,'Fatal error'
            print*,'radcl\yyyymmdd_met.r2c file not found'
            print*,'Please check the event\event.evt file'
            print*,'for correct keyword nad/or file name'
            print*
            stop 'Program aborted in sub @ 1185'
          endif

          
          call find_filetype(10)
          if(filetype.eq.'r2c')then
!         read the header:
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c            call read_rain_ef('1',conv,jz,jan) !EnSim compatible r2c file
            call read_rain_ef('1',conv,0,jan) !EnSim compatible r2c file
c          elseif(filetype(1:3).eq.'bin')then
c            call read_rain_bin('1',conv,jz,jan) 
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          elseif(filetype.eq.'nc ')then
              continue
!             no header to read
          else
            print*,'Error:'
            print*,'File type for met file =',filetype
            print*,'Old format met files not accepted'
            print*,'Please create EF _met.r2c files & rerun'
            print*
            print*,'Forward slashes / in the event files not'
            print*,'accepted on the PC (it`s ok in unix)'
            print*
            stop 'Program aborted in sub @ 830'
          endif
        endif
        
!// End Dave addition
!/////////////////////////
!//////////////////////////////////////////////

d        if(iopt.eq.2)print*,' In sub, passed location  769'
d        write(98,*)'met file header read'

!*******************
! TEMP (TEM FILE)
!*******************




        if(modelflg.eq.'n')then
          
          
          
          

        if(snwflg.eq.'y'.or.vapflg.eq.'y')then    ! added nk Jun. 28/06
!//////////////////////////////////////////////
!///////////////////////// 
!// Added by Dave
d         if(iopt.eq.2)print*,' In sub, passed location  936'
!         rev. 9.9.35  Oct.  20/14  - NK: Added keyword & file checks
          inquire(FILE=fln(15),EXIST=exists)
          if(.not.exists)then
            print*
            print*,'Fatal error <<<<<'
            print*,fln(15)(1:72)
            print*,' not found'
            print*,'Please check the event\event.evt file'
            print*,'for correct keyword and/or file name'
            print*
!           THE TEMP FILE IS MISSING WHEN IT SHOULD BE THERE -> TERMINATE
            write(*,99111)
99111     format(' This file is optional - used to input gridded temps'/
     *  ' needed if snwflg=y or vapflg=y in yymmdd.evt'/
     *  ' not needed if flgevp2 .le. 1.0'/
     *  ' OR: in config.sys have you set files=100 & buffers=50?'/)
            stop 'Program aborted in sub @ 1259'
          endif

          if(IsFileTypeR2C(fln(15)).and.
     *                  .not.FLtype(15)(1:2).eq.'nc')then
d          if(iopt.eq.2)then
d			print*,'In sub, passed location  796 before read_temp'
d             print*,'If program dies here, check the tem file'
d		 endif 
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call read_temp_ef('1',jan,jz) !EnSim  r2c file
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d          if(iopt.eq.2)print*,' In sub, passed location  798'
          elseif(.not.FLtype(15)(1:2).eq.'nc')then   !  i.e. nc is ok so go on
!           call rdtemp('1',jan,jz)
            print*,'Old format tem files not accepted'
            print*,'Please create EF _tem.r2c files & rerun'
            print*
            stop 'Program aborted in runof6 @ 857'
          endif
!// End Dave addition
!/////////////////////////
!//////////////////////////////////////////////

!     rev. 9.9.06  Jan.08/14  - NK: Add daily differences to Harfreaves ETHarg.f
          if(flgevp2.eq.4)then
            INQUIRE(FILE=fln(62),EXIST=exists)
            if(exists)then
              dlyflg=.true.
!             read the header of tempr\yyyymmdd_dif.r2c
              if(iopt99)then
                print*
                print*,'NEW  <<<<<'
                print*,'Reading the header on unit 292 ',fln(62)(1:40)
              endif
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              call read_r2c(292,62,'1',jz,newDataFlag) !EnSim  r2c file
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if(iopt99)then
                print*,'Back from reading the header ',fln(62)(1:40)
                print*,'SPL will use the daily temperature differences'
                print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
                print*   
              endif           
            else
!             in this case, monthly_climate_normals or
!             mean_dly_diff will be used depending on availability
!     rev. 9.9.36  Nov.  03/14  - NK: Revised error message for daily diff choices
              print*
              print*,'WARNING'
              print*,'flgevp is read having a value of 4'
              print*,'expecting file: ',fln(62)(1:40)
              print*,'but is is not found'
              print*,'Are you using a current evt file?????'
              print*,'If the temp data was not produced by tmp.exe,'
              print*,'the *dif.r2c file was not created.' 
              print*,'Please use diff.exe to create the *dif.r2c files.'

              stop 'Program aborted in sub @ 1455'

              print*,'You can continue using the values in the'
              print*,'basin\monthly_climate_normals.txt  file'
              print*,'for the Hargreaves & Samani (1985) equation'
              print*,'This is not ideal as daily diff`s are expected'
              print*,'You may continue by hitting Enter'
              print*,'if not running DDS'
              if(dds_flag.eq.1)then
                pause 'program will abort'
                stop  'Program aborted in sub @ 1291'
              else    
                dlyflg=.false.
                print*,'program continues with monthly_climate_normals'
                print*,'if the file is found in the basin directory'
                if(id.eq.1)pause  'Hit enter to continue'
              endif
            endif
          endif
        endif  ! added nk Jun. 28/06
        endif    !modelflg='n'

!       ***********************
!       EVAPORATION (EVP FILE) for lakes
!       ***********************

C//////////////////////////////////////////////
C///////////////////////// 
!     rev. 9.5.12  Feb.  13/08  - NK: added evaporation input file with read_r2c
d          if(iopt.eq.2)print*,' In sub, passed location  922'
        inquire(FILE=fln(51),EXIST=exists)  !gridded evaporation
        if(exists)then
          rd_evp_flg=.true.
!         if the evp file exists, it will be used for lake evaporation. 
!         there's no way not to use it if it is not there
!         it also means that we can switch during a run

          if(IsFileTypeR2C(fln(51))) then
d           if(iopt.eq.2)then
d			  print*,'In sub, passed location  929 before read_r2c'
d             print*,'If program dies here, check the evaporation file'
d	 	    endif 
!           read the header
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d           print*,'reading evaporaton r2c file'     
            call read_r2c(281,51,'1',jz,newDataFlag) !EnSim  r2c file
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if(id.eq.1)then
              if(.NOT.allocated(evap_convert))then      
              allocate(evap_convert(na),lake_evaporation(na),
     *                stat=iAllocate)
              if(iAllocate.ne.0) STOP
     *           'Error with allocation of evap_convert in sub @ 941'
                print*
              endif
!             conversion factor for m evaporation on a grid to m**3
              do n=1,naa
              evap_convert(n)=grid_area(n)*aclass(n,classcount)/1000.0
                end do
              print*,'read header fln='
              print*,fln(51)
              print*,'Time step other than 1hr to be fixed'
              print*

            endif
           
d            if(iopt.eq.2)print*,' In sub, passed location  929'
          else
            print*
            print*,'Fatal error:'
            print*,'r2c file expected'
            print*,'Wrong file type found for file fln(51)'
            print*,fln(51)(1:60)
            print*,'Delete this file if not needed or just junk'
            print*,'and try again'
            print*
            stop 'Program aborted in sub @ 942'
          endif
        else
          rd_evp_flg=.false.
        endif

!     rev. 9.5.03  Dec.  09/07  - NK: added reads for precip isotopes

        if(frcflg.eq.'y') then

!         Check for REMIiso input files. If all of the 4 required files are present, use them, otherwise check for other inputs.            
          frc_file_flg='y'
            inquire(FILE=fln(21),EXIST=exists)
          if(.not.exists)frc_file_flg='n'
            inquire(FILE=fln(47),EXIST=exists)
          if(.not.exists)frc_file_flg='n'
            inquire(FILE=fln(48),EXIST=exists)
          if(.not.exists)frc_file_flg='n'
            inquire(FILE=fln(49),EXIST=exists)
          if(.not.exists)frc_file_flg='n'

!         If all REMOiso inputs are present, then skip the next section of code (because frac_file_flg='y' will be written over).         
         if(frc_file_flg.eq.'y') GOTO 1221 
         
!         Check if the time series delta rain input is present.            
          frc_file_flg='c'  
          inquire(FILE=fln(48),EXIST=exists)
          if(.not.exists)then
            frc_file_flg='n'
            if(dds_flag.eq.0)then
              print*,'WARNING:'
              print*,'Looking for :',fln(48)(1:60)
              print*,'but file not found'
              print*,'Program continues without time series dRain data'
            endif
          else
            if(dds_flag.eq.0)then
              print*
              print*,'frc_file_flg =',frc_file_flg
              print*,'Program continues with time series dRain data'
            endif
          endif
          print*

!     rev. 10.1.01 Oct.  05/15  - NK: Isotope update: added 2H
c          this section removed by TH:
c!         if one of the needed REMOiso files is missing, or if there is no time series delta rain, use the defaults (values from isotope.par file).
c          print*,'frc_file_flg =',frc_file_flg
c          if(frc_file_flg.eq.'n')then
c            print*,'do you want to continue to run isoWATFLOOD without
c     *     time series input data?'
c            print*
c           pause 'enter to continue; cntl C to abort'
c          endif

 1221     CONTINUE
          
          if(frc_file_flg.eq.'y')then
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d         print*,'reading humidity file'      
c          call read_hum(251,21,'1',jan,jz)   ! humidity  removed by TH:
d         print*,'reading gridded snow precip'      
          call read_gsn(277,47,'1',jan,jz)   ! snow_precip
d         print*,'reading REMOiso delta rain'      
          call read_drn(278,48,'1',jan,jz)   ! delta rain
d         print*,'reading REMOiso delta snow'      
          call read_dsn(279,49,'1',jan,jz)   ! delta snow
          
          elseif(frc_file_flg.eq.'c')then
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         print*    
         print*,'reading time series delta rain'      
          call read_r2c(278,48,'1',jz,newDataFlag)   ! delta rain
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d      print*,'after reading ddrain time=',time

          endif 
          
          !TH: relative humidity from EC data based r2c
          inquire(FILE=fln(21),EXIST=exists)
          if(exists)then
            call read_r2c(251,21,'1',jz,newDataFlag)
            RH_flg=.true.
          else
            RH_flg=.false.
          endif

!         read in isotope parameter file from \basin
          inquire(FILE='basin\isotope.init',EXIST=exists)
          if(exists)then
          open(99,file='basin\isotope.init',iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file  basin\isotope.init'
            print*,'Possible cause(s):'
            print*,'file in use by another application'
            print*,'or target directory does not exist'
            stop 'Program aborted in sub.f @ 1593'
          endif
          
          read(99,9990)flg2H        ! deutirium flag 1=18O, 2=18O and 2H
          read(99,9990)ninit    ! number of initializations
          read(99,9990)h2oflg    ! classcount for isotope output

          if(.NOT.allocated(deltar))then
!         ALLOCATE ISOTOPE INPUT PARAMETERS
          allocate(isoyear(ninit),deltar(ninit),deltas(ninit),
     *    rfoffset(ninit),smoffset(ninit),
     *    deltar2H(ninit),deltas2H(ninit),
     *    rfoffset2H(ninit),smoffset2H(ninit),stat=iall)
          if(iall.ne.0) STOP 'Error allocating isotope inputs'
          
          
          read(99,9991)delta1                    ! riverwater
          read(99,9991)delta2                    ! soil water
          read(99,9991)delta3                    ! groundwater
          read(99,9991)delta4                    ! snow
          read(99,9992)(isoyear(ii),ii=1,ninit)  ! data year 
          read(99,9993)(deltar(ii),ii=1,ninit)   ! rain
          read(99,9993)(deltas(ii),ii=1,ninit)   ! snow
          read(99,9993)(rfoffset(ii),ii=1,ninit)  ! per mil depletion for refreezing
          read(99,9993)(smoffset(ii),ii=1,ninit)  ! per mil enrichment for snowmelt
          if(flg2H==2)then
           read(99,9991)delta2H1                    ! riverwater 2H
           read(99,9991)delta2H2                    ! soil water 2H
           read(99,9991)delta2H3                    ! groundwater 2H
           read(99,9991)delta2H4                    ! snow 2H
           read(99,9993)(deltar2H(ii),ii=1,ninit)   ! 2H rain
           read(99,9993)(deltas2H(ii),ii=1,ninit)   ! 2H snow
           read(99,9993)(rfoffset2H(ii),ii=1,ninit)  ! 2H per mil depletion for refreezing
           read(99,9993)(smoffset2H(ii),ii=1,ninit)  ! 2H per mil enrichment for snowmelt
           read(99,9990)isoframeflg                 ! This doohickey is to set up the framework calcs: 1: whole basin avg framework, 2: 1 framework per gauge, 3: multiple frameworks, but with gauges aggregated-> read another file

           if(isoframeflg.eq.1)then 
           nisoframe=1
           else if(isoframeflg.eq.3)then
            
            inquire(FILE='basin\isobasin_combine.txt',EXIST=exists)
            if(exists)then
             open(98,file='basin\isobasin_combine.txt',iostat=ios)
             if(ios.ne.0)    
     *         print*,'Unable to open file  basin\isobasin_combine.txt'
             allocate (isoframecom(no),stat=iall)
             if(iall.ne.0) STOP 'Error allocating isotope frame combine'
             nisoframe=1
             do n=1,no
              read(98,9994)isoframecom(n)
              if(isoframecom(n).gt.nisoframe)nisoframe=isoframecom(n)
             end do
             close(98)
            else
             STOP 
     *        'Error reading \basin\isobasin_combine.txt 
     *          - create or do not use option 3'
            endif
            
           else 
           nisoframe=no
           end if
           
          end if
          close(99)
          endif

!     REV. 10.1.21 Jan.  22/16  - NK: isotope updates
          !TH: NEW ET split pars
          inquire(FILE='basin\ET.par',EXIST=exists)
            if(exists)then
             open(98,file='basin\ET.par',iostat=ios)
             if(ios.ne.0)    
     *         print*,'Unable to open file  basin\ET.par'
             do ii=1,classcount
              read(98,9995)acg(ii)
              read(98,9995)bcg(ii)
             end do
             close(98)
            else
             do ii=1,classcount
              acg(ii)=0.3
              bcg(ii)=4.
             end do
          end if
          
          else
          STOP 
     *     'Error reading \basin\isotope.init - create file & re-run'
          
          endif  !if par exists

!     REV. 10.1.21 Jan.  22/16  - NK: isotope updates
          !TH: if 18O is coming from gridded frc file, 2H should be checked too
           !TH: added 2H flg to allow for 18O runs only with CD gridded data
          if(frc_file_flg.eq.'c'.and.flg2H==2)then
           inquire(FILE=fln(49),EXIST=exists)
              if(.not.exists)STOP 
     *     'Error reading time series 2H data. Either added gridded 2H, 
     *turn off 2H flg, or turn off gridded 18O'
           print*,'reading time series delta 2H rain'
           call read_r2c(279,49,'1',jz,newDataFlag) !TH: I stole the gridded snow file for 2H (hopefully temp)
          endif
                   
        endif

d        if(iopt.eq.2)print*, '4 in rain'
d        write(98,*)'tem file header read'

!     rev. 9.1.75  Feb.  08/05  - NK: added rdgsm (gridded soil moisture)
!     rev. 9.2.25  Dec.  13/05  - NK: ENSIM r2c gridded soil moisture 

        if(id.eq.1)then 
          INQUIRE(FILE=fln(37),EXIST=exists)  ! gridded initial soil moisture
          IF(exists)THEN
            if(IsFileTypeR2C(fln(37)))then 
!             note: gsm is read again later to overwrite soil_init.r2c
!                 when present
d             if(iopt.eq.2)then
d		  	  print*,'In sub, passed location 1083 before read_gsm'
d               print*,'If program dies here, check the gsm file'
d		    endif 
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d             print*,'reading gridded soil moisture'      
              call read_gsm_ef    !EnSim compatible r2c file
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            else
              print*,'Old format gsm files not accepted'
              print*,'Please create EF _gsm.r2c files & rerun'
              print*
              stop 'Program aborted in sub @ 841'
            endif
d            write(98,*)'gridded soil moistures read'
          else
c            if(resumflg.eq.'y')then
            if(resumflg.ne.'n')then  ! could be y of s
c             print*,'Soil moisture from soil_init.r2c file used'
            else
              print*
              print*,'WARNING:'
              print*, fln(37)
              print*, 'not found'
              print*
              print*,'Please see manual on how to create a '
              print*,'gridded soil moisture file'
              print*
              pause 'Hit enter to use the event file values'
            endif
            ssmc_firstpass='n'    
            if(.NOT.allocated(ssmc))then      
!             suggested modification by D. Watson  mar 21/06 
              allocate(ssmc(ycount,xcount),stat=iAllocate)
              if(iAllocate.ne.0) STOP
     *         'Error with allocation of ssmc in sub @ 919'
              if(iopt.eq.2)print*,
     *          'allocation done in sub',ycount,xcount
              print*
            endif
            if(.NOT.allocated(api))then  
!             suggested modification by D. Watson  mar 21/06 
              allocate(api(na,classcount),stat=iAllocate)
              if(iAllocate.ne.0) STOP
     *         'Error with allocation of api in sub @ 927'
              if(iopt.eq.2)print*,
     *         'allocation done in sub',ycount,xcount
              print*
            endif
          endif
        endif

d        if(iopt.eq.2)print*,' In sub, passed location  934'

!       THIS SIMPLIFIES THE RAINFALL INPUT SO AN ENTIRE RAINFALL 
!       SEQUENCE CAN BE INPUTTED AT ONCE BUT ONLY THE DATA FOR THE 
!       CALIBRATION PERIOD WILL BE USED FOR THE FORECAST.

d          if(iopt.eq.2)print*,' In sub - gone to read_flow_ef'

           if(IsFileTypeTB0(fln(6))) then
c!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             call read_flow_ef('1',date)  !EnSim compatible tb0 file
!            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           endif
            
d          if(iopt.eq.2)print*,' In sub - back from read_flow_ef'
d          write(98,*)'flows read'
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c        endif

!     rev. 9.1.59  Jul.  15/04  - NK: split rerout into two parts: rdresv & rerout

d          if(iopt.eq.2)print*,' In sub - gone to read_resv_ef'
c          if(IsFileTypeTB0(fln(7))) then

!     rev. 9.5.57  Apr.  13/09  - NK: added ntrlflg for natural lake flows
!     rev. 9.9.40  Nov.  19/14  - NK: Modified the 'a' option for ntrlflg
c           if(id.eq.1.and.ntrlflg.eq.'y'.or.ntrlflg.ne.'y')then
           if(id.eq.1.and.ntrlflg.eq.'a')then
!             read the rel file - running natural flows
!             For natural flows, read the rel file only once for 1st event
!             then keep the same coefficients throughout the run
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             call read_resv_ef()  !EnSim compatible tb0 file
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           elseif(ntrlflg.eq.'a')then
!            in this case we are using the coefficients from event #1           
             continue
           else
!            for ntrlflg y or n in any event file except if ntrlflg = a in event #1  
!            read the release file           
!            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             call read_resv_ef()  !EnSim compatible tb0 file
!            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            endif
            if(id.eq.1.and.noresv.ge.1)then
!             no need to do this if there are no reservoirs or lakes

!             fist event - read the rel file no matter what.
!             to run natural flows, make sure all locations have coefficients. 
!             you can rename the rel file resrl\yyyymmdd_ntrl.tb0
!             write an xyz file flow reservoir location plotting in GK
!             Added this here for FEWS  NK Jul.                
              if(iopt99.or.FLtype(6).eq.'nc')then
                if(FLtype(6).eq.'nc')then
                  open(unit=99,file='debug\reservoir_location.xyz',
     *                                           status='unknown')
                else
                  open(unit=99,file='reservoir_location.xyz',
     *                                           status='unknown')
                endif
                if(ios.ne.0)then    ! added Nov. 10/14  nk
                  print*
                  print*,'Unable to open file  reservoir_location.xyz'
                  print*,'Possible cause(s):'
                  print*,'file in use by another application'
                  stop 'Program aborted in sub.f @ 1637'
                endif
                do i=1,noresv
                  if(b1(i).gt.0.0)then
                    outlet_type='weir       '
                  else
                    outlet_type='regulated   '
                  endif
!     rev. 10.1.87 May   18/17  - NK: Added DA to reservoir_location.xyz
!     rev. 10.2.18 Mar.  12/18  - NK: Fixed array fault in read_resv_ef and sub
                  if(ires(i).le.ycount.and.jres(i).le.xcount.and.
     *                        ires(i).gt.0.and.jres(i).gt.0)then
                    n=s(ires(i),jres(i))
                    if(n.gt.0)write(99,99002)xres(i),yres(i),i,
     *                           resname(i),outlet_type,da(n)
                  endif
                end do
                close (unit=99,status='keep')
              endif   ! iopt99
              
              
!     rev. 9.5.59  Jul.  26/09  - NK: added fpet_lake for each lake in ill file
!             initialize    
              do l=1,noresv
                b6(l)=100.0     ! changed from 0.0  Dec. 3/14 nk
                b7(l)=100.0
!     rev. 9.9.38  Nov.  12/14  - NK: Added LKdepth to ill file
                LKdepth(l)=1.0  !default value if not assigned in yyyymmdd_ill.pt2
                fpet_lake(l)=-1.0
              end do 
              

!     rev. 9.5.51  Jan.  13/09  - NK: added reading yyyymmdd_ill.pt2 foa all lakes
!     rev. 9.5.59  Jul.  26/09  - NK: added fpet_lake for each lake in ill file
!             read the initial lake elevations for all lakes:
!             read the initial lake levels for all lakes:
!             test if yyyymmdd_ill.pt2 file exists. if yes, read it
!               Initial lake levels
!               read ill.pt2 file               ill   ill   ill   ill   ill   ill   ill
                INQUIRE(FILE=fln(50),EXIST=exists)
                if(exists)then
!                 read the initial lake level & datum  yyyymmdd_ill.pt2                 
!                 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  call read_pt2(280,50,nrows,ncols)
!                 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  print*
                  if(nrows.ne.noresv)then
                    write(98,*)
     *              'Warning: Number of initial lake levels does not ',
     *               'match the number of lakes & reservoirs'
c                    print*,'Check for extra blank lines in the file'
                    write(98,*)'No of reservoirs & lakes =',noresv
                    write(98,*)'No init lake levels =',nrows
                    write(98,*)'Warning: program aborted!'
                    if(dds_flag.eq.1)pause 'Pause in sub - please fix'
                    stop 'Aborted in sub @ 1347 -SEE debug\warnings.txt'
                  endif
                  write(53,*)'in sub - init lake levels & datums'
                  print*,'Finished reading ',fln(50)(1:40),nrows,ncols      

!                 Qmin set as a flag to check that qmin & safe_max
!                 are entered in the ill file
                  qmin(1)=-999.0  
                  do l=1,noresv
                    b6(l)=inarray(l,1)           ! init lake level
!     rev. 10.1.95 Sep   11/17  - NK: Fixed LKdepth bug in sub
                    lake_elv(l,1)=inarray(l,1)   ! to make sure it's assigned
                    b7(l)=inarray(l,2)           ! datum
c                    if(ncols.eq.4)then           ! 4 att's incl. name
                    if(ncols.ge.4)then
                        
                        ! 4 attributes incl. name
!                     default value = 1.0 set above                    
                      LKdepth(l)=inarray(l,3)
                      if(LKdepth(l).lt.b6(l)-b7(l))then
                          print*,
     *                 'WARNING: Lake Depth',l,' lower than Datum'
                      endif    
!                     LKdepth in the ill file is the initial lake depth below the invert
!                     like the hight of a weir - and is dead storage usually 
!                     But here it is redefined as tha actual depth of water in the lake
!                     and becomes a variable and is used in the evaporation routine
                      LKinvert(l)=b7(l)-inarray(l,3)
                      LKdepth(l)=lake_elv(l,1)-LKinvert(l)

                      write(53,*)'l=',l,'level=',b6(l),'datum=',b7(l),
     *                       'LKdepth=',LKdepth(l),'LKinve=',LKinvert(l)
                    endif  
                    if(ncols.ge.6)then
                      safe_max(l)=inarray(l,4)
                      qmin(l)=inarray(l,5)
                      DecayT(l)=inarray(l,6)
                      write(53,*)'l=',l,'level=',b6(l),'datum=',b7(l),
     *                    'LKdepth=',LKdepth(l),'LKinvert=',LKinvert(l),
     *                      'Qmin=',qmin(l),'safe_max=',safe_max(l),
     *                      'DecayT=',DecayT(l)
                    endif
c                    fpet_lake(l)=inarray(l,3)    ! fpet for lake l
c                    if(iopt.eq.1)then
c                    endif                      
                  end do
                  deallocate(inarray)
!                 NOTE:
!                 initial lake storage & outflows are calculated in flowinit.f below
                  
!                 write a file that has the datum for lakes & reservoirs
!                 so the datum can be plotted
                  if(iopt99)then
                    open(unit=99,file='results\datum.txt',
     *                               status='unknown',iostat=ios)
                    write(99,99005)0,(b7(l),l=1,noresv)
                    do i=1,ni
                      write(99,99005)i*365,(b7(l),l=1,noresv)  
                    end do
                    close(unit=99,status='keep')
                  endif
                else
                  if(dds_flag.eq.0)then
                    write(53,*)
                    write(53,*)'WARNING'
                    write(53,*)'No initial lake levels file found.'
                    write(53,*)'Looking for file:' 
                    write(53,*)fln(50)(1:40)
                    write(53,*)
     *                'Default values for datum & initial level = 0.0'
                    write(53,*)
                    print*
                    print*,'WARNING'
                    print*,'No initial lake levels file found.'
                    print*,'Looking for file:' 
                    print*,fln(50)(1:40)
                    print*,
     *                'Default values for datum & initial level = 0.0'
                    print*
                    if(lakeflg.eq.'y')then 
                      print*,'NOTE: for lake evap model, ill file with'
                      print*,'lake depths are required'
                      stop 'program aborted in sub @ 1714'
                    endif
                  endif
                endif  ! (exist)

!     rev. 9.7.27  May.  26/11  - NK: Add lake_ice_factor
!     Commented out for MRBM in 2014
!     REV. 10.1.42 Oct   20/16  - NK: Reinstated read_ice_factor.f as default if present
              icefactorfile=.false.   ! change to true if file exists
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              call read_ice_factor()    !!for lakes only!!
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            endif    ! id=1

c          else
c!           need this just in case there is an old format rel file 15/04/08 nk
c!           call rdresv()   
c            inquire(FILE=fln(7),EXIST=exists)
c            if(exists)then
c              print*,'Old format .rel files not accepted'
c              print*,'Please create yyyymmdd_rel.tb0 files & rerun'
c              stop 'Program aborted in sub @ 1052'
c            endif
c          endif
d          if(iopt.eq.2)print*,' In sub - back from read_resv_ef'
d          write(98,*)'releases read'

!     rev. 9.8.23  Aug.  03/12  - NK: Added resinid1flg to use resinflg for id=1
          if(resinflg.eq.'y')then
!         call rdresvin() 
          if(IsFileTypeTB0(fln(8))) then
d            if(iopt.eq.2)print*,' In sub - gone to read_resvin_ef'
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call read_resvin_ef()  !EnSim compatible tb0 file
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          else
!            call rdresvin()   
             print*,'Old format .tin files not accepted'
             print*,'Please create yyyymmdd_rin.tb0 files & rerun'
             stop 'Program aborted in sub @ 979'
          endif
d          if(iopt.eq.2)print*,' In sub - back from read_resvin_ef'
d          write(98,*)'reservoir inflows read'
          endif        

 
!     rev. 9.5.52  Jan.  20/09  - NK: added reading yyyymmdd_div.tb0 for diversions
!     rev. 9.8.53  Mar.  20/13  - NK: Add Lake St. Joseph diversion algorithm to REROUT.f
!                                 This added the divertflg to see whether to read the 
!                                 diversion data or generate it.
!         Diversions    Diversions    Diversions    Diversions    Diversions
          diversion=.false.
          if(divertflg.eq.'y')then
            INQUIRE(FILE=fln(52),EXIST=exists)
            if(exists)then
              diversion=.true.
                call read_divert(282,52)
              if(iopt.ge.1)print*,'Diversion file read:',fln(52)(1:40)
              if(iopt.ge.3)pause 'after read_divert in sub'
            else
!     rev. 9.9.42  Nov.  26/14  - NK: Added errer check if diversion does not exist 
              print*,'Error:'
              print*,'divertflg in the event file = `y`'
              print*,'but',fln(52)(1:40)
              print*,'i.e.  diver\yyyymmdd_div.tb0'
              print*,'is not found'
              stop 'Program aborted in sub @ 1827'
            endif
          else
           nodivert=0    ! used as a flag in rerout
          endif

!     rev. 9.8.24  Aug.  07/12  - NK: Added reading yyyymmdd_lvl.tb0 for lake levels
          lvlflg=.false.
            if(dds_flag.eq.0)then     !added Feb. 20/14 nk
c             INQUIRE(FILE=fln(53),EXIST=exists)
c             if(exists)then
                
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~          
              call read_level(283,53)
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
              
              print*,'Lake levels read:',fln(53)(1:40)
c              lvlflg=.true.
!             open levels output file results\levels.txt
              if(id.eq.1.and.lvlflg)then  
!               open 'results\levels.csv'    !lake level comparison                  
                open(unit=953,file=filename(953),status='unknown',
     *                    iostat=ios)
                if(ios.ne.0)then    ! added Nov. 10/14  nk
                  print*
                  print*,'Unable to open file',filename(953)(1:40)
                  print*,'Possible cause(s):'
                  print*,'file in use by another application'
                  print*,'or target directory does not exist'
                  stop 'Program aborted in sub.f @ 1792'
                else
                  print*,'Opened file ',filename(953)(1:50)
                endif
              endif
              print*,'~~~~~~~~~~~~~~~~~~~~~~~~'
c            endif   
          endif   ! dds_flag.eq.0

        if(icase.le.-10)then

!         TO RUN A FORECAST SIMULATION (BEYOND REAL TIME DATA)
!         AND USE ALL ENTERED (INCLUDING FORECAST) RAINFALL 
!         THE SOIL MOISTURE COMES FROM THE EVENT FILE, NOT    
!         FROM THE MET FILE.

!         WE HAVE OPTIMIZED SOIL MOISTURE OR SCALED THE RADAR
!         NOW WE WANT TO FORECAST WITH THE RAIN UP TO THE FORECAST
!         TIME.

          mhtot=nl
          nr=min(nr,mhrd)

        elseif(icase.eq.-2)then

!         FOR SCALING THE ENTIRE RAINFALL FIELD BY THE SAME AMOUNT
!         CALIBRATE TO RECORDED FLOWS ONLY UNTIL MHRD
!         SEE ALSO THE CALL TO RAIN BELOW.  

          if(nnn.le.0)write(51,6004)nr,mhrd

!         ACTIVATE THIS IF YOU WANT TO READ IN THE OPT PERIOD 
!         WITH THE STREAMFLOW INPUT MENU

!!!!!!    scale=a(1)
          mhtot=mhrd

        elseif(icase.eq.-1)then

!         FOR SMC FITTING ONLY
!         CALIBRATE TO RECORDED FLOWS ONLY FOR SMC FITTING
!         USES STREAMFLOW AND RAINFALL UNTIL MHRD
!         BUT RUNS SPL FOR THE ENTIRE PERIOD NL

          if(nnn.le.0)write(51,6004)nr,mhrd

!         ACTIVATE THIS IF YOU WANT TO READ IN THE OPT PERIOD 
!         WITH THE STREAMFLOW INPUT MENU

          mhtot=mhrd

        elseif(icase.eq.0)then

!         DO A SINGLE RUN WITH ALL AVAILABLE RAINFALL 
!         TO THE END OF THE RAINFALL RECORD
!         FOR THIS OPTION, THE AP's AT THE GAUGES WILL BE USED
!         IF AVAILABLE, OTHERWISE THE NUMBERS IN THE EVENT FILE
!         ALSO, THE SCALING FACTOR IS SET eq. TO CONV
  
          mhrd=nl
          mhtot=nl

        elseif(icase.gt.0)then

!         OR DO A PARAMETER CALIBRATION ON WHOLE HYDROGRAPH, 
!         WHERE NUMA IS THE NUMBER OF PARAMETERS TO BE OPTIMIZED
!         ALL AVAILABLE RAINFALL WILL BE USED.
!         AND FOR PARAMETER CALIBRATION USE THE SAME FOR ALL
!         ELEMENTS.

          mhtot=nl    

        endif

        if(iopt.ge.1)then
          if(id.eq.1.and.ni.gt.3.and.mhtot.gt.8000)then
          print*
          print*
          print*,' It looks like you are running more '
          print*,' than 3 years of simulation in the debug'
          print*,' mode (iopt>0) You may run out of disk space.'
          print*
          print*
c         print*,' Hit Ctrl^C to abort'
c         print*
c         pause ' Hit enter to continue with this run.'
          endif
        endif

!       VER. 9.1 - SEDIMENT/NUTRIENT COMPONENT 
!       TO CALL wqread, SEDFLG MUST BE READ IN AS 'Y' IN rdevtA
!       Read water quality data (sediments and nutrients)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(sedflg.eq.'y')call wqread
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if(iopt.ge.1)write(51,5002)nr,smc5(1),mo,conv

d        if(iopt.eq.2)print*,' In sub, passed location  234'
  
!       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        if(id.eq.1)then

d          if(iopt.eq.2)print*,'in sub, gone to flowinit'
!         make sure this call is after reading initial lake levels
!         read_pt2(280,50)....      

!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call flowinit()
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          if(routeflg.eq.'y')then
!           write a new flow_init file only if creating files for watroute  
!           first delete the old file if one exists:      
            inquire(FILE='flow_init.r2c',EXIST=exists)
            if(exists)then
              fln(99)='flow_init.r2c'
              open(99,file=fln(99),status='unknown',iostat=ios)
              close(unit=99,status='delete')
              print*
              print*,'Old ***',fln(99)(1:60),'*** deleted'
            endif
          endif

          if(initflg.eq.'y')then
            author='spl.exe (flowinit)    '
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call write_flowinit()
c            call write_lzsinit()
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            print*,'flowinit.r2c '
            print*,' written in working directory'
          endif

c          if(resumflg.ne.'y'.and.snwflg.eq.'y') then
          if(resumflg.eq.'n'.and.snwflg.eq.'y') then
!           READ THE SNOW COURSE DATA FOR INIT SNOW
!           but not if there is aresume file
!           If swe needs to be read in, set the crseflg='y'
!           read swe in the first event if snwflg = y
!           but it would be overridden by the resume file
!           values if resumflg = y
            if(IsFileTypeR2C(fln(36))) then
d             if(iopt.eq.2)then
d		   	    print*,'In sub, passed location 1430 before read_swe'
d               print*,'If program dies here, check the swe file'
d		      endif 
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d             print*,'reading swe r2c file'      
              call read_sweinit(266,36)  !EnSim compatible r2c file
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            else
              print*,'Found SWE file ',fln(36)
              print*,'Old format swe files not accepted'
              print*,'Please create EF ????_swe.r2c files & rerun'
              stop 'Program aborted in sub @ 1216'
              if(iopt.eq.2)print*,'in sub, passed location 1114'
            endif
d            write(98,*)'swe file read'
d            if(iopt.eq.2)print*,'in sub, back from read_sweinit'
          else
!           ADDED DEC. 13/98  NK
!           default swe = 0.0 - replaced be rdresume or rdswe 
!           if files are present & read.
            do n=1,naa
              do ii=1,classcount
              snowc(n,ii)=0.0
              end do
            end do
          endif

d          if(iopt.eq.2)print*,'in sub, passed location 237'
          do n=1,naa
            if(slope(n).gt.0.0)then  
              sump(n)=0.0
              sumrff(n)=0.0
              do ii=1,classcount
                ssumr(n,ii)=0.0
                intev(n,ii)=0.0
                ev(n,ii)=0.0
              end do
            endif
          end do
d          if(iopt.eq.2)print*,'in sub, passed location 238'

!     rev. 9.2.30  Feb.  07/06  - NK: Added class_distribution.txt to output

          if(.NOT.allocated(areaclass))then     
c           allocate(areaclass(no,classcount),stat=iAllocate)

!            fixed april 1/10 - has to be max # gauges  # grids  nk
c            allocate(areaclass(na,classcount),areasum(na),stat=iAllocate)

            n=max(no,na)
            allocate(areaclass(n,classcount),areasum(n),stat=iAllocate)
            if(iAllocate.ne.0) STOP
     *         'Error with allocation of areaclass in sub @ 1351'
            print*
          endif

!     rev. 9.5.65  Sep.  23/09  - NK: change class frac to whole basin values
c          if(numa.eq.0.and.iopt.ge.1)then
          if(numa.eq.0)then
            open(unit=99,file='class_distribution.txt',
     *            status='unknown',iostat=ios)
            if(ios.ne.0)then    ! added Nov. 10/14  nk
              print*
              print*,'Unable to open file  class_distribution.txt'
              print*,'Possible cause(s):'
              print*,'file in use by another application'
              stop 'Program aborted in sub.f @ 2004'
            endif
            write(51,6016)
            write(51,6017)classcount
            write(99,6017)classcount

            do n=1,na
!             initialize values
              areasum(n)=0.0
              do ii=1,classcount
                areaclass(n,ii)=0.0
              end do
            end do

            do n=1,naa
c             write(51,6014)yyy(n),xxx(n),
c    *             frac(n),(aclass(n,ii),ii=1,classcount)
c             write(*,6014)yyy(n),xxx(n),
c     *            frac(n),(aclass(n,ii),ii=1,classcount)
     *            
              if(n.gt.0.and.next(n).gt.0)then        ! ??????????????????????
                areasum(next(n))=areasum(next(n))+areasum(n)+frac(n)
                areasum(n)=areasum(n)+frac(n)
                do ii=1,classcount
                  areaclass(next(n),ii)=areaclass(next(n),ii)+
     *              areaclass(n,ii)+frac(n)*aclass(n,ii)
                  areaclass(n,ii)=areaclass(n,ii)+frac(n)*aclass(n,ii)
                end do
c               write(*,6018)n,areasum(n),(areaclass(n,ii),ii=1,classcount)
              endif
            end do

            do l=1,no
!     rev. 9.5.77  Oct.  26/09  - NK: fixed some inits for out of basin gauges
              if(inbsnflg(l).eq.1)then
                i=(ystr(l)-yorigin)/ydelta+1
                j=(xstr(l)-xorigin)/xdelta+1
                n=s(i,j)
                write(99,6013)xstr(l),ystr(l),l,gage(l),
     *            areasum(l),
     *            ((areaclass(n,ii)/areasum(n)),ii=1,classcount)
                write(51,6014)
              endif
            end do
            close(unit=99,status='keep')
          endif
c
c!     rev. 9.2.07  Jul.  29/05  - NK: soilinit moved from runoff to sub 
c!         tdum converts mm/hour to m^3/sec for a grid with frac=1.0
c           tdum=1000.*step2/3600.

          if(modelflg.eq.'r'.or.modelflg.eq.'l'
     *          .or.modelflg.eq.'i')then
!           skip soilinit for watroute
            continue
          else

!           ~~~~~~~~~~~~~~~~~~~
            call soilinit()
!           ~~~~~~~~~~~~~~~~~~~
            

            if(initflg.eq.'y')then
!             ~~~~~~~~~~~~~~~~~~~~~
              call write_soilinit()
!             ~~~~~~~~~~~~~~~~~~~~~
            endif

            endif
            
            
            
            
            
            
            
            
!     rev. 9.2.05  Jul.  15/05  - NK: reversed order of reading resume file 
!         RESUME  RESUME  RESUME  RESUME  RESUME  RESUME  RESUME  RESUME  RESUME  RESUME              
!         RESUME  RESUME  RESUME  RESUME  RESUME  RESUME  RESUME  RESUME  RESUME  RESUME              
!         RESUME  RESUME  RESUME  RESUME  RESUME  RESUME  RESUME  RESUME  RESUME  RESUME              
!         RESUME  RESUME  RESUME  RESUME  RESUME  RESUME  RESUME  RESUME  RESUME  RESUME              
!         RESUME  RESUME  RESUME  RESUME  RESUME  RESUME  RESUME  RESUME  RESUME  RESUME              
          if(resumflg.ne.'n')then   !NOTE THE .NE. 
!         .NE. IS USED BECAuSE IT CAN BE OTHER THAN 'n' OR 'y'
!         THIS SECTION IS USED TO READ IN ALL STATE VARIABLES
!         SO PROGRAM CAN BE CONTINUED WHERE IT LEFT OFF
!         IN RUNOF6, SOILINIT.F IS not SKIPPED but over written
!         ~~~~~~~~~~~~~~~~~~~~~~~~

!         NOTE: these 3 calls go together for a resumed run:
!     rev. 9.5.79  Nov.  04/09  - NK: added resumflg='s' for read_soilinit ONLY
            if(resumflg.ne.'s')then
!             these are skipped for the special case when we want to 
!             read the soilinit file (below) but not the flow init and other stuff
!             i.e. we skip this is the resumflg='s'
!             ~~~~~~~~~~~~~~~~~~~~~
              call rdresume()
!             ~~~~~~~~~~~~~~~~~~~~~
              flnNum=99
              fln(99)='flow_init.r2c'
!             There are 2 flowinit possibilities:
!             One with a resume: the flow_init.r2c file will be in the working directory
!             the other when is it done on the fly as an updating tool with
!                       strfw\yyyymmdd_fli.r2c
!             This is with a resume:
!             ~~~~~~~~~~~~~~~~~~~~~
              call read_flowinit(flnNum)
!             ~~~~~~~~~~~~~~~~~~~~~
!     rev. 9.8.91  Oct.  30/13  - NK: Got rid of lzs_init.r2c - data is in flow_init.r2c already
!             This data is in the flow_init.r2c file so it's not needed 
c              read the gridded lzs  -- but why bother?????
c              call read_r2c(268,38,'1',jz,newDataFlag)
!             ~~~~~~~~~~~~~~~~~~~~~
            endif
!           For resumflg = ‘s’, only the soil_init.r2c file will be read 
!           but the lzs and all flow variables will be initialized with 
!           streamflow.
!           Soilinit is read every time we use the resume files
            fln(99)='soil_init.r2c'
d           print*,'reading soil_init.r2c file'     
!           ~~~~~~~~~~~~~~~~~~~~~
            call read_soilinit()
!           ~~~~~~~~~~~~~~~~~~~~~
!           NOTE: if the _gsm.r2c and _swe.r2c files exist, values therein
!           will replace values from the soil_init.r2c file 
            if(resumflg.ne.'s')then           
!             i.e. we skip this if the resumflg='s'
d             print*,'reading gridded soil moisture file'      
!             ~~~~~~~~~~~~~~~~~~~~~
              call read_gsm_ef()
!             ~~~~~~~~~~~~~~~~~~~~~
d             print*,'reading swe r2c file' 
!     rev. 9.9.10  Mar.  20/14  - NK: Update swe anytime a file is found
!             but don't call if the resumflg = 'y' 
!             It will update later if the sweflag = 'u'
              if(.not.resumflg.eq.'y')then     
!               ~~~~~~~~~~~~~~~~~~~~~
                call read_sweinit(266,36)  
!               ~~~~~~~~~~~~~~~~~~~
              endif
              endif
            
!     rev. 9.9.11  Mar.  20/14  - NK: Added lake_level_init.pt2 file for a resume
!           >>>>>resume with updated lake storage
!           No matter what, we need coefficients (safe max etc.) from this file
            fln(99)='lake_level_init.pt2'   ! for resume only - generic name
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call read_pt2(99,99,nrows,ncols)
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d           print*
d           print*,'Back from read_pt2 to read initial lake levels'
d           print*
            if(nrows.ne.noresv)then
              print*
              print*,'Number of initial lake levels does not'
              print*,'match the number of lakes & reservoirs'
              print*,'No of reservoirs & lakes =',noresv
              print*,'No init lake levels =',nrows
              print*
              if(dds_flag.eq.1)pause 'Pause in sub - please fix'
              stop 'Program aborted in sub @ 1966'
            endif
            write(53,*)'in sub @ 1951 - init lake levels & datums'
            print*,'Finished reading ',fln(99)(1:40)


            
            
c take this out - it's all in the resume files <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<            
            
            
            do l=1,noresv                                                            
               b6(l)=inarray(l,1)           ! init lake level                         
!         rev. 10.1.95 Sep   11/17  - NK: Fixed LKdepth bug in sub                    
               lake_elv(l,1)=inarray(l,1)   ! to make sure it's assigned              
               b7(l)=inarray(l,2)           ! datum                                   
c               if(ncols.eq.4)then           ! 4 att's incl. name                     
               if(ncols.ge.4)then                                                     
                                                                                      
                   ! 4 attributes incl. name                                          
!                default value = 1.0 set above                                        
                 LKdepth(l)=inarray(l,3)                                              
                 if(LKdepth(l).lt.b6(l)-b7(l))then                                    
                     print*,                                                          
     *            'WARNING: Lake Depth',l,' lower than Datum'                         
                 endif                                                                
                 LKinvert(l)=b7(l)-inarray(l,3)                                       
                 LKdepth(l)=lake_elv(l,1)-LKinvert(l)                                 
                                                                                      
                 write(53,*)'l=',l,'level=',b6(l),'datum=',b7(l),                     
     *                  'LKdepth=',LKdepth(l),'LKinve=',LKinvert(l)                   
               endif                                                                  
               if(ncols.ge.6)then                                                     
                 safe_max(l)=inarray(l,4)                                             
                 qmin(l)=inarray(l,5)                                                 
                 DecayT(l)=inarray(l,6)                                               
                 write(53,*)'l=',l,'level=',b6(l),'datum=',b7(l),                     
     *               'LKdepth=',LKdepth(l),'LKinvert=',LKinvert(l),                   
     *                 'Qmin=',qmin(l),'safe_max=',safe_max(l)                        
               endif                                                                  
c               fpet_lake(l)=inarray(l,3)    ! fpet for lake l                        
c               if(iopt.eq.1)then                                                     
c               endif                                                                 
            end do                                                                  

               
            do l=1,noresv
c              b6(l)=inarray(l,1)           ! init lake level
c              b7(l)=inarray(l,2)           ! datum
              LKinvert(l)=b7(l)-inarray(l,3)
              lake_elv(l,1)=inarray(l,1)   ! to make sure it's assigned
c             fpet_lake(l)=inarray(l,3)    ! fpet for lake l
c             if(iopt.eq.1)then
c                write(53,*)'l=',l,'b6=',b6(l),'b7=',b7(l),
c     *                              'fpet=',fpet_lake(l)
c             endif        
!             lake storage must be recalculated 
!             as lake levels could be reset
!             this will overwrite the storaged read from flow_init.r2c

!     rev. 9.9.34  Oct.  17/14  - NK: Added re-compute of lake storage re: new lake levels
	        if(b6(l).ne.0.0.or.b7(l).ne.0.0)then
!               initial lake levels can be lower than the datum	        
c	          if(b7(l).lt.b6(l))then
                write(53,*)
                write(53,*)'Override flow based initial storage:'
c                store1(n)=(b6(l)-b7(l))*lake_area(l)
!               store1 & store2 are now total lake storage
 !     rev. 10.2.45 Jan.  21/19  - NK: Fixed bug in reservoir initialization in sub - n was undefined
                i=ires(l)
                j=jres(l)
                n=s(i,j)
                store1(n)=(b6(l)-b7(l))*lake_area(l)+store_dead(l)
                store2(n)=store1(n)
                if(store2(n).ge.0.0)then
                  if(b3(l).eq.0.0)then
c                    qo2(n)=b1(l)*store2(n)**b2(l)
                    qo2(n)=b1(l)*(store2(n)-store_dead(l))**b2(l)
                  else
c                    qo1(n)=store2(n)*(b1(l)+store2(n)*(b2(l)+store2(n)*
c     *                (b3(l)+store2(n)*(b4(l)+b5(l)*store2(n)))))
                    store_live=store2(l)-store_dead(l)
                    qo1(n)=store_live*
     *                (b1(l)+store_live*
     *                (b2(l)+store_live*
     *                (b3(l)+store_live*
     *                (b4(l)+b5(l)*store_live))))
                  endif   
                  qo2(n)=qo1(n)
                else
                  qo2(n)=0.0
                endif
              endif
!     end rev. 9.9.34  Oct.  17/14  - NK: Added re-compute of lake storage re: new lake levels
            end do


                
                
                
                
            print*
            print*,'NEW        Rev. 9.9.34        Oct.  17/14'
            print*,'Recalculated lake storage & outflow based'
            print*,'on the file lake_level_init.pt2'
            print*,'in the working directory'
            print*
            deallocate(inarray)
              endif   ! resume
!         END RESUME  END RESUME  END RESUME  END RESUME  END RESUME  END RESUME  END RESUME                 
!         END RESUME  END RESUME  END RESUME  END RESUME  END RESUME  END RESUME  END RESUME                 
!         END RESUME  END RESUME  END RESUME  END RESUME  END RESUME  END RESUME  END RESUME                 
!         END RESUME  END RESUME  END RESUME  END RESUME  END RESUME  END RESUME  END RESUME                 
!         END RESUME  END RESUME  END RESUME  END RESUME  END RESUME  END RESUME  END RESUME                 
!         END RESUME  END RESUME  END RESUME  END RESUME  END RESUME  END RESUME  END RESUME            


          if(iopt.ge.1)print*,'In SUB: initialization completed '
          
!     rev. 9.8.69  Jun   17/13  - NK: Fixed bug in allocating clumnunits in SUB.f
!         moved from below to make it general
          if(allocated(column_type))then
            deallocate(column_type,column_units,stat=iDeallocate)
            if(iDeallocate.ne.0)then
              print*,'Error with deallocation of resv stuff'
              print*
              stop 'Program aborted in subf @ 1575'
            endif
          endif
!         allocate for the largest number needed              
          n=max(noresv,2*no)
          allocate(column_type(n),column_units(n),stat=iAllocate)
          if(iAllocate.ne.0) STOP
     *          'Allocation Error:  arrays in sub @1558'

!         Net Basin Supply header for Great Lakes and Mackenzie
          if(noresv.gt.0.and.dds_flag.eq.0)then
!           write the header for the nbs.tb0 file
            if(resname(1).eq.'Superior     '.or.
     *        routeflg.eq.'q'     )then
!             create a nbs.tb0 file for 1D
              author='watflood                                '
              if(resname(1).eq.'Superior     ')then
                name='Net_Basin_Supply_GLAKE                  '
              elseif(routeflg.eq.'q')then                                
                name='MRBB_MASTER_INFLOWS                     '
              else
                name='Net_Basin_Supply                        '
              endif
              coordsys_temp=coordsys1
!             GreenKenue uses LatLong - code below uses LATLONG
             if(coordsys_temp.eq.'LatLong   ')coordsys_temp='LATLONG   '
             if(coordsys_temp.eq.'Cartesian ')coordsys_temp='CARTESIAN '
              zone_temp=zone1
              datum_temp=datum1
              xorigin_temp=xorigin
              yorigin_temp=yorigin
              author='WATROUTE'

              if(mo1.lt.10.and.day1.lt.10)then
                write(line,77661)year1,mo1,day1
77661           format(i4,'/',i1'/',i1)
              elseif(mo1.lt.10.and.day1.ge.10)then
                write(line,77662)year1,mo1,day1
77662           format(i4,'/',i1'/',i2)
              elseif(mo1.ge.10.and.day1.lt.10)then
                write(line,77663)year1,mo1,day1
77663           format(i4,'/',i2'/',i1)
              else
                write(line,77664)year1,mo1,day1
77664           format(i4,'/',i2'/',i2)
              endif
              startdate=line
              write(line,77665)hour1
77665         format(i2,':00:00')
              starttime=line
c              startdate='0000/01/01'
c              starttime='0:00      '
              unit_conversion=1.0
!             attribute_name='sum precipitation                       '
              attribute_units='m3/s' 
!              attribute_type='Runoff                                  '  
              source_file_name='last spl run' 
  
              if(allocated(column_type))then
                deallocate(column_type,stat=iDeallocate)
                deallocate(column_units,stat=iDeallocate)
                if(iDeallocate.ne.0)then
                  print*,'Error with deallocation of resv stuff'
                 print*
                  stop 'Program aborted in subf @ 1575'
                endif
              endif
!             allocate for the largest number needed              
              n=MAX0(noresv,2*no)
              allocate(column_type(n),stat=iAllocate)
              allocate(column_units(n),stat=iAllocate)
              if(iAllocate.ne.0) STOP
     *         'Allocation Error:  arrays in sub @1558'
              do l=1,n
                column_type(l)='float' 
                column_units(l)='m**3/s'   
              end do
!             write the header for the nbs.tb0 file
!             
!     REV. 10.1.15 Jan.  08/16  - NK: Custom coding for Mackenzie River Basin Hydraulic Model
!             write the header for the mrb_master_inflows.tb0 file
!             the last 20 lakes are not to be writted to this file but 
!             passed through rerout for natural lake routing
!                                   un,fn,nfg,ng,no_signf
!             # of MRBHM nodes = 126
              call write_flow1d_tb0(70,70,nfg,126,-1)
            endif
          endif

        endif    !id.eq.1
!       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        
!     rev. 9.8.09  Nov.  22/11  - NK: nopt(l)=0 for area_error(l) > 10%
c          if(dds_flag.eq.1)then
c            i=0
c            do l=1,no
c              if(abs(area_error(l)).gt.10.0)then
c                nopt(l)=0
c                i=i+1
c              endif
c            end do
c            if(i.gt.0)then
c              print*,'no error calcs for',i,'stations w/area errors'
c            endif
c          endif  


!     rev. 9.7.00  May.  26/10  - NK: dds with pre-emption
!       this section of code needs to stay here because we need to
!       know the no of flow & reservoir locations
!       pre-emption
        if(id.eq.1)then
        
          sum_sq_error=0.0
          dds_error=0.0

          if(.not.allocated(mean_observed))then
            allocate(mean_observed(no+noresvi),pre_sse(no+noresvi),
     *           stat=iAllocate)
            if(iAllocate.ne.0) STOP 
     *       'Error with allocation in sub @ 1598'
	      allocate(qhyd_sum(no+noresvi),qsyn_sum(no+noresvi),
     *           qhyd_mean(no+noresvi),qsyn_mean(no+noresvi),
     *           qhyd_mean_evt(no+noresvi,ni),
     *           qsyn_mean_evt(no+noresvi,ni),
     *           sum_den(no+noresvi),sum_num(no+noresvi),
     *           num_obs(no+noresvi),
     *           stat=iAllocate)
            if(iAllocate.ne.0) STOP 'Problem allocating qhyd_sum etc'
	    endif

!         initialize the sums
          do l=1,no+noresvi
            qhyd_sum(l)=0.0
            qsyn_sum(l)=0.0
            num_obs(l)=0
	      sum_num(l)=0.0
	      sum_den(l)=0.0
	      do i=1,ni
              qhyd_mean_evt(l,i)=0.0
              qsyn_mean_evt(l,i)=0.0
            end do
          end do

          if(dds_flag.eq.1)then
!           open the dds_log.txt file:
            open(unit=30,file=filename(30),status='unknown',iostat=ios)
            if(ios.ne.0)then    ! added Nov. 10/14  nk
              print*
              print*,'Unable to open file',filename(30)(1:40)
              print*,'Possible cause(s):'
              print*,'file in use by another application'
              print*,'or target directory does not exist'
              print*
              pause 'For single run, to continue - hit Enter'
              dds_flag=0
c              stop 'Program aborted in sub.f @ 2351'
            endif
!           go to the end of the file so new data is appended  
            DO WHILE(.NOT.EOF(30))
              read(30,*,iostat=ios)
            end do
          endif

!     rev. 9.7.08  Sep.  21/10  - NK: revised mean squared error weighting for DDS
          if(.not.allocated(sta_weight))then
	      allocate(sta_weight(no+noresvi),mse(no+noresvi),
     *		 stat=iAllocate)
            if(iAllocate.ne.0) STOP
     *      'Error with allocation of sta_weight in sub @491'   
          endif 

!     rev. 9.7.00  May.  26/10  - NK: dds with pre-emption
c          if(abs(dds_flag).eq.1.or.dds_flag.eq.3.or.
c     *           dds_flag.eq.6.or.dds_flag.eq.7.or.dds_flag.eq.8)then
          if(abs(dds_flag).eq.1)then
          inquire(FILE='mean_observed_flows.txt',EXIST=exists)
          if(exists)then
            open(unit=99,file='mean_observed_flows.txt',
     *                    status='unknown',iostat=ios)
            if(ios.ne.0)then    ! added Nov. 10/14  nk
              print*
              print*,'Unable to open file  mean_observed_flows.txt'
              print*,'Possible cause(s):'
              print*,'file in use by another application'
              pause 'Program will abort with hit enter'
              stop 'Program aborted in sub.f @ 2380'
            endif
d            print*,'reading the mean flows in sub'

!     rev. 9.8.03  Aug.  08/11  - NK: check no of mean observed flows in file are ok
!     rev. 10.1.71 Mar.  14/17  - NK: Revised reading mean_observed_flows in sub 
              i=0
              ios=0
              read(99,*,iostat=ios)
              DO WHILE(.NOT.EOF(99))
                read(99,*,iostat=ios)i,mean_observed(i)
d               write(*,*)i,mean_observed(i)
              end do
              ios=0
              if(iopt.ge.1)then
                print*,'no of entries in mean_observed_flows.txt =',i
                print*
              endif
              if(errflg.eq.1.or.errflg.eq.3.or.errflg.eq.5.or.
     *           errflg.eq.6.or.errflg.eq.7.or.errflg.eq.8.or.
     *           errflg.eq.100)then

                if(i.ne.no+noresvi)then
                  print*,'no of mean_observed_flows.txt =',i
                  print*,'no of flow & resv inflow statns =',no+noresvi
                  print*,'These numbers must be equal'
                  print*,'Remove file & run splx to get a new file'
                  print*,'It is needed only for DDS'
                  print*,'program paused in sub @ 1830'
                  if(dds_flag.eq.1)pause 'Hit enter to abort splx'
                  stop 'Program aborted in sub @ 1904'
                endif
              endif
c             rewind(unit=99)
            close(unit=99,status='keep')
!     rev. 9.8.03 -- end              

!     rev. 9.8.16  Mar.  21/11  - NK: reinstate reservoir inflow error for dds
c              read(99,*,iostat=ios)
c            if(resinflg.ne.'y')noresvi=0
c            do l=1,no+noresvi
c              read(99,*,iostat=ios)i,mean_observed(l)
cd             write(*,*)i,mean_observed(l)
c            end do
            if(i.ne.no+noresvi)then
                print*,'File found but'
                print*,'Problems reading mean_observed_flows.txt'
                print*,'# entries found =',i
                print*,'# entries needed =',no+noresvi
                print*,'for # flow stations =',no
                print*,'and # reservoir inflows =',noresvi
                print*,'Have you changed the number of stations?'
                print*,'Do you have a reservoir or lake inflow file?'
                print*,'Please fix or remove the file'
                print*,'Remove file & run splx to get a new file'
                print*,'It is needed only for DDS'
                stop 'Program aborted in sub @ 1792'
            endif
c            close(unit=99,status='keep')
!     rev. 9.7.08  Sep.  21/10  - NK: revised error weighting for DDS
            if(errflg.eq.1.and.dds_flag.ge.0)then
              print*
              print*,'mean_observed_flows.txt file  found'
              print*,'station mean squared errors will be weighted'
              print*,'with station weights based on'
              Print*,'Weight(i)=sta_MSE(i)/sum(sta_MSE(n))'
	        print*,'where n is the number of stations used for'
	        print*,'calibration'
              print*
	        print*,'Calculating station weights for DDS'
	        print*,'Table in spl.txt'
	        print*
            endif

!     rev. 9.7.08  Sep.  21/10  - NK: revised mean squared error weighting for DDS
            i=0
	      sum_mean_flow=0.0
            do l=1,no
              if(nopt(l).eq.1.and.mean_observed(l).gt.0.0)then
                sum_mean_flow=sum_mean_flow+mean_observed(l)
	          i=i+1
	        endif
	      end do
!     rev. 9.8.16  Mar.  21/11  - NK: reinstate reservoir inflow error for dds
            do l=1,noresvi
              if(nopti(l).eq.1.and.mean_observed(no+l).gt.0.0)then
                sum_mean_flow=sum_mean_flow+mean_observed(no+l)
	          i=i+1
	        endif
	      end do

	      write(51,*)
	      write(51,*)'Station weights for DDS error calculation:'
	      do l=1,no
	        if(nopt(l).eq.1.and.mean_observed(l).gt.0.0)then
	          sta_weight(l)=mean_observed(l)/sum_mean_flow
	        else
	          sta_weight(l)=0.0
	        endif
	        if(nopt(l).eq.1)then
	          write(51,*)l,mean_observed(l),sum_mean_flow,sta_weight(l)
	        endif
	        num_obs(l)=0          !initialize number of observations
	      end do
!     rev. 9.8.16  Mar.  21/11  - NK: reinstate reservoir inflow error for dds
	      do l=1,noresvi
	        if(nopti(l).eq.1.and.mean_observed(no+l).gt.0.0)then
	          sta_weight(no+l)=mean_observed(no+l)/sum_mean_flow
	        else
	          sta_weight(no+l)=0.0
	        endif
	        if(nopti(l).eq.1)then
	          write(51,*)l,mean_observed(no+l),
     *	                	  sum_mean_flow,sta_weight(no+l)
	        endif
	        num_obs(no+l)=0          !initialize number of observations
	      end do

          else
d            print*,'no=',no
            do l=1,no+noresvi
              mean_observed(l)=1.0
            end do

            i=0
	      do l=1,no
	        if(nopt(l).eq.1)then
	          i=i+1
	        endif
	      end do
	      do l=1,no
	        if(nopt(l).eq.1)then
	          sta_weight(l)=1.0/float(i)
	        else
	          sta_weight(l)=0.0
	        endif
	        if(nopt(l).eq.1)then
	          write(51,*)l,mean_observed(l),sum_mean_flow,sta_weight(l)
	        endif
	      end do
!     rev. 9.8.16  Mar.  21/11  - NK: reinstate reservoir inflow error for dds
	      do l=1,noresvi
	        if(nopti(l).eq.1)then
	          i=i+1
	        endif
	      end do
	      do l=1,noresvi
	        if(nopt(no+l).eq.1)then
	          sta_weight(no+l)=1.0/float(i)
	        else
	          sta_weight(no+l)=0.0
	        endif
	        if(nopti(no+l).eq.1)then
	          write(51,*)l,mean_observed(no+l),
     *			  sum_mean_flow,sta_weight(no+l)
	        endif
	      end do
            print*,'WARNING:'
            print*,'mean_observed_flows.txt file not found'
            print*,'all flagged stations will be equally weighted'
	      print*,'--including stations without any flows if flagged--'
            print*
            if(errflg.eq.1.or.errflg.eq.7.or.
     *         errflg.eq.8.or.errflg.eq.100)then
              print*,'For errflg=',errflg,' mean flow file is required'
              print*
            endif
          endif
        endif
        endif  ! dds_flag.eq.1

!       read in new snow course data and replace computed swe
!       moved  April 12/01  nk
!       moved here June 09/06  nk
!       moved here so swe can be read in even if run starts with
!       a resume file. SWE from here overrides the resume file values
!       I.e. to use swe from the swe file, it has to be explicitly 
!       asked for with the crseflg. Just the snwflg is not enough!

        if(crseflg.eq.'y'.and.snwflg.eq.'y')then
!         this is where swe is updated at the start of an new event        
!         in a sequence to update the swe
          if(IsFileTypeR2C(fln(36))) then
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d           print*,'reading swe r2c file @ 1905'      
            call read_sweinit(266,36)    !EnSim compatible r2c file
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          else
            print*,'Old format swe files not accepted'
            print*,'Please create EF ????_swe.r2c files & rerun'
            stop 'Program aborted in sub @ 1321'
          endif
          write(98,*)'swe file read'
        endif
        
        if(resumflg.eq.'y'.and.crseflg.eq.'y')then
          print*
          print*,'WARNING  WARNING    WARNING  WARNING    WARNING  '
          print*,'Snow course swe used - not swe in the resume file'
          print*
          write(51,*)
          write(51,*)'WARNING'
          write(51,*)'Snow course swe used - not swe in the resume file'
          write(51,*)
        endif

!       END OF INITIALIZATION SECTION

! ^^^^^^^^^^^^^^^^^^^^^^^^^



!       WATROUTE START    WATROUTE START  WATROUTE START  WATROUTE START

        if(modelflg.ne.'n')then
          if(id.gt.1)close(unit=261,status='keep')
!         read the header in the runoff file:
d         print*,'reading runoff file'      
          call read_r2c(261,31,'1',jz,newDataFlag)
          if(modelflg.eq.'l')then
          if(id.gt.1)close(unit=263,status='keep')
!           read the header in the leakage (baseflow) file:
d           print*,'reading leakage file'      
            call read_r2c(263,33,'1',jz,newDataFlag)
          elseif(modelflg.eq.'r')then
            if(id.gt.1)close(unit=262,status='keep')
!           read the header in the recharge file:
d           print*,'reading recharge file'      
            call read_r2c(262,32,'1',jz,newDataFlag)
          endif
        endif

!       WATROUTE END   WATROUTE END    WATROUTE END    WATROUTE END 



        if(iopt.eq.2)print*,' In sub, passed location  244'

        a66=a6

c!         TIMER SETS ALL THE CLOCKS i.e. SUM HOURS, SUM SECONDS, ETC.
c
c        time=0.0

        if(iopt.eq.99)then
!         THIS OPTION IS TO CHECK ALL INPUT FILES
!         see above for more
          mhtot=kt*2
        endif

!       rev. 9.8.21  Jun.  18/12  - NK: Assed swe observed date & report
!       check & see if there is time series snow course data
!       this subroutine will check if the file yyyymmdd_swe.tb0 exists
!       If it exists, compute swe will be compared with observed values
!       unit number is 951 & filename(951) is used
!     rev. 10.2.63 Sep.  09/19  - NK Fixed check for file exists for fln(54) - swe time series
        call find_filetype(54)  
        inquire(file=fln(54),EXIST=exists)
        if(.not.exists)courseflg=.false.
        if(courseflg)then
          if(filetype(1:3).eq.'tb0')then
d           if(iopt.eq.2)print*,' Gone to read_swe'
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call read_swe()
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          elseif(filetype(1:3).eq.'ts5')then
!     rev. 10.2.51 Apr.  03/19  - NK: New section to read ts5 format file for swe
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call read_ts5(284,54)        !( not an .nc read)
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            print*,'Finished reading ',fln(54)(1:30)
		if(id.eq.1)then
			n_max=xcount_temp
			nold=xcount_temp
			nlines=ycount_temp
			nlines_old=ycount_temp
              nswe=xcount_temp
d              print*,'nswe',nswe
	      allocate(course_obs(nswe,max0(366,nlines)),
     *              course_calc(nswe,max0(366,nlines)),
     *             	      indomainflg(nswe),stat=iAllocate)
!           corrected bug above - replaced na by nswe   Nov. 19/14  nk     
            If(iopt.ge.1)print*,'course_obs dimensioned for:',
     *                        max0(366,nlines)
			if(iAllocate.ne.0) STOP
     *		'Error with allocation of  arrays in read_swe @ 147'
		else                    !  firstpass
!         check to see memory allocated is adequate      

			if(n.ne.nold)then
				print*,'No of swe stations has been changed in'
				print*,'nold=',nold,' n=',n
				print*,'in file ',fln(99)(1:60)
				print*,'This is not permitted'
				print*
				stop 'Program aborted in read_swe @ 156'
			endif
			if(nlines.gt.nlines_old)then
				print*,'No of data lines have changed in'
				print*,'nlines_old=',nlines_old,' nlines=',nlines
				print*,'in file ',fln(99)(1:60)
				print*
c			if(n.gt.n_max)then
				n_max=n

!           the file length is longer than any of the previous events so 
!           more memory has to be allocated

!           DEALLOCATION OF ARRAYS FROM AREA10A:
				deallocate(course_obs,course_calc,stat=iDeallocate)
				if(iDeallocate.ne.0)then
					print*,'Error with deallocation ofswe)obs'
					print*
					stop 'Program aborted in read_swe @ 170'
				endif

                print*,'reallocation for more swe locations',nswe

				allocate(course_obs(nswe,max(366,nlines)),
     *                    course_calc(nswe,max(366,nlines)),
     *         				stat=iAllocate)
				if(iAllocate.ne.0) STOP
     *			'Allocation Error: arrays in read_swe @177'
                nlines_old=nlines
d               print*,'course_obs re-dimensioned for:',max(366,nlines)
			endif
		endif                   !  firstpass
		
          if(.not.allocated(xswe))then
              allocate(xswe(nswe),yswe(nswe),gname_swe(nswe),
     *                   ewg_swe(nswe),sng_swe(nswe),stat=iAllocate)
              if(iAllocate.ne.0)then
                  print*,'Error with allocation of station descripters'
                  print*,' in read_swe'
                  STOP 'Program aborted in read_swe @ 213'
              endif
          endif
!     rev. 9.1.68  Dec.  19/04  - NK: rewrote read_tbo c/w memory allocation 
!       turn into local coordinates
          do n=1,nswe
            ewg_swe(n)=int((x_temp(n)-xorigin)/xdelta+1.0)
            sng_swe(n)=int((y_temp(n)-yorigin)/ydelta+1.0)
            gname_swe(n)=gname_temp(n)
          end do
d          print*,(ewg_swe(n),n=1,nswe)
d          print*,(sng_swe(n),n=1,nswe)
d          print*,(gname_swe(n),n=1,nswe)
            
          do l=1,xcount_temp
!           convert to local coordinate unit system 
!           and check that the stations are in the watershed
            j=int((x_temp(l)-xorigin)/xdelta)+1
            i=int((y_temp(l)-yorigin)/ydelta)+1
!           find the rank for this data point
            if(s(i,j).le.0)then
              print*,'coordinates are outside watershed'
              print*,'for location No.',i,j 
            endif
          end do
          do j=1,ycount_temp
            do i=1,xcount_temp
              course_obs(i,j)=inarray(j,i)
            end do
d           print*,(course_obs(i,j),i=1,xcount_temp)
          end do  

!         Write the header in the results\swe.csv file          
          write(951,91101)(gname_swe(n),gname_swe(n),n=1,nswe)
91101   format(<2*nswe>(a9,','))      

!         check to see if stations are in the model domain
          do n=1,nswe
            if(ewg_swe(n).ge.1.and.ewg_swe(n).le.xcount.and.
     *        sng_swe(n).ge.1.and.sng_swe(n).le.ycount)then
              if(s(sng_swe(n),ewg_swe(n)).gt.0)then
              indomainflg(n)=.true.
d             print*,n,xswe(n),yswe(n),s(sng_swe(n),ewg_swe(n)),
d    *         indomainflg(n)
              else
                  indomainflg(n)=.false.
              endif
            else
              indomainflg(n)=.false.
            endif
          end do
!     END rev. 10.2.51 Apr.  03/19  - NK: New section to read ts5 format file for swe
          
          
            open(unit=99,file='swe_location.xyz',status='unknown')
            do l=1,xcount_temp
                write(99,*)x_temp(l),y_temp(l),l,gname_swe(l)
                write(*,*)x_temp(l),y_temp(l),l,gname_swe(l)
            end do
            close(unit=99,status='keep')
            
          else
            if(dds_flag.eq.0)then
              print*,fln(54)(1:60),'not found'
	        print*,'program continues without swe analysis'
	      endif
	      courseflg=.false.
          endif
        endif
        
!     rev. 9.8.49  Feb.  20/13  - NK: Added n=municipla & irrigation withdrawals
d       if(iopt.eq.2)print*,' Gone to withdraw'
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call withdraw()
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     rev. 10.1.11 Dec.  11/15  - NK: Revised ice factor initialization and calculation   
!       re-initialize ice factors Jan 1 each year.
c        if(jul_day_now.eq.1)then
c          do n=1,naa
c            ice_fctr_min(n)=0.5
c            ice_fctr_max(n)=0.5
c          end do
c        end if

        if(numa.eq.0.and.dds_flag.eq.0)print*

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!         THIS IS THE MAIN TIME LOOP, EXECUTED FOR EACH TIME STEP
!         --------------------------------------------------
!         --------------------------------------------------
!         --------------------------------------------------
!         --------------------------------------------------
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!       start time loop
!       start time loop
!       start time loop
!       start time loop
!       start time loop

        if(.NOT.allocated(inarray))then    
!         inarray not previously allocated
          allocate(inarray(ycount,xcount),stat=iAllocate)
          if(iAllocate.ne.0)then
            STOP 'Error with allocation of inarray in read_sub @ 1623'      
          end if
        endif

!        do while(time.le.float(mhtot-1))
!        do while(time.le.float(mhtot))
!     rev. 9.2.43  Jun.  21/06  - NK: fixed spikes in route
!        do while(time.lt.float(mhtot+kt))

!      cant do this when running opt
cd     mhtot=2
c      print*,'in sub @ 1553 mhtot set to 24 hours'

        do while(time.lt.float(mhtot))
          time=time+1.000
          totaltime=totaltime+1.0

          if(iopt99)then
            call date_and_time(cday,ctime)
            read(ctime,11111)hhhhh
            read(ctime,11112)mmmmm
            read(ctime,11113)sssss
11111       format(f2.0)
11112       format(2x,f2.0)
11113       format(4x,f6.0)      
          endif

d          if(iopt.eq.2)print*,' Gone to timer'
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call timer(iz,jz,mz,clock,time,t,thr,dtmin,dtmax,div,m,
     *              ju,a66)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d          if(iopt.eq.2)print*,' Back from timer'

!     rev. 10.1.54 Nov.  25/16  - NK: Moved tdum under call timer in sub
      tdum=1000.*step2/t

!     rev. 9.8.22  Jul.  17/12  - NK: Added resetflg to reset cumm. precip Sept.1
        if(resetflg)then
          if(jul_day_now.gt.265.and.jul_day_now.lt.275)then !Oct. 1 approx
            do n=1,naa
              sump(n)=0.0
              wfo_cum_p(n)=0.0
              do ii=1,classcount
                intevt(n,ii)=0.0
                evt(n,ii)=0.0
                ssumr(n,ii)=0.0
                sum_sublim(n,ii)=0.0
              end do
            end do
          endif
        endif
!     rev. 9.8.22  Jul.  17/12  - NK: Added resetflg to reset cumm. precip Sept.1
        if(resetflg)then
!         reset the cummulative precip to 0 at the beginning of Sept.        
          do n=1,naa
            if(jul_day_now.gt.265.and.jul_day_now.lt.275)wfo_cum_p(n)=0.0
          end do
        endif


        if(iz.lt.jz)then 
!         STATUS LINE:
          if(iopt.lt.99.and.abs(dds_flag).ne.1)then
!           THIS IS TO PREVENT DIV BY 0
            if(mhtot.gt.10)then       ! added Mar. 15/14 NK to prevent div/0
              if(mod(jz,mhtot/10).eq.0.and.mhtot.gt.10)then
                if(icase.le.0)then
                  if(icase.eq.-1.or.icase.eq.-11)then
                    write(6,5001)id,ni,jz,mhtot,nnn,maxn,smc5,
     *                        optlow
                  elseif(numa.eq.0.and.na.gt.100)then 
                    write(6,5003)id,ni,jz,mhtot,nnn,maxn,scale,
     *                    optlow
                  endif
                endif
              endif
            endif
          endif
!         FIRST AND SUBSEQUENT PASSES  
!         IF FLGEVP = 'n' IN THE EVENT FILE, FLGEVP2 IS SET TO 0
!         IN SPL
!         DON'T NEED TEMPERATURES FOR EVAP OPTION .LE. 1.0

!     rev. 9.9.10  Mar.  20/14  - NK: Update swe anytime a file is found
!         create a file name for today's snow course
!         use the 'u' for update - this will do the first event also
!         for crseflg = 'y', swe will only be used at the beginning of an event
!         and an error will result if the file is not there.
          if(crseflg.eq.'u')then
            if(hour_now.eq.12)then
              if(month_now.lt.10.and.day_now.lt.10)then
                write(line,77671)year_now,month_now,day_now
77671           format('snow1\',i4,'0'i1,'0',i1,'_swe.r2c')
              elseif(month_now.lt.10.and.day_now.ge.10)then
                write(line,77672)year_now,month_now,day_now
77672           format('snow1\',i4,'0',i1,i2,'_swe.r2c')
              elseif(month.ge.10.and.day_now.lt.10)then
                write(line,77673)year_now,month_now,day_now
77673           format('snow1\',i4,i2,'0',i1,'_swe.r2c')
              else
                write(line,77674)year_now,month_now,day_now
77674           format('snow1\'i4,i2,i2,'_swe.r2c')
              endif
              fln(99)=line
              INQUIRE(FILE=fln(99),EXIST=exists)
              if(exists)then
                print*,'Updating the swe with fln',fln(99)(1:40)
                call read_sweinit(99,99)
                print*,'NEW <<<<<'
                print*,'swe updated with file:'
                print*,fln(99)(1:72)
              endif
            endif
          endif
          
!     rev. 10.2.55 June  09/19  - NK Added read_swe_update.f90 for read in swe adjustment factors
          if(swe_update)then 
!             This order is with the least # checks - but fix so we find hour from sart so there needs to be just one check              
              if(year_now.eq.yy_swe)then
              if(month_now.eq.mm_swe)then
              if(day_now.eq.dd_swe)then
              if(hour_now.eq.hh_swe)then
              if(swe_mult.gt.-998.)then
                  do n=1,naa
                      do ii=1,classcount
                          snowc(n,ii)=snowc(n,ii)*swe_mult
                      end do
                  end do
                  if(swe_mult.gt.-998..and.swe_add.gt.-998.)then
                      write(98,*)'Warning:'
                      write(98,*)'Both swe mult & add factors > -999'
                      write(98,*)'Mult applied on '
                      write(98,*)year_now,month_now,day_now,hour_now
                      write(98,*)
                  endif
             elseif(swe_add.gt.-998.)then
                  do n=1,naa
                      do ii=1,classcount
                          snowc(n,ii)=snowc(n,ii)+swe_add
                          snowc(n,ii)=amax1(0.0,snowc(n,ii))
                      end do
                  end do
              endif
              endif
              endif
              endif
              endif
          endif

!     rev. 10.2.62 Sep.  09/19  - NK Added read_sm_update.f90 for read in sm adjustment factors
          if(uzs_update)then 
!             This order is with the least # checks - but fix so we find hour from sart so there needs to be just one check              
              if(day_now.eq.dd_uzs)then            ! check once a day
              if(hour_now.eq.hh_uzs)then           ! then check once an hour
              if(year_now.eq.yy_uzs)then
              if(month_now.eq.mm_uzs)then
              if(uzs_mult.gt.-998.)then
!                  print*,uzs_mult
                  do n=1,naa
                      do ii=1,classcount
!                         set upper limit for uzs adjustment                          
                          uzs(n,ii)=amin1(retn(ii),uzs(n,ii)*uzs_mult)
                      end do
                  end do
                  if(uzs_mult.gt.-998..and.uzs_add.gt.-998.)then
                      write(98,*)
     *                     'Warning: Both uzs mult & add factors > -999'
                      write(98,*)'INFO: Mult applied on '
     *                      ,year_now,month_now,day_now,hour_now
                      write(98,*)
                  endif
             elseif(uzs_add.gt.-998.)then
                  do n=1,naa
                      do ii=1,classcount
!                         set limits for uzs adjustment                          
                          uzs(n,ii)=amin1(retn(ii),uzs(n,ii)+uzs_add)
                          uzs(n,ii)=amax1(0.0,uzs(n,ii))
                      write(98,*)'INFO: UZS addition applied on '
     *                      ,year_now,month_now,day_now,hour_now
                      end do
                  end do
              endif
              endif
              endif
              endif
              endif
          endif

!     rev. 10.2.21 Apr.  14/18  - NK: Added Lake Level update 
          if(ruleflg)then
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              call lake_lvl_update(jz,date,time)
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          endif
          
!       WATROUTE START    WATROUTE START  WATROUTE START  WATROUTE START

          if(modelflg.ne.'n')then     ! i.e. modelflg.eq.i or l or r

!           WATROUTE only <<<<<< !!!!!!
!           the headers have been read above.
            newDataFlag=.false.
            call read_r2c(261,31,'0',jz,newDataFlag)
            if(xcount.eq.xcount_temp.and.ycount.eq.ycount_temp)then
!             vectorize & convert mm to flow

      
              do n=1,naa
              i=yyy(n)
              j=xxx(n)
!     rev. 10.1.58 Dec.  06/16  - NK: corrected tdum > tdum1 for modelflg=i
!             need tdum1 here for hour time step conversion              
              qr(n)=inarray(i,j)*tdum1*frac(n)
              end do
            else
              print*,'runoff grid size does not match the shed grid'
              print*
              stop 'Program aborted in sub @ 1406'
            endif
            if(modelflg.eq.'r')then
!             read the recharge and route through the lz
              call read_r2c(262,32,'0',jz,newDataFlag)
              if(xcount.eq.xcount_temp.and.ycount.eq.ycount_temp)then
!             vectorize & convert mm to flow
              do n=1,naa
                i=yyy(n)
                j=xxx(n)
                lzs(n)=lzs(n)+inarray(i,j)
!     rev. 9.4.08  May.  29/07  - NK: changed baseflow argument list
!               call baseflow(n,ibn(n),dlz,sdlz,tdum)
                call baseflow(n,dlz,sdlz,tdum)
                qr(n)=qr(n)+qlz(n)
              end do
              else
              print*,'runoff grid size does not match the shed grid'
              print*
              stop 'Program aborted in sub @ 11420'
              endif
            endif
            if(modelflg.eq.'l')then
!             read qlz = groundwater flow (leakage/baseflow)       
              call read_r2c(263,33,'0',jz,newDataFlag)
              if(xcount.eq.xcount_temp.and.ycount.eq.ycount_temp)then
!             vectorize & convert mm to flow
              do n=1,naa
                i=yyy(n)
                j=xxx(n)
!     rev. 10.1.58 Dec.  06/16  - NK: corrected tdum > tdum1 for modelflg=i
!               need tdum1 here for hour time step conversion              
                qr(n)=qr(n)+inarray(i,j)*tdum1
              end do
              else
              print*,'runoff grid size does not match the shed grid'
              print*
              stop 'Program aborted in sub @ 1434'
              endif
            endif

!     rev. 9.5.75  Oct.  26/09  - NK: commented "deallocate in sub for watroute reads
!     this deallocate killed the watroute run. Used only for watroute 
!     so we can probably leave it out at theis point as we need to use 
!     inarray to read the files for each following time step.
c               deallocate(inarray)

!       WATROUTE END   WATROUTE END    WATROUTE END    WATROUTE END 

          else         ! modelflg.eq.........

!           DO THE MODELLING!!!!!!!!!

!           REV. 8.90 - Dec.  04/98 - INPUT TO MEMORY FOR OPT RUNS
!           RAIN PUTS PROPER AMOUNT OF PRECIPITATION ON EACH SQUARE, 
!           ONCE, ONLY AT THE BEGINNING OF EACH HOUR
!           AVERP IS THE AVERAGE OF THE PRECIP READ IN AT THE 
!           PRECIP STATIONS.      


!           Add the binary file option here - unit 91
d            if(iopt.eq.2)print*,' Gone to rain'
!     rev. 9.5.30  May.  26/08  - NK: conv back in read_rain & process_rain arg. list
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d     if(iopt.eq.4)pause 5000

!     rev. 10.2.46 Jan.  21/19  - NK: Added netCDF capability
          if(.not.netCDFflg)then
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              call read_rain_ef('0',conv,jz,jan) !EnSim compatible r2c file
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          elseif(FLtype(10)(1:2).eq.'nc')then
              do i = 1,ycount
                  do j=1,xcount
                      p(i,j)=0.0
                  end do    
              end do
              if(mod(jz+5,6).eq.0)then
!                 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  call read_2D_pcp_nc((jz+5)/6,10,1)                ! netCDF
!                 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  endif
!              ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              call disaggregate(jz)
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          endif
d     if(iopt.eq.4)pause 5100

c            else
c              call read_rain_bin('0',conv,jz,jan) !EnSim compatible r2c file
c            endif
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call process_rain(conv,scale)
!          call process_rain_gcm(conv,scale)
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     rev. 9.7.15  Dec.  14/10  - NK: Create reduced precip & temp files for sub-basins

            if(new_precip_grid_flg)then
            if(xcount.ne.xcount2.or.ycount.ne.ycount2)then
!             create a new name for the met file       
              if(jz.eq.1)then
                write(line,66001)fln(10)(1:6),'new_grid\',fln(10)(7:73)
c                write(*,66001)fln(10)(1:6),'new_grid\',fln(10)(7:73)
66001           format(a6,a9,a66)
                read(line,*)fln(12)
c                print*,fln(12)              
!               write the header for the new condensed met file    
                author='watflood                    '
                name='Precipitation                       '
                coordsys_temp=coordsys1
!               GreenKenue uses LatLong - code below uses LATLONG
             if(coordsys_temp.eq.'LatLong   ')coordsys_temp='LATLONG   '
             if(coordsys_temp.eq.'Cartesian ')coordsys_temp='CARTESIAN '
                zone_temp=zone1
                datum_temp=datum1
                xorigin_temp=xorigin
                yorigin_temp=yorigin
                xcount_temp=xcount
                ycount_temp=ycount
                xdelta_temp=xdelta
                ydelta_temp=ydelta
                attribute_name='Precipitation                     '
                attribute_units='mm                            ' 
                attribute_type='                              '
                unit_conversion=1.0   
                source_file_name=fln(10)     
                no_frames=mhtot*ni+1
!               write the header for new precipitation file
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                call write_r2c(42,12,no_frames,1,0,1,1)   
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              endif
              if(new_precip_flg)then
                dataflg=.false.
                do i=1,ycount
                  do j=1,xcount
                    if(p(i,j).gt.0.0)dataflg=.true.
                  end do
                end do
                if(dataflg)then      
                  write(42,66003)precip_frame_header
66003             format(a60)                
                  do i=1,ycount
                    write(42,66002)(p(i,j),j=1,xcount)
66002               format(9999f6.1)                  
                  end do
                endif
              endif  
            endif
            endif
        
            do n=1,naa
              i=yyy(n)
              j=xxx(n)
c              print*,p(i,j)
c              domain_precip=domain_precip+p(i,j)*grid_area(n)
            end do

!           adjust the precip
d            if(iopt.eq.2)print*,' Back from rain, gone to rdtemp'

!            if(flgevp2.ge.2.0.or.snwflg.eq.'y')then
            if(vapflg.eq.'y'.or.snwflg.eq.'y')then
              if(time.le.float(mhtot))then

c              if(snwflg.eq.'y'.or.vapflg.eq.'y')then  ! added nk Jun. 28/06
c                if(IsFileTypeR2C(fln(15))) then

                if(netCDFflg)then
!     rev. 10.2.46 Jan.  21/19  - NK: Added netCDF capability
c              if(mod(jz,deltaT3).eq.0)then
                  if(mod(jz+5,6).eq.0)then
!                     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      call read_2D_tmp_nc((jz+5)/6,15,1)    
 !                    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  endif
c             endif
c                  endif
                 else
 !                ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  call read_temp_ef('0',jan,jz)     !EnSim compatible r2c file
  !                ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                endif

!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  call process_temp(jz)
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                 call process_temp_gcm(jz)
                  if(new_temp_grid_flg)then
                    if(xcount.ne.xcount3.or.ycount.ne.ycount3)then
!                     create a new name for the met file       
                      if(jz.eq.1)then
                        write(line,66001)fln(15)(1:6),
     *                             'new_grid\',fln(15)(7:73)
c                        write(*,66011)fln(15)(1:6),
c     *                             'new_grid\',fln(15)(7:73)
66011                   format(a6,a9,a66)
                        read(line,*)fln(9)
c                        print*,fln(9)              
!                       write the header for the new condensed met file    
                        author='watflood                    '
                        name='temperature                         '
                        coordsys_temp=coordsys1
!            GreenKenue uses LatLong - code below uses LATLONG
             if(coordsys_temp.eq.'LatLong   ')coordsys_temp='LATLONG   '
             if(coordsys_temp.eq.'Cartesian ')coordsys_temp='CARTESIAN '
                        zone_temp=zone1
                        datum_temp=datum1
                        xorigin_temp=xorigin
                        yorigin_temp=yorigin
                        xcount_temp=xcount
                        ycount_temp=ycount
                        xdelta_temp=xdelta
                        ydelta_temp=ydelta
                        attribute_name='Temperature                    '
                        attribute_units='dC                            ' 
                        attribute_type='                              '
                        unit_conversion=0.0   
                        source_file_name=fln(10)     
                        no_frames=mhtot*ni+1
!                       write the header for new precipitation file
!                       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        call write_r2c(39,9,no_frames,1,0,1,1)   
!                       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      endif
                      if(new_temp_flg)then
                        dataflg=.false.
                        do i=1,ycount
                          do j=1,xcount
                            if(ttemp(i,j).gt.-99.0)dataflg=.true.
                          end do
                        end do
                        if(dataflg)then      
                          write(39,66013)temp_frame_header
66013                     format(a60)                
                          do i=1,ycount
                            write(39,66012)(ttemp(i,j),j=1,xcount)
66012                       format(9999f7.1)                  
                          end do
                        endif
                      endif  
                    endif
                  endif
c                else
c                  print*,'Old format .tem files not accepted'
c                  print*,'Create yyyymmdd_tem.r2c files & rerun'
c                  stop 'Program aborted in sub @ 192'
c                endif
               
c              endif  ! added nk Jun. 28/06

              endif  
            end if

d            if(iopt.eq.2)print*,' Back from rdtemp'





c!     rev. 9.9.00  Dec.  08/13  - NK: Added Lake Evaporation model
c            if(lakeflg.eq.'y')then
c!             read the data          
cd              print*,'reading wind speed on unit',288,'hour',jz
c                newDataFlag=.false.
cc              call read_r2c(288,58,'0',jz,newDataFlag) !EnSim compatible r2c file
cd              print*,'back from reading',fln(58)(1:40)
c              if(newDataFlag)then
c!               if the hdrflg = 'y' there is new data              
c!               vectorize 
c                do n=1,naa
c                  i=yyy(n)
c                  j=xxx(n)
cc                  windspd(n)=inarray(i,j)
c                end do
c              endif
cd                print*
c              
cd              print*,'reading wind direction on unit',289,'hour',jz
c              newDataFlag=.false.
cc              call read_r2c(289,59,'0',jz,newDataFlag) !EnSim compatible r2c file
cd              print*,'back from reading',fln(59)(1:40)
c              if(newDataFlag)then
c!               if the hdrflg = 'y' there is new data              
c!               vectorize 
c                do n=1,naa
c                  i=yyy(n)
c                  j=xxx(n)
cc                  winddir(n)=int(inarray(i,j))
c                end do
c              endif
c            endif   !lakeflg





!     rev. 9.9.06  Jan.  08/14  - NK: Add daily differences to Harfreaves ETHarg.f
            if(dlyflg)then

                
                
      deallocate(inarray)     !fix <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      allocate(inarray(ycount,xcount))
      
      
      
c              if(ju.ne.julast)then
!               read the data in tempr\yyyymmdd_dif.r2c
!               it is assumed that this file will have just one grid per day
!               but since it is written at the end of the day, we can not use the 
!               timestamp to 
                newDataFlag=.false.
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                call read_r2c(292,62,'0',jz,newDataFlag) !EnSim  r2c file
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if(newDataFlag)then
                  do n=1,na
                    i=yyy(n)
                    j=xxx(n)
                    dly_diff(n)=inarray(i,j)
                  end do
d                 if(iopt.ge.2)then
d                  print*
d                  print*
d                  print*,'new dly_diff found'
d                  print*
d                  print*
d                  print*
d                 endif
                endif  
c                julast=ju
c              endif
            endif

!     rev. 9.5.03  Dec.  09/07  - NK: added reads for precip isotopes
!     rev  9.9.12  Oct.  15/14  - CJD: added in frc_file_flg='c' for time series dPPT inputs.
!     rev. 10.1.01 Oct.  05/15  - NK: Isotope update: added 2H

            if(frcflg.eq.'y')then
            !TH: reading in relative humidity from EC data based r2c
              newDataFlag=.false.
              if(RH_flg.eq..true.)then
               call read_r2c(251,21,'0',jz,newDataFlag)
               if(newDataFlag)then
                do n=1,naa                
                i=yyy(n)
                j=xxx(n)           
                relh(n)=inarray(i,j)/100.
                if(relh(n).gt.0.98)relh(n)=0.98
                end do
               endif
              else
               do n=1,naa
                 relh(n)=(10/1000.)/(((1e-06*tempv(n)**4.-3e-05*
     *                 tempv(n)**3.+0.0036*tempv(n)**2.+0.0182*tempv(n)
     *                 +0.6607)*1000.)/101325.*0.625)
                  if(relh(n).gt.0.98) relh(n)=0.98  ! upper limit due to wind turbulence
	            if(relh(n).lt.0.30) relh(n)=0.30  ! can't be less than 30%
	         end do
	        endif

              if(frc_file_flg.eq.'y')then

!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!              call read_hum(251,21,'0',jan,jz)   ! humidity for REMOiso
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!             vectorize 
!              do n=1,naa
!              i=yyy(n)
!              j=xxx(n)
!              spech(n)=inarray(i,j)
!              call process_hum(time,n)
!              end do

!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              call read_gsn(277,47,'0',jan,jz)   ! snow_precip for REMIiso
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!             vectorize 
              do n=1,naa
              i=yyy(n)
              j=xxx(n)
              do ii=1,classcount
                snowf(n,ii)=inarray(i,j)*0.01
              end do
              sum_precip(i,j)=sum_precip(i,j)+inarray(i,j)
              end do
              
c!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c              call read_drn(278,48,'0',jan,jz)   ! REMIiso delta rain
c!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c!             vectorize 
c              do n=1,naa
c              i=yyy(n)
c              j=xxx(n)
c              dlt_rain(n)=inarray(i,j)
c              end do

c!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c              call read_dsn(279,49,'0',jan,jz)   ! REMOiso delta snow
cc!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c!             vectorize 
c              do n=1,naa
c              i=yyy(n)
c              j=xxx(n)
c              dlt_snow(n)=inarray(i,j)
c              end do
c!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!             Check if times series delta rain file exists             
              elseif(frc_file_flg.eq.'c')then
              newDataFlag=.false.
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              call read_r2c(278,48,'0',jz,newDataFlag)   ! time series delta rain
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d             if(iopt.eq.2)print*,'back from reading',fln(278)             
!              vectorize 
!     rev. 10.1.01 Oct.  05/15  - NK: Isotope update: added 2H
              if(newDataFlag)then
                do n=1,naa                
                i=yyy(n)
                j=xxx(n)           
                dlt_rain(n)=inarray(i,j)
                dlt_snow(n)=inarray(i,j)
                spech(n)=-999.9
  !              relh(n)=0.70 
                
                  do ii=1,classcount
                    snowf(n,ii)=0.0
                  end do
        
                  if(n.eq.nnprint)write(866,*)
     *                     totaltime,newDataFlag,dlt_rain(n)
                end do
              endif
!     rev. 10.1.01 Oct.  05/15  - NK: Isotope update: added 2H
              if(flg2H.eq.2)then
                newDataFlag=.false.
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                call read_r2c(279,49,'0',jz,newDataFlag)   ! time series delta rain
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~          
!               vectorize 
                if(newDataFlag)then
                  do n=1,naa           
                  i=yyy(n)
                  j=xxx(n)        
                  dlt2H_rain(n)=inarray(i,j)
                  dlt2H_snow(n)=inarray(i,j)
                  spech(n)=-999.9
                  end do
                endif
              endif

!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              deallocate(inarray)
                         
!             if one of the needed files is mising in this event, use the defaults
 !            else  ! use default values from isotope.par

!               do n=1,naa
!               dlt_rain(n)=deltar(1)      ! avg precip @ start of simulation
!               dlt_snow(n)=deltas(1)      ! avg snow conc @ start of simulation
!               spech(n)=-999.9
    !           relh(n)=0.70        ! default value: calculate in craig_gordon.f
!               do ii=1,classcount
!                 snowf(n,ii)=0.0
!               end do
!               end do

             endif  ! frcfileflg
            endif

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if(time.le.float(mhtot))then
!              if(flgevp2.eq.3.0)call rdrad(jan)
            endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!!!             mon=mo         !changed it back nk may 13/02
!           Changed this line to mon=mo1  AB May 10, 2002 again Aug. 1/02 nk
            mon=mo1
            dtmax=amin1(float(kt)*3600.0,float(itogo)*3600.0)
            
            aintvl=time-clock
            sintvl=aintvl*3600.0

            clock=time    ! used for r/w resume files
            tot2=tot1*conv*scale/float(naa)
c              arain=tot2  ! not used

!     rev. 9.5.56  Mar.  26/09  - NK: Fix bug with month in yearly events
            mon=mo

d            if(iopt.eq.2)print*,' Gone to runof6'
          
!     REV. 10.1.16 Jan.  11/16  - NK: Added subroutine ice_factor.f
!     rev. 10.1.44 Dec.  02/15  - NK: Reworked icerivflg & icelakeflg
!           note: strloss is set to zero in rinof6 when ice_fctr < 0
!           note: icerivflg = n for watroute
!           note: ice_factor is called if icerivflg and/or icelakeflg = y
            if(iceflg.eq.'y'.and.mod(int(time),24).eq.0)then
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              call ice_factor 
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            endif

            if(iopt99)then
              call date_and_time(cday,ctime)
              read(ctime,11111)hhhhh
              read(ctime,11112)mmmmm
              read(ctime,11113)sssss
              now=hhhhh*3600.0+mmmmm*60.0+sssss
              input=input+now-last
              last=hhhhh*3600.0+mmmmm*60.0+sssss
            endif
          
!           RUNOF5 CALC'S RUNOFF FOR EACH HOUR (OR LARGER INTERVALS)
!           fixed the time step 3600.0 sec & 1.0 min rather than 
!           basin it on aintvl & sintvl which ends up with -ve 
!           values at new events which creates problems with ak & akfs.      
!     rev. 10.1.75 Apr.  03/17  - NK: Fixed time & thr in runof6 arg list
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call runof6(jan,time,3600.0,1.0,mon,e1,mz,ju)
c            call runof6(jan,time,sintvl,aintvl,mon,e1,mz,ju)
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d            if(iopt.eq.2)print*,' Back from runof5_jan'

            if(iopt99)then 
              call date_and_time(cday,ctime)
                  read(ctime,11111)hhhhh
                  read(ctime,11112)mmmmm
                  read(ctime,11113)sssss
                  now=hhhhh*3600.0+mmmmm*60.0+sssss
                  runoff=runoff+now-last
                  last=hhhhh*3600.0+mmmmm*60.0+sssss
            endif

!           Used for Great Lakes model only
!           For future lake evap model 
!           This calculates strloss - replaces strloss from runof6
c            if(na.eq.4029.and.
c     *         int(xorigin).eq.-93.and.int(yorigin).eq.40.)then

            if(noresv.gt.0)then
              if(resname(1)(1:8).eq.'Superior'.and.rd_evp_flg)then
!     rev. 9.5.12  Feb.  13/08  - NK: added evaporation input file with read_r2c
!               this data will be used in lake_evap.f
c                if(rd_evp_flg)then

!                 read the lake evaporation from GEM output in r2c format
!                 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d                 print*,'reading the lake evap file @ 2396'
                  call read_r2c(281,51,'0',jz,newDataFlag) !EnSim  r2c file
!                 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  do n=1,naa
!     rev. 9.5.50  Jan.  05/09  - NK: read evap data for reaches only
                    if(ireach(n).gt.0)then
                      i=yyy(n)
                      j=xxx(n)
!                     inarray is evaporation in mm for the data timestep
!     replace /t by /data_time_step <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!     this replaces the strloss calculated by runof6
                      strloss(n)=inarray(i,j)*evap_convert(n)/3600.    !splx
                      sum_et(n,classcount)=
     *                        sum_et(n,classcount)+inarray(i,j) 
                    endif
                  end do
c                endif    !  rd_evp_flg
c              elseif(lakeflg.eq.'y')then
c              else
c!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c                call lake_evap
c!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              endif
            endif    !(noresv.gt.0)
            
!     rev. 10.1.83 May   09/17  - NK: Fixed lake evap bug - moved it outside lake-only loop
!     Bug created 10.1.45 Aug. 2016
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call lake_evap
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!           This is THE place to do this addition so runof6 and 
!           route and rerout can be common for the seperate programs. 
!           This way precip & loss is included in the gridflow.r2c file

!           But for wetland routing, it is taken off again in route
!           because qr is the inflow to the wetland!!

!     rev. 9.5.07  Feb.  05/08  - NK: fixed double counting of strloss & qstream

d            if(iopt.eq.2)print*,' Back from lake_evap'
 
            if(sedflg.eq.'y')then
!             send to sediment calculations >>>>>>>>>
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              call wqsed(jan,aintvl)
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!             send to nutrient calculations >>>>>>>>>
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              call wqnut
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            endif
!           REV. 8.82 - July 10/98 - ADDED RUNOFF OUTPUT OPTION: 
!           rev. 9.1.45  Jun.  11/03  - runoff, recharge and leakage files added 

!           ROUTEFLG
!           write the output for WATROUTE:
            if(routeflg.eq.'y'.and.numa.eq.0)then
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!              call write_flux(jz,dtmin,ju,tdum)
!            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     rev. 9.3.02  Jul.  18/06  - NK: converted runof, rchrg & lkage to r2c
!            runoff in mm per hour (does not include ground water)
             do i=1,ycount
               do j=1,xcount
                 outarray(i,j)=0.0
               end do
             end do
             do n=1,naa
              i=yyy(n)
              j=xxx(n)
!     rev. 10.1.58 Dec.  06/16  - NK: corrected tdum > tdum1 for modelflg=i
!             need tdum1 here for hour time step conversion              
              outarray(i,j)=
     *         (qr(n)+qstream(n)-strloss(n)-qlz(n))/tdum1   
!     *       (qr(n)-qlz(n))/tdum/frac(n) 
!     rev. 9.5.07  Feb.  05/08  - NK: fixed double counting of strloss & qstream
!             qstream(n) & strloss are not added in runof6 so it's done here
!             and in route
             end do   
             frame_no3=frame_no3+1
!     rev. 9.9.23  Aug.  10/14  - NK: Uniform arguments for write_r2c
!            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c            call write_r2c(261,31,mhtot,1,frame_no3,0,3)    
             call write_r2c(261,31,mhtot,0,frame_no3,0,3)    
!            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!            recharge in mm
             do n=1,naa
              i=yyy(n)
              j=xxx(n)
               outarray(i,j)=rechrg(n)
             end do   
!     rev. 9.9.23  Aug.  10/14  - NK: Uniform arguments for write_r2c
!            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c            call write_r2c(262,32,mhtot,1,frame_no3,0,3)    
             call write_r2c(262,32,mhtot,0,frame_no3,0,3)    
!            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!            leakage in mm/hr
             do n=1,naa
              i=yyy(n)
              j=xxx(n)
!     rev. 10.1.58 Dec.  06/16  - NK: corrected tdum > tdum1 for modelflg=i
!             need tdum1 here for hour time step conversion              
              outarray(i,j)=qlz(n)/tdum1
             end do     
!     rev. 9.9.23  Aug.  10/14  - NK: Uniform arguments for write_r2c
!            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c            call write_r2c(263,33,mhtot,1,frame_no3,0,3)    
             call write_r2c(263,33,mhtot,0,frame_no3,0,3)    
!            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            elseif(routeflg.eq.'m')then
!             write output for MODFLOW:  added May, 2000 NK
!             recharge is accumulated over 24 hour in this case
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              call write_modflow(ju,juold,jz)
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            endif    !  if(routeflg......

d            if(iopt.eq.2)print*,'  1003'
            juold=ju
            jan=2
!           JAN HAS TO STAY HERE BECAUSE ROUTE MAY BE CALLED MORE
!           THAN FOR EACH TIME RUNOF5 IS CALLED - SO WE CAN'T USE 
!           JAN IN ROUTE.

          endif     !  modelflg.eq......       
          endif       ! CALL RUNOFF ONLY IN NEW HOUR 

!         ROUTE ROUTES WATER THROUGH EACH SQUARE USING STORAGE ROUTING
!         RESET QDWPR=0 INCASE THIS IS SECOND DT DURING FLOW INCREMENT

!         >>>>>>>IS THIS LOOK OK?? CHANGED 6 TO NORESVO   
 
!     rev. 9.5.52  Jan.  20/09  - NK: added reading yyyymmdd_div.pt2 for diversions
!         Diversions    Diversions    Diversions    Diversions    Diversions
!     rev. 9.7.17  Jan.  05/11  - NK: Fixed diversions outside sub-basin
!     rev. 9.7.22. Mar.  07/11  - NK: Changed diversion code: give/route take/rerout
!     rev. 9.9.41  Nov.  20/14  - NK: Added check if diversion = in-basin 
          
          if(diversion)then
            do nd=1,nodivert
              if(divert_inbasin_flg(nd))then
                qr(gridtake(nd))=qr(gridtake(nd))-qdivert(nd,jz)
                qr(gridgive(nd))=qr(gridgive(nd))+qdivert(nd,jz)
              endif 
            end do
          endif
          
          if(irdt.le.kt)then
!           old way of doing things
            if(nastart.eq.1)then
!     rev. 9.3.11  Feb.  17/07  - NK: force hourly time steps for runoff
!     rev. 9.3.12  Feb.  20/07  - NK: changed dtmin & call to route
              no_dt=max(int(3599./dtmin)+1,1)
              route_dt=3600.0/float(no_dt)
              sec_div=route_dt/2.0
              hr_div=sec_div/3600.
!             write(55,*)dtmin,ndt,rdt,sec_dt,hr_dt
!             The value of dtmin has been used to determine how many
!             times route is called. Route will determine a new dtmin
!             for the next time step.
              dtmin=3600.0
              time2=time  ! so time incrementing doesn't get affected

              do n=1,no_dt

!     TS: added time update since routing timestep is less than increment of +1.00: Jan 4/08
                if(n.gt.1)time2=time+route_dt/3600.*float(n-1)


              if(iopt99)then
                  call date_and_time(cday,ctime)
                  read(ctime,11111)hhhhh
                  read(ctime,11112)mmmmm
                  read(ctime,11113)sssss
                  now=hhhhh*3600.0+mmmmm*60.0+sssss
                  input=input+now-last
                  last=hhhhh*3600.0+mmmmm*60.0+sssss
              endif
      
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c               SUBROUTINE route(div,    thr   ,dtmin,jz,iz,time,date)
                call route(sec_div,hr_div,dtmin,jz,iz,time2,date,tdum)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

              if(iopt99)then
                  call date_and_time(cday,ctime)
                  read(ctime,11111)hhhhh
                  read(ctime,11112)mmmmm
                  read(ctime,11113)sssss
                  now=hhhhh*3600.0+mmmmm*60.0+sssss
                  router=router+now-last
                  last=hhhhh*3600.0+mmmmm*60.0+sssss
               endif

!         rev. 9.1.47  Oct   30/03    - Tracer module parameter list changed
                if(numa.eq.0)then
!                 tracer turned off for optimization
                  if(trcflg.eq.'y')then
                    if(n.eq.no_dt)then
d                     if(iopt.eq.2)print*,'gone to tracer    (1)'                  
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      call tracer(iz,-iz,time2,route_dt,jan,tdum,1)
d                     if(iopt.eq.2)print*,'back from tracer  (1)'                  
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    else
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d                     if(iopt.ge.2)print*,'gone to tracer    (2)'                   
                      call tracer(iz,iz,time2,route_dt,jan,tdum,0)
d                     if(iopt.eq.2)print*,'back from tracer  (2)'                  
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    endif
                  endif
!     rev. 10.1.01 Oct.  05/15  - NK: Isotope update: added 2H
                  if(frcflg.eq.'y')then
                   if(n.eq.no_dt)then
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                     call isotopes(iz,jz,time2,route_dt,1)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   else
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      call isotopes(iz,jz,time2,route_dt,0)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   endif
!     end rev. 10.1.01 Oct.  05/15  - NK: Isotope update: added 2H
                  endif
                endif
d                if(iopt.eq.2)print*,' Back from tracer'
              end do
!             This little ditty is here if we want routing dt's
!             longer than 1 hr. This may be needed for very large grids.
!             do ijk=1,2
!!              call route(div,thr,dtmin,jz,iz,time,date)
!               call route(1800.,.5,dtmin,jz,iz,time,date)
!               dtmin=3600.
!             end do
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
            endif 
          else
!           rev. 9.1.38  Mar.  31/03    - revised str header and routing dt selectable
!           For longer time steps, route is called at the desired intervals. 
!           This was done just for the DMIP paper    
            if(mod(jz,irdt).eq.0)then
              if(nastart.eq.1)then
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                call route(rdt*1800.0,thr,dtmin,jz,iz,time,date,tdum)
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              endif 
            endif
          endif

c          if(wetflg.eq.'y')then
c!              call decay(time,n,l,t,ii,ha,fpw,kdn,nratio)
c          endif

d          if(iopt.eq.2)print*,' Back from route'
    
!         REV. 9.1 - SEDIMENT COMPONENT
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(sedflg.eq.'y')call wqroute(t)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d          if(iopt.eq.2)print*,' Gone to tracer'

!         rev. 9.1.47  Oct   30/03    - Tracer module parameter list changed
!          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!          if(trcflg.eq.'y')call tracer(iz,jz,time,t,jan,tdum)
!          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c              if(iopt.eq.2)print*,' Back from tracer'

!         INTERPOLATE TO GET FLOWS AT reporting TIMES:
          if(irdt.le.kt)then
!         old way of doing things
          if(nastart.eq.1)then
            if(iz.lt.jz.and.jz.ge.1)then
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              call synflw(time,thr)
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            endif      
          endif 
          else
!         rev. 9.1.38  Mar.  31/03  - revised str header and routing dt selectable
          if(mod(jz,irdt).eq.0)then
            if(nastart.eq.1)then
              if(iz.lt.jz.and.jz.ge.1)then
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              call synflw(time,thr)
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              endif        
            endif 
          endif
          endif

!         print*,' ireport=',ireport,iz,jz,mod(jz,ireport)

d          if(iopt.eq.2)print*,' Back from synflw'

          if(jz.ne.iz.and.mod(jz,deltat_report).eq.0)then
!     rev. 9.5.22  Mar.  12/08  - NK: added grdflg to print gridded flow, swe & evap
 
          if(id.lt.ni)then
!           n_trick=mhtot+100  
            n_trick=mhtot*ni+1 ! tricks write_r2c into NOT closing the file
          else
            n_trick=mhtot
          endif
          
!     rev. 10.2.52 Apr.  15/19  - NK: Added total UZS for reporting in FEWS
          if(netCDFflg)then
            do n=1,naa
              totd1(n)=0.0
              totuzs(n)=0.0
              do ii=1,classcount-3
                  if(aclass(n,ii).gt.0.0)then
c                      totd1(n)=totd1(n)
c     *                 +d1(n,ii)*aclass(n,ii)*(1.0-sca(n,ii))
c     *                 +d1fs(n,ii)*aclass(n,ii)*sca(n,ii)
                     totuzs(n)=totuzs(n)
     *                +(uzs(n,ii)/retn(ii))*aclass(n,ii)*(1.0-sca(n,ii))
     *                +(uzsfs(n,ii)/retn(ii))*aclass(n,ii)*sca(n,ii)
                  endif
              end do
c              totuzs(n)=totuzs(n)+totd1(n)    ! include d1 with uzs for reporting in FEWS
            end do
c      write(401,*)hour_now,',',totuzs(na/2)
          endif
          
!     rev. 9.5.31  May.  27/08  - NK: moved totsnw(n) computation in sub
          if(numa.eq.0.and.ensimflg.eq.'y'.or.netCDFflg)then
            do n=1,naa
              totsnw(n)=0.0
              do ii=1,classcount
              if(aclass(n,ii).gt.0.0)then
                totsnw(n)=totsnw(n)
     *              +snowc(n,ii)*aclass(n,ii)*sca(n,ii)
              endif
              end do
            end do     
          endif

!     rev. 9.5.80  Dec.  20/09  - NK: added swe_locations.txt file for swe input
!     write swe at designated snow courses
!     DDS_flag added April 4, 2010  nk
      if(.not.netCDFflg)then
      if(id.eq.1.and.no_swe_l.eq.0.and.dds_flag.ne.1)then
        inquire(FILE='basin\swe_locations.txt',EXIST=exists)
        if(exists)then
          open(99,file='basin\swe_locations.txt',
     *            status='unknown',iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file  basin\swe_locations.txt'
            print*,'Possible cause(s):'
            print*,'file in use by another application'
            print*,'or target directory does not exist'
            stop 'Program aborted in sub.f @ 3594'
          endif
          read(99,*)no_swe_l
          read(99,*)(swe_l(i),i=1,no_swe_l)
          close(unit=99,status='keep')
          write(*,64000)
          write(*,64001)no_swe_l
          write(*,64001)(swe_l(i),i=1,no_swe_l)
          if(writeflg(64))write(64,64002)(swe_l(i),i=1,no_swe_l)
          if(ensimflg.ne.'y')then
            print*,'ensimflg=',ensimflg
            print*,'ensimflg needs to be `y` or `a` to get output for'
            print*,'swe at snow courses whe specified in'
            print*,'basin\swe_locations.txt'
            print*
            print*,'If DDS_flag = 1, ensimflg is set to n'
            print*
c           stop 'Program aborted in sub @ 2204'
          endif
          print*
        endif
      endif
          
      if(ensimflg.eq.'y')then
        if(ju.eq.60.or.ju.eq.91.or.ju.eq.121.and.no_swe_l.ge.0)then
          if(writeflg(64))write(64,64000)year_now,ju,(totsnw(swe_l(i)),
     *                             i=1,no_swe_l)
64000     format(2i5,99f10.2)
64001     format(99i5)
64002     format(99i10)
        endif
      endif
      endif   !netCDFflg

!     rev. 9.5.22  Mar.  12/08  - NK: added grdflg to print gridded flow, swe & evap
          if(grdflg.eq.'y')then
!     rev. 9.9.21  Jul.  27/14  - NK: Added allocation for outarray in sub
            if(allocated(outarray))deallocate(outarray)
            allocate(outarray(ycount,xcount),stat=iAllocate)
            if(iAllocate.ne.0)then
              print*,'Error with allocation of outarray in sub @ 594'
              STOP 'Program aborted in sub @ 595'
            endif
            do i=1,ycount
              do j=1,xcount
                outarray(i,j)=0.0
              end do
            end do
            xcount_temp=xcount
            ycount_temp=ycount

            if(iopt99)then
            if(numa.eq.0.and.ensimflg.eq.'y')then
              do n=1,naa
              i=yyy(n)
              j=xxx(n)
              outarray(i,j)=qo2(n)
              end do     
!             unit=72  fln(72)=gridflow.r2c  - gridded flow - r2c file
              frame_no2=frame_no2+1
!     rev. 9.9.23  Aug.  10/14  - NK: Uniform arguments for write_r2c
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c              call write_r2c(72,72,n_trick,1,int(totaltime),1,8)  
              call write_r2c(72,72,mhtot+8784,0,frame_no2,0,8)  
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!             use 8784 as the total # of frames so the file only closes when the run is done.
            endif

!           this added for Vincent Fortin:
            if(numa.eq.0.and.ensimflg.eq.'y')then
              do n=1,naa
              i=yyy(n)
              j=xxx(n)
              outarray(i,j)=totsnw(n)
              end do     
              frame_no4=frame_no4+1
!             unit=61  fln(61)
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              call write_r2c(61,61,mhtot+8784,1,frame_no4,0,0)  
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            endif

!     rev. 9.5.09  Feb.  12/08  - NK: added evap.r2c to the output files
!           this added for Vincent Fortin:
            if(numa.eq.0.and.ensimflg.eq.'y')then
            do n=1,naa
              i=yyy(n)
              j=xxx(n)
              outarray(i,j)=eloss(n)
            end do     
            frame_no5=frame_no5+1
!           unit=100  fln(100)=evap.r2c  - gridded evap - r2c file
!     rev. 9.9.23  Aug.  10/14  - NK: Uniform arguments for write_r2c
!     rev. 9.9.60  Mar.  06/15  - NK: In sub: fixed call write_r2c for close condition 
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c            call write_r2c(100,100,mhtot+8784,0,frame_no5,0,8)  
            call write_r2c(100,100,int(totaltime),0,frame_no5,0,8)  
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            endif
            endif    ! iopt99
          endif   !  if(grdflg.eq.'y'
          
          if(iopt99)then
          if(frcflg.eq.'y'.and.numa.eq.0)then 
            if(allocated(outarray))deallocate(outarray)
            allocate(outarray(ycount,xcount),stat=iAllocate)
            if(iAllocate.ne.0)then
              print*,'Error with allocation of outarray in sub @ 594'
              STOP 'Program aborted in sub @ 595'
            endif
d            print*,'outarray allocated for '
d            print*,ycount,' rows'
d            print*,xcount,' cols'
            do n=1,naa
              i=yyy(n)
              j=xxx(n)
              l=nbasin(i,j)
              if(l.ne.0)then
                outarray(i,j)=dSTRconc2(n) !isoSTRconc2(n,l)
              else 
                outarray(i,j)=0.0
              endif
            end do     
            ycount_temp=ycount
            xcount_temp=xcount
            frame_no6=frame_no6+1
 !          unit=72  fln(72)- gridded flow - r2c file
!     rev. 9.9.60  Mar.  06/15  - NK: In sub: fixed call write_r2c for close condition 
!           changed mhtot+8784 to totaltime to make it work for monthly events
!           longer than one year
 !          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c           call write_r2c(200,200,mhtot+8784,0,frame_no6,0,8) 
            call write_r2c(200,200,8784000,0,frame_no6,0,8)  
 !          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          endif
          endif   !iopt99

          endif   !  if(jz.ne.iz.and.mod(jz,deltat_report).eq.0

          if(ensimflg.eq.'y'.and.jz.ne.iz)then
          if(ireport.ge.1.and.wfo_open_flg.eq.'y')then
!     rev. 9.2.20  Oct.  28/05  - NK: WFO_SPEC - reporting start & finish times 
            if(int(totaltime).gt.ireport_start
     *        .and.int(totaltime).le.ireport_end)then
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              call write_wfo(t,jz,iz)
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            endif
          endif
          endif

!>>>>>>>>>>>>>>>>>>>>





!>>>>>>>>>>>>>>>>>>>>
 
          do n=1,naa
!           SELECT THE MAXIMUM FLOW FOR EACH GRID DURING THE EVENT 
            qmax(n)=amax1(qmax(n),qo2(n))   
!           FIND THE MINIMUM SNOWC FOR ERROR MSG 
            do ii=1,classcount
              snowcmin(n,ii)=amin1(snowc(n,ii),snowcmin(n,ii))
            end do
          end do
   82     m=m+1
 
 
!     rev. 9.8.21  Jun.  18/12  - NK: Added swe observed date & report
!         extract the computed swe at snow course locations
          if(mod(jz,24).eq.0)then
            do i=1,nswe
              if(indomainflg(i))then
                n=s(sng_swe(i),ewg_swe(i))
c      print*,i,sng_swe(i),ewg_swe(i),n
                sum_swe=0.0
                class_sum=0.0
                do ii=1,classcount
!                 we have to not include the water class as it has 0 swe
                  if(ii.ne.ii_water)then
                    class_sum=class_sum+aclass(n,ii)
                    sum_swe=sum_swe+snowc(n,ii)*aclass(n,ii)*sca(n,ii)
                  endif
                end do
c                course_calc(i,ju)=sum_swe/class_sum
                course_calc(i,jz/24)=sum_swe/class_sum
              else
c                course_calc(i,ju)=-1.0            
                course_calc(i,jz/24)=-1.0            
              endif
              end do
c              write(777,*)ju,(course_calc(i,jz/24),i=1,nswe)
          endif


c          call swe_extract(ju)
 
!      May 1/11  NK 
!      save this little ditty. Used to create WSC file format for 
!      missing flows at any station
c          if(mod(jz,24).eq.0)then
c            if(qhyd(68,jz).lt.0.0)then
c              write(900,828)
c     *         '05PF999',year1,month_now,day_now,qsyn(68,jz)
c            else
c              write(900,829)
c     *         '02PF999',year1,month_now,day_now,qhyd(68,jz)
c            endif
c          endif
c828       format(a7,',1,',i4,',',i2,',',i2,',',f10.3,',c,')     
c829       format(a7,',1,',i4,',',i2,',',i2,',',f10.3,',o,')     

        if(iopt99)then
          call date_and_time(cday,ctime)
          read(ctime,11111)hhhhh
          read(ctime,11112)mmmmm
          read(ctime,11113)sssss
          now=hhhhh*3600.0+mmmmm*60.0+sssss
          output=output+now-last
          last=hhhhh*3600.0+mmmmm*60.0+sssss
        endif
          
!     rev. 10.1.78 Apr.  17/17  - NK: New s/r dds_UZS to calculate low flow penalty
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          if(dds_flag.gt.0)then
              if(errflg.eq.10.and.id.gt.idskip)then
                  call dds_uzs(score)
              endif
          endif
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          
!     rev. 10.2.46 Jan.  21/19  - NK: Added netCDF capability
!         write the CHARM_output.nc file for one time step
!         Assume this file is wanted if one of str, met or tem file is .nc       
!         Common 6 hr output time step - see Ben's e-mail Apr. 8/19      
!         Note: for netCDF (FEWS) applications, we use 6 hour time steps 
!               so qsum is set = 0 for that application          
!               For other uses, flow is summed for the whole run so we can
!               calculate runoff in mm for the whole upstream area at every grid.
!               This can help find discontinuities.
          do n=1,naa
                i=yyy(n)
                j=xxx(n)
                QQsum(n)=QQsum(n)+qo2(n)
          end do  
          if(netCDFflg)then
!             calculate the mean flow for the time step              
              deltat=6
              if(mod(int(totaltime),6).eq.0)then
                do n=1,naa
                  i=yyy(n)
                  j=xxx(n)
                  outarray(i,j)=QQsum(n)/6.0
                  QQsum(n)=0.0
                end do     
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
              call write_2D_flow(int(totaltime/deltat),mhtot/deltat,201)
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                do n=1,naa
                  i=yyy(n)
                  j=xxx(n)
                  outarray(i,j)=totsnw(n)
                end do     
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
               call write_2D_swe(int(totaltime/deltat),mhtot/deltat,202)
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                do n=1,naa
                  i=yyy(n)
                  j=xxx(n)
                  outarray(i,j)=totuzs(n)
                end do     
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
               call write_2D_uzs(int(totaltime/deltat),mhtot/deltat,203)
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!     rev. 10.2.56 June  13/19  - NK Added .nc output for grid_runoff & cumm ET
                do n=1,naa                  ! added June 13/19 NK
                  i=yyy(n)
                  j=xxx(n)
                  outarray(i,j)=qr(n)
                end do   
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
               call write_2D_grid_runoff(int(totaltime/deltat),
     *                                                mhtot/deltat,204)
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                do n=1,naa                  ! added June 13/19 NK
                  i=yyy(n)
                  j=xxx(n)
                  outarray(i,j)=eloss(n)
                end do    
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
               call write_2D_cumm_ET(int(totaltime/deltat),
     *                                                mhtot/deltat,205)
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
            endif
          endif
          
          

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
        end do     !while(time.lt.float(mhtot)
!         END TIME LOOP
!         --------------------------------------------------
!         --------------------------------------------------
!         --------------------------------------------------
!         --------------------------------------------------
!         --------------------------------------------------
!         --------------------------------------------------
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

        if(numa.eq.0.and.dds_flag.eq.0)print*

   83   CONTINUE       ! FROM STOP COMMAND

!       rev. 9.8.21  Jun.  18/12  - NK: Assess swe observed date & report
!       rev. 9.8.44  Jan.  31/13  - NK: fixed bug in sub.f : uninitialized course_calc(n,j)
!       rev. 9.8.92  Nov.  06/13  - NK: Changed output file swe.txt to swe.csv
!       write the file only if snow course data is availabel in file 951
        if(courseflg)then
          do j=1,mhtot/24
          if(iopt99)write(951,91100)
     *                (course_obs(n,j),course_calc(n,j),n=1,nswe)
91100       format(f10.2,<2*nswe-1>(',',f10.2))
          end do
!         calculate the error
          do n=1,nswe
            if(indomainflg(n))then
              do j=1,mhtot/24
                if(course_obs(n,j).ge.0.0)then
                  swe_error=swe_error
     *                  +(course_obs(n,j)-course_calc(n,j))**2
                  no_swe_obs=no_swe_obs+1
                endif
              end do    
            endif
          end do    
          swe_penalty=swe_error/no_swe_obs 
          print*,'>>NEW<< swe_penalty =',swe_penalty 
        endif

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(iopt.eq.2)then
          print*,' In sub before call lst'
          print*
        endif

!       Although it would be nice to write the hydrographs whenever a new.par
!       file is written, the data for the first  
!       ni-1 events is not available because the spl.csv file is
!       appended with each event and the data is lost. So it would have to be 
!       put into memory (which is also not attractive) so it could be written
!       at the time each new.par is written.

        mhtot=mhrd    
        if(.not.netCDFflg.or.dds.eq.0.or.iopt99)then
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call lst(scale,igrdshft)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        else
          print*,'Note'
          print*,'Supressing normal output in lst.f'
          print*,'---------------------------------'
        endif

!       STATUS LINE:
        if(icase.gt.0.and.numa.eq.0.and.dds_flag.eq.0)then
          write(6,5003)id,ni,mz,mhtot,nnn,maxn,scale,optlow
        endif

        if(mhtot.gt.8000)then
!         ONLY FOR ANNUAL RUNS
!         This is to apply a penaly if the snow does not melt each year. 
!         In the high mountains this can happen if the base temp is too high
!         from dds - the program will just accumulate swe if there's too much
!         water. Instead, this forces more evaporation to get rid of it.            
          minsnwflg=0
          do n=1,naa
          do ii=1,classcount
            if(aclass(n,ii).gt.0.0)then
              if(snowcmin(n,ii).gt.0.00001)then
              if(numa.le.0)then
                write(98,6400)n,ii,snowcmin(n,ii),id
              endif
              minsnwflg=1
              endif
            endif
          end do
          end do
        endif   !  mhtot.gt.8000
        if(minsnwflg.eq.1)then
          print*,'Warning: check spl.err file for snow problems'
          dds_penalty=10.0
        endif
        if(mod(id,12).eq.0.and.iopt99)then
           write(29,*)id,domain_precip/domain_area
           domain_precip=0.0
        endif
       
!     rev. 10.1.77 Apr.  17/17  - NK: Moved DDS err calcs to new dds_code s/r's
        if(abs(dds_flag).eq.1.and.id.gt.idskip)then
c            print*,trcflg,penalty,score
c            pause 6565
          penalty=penalty+score
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          call dds_options(dds_error)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          if(errflg.eq.10)dds_error=dds_error+score
        endif           ! if(abs(dds_flag).eq.1.and.......)
        
!       for insurance:	
	  close(unit=40,status='keep',iostat=ios)
	  close(unit=45,status='keep',iostat=ios)
        
        
        if(tbcflg.eq.'y')then
!         THIS SECTION IS USED TO WRITE ALL STATE VARIABLES
!         SO PRGRAM CAN BE CONTINUED WHERE IT LEFT OFF
!         Called only after the last event is complete.
!         REV. 8.1  - Feb.    15/96 -  TBC & RSM (TO BE CONTINUED & RESUME)
!         rev. 10.2.47 Feb.  10/19  - NK: Revised write_resume
!         soil_init, flow_init & lake_level_init moved there          
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call write_resume(jz)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        endif    !   if(tbcflg.eq.'y')....

      end do     ! EVENT LOOP END
 
!     * * * * * * * *  EVENT LOOP END  * * * * * * * * * *
!     * * * * * * * *  EVENT LOOP END  * * * * * * * * * *
!     * * * * * * * *  EVENT LOOP END  * * * * * * * * * *
!     **************************************************** 


!     rev. 10.2.46 Jan.  21/19  - NK: Added netCDF capability
!     Close the CHARM_output.nc file with nrecs = -1
      if(netCDFflg)then
!         unit=???  fln(201)-Gridded output for FEWS       CHARM_output.nc
!         close the files          
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          call write_2D_flow(int(totaltime),-int(totaltime),201)
          call write_2D_swe(int(totaltime),-int(totaltime),202)
          call write_2D_uzs(int(totaltime),-int(totaltime),203)
          call write_2D_grid_runoff(int(totaltime),-int(totaltime),204)
          call write_2D_cumm_ET(int(totaltime),-int(totaltime),205)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
          deltat=6
          fln(99)='results\CHARM_flow_vector.nc'
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call write_ts_flow(no,mhtot/deltat,99) 
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          
          fln(99)='results\CHARM_lake_levels.nc'
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call write_ts_lake_levels(noresv,mhtot/deltat,99) 
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          
          fln(99)='results\CHARM_lake_inflow.nc'
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call write_ts_lake_inflow(noresv,mhtot/deltat,99) 
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          
          fln(99)='results\CHARM_lake_outflow.nc'
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call write_ts_lake_outflow(noresv,mhtot/deltat,99) 
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      endif


!     rev. 9.8.08  Nov.  18/11  - NK: added fratio for interception hight optimization
!     set interception capacity to values read in originally in read_par_parser
      if(fratioflg)then
        do i=1,12
          do ii=1,classcount
            h(i,ii)=h(i,ii)/fratio(ii)
          end do
        end do
      endif

!     rev. 9.7.19  Jan.  18/11  - NK: Added sensitivity analysis
!       for sensitivity analysis and pattern search
c      if(dds_flag.eq.-1)then

      if(dds_flag.lt.0)then
	      optim=dds_error
	      if(optim.le.0.0)then
	        print*,'For sensitivity run:'
	        print*,'calculated error = zero'
	        print*,'Please make sure there is streamflow and/or'
	        print*,'value1 in the str files are set to 1'
	        stop 'Program aborted in sub @ line 3696'
	      endif
	      return                     !<<<<<<<<<<<<<<<<<<<<<return
      endif

!     rev. 10.2.46 Jan.  21/19  - NK: Added netCDF capability
!     stats not caculated for FEWS      
      if(numa.eq.0.and.iopt.lt.99.and.FLtype(6).eq.'tb0')then
!       call stats only when there's no pre-emption:
!       but not for the pattern search as it get here with every trial
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call stats(28,28)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      endif

!     rev. 9.7.00  May.  26/10  - NK: dds with pre-emption
!     rev. 9.6.01  Mar.  01/10  - NK: DDS capability added

c      if(errflg.eq.7)dds_error=optim  !  for nash obj fn only

!     write the function value for DDS
      if(dds_flag.eq.1)then
        open(unit=99,file='dds\function_out.txt',
     *          status='unknown',iostat=ios)
        if(ios.ne.0)then    ! added Nov. 10/14  nk
          print*
          print*,'Unable to open file  dds\function_out.txt'
          print*,'Possible cause(s):'
          print*,'file in use by another application'
          print*,'or target directory does not exist'
          stop 'Program aborted in sub.f @ 4299'
        endif
c       write(99,*)temp_ashnum
        write(99,*)dds_error
        close(unit=99,status='keep')
        write(30,*)  !blank line between trials in dds_log.txt
         
!       this copied here from lst - error is needed each time an improvement 
!       is found for the dds\best directory
!       fix: this could be eliminated here by moving the call lst after the ddserror
!       calculation above        
c        if(dds_error.lt.pre_emption_value)then
!     rev. 9.3.07  Dec.  29/06  - NK: added error field for whole domain
!     rev. 9.8.15  Mar.  12/11  - NK: write error.txt for every dds evaluation
!         (because it's needed whe doing on-the-fly validation runs)
          author='watflood                                '
          name='Volume errors in %                      '
          coordsys_temp=coordsys1
!         GreenKenue uses LatLong - code below uses LATLONG
          if(coordsys_temp.eq.'LatLong   ')coordsys_temp='LATLONG   '
          if(coordsys_temp.eq.'Cartesian ')coordsys_temp='CARTESIAN '
          zone_temp=zone1
	    datum_temp=datum1
	    xorigin_temp=xorigin
	    yorigin_temp=yorigin
	    xcount_temp=xcount
	    ycount_temp=ycount
	    xdelta_temp=xdelta
	    ydelta_temp=ydelta
	    startdate='unknown   '
	    starttime='unknown   '
          unit_conversion=1.0
          attribute_count=1
	    attribute_name='error                                    '
	    attribute_units='percent                                 ' 
          source_file_name='last spl run'     
          do j=1,xcount
	      do i=1,ycount
              outarray(i,j)=basinerr(i,j)
	      end do
  	    end do
        
!         write the header
!         written in lst for ddsflg=0
!     rev. 9.9.23  Aug.  10/14  - NK: Uniform arguments for write_r2c
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c        call write_r2c(67,67,1,0,0,0,11)   
        call write_r2c(67,67,0,0,0,0,11)   
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       write the data
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c        call write_r2c(67,67,1,1,1,1,11)   
        call write_r2c(67,67,0,0,0,1,11)   
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c        endif
      endif

!     STATUS LINE:   PS optimization only
      if(numa.gt.0.and.dds_flag.eq.0)then
        if(nnn.eq.0)optim=optim*4.0
        write(6,5006)id-1,ni,mz,mhtot,nnn,maxn,optlow
      endif

      close(unit=271,status='keep')   ! gridflow.r2c

      id=ni

!     rev. 9.5.01  Oct.  15/07  - NK: added wetland continuity check
      
      if(iopt.ge.1)then
        write(51,*)'Continuity check for wetlands:'
        write(51,*)
     *'grid #,    net_inflow,    delta_storage,      difference,     %'
      if(wetflg.eq.'y')then
        do n=1,naa
          xtemp=qiwetsum(n)-qowetsum(n)
          ytemp=wstore2(n)-wstoreinit(n)
          ztemp=(xtemp-ytemp)/xtemp
          if(ztemp.gt.0.0001)then
!            write(51,5103)n,xtemp,ytemp,ztemp*100.0
            write(51,*)n,xtemp,ytemp,xtemp-ytemp,ztemp*100.0
          endif
          end do
        endif
5103    format(i10,2f20.3)
      endif

!     rev. 9.6.05  Apr.  06/10  - NK: added store_error_flag for -ve storage grids
      errflg_store=.false.
c      if(numa.eq.0.and.dds_flag.eq.0)then
        do n=1,naa
          if(store_error_flag(n))then
            Print*,'WARNING: -ve storage(s) in grid no ',n
            errflg_store='true'
          endif
        end do
c     endif
      if(errflg_store)then
        print*,'-ve storage errors found in the grids listed above'
        print*,'Probable cause:'
        print*,'1. steep slopes with low Manning n'
        print*,'2. grids with 100% water not marked as a reach'
        print*,'3. cummulative evaporation exceeds inflow + precip'
        print*,'Storages are set to 0.0 when -ve'
        print*,'and you get this annoying message'
        print*
      endif


c      if(tbcflg.eq.'y')then
!       THIS SECTION IS USED TO WRITE ALL STATE VARIABLES
!       SO PRGRAM CAN BE CONTINUED WHERE IT LEFT OFF
!       Called only after the last event is complete.
!     REV. 8.1  - Feb.    15/96 -  TBC & RSM (TO BE CONTINUED & RESUME)
!     rev. 10.2.47 Feb.  10/19  - NK: Revised write_resume
!     soil_init, flow_init & lake_level_init moved there          
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c        call write_resume(jz)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c      endif    !   if(tbcflg.eq.'y')....

!     end writing the resume file

!     CALL RUNOF5 TO PRINT SUMMARY WATERBALANCE FILE WATBAL.TXT

      oldjan=jan
      if(iopt.eq.2)then
        print*,' In sub before call runof6(3,...'
        print*
      endif

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call runof6(3,time,sintvl,aintvl,mon,e1,mz,ju)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

      jan=oldjan
      
c      if(iopt.ge.1)then
      if(writeFlg(1))then
      open(unit=99,file='rff_grid_specs.txt',status='unknown',
     *        iostat=ios)
          write(99,*)'Frac ,',frac(nnprint)
          write(99,*)'GridArea,',grid_area(nnprint)
       
      
        do ii=1,classcount
          if(aclass(nnprint,ii).gt.0.0)then
            write(99,*)
            write(99,*)'Class ',ii,nclass(ii)
            write(99,*)'ClassRatio,',aclass(nnprint,ii)
            write(99,*)'ClassArea,',
     *                grid_area(nnprint)*aclass(nnprint,ii)
            write(99,*)'ClassFraction,',
     *                frac(nnprint)*aclass(nnprint,ii)
          endif
        end do
      endif
      
!     if(iopt.eq.2)print*,' In sub before call snout3'
!
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     if(iopt.ge.1)call snout3()
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     rev. 9.9.39  Nov.  14/14  - NK: Modifications for watroute
      if(ittoflg.eq.0.and.vapflg.eq.'y'.and.modelflg.eq.'n')then
!       only needed if we have evaporation turned on
        write(98,*)'WARNING: In .par file, temp3 set too low'
        write(98,*)'Results in underestimated evaporation'
        write(98,*)'Please see manual section 2.4.2'
      endif

      if(climate_delta_flg)then
        write(98,*)'Reminder:'
	  write(98,*)'Montly climate deltas from file:'
	  write(98,*)'basin\monthly_climate_deltas.txt'
	  write(98,*)'Temperatures were adjusted by dC:'
	  write(98,*)monthly_temperature_delta
	  write(98,*)'Precipition were multiplied by:'
	  write(98,*)monthly_precipitation_delta
      endif

      firstpass=.false.

      if(msgflg.and..not.netCDFflg)then
        write(98,*)'WARNING:'
        write(98,*)'Tracer model data may not reflect'
        write(98,*)'the actual flow components at locations downstream'
        write(98,*)'of stations with nudging'
      endif

      RETURN

! FORMATS

  201 format(999f5.2)
 2000 format(5i5)
 2001 format(60i1)
 2002 format(a1)
 2003 format(a20,i50)
 3004 format(a20,a10,2x,a10)
 5001 format('+',1x,'id=',i3,'/',i3,' mz=',i3,'/',i3,' nnn=',i5,'/',i5,
     *         ' smc=',5f4.2,' optlow=',e16.10)
 5002 format(i5,f5.2,i5,f5.1)
 5003 format('+',1x,'id=',i3,'/',i3,' mz=',i5,'/',i5,' nnn=',i5,'/',i5,
     *         ' scale=',f4.2,' optlow=',e16.10)

 5004 format(a20,a10)
 5005 format(a20,i50)
 5006 format('+',1x,'id=',i3,'/',i3,' mz=',i5,'/',i5,' nnn=',i5,'/',i5,
     *         ' optlow=',e16.10)

 5010 format(' id,l,mhtot,kt,nu,delta,opt/',5i5,2e10.3)
 5101 format('new shed file read in')
 5102 format(' id = ',i5,' filename = ',a30)
 5558 format(' n,nnx,qda(n),qda(nnx),da(n),da(nnx)/',2i3,4f10.1)
 6001 format(5x,' no of storms =',i5)
 6002 format(35x,' * * * time =',f7.2,' hours * * *',f12.3)
 6004 format('+',i5,' hrs rainfall available, ',i5,' hrs used')
 6005 format(' Warning: no streamflow stations selected for '/
     *' error calculation. Enter data in 1st line of .str file.')
 6007 format(' no,nl,mhtot,kt/',4i5)
 6008 format(' id,nr,tj1,mo,conv/',2i5,f5.2,i5,f5.2)
c 6012 format(' ',2f12.3,i5,1x,a12,3x,2a1,3x,2a1)
 6013 format(' ',2f12.3,i5,1x,a12,99f5.2)
 6014 format(' ',2i5,99f5.2)
 6015 format(' ','sub basin file'/
     *   '      sub-basin percent of land covers')
 6016 format(' ','sub basin file'/
     *   '      can be used to find sub-basin percent of land covers')
 6017 format('           xx          yy    # name         sum   ',
     *                        ' classes 1-',i5)
 6018 format(' ',i7,99f7.2)
 6100 format(' [in sub]'/
     *'   n    i    j      da     qda  qbase')
 6102 format(' ',3i5,3f10.2)
 6103 format(' element #',i5,' base flow = 1 cms assumed')
 6104 format(25i3)
 6105 format(' initial flows modified by reservoir releases') 
 6106 format(' i,qinit(i)/',i5,f10.3)
 6107 format(' n,yyy(n),xxx(n),i,ires(i),jres(i)/',6i3)
 6108 format(' initial base flows prorated upstream')
 6109 format(' reservoir releases added back in') 
 6110 format(' initial flow at streamflow stations')
 6120 format(a1)
 6121 format(5i10,2f10.2)
 6226 format(' error encountered opening unit=37 fln=',a31/)
 6227 format(' error encountered on unit=37 fln=',a31/)
 6229 format(' error encountered on unit=38 fln=',a31/)
 6231 format(' error encountered on unit=40 fln=',a31,'on line 5'/)
 6233 format(' error: problems opening unit=45 fln=',a31/
     *    '     OK if vapflg=y and flgevp.le.1')
 6235 format(' end of file encountered on unit=36 fln=',a31/) 
 6237 format(' end of file encountered on unit=40 fln=',a31/)
 6250 format(' error encountered on unit=99 fln=',a31/)
 6263 format(' error encountered on unit=49 fln=',a31/)
 6273 format(' error encountered on unit=50 fln=',a31/)
 6300 format(' ','Mode is set to check all input files'/
     *' iopt = ',i5, ' To run the event, change debug level to'/
     *' a value less than 99'/)
 6400 format(' Warning: min swe on (grid,class)',i4,' ',i2,
     *    ' =',f8.2,'mm in event #',i3)
 6501 format(' number of grids,',i5)
 6502 format(' number of rows,',i5)
 6503 format(' number of colums,',i5)
 6504 format(' number of land cover classes = classcount,',i5) 
 6505 format(' The reference grid is the bottom left hand corner'/
     *       ' there are blanks around the watershed'/
     *       ' All units are mm or cms')
 6508 format(' time=',f14.2)
 6509 format(/' Colunm Key:'/ 
     *'      grid no, row no, column no, precip'/ 
     *i5,' liquid surface storage bare/snow (2*classcount)'/
     *i5,' snow water equivalent (classcount)'/
     *i5,' snow covered area (classcount)'/
     *i5,' upper zone storage bare/snow (2*classcount)'/
     *' lower zone storage'/
     *' grid flow, channel flow') 
 6510 format(3i5,60e10.2)

 9001 format(f25.0)
 9700 format(19e14.7)
 9701 format(a3)
 9702 format(' resume.txt file length does not match current'/
     *       ' memory array spec. Regenerate the file'/)
 9703 format('block no.',i2)
 9704 format(9x,i2)
 9706 format(16i5)
 9999 format(15x,f16.6)
 9990 format(i5)
 9991 format(12x,f9.6)
 9992 format(12x,<ninit>(i9))
 9993 format(12x,<ninit>(f9.6))
 9994 format(5x,i9)
 9995 format(7x,f9.6)
99000 format(f5.1)
99001 format(f25.0)
99002 format(2f15.3,i6,4x,2a12,f12.1)
99003 format(i1,5x,a50)
99004 format(i5)
99005 format(i10,<noresv>g12.3)
99182   format(' Warning: Error opening or reading fln:',a30/
     *  ' Probable cause: missing strfw/yymmdd.str input file'/
     *  ' OR: in config.sys have you set files=100 & buffers=50?'/)
98001 format(' Warning: Error reading resume.txt'/
     *' Probable cause: wrong format or end of file reached'/
     *' This could be due to a change in grid characteristics'/
     *' since creating the resume file'/
     *' Solution: create a new resume file with current executable'/
     *' block= ',i5)

      end  SUBROUTINE sub

