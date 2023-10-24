      SUBROUTINE read_par_parser(unitNum,flnNum)

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

!    You should have recieved a copy of the GNU Lesser General Public License
!    along with WATFLOOD.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!***********************************************************************
!
!     rev. 9.8.05  Oct.  18/11  - NK: New read_par_parser subroutine
!     rev. 9.8.08  Nov.  18/11  - NK: added fratio for interception hight optimization
!     rev. 9.8.68  Jun   17/13  - NK: Added dds_override file
!     rev. 9.8.84  Sep.  15/13  - NK: Added fratio to list of equal values for bog & fen
!     rev. 9.9.03  Dec.  15/13  - NK: Change to gridded latitude for etharg
!     rev. 10.5.01 Jan.  06/23  = NK Added check on rlake pars to avoid instability
!     rev. 10.5.09 May   11/23  = NK Added fpFactor - flood plain width/depth .17 = 100/1
!     rev. 10.5.15 June  11/23  = NK Added sdcd, r1n, aa2, aa3, widep, pool to Ostin.txt list
!     rev. 10.5.15 June  11/23  = NK Added sdcd, r1n, aa2, aa3, widep, pool to Ostin.txt list
!
!  version numbers are added for version 7.0 and later. this allows
!  parameter files to backward compatible when more parameters are 
!  added to the file in the future. in version 7.0 some headers are 
!  added as well for readability.
!
!
!  iprtflg - if eq 2 write a new parameter file after optimization
!  r1      - the roughness of the flood plain
!  r2      - the roughness of the largest part of the channel
!  r1      - factor for raising r2 ie if r1=2 then f.p. roughness
!            is 2 times channel roughness
!  zz      - is an exponent in calculating the mannings  n
!  h       - crop height
!  ix      - exponent in defining porosity
!  e1       - void ratio
!  ak      - permeability
!  sm      -soil moisture(average for month)
!  nbsn    - no of basins with different r2 to be optimized
!          - must be smaller than nrvr
!  classcount_local  - no of classes to be optimized max=4
!  errflg  - picks the type of error calculation 
!  a1 .. a12 parameters variously used mostly in runoff
!  a7      - weighting factor for old vs. new sca value default=.90
!  a8      - time offset to check for min temp on rdtemp
!  a9      - heat deficit to swe ratio  default=0.333
!  a10     - power on interflow discharge function default = 1.0
!  a11     - equivalent vegetation height for bare ground evaporation
!  a12     - min precip rate for smearing used in rain.for default=0.0
!
!***********************************************************************

      use area_watflood
      use area_debug
      implicit none

      DIMENSION     :: ajunk(17)

      CHARACTER(10) :: time
      CHARACTER(8)  :: cday
      CHARACTER(128):: qstr
      CHARACTER(60) :: junk
      CHARACTER(30) :: filename1
      CHARACTER(1)  :: errorflg,answer,linetype,tempflg
      INTEGER(kind=2) :: result1,nnts
      INTEGER       :: n,classcount_local
      INTEGER       :: nrvr1,nchr,iprtflg,linecount,i,unitNum,flnNum,&
                       ios,iverflg,ix,iallocate,ijunko,ii,j,numb
      INTEGER :: GlobalParLine ,EndGlobalParLine
      INTEGER :: RoutingParLine,EndRoutingParLine
      INTEGER :: HydrologicalParLine,EndHydrologicalParLine
      INTEGER :: SnowParLine,EndSnowParLine
      INTEGER :: InterceptionParLine,EndInterceptionParLine
      INTEGER :: MonthlyEvapParLine,EndMonthlyEvapParLine
      INTEGER :: OptimizationParLine,EndOptimizationParLine
      INTEGER :: APILimitsLine,EndAPILimitsLine
      INTEGER :: HydrologicalParLimitsLine,EndHydrologicalParLimitsLine
      INTEGER :: GlobalSnowParLimitsLine,EndGlobalSnowParLimitsLine
      INTEGER :: SnowParLimitsLine,EndSnowParLimitsLine
      INTEGER :: RoutingParLimitsLine,EndRoutingParLimitsLine
      INTEGER :: GlobalParLimitsLine,EndGlobalParLimitsLine
      real(4)  :: hmax,e1,ajunk,kcondflg,xxx1
      real(4)  :: lat_origin,lat_delta,lat_max
    character(1)  :: firstpass
    character(512):: line(1000)
    logical       :: exists,fratioflg1


!  DATA ntest/-10906/qstr/'param'/nchr/5/


      data firstpass/'y'/

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

!  i=max(nrvr,classcount)
!  allocate(ajunk(i),stat=iAllocate)
!  if(iAllocate.ne.0) STOP 'Error with allocation of ajunk in rdpar'

      if(firstpass.eq.'y')then

!     TS - ALLOCATIONS OF AREA4A ARRAYS
!     classcount for the number of land cover classes
!     nrivertype for the number of channel or basin types
!     moved here from spl9  nk 06/07/00
!     then moved from rdpar nk 28/12/04
!     moved back from rdshed 27/07/06 because needed for bsn.for  nk

      if(.NOT.allocated(rivtype))then
!     these are not needed for watroute 
        allocate(rivtype(nrvr),stat=iAllocate)
        if(iAllocate.ne.0)then 
            write(98,*)'Error: with allocation of rivtype in rdpar'
            write(98,*)'Error: program aborted in read_par_parser @ 122'
            STOP 'Error with allocation of rivtype in rdpar'
        endif
      endif

!     TS - ALLOCATIONS OF AREA4A ARRAYS
!     classcount for the number of land cover classes
!     nbsn for the number of channel or basin types
!     moved here from rdshed 27/07/06 because needed for bsn.for  nk
!     rev. 9.5.27  Apr.  15/08  - NK: fixed allocation for chnl in rdpar
      if(.NOT.allocated(ds))then
      allocate(ds(classcount),dsfs(classcount),chnl(5),&  ! note: chnl = always 5
       r3(classcount),r4(classcount),r3fs(classcount),rec(classcount),&
       ak(classcount),akfs(classcount),&
       r3low(classcount),r3fslow(classcount),reclow(classcount),&
       aklow(classcount),akfslow(classcount),ak2fslow(classcount),&
       r3hgh(classcount),r3fshgh(classcount),rechgh(classcount),&
       akhgh(classcount),&
       akfshgh(classcount),ak2fshgh(classcount),r3dlt(classcount),&
       r3fsdlt(classcount),recdlt(classcount),akdlt(classcount),&
       akfsdlt(classcount),ak2fsdlt(classcount),retn(classcount),&
        retnlow(classcount),ak2low(classcount),&
       retnhgh(classcount),ak2hgh(classcount),&
       retndlt(classcount),ak2dlt(classcount),&
       retfs(classcount),ak2fs(classcount),fpet(classcount),&
       fpetdlt(classcount),fpetlow(classcount),fpethgh(classcount),&
       ftall(classcount),ftalldlt(classcount),&
       ftalllow(classcount),ftallhgh(classcount),&
       fratio(classcount),fratiodlt(classcount),&
       fratiolow(classcount),fratiohgh(classcount),&
       ak2(classcount),stat=iAllocate)      
        if(iAllocate.ne.0)then
            write(*,*)'Error with allocation hydrological arrays @ 145'
            stop
        endif
      allocate(&
!       pot(classcount),potfs(classcount),&   ! moved from read_shed_ef may 15/07 nk
       iiclass(classcount*2),h(12,classcount),fpetmo(12,classcount),&
       nclass(classcount),stat=iAllocate)
        if(iAllocate.ne.0)then
            write(98,*)'Error: with allocation of hydrological arrays @ 157'
            write(98,*)'Error: program aborted in read_par_parser @ 158'
            write(*,*)'Error with allocation hydrological arrays @ 159'
            STOP 'Error: program aborted in read_par_parser @ 159'
        endif

        ParFileCommentNumber=0

    endif
!     TS - ALLOCATIONS OF AREAOPTSA ARRAYS
      if(.NOT.allocated(fmdlt))then
        allocate(fmdlt(classcount),fmlow(classcount),fmhgh(classcount),fmndlt(classcount),&
       fmnlow(classcount),fmnhgh(classcount),uajdlt(classcount),uajlow(classcount),&
       uajhgh(classcount),mbsdlt(classcount),mbslow(classcount),mbshgh(classcount),&
       basdlt(classcount),baslow(classcount),bashgh(classcount),&
       subdlt(classcount),sublow(classcount),subhgh(classcount),&
       sdcddlt(classcount),sdcdlow(classcount),sdcdhgh(classcount),stat=iAllocate)
        if(iAllocate.ne.0)then
            write(98,*)'Error:with allocation of areaoptsa arrays @ 174'
            write(98,*)'Error: program aborted in read_par_parser @ 175'
            STOP  'error with allocation of areaoptsa arrays @ 196'
        endif
      endif

!     TS - ALLOCATIONS OF AREAMELTA ARRAYS (PARTIAL)
!     SNW,DSN,TTEMP,TMX,TMN,EL ALLOCATED IN SHEDA.FOR
!     SDCD,SDCSCA ALLOCATED IN RDSDCA.FOR
!     allocate(snowc(na,classcount),dsnow(na),tempv(na),tempvmin(na),
!     moved here from rdshed 27/07/06 because needed for bsn.for  nk
      if(.NOT.allocated(base))then
        allocate(base(classcount),fm(classcount),fmn(classcount),&
       whcl(classcount),tipm(classcount),uadj(classcount),&
       rho(classcount),stat=iAllocate)
        if(iAllocate.ne.0)then
            write(98,*)'Error: with allocation of areamelta arrays in read_par'
            write(98,*)'Error: program aborted in read_par_parser @ 191'
            STOP  'Warning: error with allocation of areamelta arrays in read_par'
        endif
      endif

!     TS - ALLOCATIONS OF AREAETA ARRAYS (PARTIAL)
!     RAD ALLOCATED IN SHEDA.FOR
!     TS - ADDED ALLOCATIONS FOR EVAP-SEPARATION PARAMS (22/03/06)
!     TS: CHANGED ALLOCATIONS OF alb,pet, evap TO classcount (27/03/06)
!     rev. 9.1.80  Mar.  31/05  - NK: added sublimation   (sublim)
!     moved here from rdshed 27/07/06 because needed for bsn.for  nk
      if(.NOT.allocated(evap))then
        allocate(evap(classcount,12),sublim_factor(classcount),&
       sublim_rate(classcount),flint(classcount),&
       diff(12),hu(12),ffcap(classcount),fcap(classcount),&
       spore(classcount),alb(classcount),pres(12),stat=iAllocate)
 !      *rh(na),stat=iAllocate)
         if(iAllocate.ne.0)then
             write(98,*)'Error: error with allocation of evap arrays in read_par'
             write(98,*)'Error: program aborted in read_par_parser @ 209'
             STOP 'Warning: error with allocation of evap arrays in read_par'
         endif
      endif



!    allocate(nsdc(classcount),snocap(classcount),idump(classcount),
! *        stat=iAllocate)
!    if(iAllocate.ne.0) STOP
! *   'Error with allocation of areamelta arrays in spl9 @ 158'

      print*,allocated(r2nlow)
      
      if(.not.allocated(r2nlow))then
        allocate(r2nlow(nrvr),r2nhgh(nrvr),r2ndlt(nrvr),    &    
              thetadlt(nrvr),thetalow(nrvr),thetahgh(nrvr),&
              widepdlt(nrvr),wideplow(nrvr),widephgh(nrvr),&
              kconddlt(nrvr),kcondlow(nrvr),kcondhgh(nrvr),&
              flzlow(nrvr),flzhgh(nrvr),flzdlt(nrvr),&
              pwrlow(nrvr),pwrhgh(nrvr),pwrdlt(nrvr),&
              rlakelow(nrvr),rlakehgh(nrvr),rlakedlt(nrvr),&
              r1nlow(nrvr),r1nhgh(nrvr),r1ndlt(nrvr),    &
              aa2low(nrvr),aa2hgh(nrvr),aa2dlt(nrvr),    &
              aa3low(nrvr),aa3hgh(nrvr),aa3dlt(nrvr),    &
              poollow(nrvr),poolhgh(nrvr),pooldlt(nrvr),    &
              stat=iAllocate)
        if(iAllocate.ne.0)then
          write(98,*)'Error: with allocation of area4 arrays in rdpar @ 524'
          write(98,*)'Error: Try setting numa = 0 in the par file'
          write(98,*)'Error: program aborted in rdpar @ 626'
          print*,'Error with allocation of area4 arrays in rdpar @ 524'
          print*,'Try setting numa = 0 in the par file'
          print*
        if(iAllocate.ne.0) write(98,*)'Error: program aborted in rdpar @ 626'
          write(98,*)'Error: program aborted in read_par_parser @ 240'
          stop 'program aborted in rdpar @ 626'
        endif
    endif

!  endif

      if(.not.allocated(flz_o))then
        allocate(flz_o(nrvr),pwr_o(nrvr),r1n_o(nrvr),&        
              r2n_o(nrvr),mndr_o(nrvr),aa2_o(nrvr),&
              aa3_o(nrvr),aa4_o(nrvr),theta_o(nrvr),&
              widep_o(nrvr),kcond_o(nrvr),pool_o(nrvr),&
              rlake_o(nrvr),fpFactor_o(nrvr),stat=iAllocate)
        if(iAllocate.ne.0)then
            write(98,*)'Error: with allocation of area4 arrays in rdpar -2'
            write(98,*)'Error: program aborted in read_par_parser @ 255'
            STOP 'Error with allocation of area4 arrays in rdpar -2'
        endif
    endif

      if(.not.allocated(sdcd))then
!       nnts=nsdc(classcount)
        nnts=2
        allocate(sdcd(nnts,classcount),sdcsca(nnts,classcount),stat=iAllocate)
        if(iAllocate.ne.0)then
             write(98,*)'Error: with allocation of areamelta arrays in read_par'
             write(98,*)'Error: program aborted in read_par_parser @ 266'
             STOP 'Error with allocation of areamelta arrays in read_par' 
        endif
      endif

        if(debug_output)write(63,*)'Allocations done in rdpar',classcount,nrvr

!     rev. 9.8.08  Nov.  18/11  - NK: added fratio for interception hight optimization
      do ii=1,classcount
        fratio(ii)=1.0
        fratiodlt(ii)=-1.0
        fratiolow(ii)=0.1
        fratiohgh(ii)=10.0
      end do
!     rev. 10.4.30 Dec.  18/20  = NK Added rainsnowtemp parameter for rain/snow tmp > base tmp
       rainsnowtemp=-999.00 ! Default value for old par files use base(ii)

      do i=1,nrvr
!          fpFactor_o(i)=.17   ! default value for floodplain w/d = 100
          fpFactor_o(i)=100.   ! default value for floodplain w/d = 100
      end do
      
      sdcddlt=-1.0
      sdcdlow=0.0
      sdcdhgh=100.0
      r1ndlt=-1.0
      r1nlow=0.013
      r1nhgh=1.0
      aa2dlt=-1.0
      aa2low=1.0
      aa2hgh=2.0
      aa3dlt=-1.0
      aa3low=0.0
      aa3hgh=0.1
      widepdlt=-1.0
      wideplow=1.0
      widephgh=100.0
      pooldlt=-1.0
      poollow=0.0
      poolhgh=1.0
      
       
      endif  !firtspass

!     ADDED THIS FROM TODD'S ETPAR.FOR - FRANK S: NOV/97 
      alamb=2.478
      den=999.
      alpha=1.35
!      dds_flag=0 ! default for old par files
      
      
      type1=float(itype)
!  itype=int(type1)

! READ SNOW COVER PARAMETERS:   
!  - CALLED ON FIRST TIME STEP FOR SIMULATION
!  - CALLED TO WRITE TO newsnow.par DURING OPTIMIZATION 
!    IF OPTIMIZATION PARAMETERS ARE UPDATED
!               -- first the snow output flags
!               -- melt factors
!               -- now mbase
!               -- negative melt factors
!               -- wind function
!               -- ati decay
!               -- conversion snow density depth to we
!               -- liquid water holding capacity
!               -- daily ground melt
!---------------------------------------------------------

!     rev. 9.5.06  Feb.  05/08  - NK: added pool and pool_o in rdpar & route

!     rev. 9.9.35  Oct.  20/14  - NK: Added keyword & file checks
      inquire(FILE=fln(flnNum),EXIST=exists)
      if(.not.exists)then
        print*
        print*,'Fatal error'
        print*,'Attempting to open the par file'
        print*,'file name ',fln(flnNum)(1:40)
        print*,'but it is not found (in this location)'
        print*,'check for correct keyword and/or file name'
        print*
        write(98,*)'Error: Attempting to open the par file'
        write(98,*)'Error: file name ',fln(flnNum)(1:40)
        write(98,*)'Error: but it is not found (in this location)'
        write(98,*)'Error: check for correct keyword and/or file name'
        write(98,*)'Error: Program aborted in read_par_parser @ 326'
        stop 'Program aborted in read_par_parser @ 280'
      endif

      open(unitNum,file=fln(flnNum),status='old',iostat=ios)
      if(ios.ne.0)then
        print*,' Error opening',fln(flnNum)(1:40)
        print*,'on unit=',unitNum
        print*
        write(98,*)'Error: opening ',fln(flnNum)(1:50),'on unit=',unitNum
        write(98,*)'on unit=',unitNum	
        write(98,*)'Error: Program aborted in read_par_parser @ 289'
        stop ' Program aborted in read_par_parser @ 289'
      else
        write(98,*)'Info: Opened unit',unitNum,'fln=',fln(flnNum)(1:30)	
    endif

    numb=0  ! so new.par can be copied to bsnm.par for batch job
      
      
!     initialize the line counts
      linecount=0
          GlobalParLine=linecount
          EndGlobalParLine=linecount
          RoutingParLine=linecount
          EndRoutingParLine=linecount
          HydrologicalParLine=linecount
          EndHydrologicalParLine=linecount
          SnowParLine=linecount
          EndSnowParLine=linecount
          InterceptionParLine=linecount
          EndInterceptionParLine=linecount
          MonthlyEvapParLine=linecount
          EndMonthlyEvapParLine=linecount
          OptimizationParLine=linecount
          EndOptimizationParLine=linecount
          APILimitsLine=linecount
          EndAPILimitsLine=linecount
          HydrologicalParLimitsLine=linecount
          EndHydrologicalParLimitsLine=linecount
          GlobalSnowParLimitsLine=linecount
          EndGlobalSnowParLimitsLine=linecount
          SnowParLimitsLine=linecount
          EndSnowParLimitsLine=linecount
          RoutingParLimitsLine=linecount
          EndRoutingParLimitsLine=linecount
          GlobalParLimitsLine=linecount
          EndGlobalParLimitsLine=linecount


      do while(.not.eof(unitNum))
      ios=0
!      do while(ios.eq.0)
        linecount=linecount+1
        read(unitNum,51200,iostat=ios)line(linecount)
        if(debug_output)write(63,*)linecount,line(linecount)(1:50)
51200   format(a512)        
        if(line(linecount)(1:19).eq.':GlobalParameters')then
          GlobalParLine=linecount
        elseif(line(linecount)(1:20).eq.':EndGlobalParameters')then
          EndGlobalParLine=linecount
        elseif(line(linecount)(1:18).eq.':RoutingParameters')then
          RoutingParLine=linecount
        elseif(line(linecount)(1:21).eq.':EndRoutingParameters')then
          EndRoutingParLine=linecount
        elseif(line(linecount)(1:18).eq.':HydrologicalParam')then
          HydrologicalParLine=linecount
        elseif(line(linecount)(1:21).eq.':EndHydrologicalParam')then
          EndHydrologicalParLine=linecount
        elseif(line(linecount)(1:15).eq.':SnowParameters')then
          SnowParLine=linecount
        elseif(line(linecount)(1:18).eq.':EndSnowParameters')then
          EndSnowParLine=linecount
        elseif(line(linecount)(1:21).eq.':InterceptionCapacity')then
          InterceptionParLine=linecount
        elseif(line(linecount)(1:24).eq.':EndInterceptionCapacity')then
          EndInterceptionParLine=linecount
        elseif(line(linecount)(1:12).eq.':MonthlyEvap')then
          MonthlyEvapParLine=linecount
        elseif(line(linecount)(1:15).eq.':EndMonthlyEvap')then
          EndMonthlyEvapParLine=linecount
        elseif(line(linecount)(1:21).eq.':OptimizationSwitches')then
          OptimizationParLine=linecount
        elseif(line(linecount)(1:24).eq.':EndOptimizationSwitches')then
          EndOptimizationParLine=linecount
        elseif(line(linecount)(1:10).eq.':APILimits')then
          APILimitsLine=linecount
        elseif(line(linecount)(1:13).eq.':EndAPILimits')then
          EndAPILimitsLine=linecount
        elseif(line(linecount)(1:22).eq.':HydrologicalParLimits')then
          HydrologicalParLimitsLine=linecount
        elseif(line(linecount)(1:25).eq.':EndHydrologicalParLimits')then
          EndHydrologicalParLimitsLine=linecount
        elseif(line(linecount)(1:20).eq.':GlobalSnowParLimits')then
          GlobalSnowParLimitsLine=linecount
        elseif(line(linecount)(1:23).eq.':EndGlobalSnowParLimits')then
          EndGlobalSnowParLimitsLine=linecount
        elseif(line(linecount)(1:14).eq.':SnowParLimits')then
          SnowParLimitsLine=linecount
        elseif(line(linecount)(1:17).eq.':EndSnowParLimits')then
          EndSnowParLimitsLine=linecount
        elseif(line(linecount)(1:17).eq.':RoutingParLimits')then
          RoutingParLimitsLine=linecount
        elseif(line(linecount)(1:20).eq.':EndRoutingParLimits')then
          EndRoutingParLimitsLine=linecount
        elseif(line(linecount)(1:16).eq.':GlobalParLimits')then
          GlobalParLimitsLine=linecount
        elseif(line(linecount)(1:19).eq.':EndGlobalParLimits')then
          EndGlobalParLimitsLine=linecount
        endif
      end do

    rewind (unit=unitNum)
      if(debug_output)write(63,*)'GlobalParLine - line',GlobalParLine
      do i=1,GlobalParLine
        read(line(i),51200)line(i)
      if(debug_output)write(63,*)line(i)(1:50)
!     REV. 10.1.29 May   04/16  - NK: Added parfile comments
        if(line(i)(1:1).eq.'#')then
          if(ParFileCommentNumber.le.99)then
!           read a commnet line and save as numberd comment          
            ParFileCommentNumber=ParFileCommentNumber+1
            parfilecomment(ParFileCommentNumber)=line(i)(1:72)
      if(debug_output)write(63,*)parfilecomment(ParFileCommentNumber)(1:72)
          else
            print*,'max number of comments in the parfile '
            print*,'GlobalParameter section exceeds 100'
            print*,'Rest of comments are not passed on'
            write(*,'(''paused, type [enter] to continue'')') 
            read (*,*) 
          endif
        endif
      end do

      if(debug_output)write(63,*)!     loop through all the lines:
      do i=GlobalParLine+1,EndGlobalParLine-1
          if(line(i)(1:5).eq.':iopt')then
            read(line(i),*,iostat=ios)junk,iopt         !# debug level
      if(debug_output)write(63,*)'iopt=',iopt
          elseif(line(i)(1:6).eq.':itype')then
            read(line(i),*)junk,itype          !# channel type - floodplain/no 
              if(debug_output)write(63,*)'itype=',itype
          elseif(line(i)(1:7).eq.':itrace')then
            read(line(i),*)junk,itrace       !# Tracer choice
              if(debug_output)write(63,*)'itrace=',itrace
          elseif(line(i)(1:4).eq.':a1,')then
            read(line(i),*)junk,a1            !# ice cover weighting factor
              if(debug_output)write(63,*)'a1=',a1
          elseif(line(i)(1:4).eq.':a2,')then
!     rev. 10.4.23 Apr.  30/20  = NK repurposed a2 for swe threshold"
            read(line(i),*)junk,a2            !# swe correction threshold
              if(debug_output)write(63,*)'a2=',a2
          elseif(line(i)(1:4).eq.':a3,')then
            read(line(i),*)junk,a3            !# error penalty coefficient
              if(debug_output)write(63,*)a3
          elseif(line(i)(1:4).eq.':a4,')then
            read(line(i),*)junk,a4            !# error penalty threshold
              if(debug_output)write(63,*)a4
          elseif(line(i)(1:4).eq.':a5,')then
            read(line(i),*)junk,a5            !# API coefficien
              if(debug_output)write(63,*)a5
          elseif(line(i)(1:4).eq.':a6,')then
            read(line(i),*)junk,a6           !# Minimum routing time step in seconds
              if(debug_output)write(63,*)a6
          elseif(line(i)(1:4).eq.':a7,')then
            read(line(i),*)junk,a7            !# weighting - old vs. new sca value
              if(debug_output)write(63,*)a7
          elseif(line(i)(1:4).eq.':a8,')then
            read(line(i),*)junk,a8            !# min temperature time offset
              if(debug_output)write(63,*)a8
          elseif(line(i)(1:4).eq.':a9,')then
            read(line(i),*)junk,a9            !# max heat deficit /swe ratio')
              if(debug_output)write(63,*)a9
          elseif(line(i)(1:4).eq.':a10')then
            read(line(i),*)junk,a10           !# exponent on uz discharce function')
              if(debug_output)write(63,*)'a10=',a10
          elseif(line(i)(1:4).eq.':a11')then
            read(line(i),*)junk,a11           !# bare ground equiv. veg height for ev')
              if(debug_output)write(63,*)a11
          elseif(line(i)(1:4).eq.':a12')then
            read(line(i),*)junk,a12           !# min precip rate for smearing')
              if(debug_output)write(63,*)'a12=',a12
          elseif(line(i)(1:9).eq.':fmadjust')then
            read(line(i),*)junk,fmadjust      !# snowmelt ripening rate')
              if(debug_output)write(63,*)fmadjust
          elseif(line(i)(1:7).eq.':fmalow')then
            read(line(i),*)junk,fmalow        !# min melt factor multiplier')
              if(debug_output)write(63,*)fmalow
          elseif(line(i)(1:8).eq.':fmahigh')then
            read(line(i),*)junk,fmahigh       !# max melt factor multiplier')
              if(debug_output)write(63,*)fmahigh
          elseif(line(i)(1:9).eq.':gladjust')then
            read(line(i),*)junk,gladjust     !# glacier melt factor multiplier')
              if(debug_output)write(63,*)gladjust
          elseif(line(i)(1:7).eq.':rlapse')then
            read(line(i),*)junk,rlapse        !# precip lapse rate fraction/km')
              if(debug_output)write(63,*)rlapse
          elseif(line(i)(1:7).eq.':tlapse')then
            read(line(i),*)junk,tlapse        !# temperature lapse rate dC/m')
              if(debug_output)write(63,*)tlapse
          elseif(line(i)(1:7).eq.':elvref')then
            read(line(i),*)junk,elvref        !# reference elevation')
              if(debug_output)write(63,*)'elvref=',elvref
          elseif(line(i)(1:13).eq.':rainsnowtemp')then
            read(line(i),*)junk,rainsnowtemp  !# rain/snow temperature')
              if(debug_output)write(63,*)rainsnowtemp
          elseif(line(i)(1:13).eq.':radiusinflce')then
            read(line(i),*)junk,radinfl       !# radius of influence km')
              if(debug_output)write(63,*)radinfl
          elseif(line(i)(1:11).eq.':smoothdist')then
            read(line(i),*)junk,smoothdist    !# smoothing diatance km')
              if(debug_output)write(63,*)smoothdist
          elseif(line(i)(1:8).eq.':flgevp2')then
            read(line(i),*)junk,flgevp2       !# 1=pan;2=Hargreaves;3=Priestley-Taylor')
              if(debug_output)write(63,*)flgevp2
          elseif(line(i)(1:5).eq.':albe')then
            read(line(i),*)junk,albe          !# albedo????')
              if(debug_output)write(63,*)albe
          elseif(line(i)(1:7).eq.':tempa2')then
            read(line(i),*)junk,tempa2
              if(debug_output)write(63,*)tempa2
          elseif(line(i)(1:7).eq.':tempa3')then
            read(line(i),*)junk,tempa3
              if(debug_output)write(63,*)tempa3
          elseif(line(i)(1:5).eq.':tton')then
            read(line(i),*)junk,tton
              if(debug_output)write(63,*)tton
          elseif(line(i)(1:4).eq.':lat')then
            read(line(i),*)junk,lat
              if(debug_output)write(63,*)lat
          elseif(line(i)(1:8).eq.':chnl(1)')then
            read(line(i),*)junk,chnl(1)       !# manning`s n multiplier')
              if(debug_output)write(63,*)'chnl(1)=',chnl(1)
          elseif(line(i)(1:8).eq.':chnl(2)')then
            read(line(i),*)junk,chnl(2)       !# manning`s n multiplier')
              if(debug_output)write(63,*)chnl(2)
          elseif(line(i)(1:8).eq.':chnl(3)')then
            read(line(i),*)junk,chnl(3)       !# manning`s n multiplier')
              if(debug_output)write(63,*)chnl(3)
          elseif(line(i)(1:8).eq.':chnl(4)')then
            read(line(i),*)junk,chnl(4)       !# manning`s n multiplier')
              if(debug_output)write(63,*)chnl(4)
          elseif(line(i)(1:8).eq.':chnl(5)')then
            read(line(i),*)junk,chnl(5)       !# manning`s n multiplier')
              if(debug_output)write(63,*)'chnl(5)=',chnl(5)
          endif
      end do
      if(ios.ne.0)then
            print*,'error in :GlobalParameters'
            
            write(98,*)'Error: in :GlobalParameters'
            write(98,*)'Error: program aborted in read_par_parser @ 572'
            stop 'Program aborted in read_par_parser @ 322'
      endif


!     if debugflg = T, iopt is not changed - left as read in thepar file      
!	if(.not.debugflg)iopt=0                      !# needed for sensitivity
    
               if(debug_output)write(63,*)'RoutingParLine - line',RoutingParLine
!     loop through all the routing lines:
      do j=RoutingParLine+1,EndRoutingParLine-1
        if(line(j)(1:13).eq.':RiverClasses')then
          read(line(j),*)junk(1:13),nbsn  
               if(debug_output)write(63,*)junk(1:13),nbsn
          if(nbsn.ne.nrvr)then
            print*
            print*,'No of river classes in the par file does not match'
            print*,'the number in the shd file'
            print*,'shd file number =',nrvr
            print*,'par file number =',nbsn
            print*
            write(98,*)'Error: No of river classes in the par file does not match'
            write(98,*)'Error: the number in the shd file'
            write(98,*)'Error: shd file number =',nrvr
            write(98,*)'Error: par file number =',nbsn
            write(98,*)'Error: program aborted in read_par_parser @ 598'
            stop 'Fatal error: program aborted in read_par_parser @ 599'
        endif
        elseif(line(j)(1:15).eq.':RiverClassName')then
          read(line(j),*,iostat=ios)junk,(rivtype(i),i=1,nrvr) !#RiverClassName, 
        elseif(line(j)(1:4).eq.':flz')then
          read(line(j),*,iostat=ios)junk,(flz_o(i),i=1,nrvr)   !# lower zone oefficient')	  
        elseif(line(j)(1:4).eq.':pwr')then
          read(line(j),*,iostat=ios)junk,(pwr_o(i),i=1,nrvr)   !# lower zone exponent')	  
        elseif(line(j)(1:4).eq.':r1n')then
          read(line(j),*,iostat=ios)junk,(r1n_o(i),i=1,nrvr)   !# overbank Manning`s n')
        elseif(line(j)(1:4).eq.':r2n')then
          read(line(j),*,iostat=ios)junk,(r2n_o(i),i=1,nrvr)   !# channel Manning`s n')	  
        elseif(line(j)(1:5).eq.':mndr')then
          read(line(j),*,iostat=ios)junk,(mndr_o(i),i=1,nrvr)  !# meander channel length multiplier')	  
        elseif(line(j)(1:4).eq.':aa2')then
          read(line(j),*,iostat=ios)junk,(aa2_o(i),i=1,nrvr)   !# channel area intercept = min channel xsect area')	  
        elseif(line(j)(1:4).eq.':aa3')then
          read(line(j),*,iostat=ios)junk,(aa3_o(i),i=1,nrvr)   !# channel area coefficient')	  
        elseif(line(j)(1:4).eq.':aa4')then
          read(line(j),*,iostat=ios)junk,(aa4_o(i),i=1,nrvr)   !# channel area exponent')	  
        elseif(line(j)(1:6).eq.':theta')then
          read(line(j),*,iostat=ios)junk,(theta_o(i),i=1,nrvr) !# wetland or bank porosity')	  
        elseif(line(j)(1:6).eq.':widep')then
          read(line(j),*,iostat=ios)junk,(widep_o(i),i=1,nrvr) !# channel width to depth ratio')	  
        elseif(line(j)(1:6).eq.':kcond')then
          read(line(j),*,iostat=ios)junk,(kcond_o(i),i=1,nrvr) !# wetland/bank lateral conductivity')	  
        elseif(line(j)(1:5).eq.':pool')then
          read(line(j),*,iostat=ios)junk,(pool_o(i),i=1,nrvr)  !# average area of zero flow pools')	  
        elseif(line(j)(1:6).eq.':rlake')then
          read(line(j),*,iostat=ios)junk,(rlake_o(i),i=1,nrvr) !# in channel lake retardation coefficient')	  
        elseif(line(j)(1:9).eq.':fpFactor')then
          read(line(j),*,iostat=ios)junk,(fpFactor_o(i),i=1,nrvr) !# in channel lake retardation coefficient')	  
        endif
      end do
               if(debug_output)write(63,*)'EndRoutingParLine - line',EndRoutingParLine

      if(ios.ne.0)then
          print*,'read error in :RoutingParameters'
          write(98,*)'Error: read error in :RoutingParameters'
          write(98,*)'Error: program aborted in read_par_parser @ 636'
          stop 'Program aborted in read_par_parser @ 637'
      endif

!     loop through all the hydrology lines
               if(debug_output)write(63,*)'HydrologicalParLine - line ',HydrologicalParLine
      do j=HydrologicalParLine+1,EndHydrologicalParLine-1
        if(line(j)(1:17).eq.':LandCoverClasses')then
          read(line(j),*,iostat=ios)junk,classcount_local  !# No land cover classes'
!         check added Mar. 08/11 nk
        if(classcount_local.ne.classcount)then
          print*,'WARNING:'
          print*,'no of land cover classes in this par file'
          print*,'classcount parfile=', classcount_local 
          print*,'does not match'
          print*,'the number in the shd file'
          print*,'classcount shd file=',classcount
          print*,'please check'
          print*
          write(98,*)'Error: no of land cover classes in this par file'
          write(98,*)'Error: classcount parfile=', classcount_local 
          write(98,*)'Error: does not match'
          write(98,*)'Error: the number in the shd file'
          write(98,*)'Error: classcount shd file=',classcount
          write(98,*)'Error: please check'
          write(98,*)'Error: program aborted in read_par_parser @ 661'
          stop 'progran aborted in read_par_parser @ 662'
        endif
        elseif(line(j)(1:10).eq.':ClassName')then
          read(line(j),*,iostat=ios)junk,(nclass(ii),ii=1,classcount)       !# class name')	  
               if(debug_output)write(63,*)nclass
        elseif(line(j)(1:4).eq.':ds,')then
          read(line(j),*,iostat=ios)junk,(ds(ii),ii=1,classcount)           !# depression storage bare ground mm')	  
        elseif(line(j)(1:5).eq.':dsfs')then
          read(line(j),*,iostat=ios)junk,(dsfs(ii),ii=1,classcount)         !# depression storage snow covered area mm')	  
        elseif(line(j)(1:4).eq.':rec')then
          read(line(j),*,iostat=ios)junk,(rec(ii),ii=1,classcount)          !# interflow coefficient')	  
        elseif(line(j)(1:4).eq.':ak,')then
          read(line(j),*,iostat=ios)junk,(ak(ii),ii=1,classcount)           !# infiltration coefficient bare ground')	  
        elseif(line(j)(1:5).eq.':akfs')then
          read(line(j),*,iostat=ios)junk,(akfs(ii),ii=1,classcount)         !# infiltration coefficient snow covered ground')	  
        elseif(line(j)(1:5).eq.':retn')then
          read(line(j),*,iostat=ios)junk,(retn(ii),ii=1,classcount)         !# upper zone retention mm')	  
        elseif(line(j)(1:5).eq.':ak2,')then
          read(line(j),*,iostat=ios)junk,(ak2(ii),ii=1,classcount)          !# recharge coefficient bare ground')	  
        elseif(line(j)(1:6).eq.':ak2fs')then
          read(line(j),*,iostat=ios)junk,(ak2fs(ii),ii=1,classcount)        !# recharge coefficient snow covered ground')	  
        elseif(line(j)(1:4).eq.':r3,')then
          read(line(j),*,iostat=ios)junk,(r3(ii),ii=1,classcount)           !# overland flow roughness coefficient bare ground')	  
        elseif(line(j)(1:5).eq.':r3fs')then
          read(line(j),*,iostat=ios)junk,(r3fs(ii),ii=1,classcount)   !# overland flow roughness coefficient snow covered grnd')	  
        elseif(line(j)(1:3).eq.':r4')then
          read(line(j),*,iostat=ios)junk,(r4(ii),ii=1,classcount) !# overland flow roughness coefficient impervious area')	  
        elseif(line(j)(1:5).eq.':fpet')then
          read(line(j),*,iostat=ios)junk,(fpet(ii),ii=1,classcount)         !# interception evaporation factor * pet')	  
        elseif(line(j)(1:6).eq.':ftall')then
          read(line(j),*,iostat=ios)junk,(ftall(ii),ii=1,classcount)        !# reduction in PET for tall vegetation')	  
        elseif(line(j)(1:6).eq.':flint')then
          read(line(j),*,iostat=ios)junk,(flint(ii),ii=1,classcount)        !# interception flag  1=on  <1=off')	  
        elseif(line(j)(1:5).eq.':fcap')then
          read(line(j),*,iostat=ios)junk,(fcap(ii),ii=1,classcount)         !# not used - replaced by retn (retention)')	  
        elseif(line(j)(1:6).eq.':ffcap')then
          read(line(j),*,iostat=ios)junk,(ffcap(ii),ii=1,classcount)        !# wilting point - mm of water in uzs')	  
        elseif(line(j)(1:6).eq.':spore')then
          read(line(j),*,iostat=ios)junk,(spore(ii),ii=1,classcount)        !# soil porosity')	  
        elseif(line(j)(1:7).eq.':fratio')then
          read(line(j),*,iostat=ios)junk,(fratio(ii),ii=1,classcount)        !# interception capacity ratio')	  
        endif
      end do
               if(debug_output)write(63,*)'EndHydrologicalParLine - line ',EndHydrologicalParLine

      if(ios.ne.0)then
         print*,'read error in :HydrologicalParameters'
         print*
         write(98,*)'Error: read error in :HydrologicalParameters'
         write(98,*)'Error: program aborted in read_par_parser @ 711'
         stop 'Program aborted in read_par_parser @ 711'
      endif

!     loop through all the hydrology lines
      do j=SnowParLine+1,EndSnowParLine-1
        if(line(j)(1:4).eq.':fm,')then
          read(line(j),*,iostat=ios)junk,(fm(ii),ii=1,classcount)         !# melt factor mm/dC/hour')	  
        elseif(line(j)(1:5).eq.':base')then
           read(line(j),*,iostat=ios)junk,(base(ii),ii=1,classcount)       !# base temperature dC')	  
        elseif(line(j)(1:4).eq.':fmn')then
          read(line(j),*,iostat=ios)junk,(fmn(ii),ii=1,classcount)        !# -ve melt factor')	  
        elseif(line(j)(1:5).eq.':uadj')then
          read(line(j),*,iostat=ios)junk,(uadj(ii),ii=1,classcount)       !# not used')	  
        elseif(line(j)(1:5).eq.':tipm')then
          read(line(j),*,iostat=ios)junk,(tipm(ii),ii=1,classcount)       !# coefficient for ati')	  
        elseif(line(j)(1:4).eq.':rho')then
          read(line(j),*,iostat=ios)junk,(rho(ii),ii=1,classcount)        !# snow density')	  
        elseif(line(j)(1:5).eq.':whcl')then
          read(line(j),*,iostat=ios)junk,(whcl(ii),ii=1,classcount)       !# fraction of swe as water in ripe snow')	  
        elseif(line(j)(1:4).eq.':alb')then
          read(line(j),*,iostat=ios)junk,(alb(ii),ii=1,classcount)          !# albedo')	  
        elseif(line(j)(1:14).eq.':sublim_factor')then
!         rev. 9.7.29  Jul.  07/11  - NK: Add sublim_rate to set sublimation rate/day to pat file
          read(line(j),*,iostat=ios)junk,(sublim_factor(ii),ii=1,classcount)!# sublimation factor ratio')	  
          sublimflg=1  ! default for sublim_factor
        elseif(line(j)(1:12).eq.':sublim_rate')then
          read(line(j),*,iostat=ios)junk,(sublim_rate(ii),ii=1,classcount)!# sublimation factor ratio')	  

          sublimflg=2
        do ii=1,classcount
          sublim_rate(ii)=sublim_rate(ii)/24.0
          sublim_factor(ii)=0.0
        end do
        elseif(line(j)(1:6).eq.':idump')then
          read(line(j),*,iostat=ios)junk,(idump(ii),ii=1,classcount)   !# receiving class for snow redistribution')	  
        elseif(line(j)(1:7).eq.':snocap')then
          read(line(j),*,iostat=ios)junk,(snocap(ii),ii=1,classcount)   !# max swe before redistribution')	  
        elseif(line(j)(1:5).eq.':nsdc')then
          read(line(j),*,iostat=ios)junk,(nsdc(ii),ii=1,classcount)   !# no of points on scd curve - only 1 allowed')	  
        elseif(line(j)(1:7).eq.':sdcsca')then
          read(line(j),*,iostat=ios)junk,(sdcsca(2,ii),ii=1,classcount)   !# snow covered area - ratio=1.0')	  
        elseif(line(j)(1:5).eq.':sdcd')then
          read(line(j),*,iostat=ios)junk,(sdcd(2,ii),ii=1,classcount)   !# swe for 100% snow covered area')	  
        endif
      end do
      if(debug_output)write(63,*)'Passed EndSnowParLine - line ',EndSnowParLine

      !     rev. 10.4.30 Dec.  18/20  = NK Added rainsnowtemp parameter for rain/snow tmp > base tmp
      do ii=1,classcount
          if(rainsnowtemp.lt.base(ii))rainsnowtemp=base(ii)  !old way if rainsnowtemp not in the par file     
      end do
               
      if(ios.ne.0)then
        print*,'error in :SnowParameters'
        print*
        write(98,*)'Error: read error in :SnowParameters'
        write(98,*)'Error: program aborted in read_par_parser @ 763'
        stop 'Program aborted in read_par_parser @ 764'
      endif

      do ii=1,classcount
        sdcsca(1,ii)=0.0
        sdcd(1,ii)=0.0
      end do

      do j=InterceptionParLine+1,EndInterceptionParLine-1
        if(line(j)(1:11).eq.':IntCap_Jan')then
          read(line(j),*,iostat=ios)junk,(h(1,ii),ii=1,classcount)    !# interception capacity jan mm')	  
        elseif(line(j)(1:11).eq.':IntCap_Feb')then
          read(line(j),*,iostat=ios)junk,(h(2,ii),ii=1,classcount)    !# interception capacity feb mm')	  
        elseif(line(j)(1:11).eq.':IntCap_Mar')then
          read(line(j),*,iostat=ios)junk,(h(3,ii),ii=1,classcount)    !# interception capacity mar mm')	  
        elseif(line(j)(1:11).eq.':IntCap_Apr')then
          read(line(j),*,iostat=ios)junk,(h(4,ii),ii=1,classcount)    !# interception capacity apr mm')	  
        elseif(line(j)(1:11).eq.':IntCap_May')then
          read(line(j),*,iostat=ios)junk,(h(5,ii),ii=1,classcount)    !# interception capacity may mm')	  
        elseif(line(j)(1:11).eq.':IntCap_Jun')then
          read(line(j),*,iostat=ios)junk,(h(6,ii),ii=1,classcount)    !# interception capacity jun mm')	  
        elseif(line(j)(1:11).eq.':IntCap_Jul')then
          read(line(j),*,iostat=ios)junk,(h(7,ii),ii=1,classcount)    !# interception capacity jul mm')	  
        elseif(line(j)(1:11).eq.':IntCap_Aug')then
          read(line(j),*,iostat=ios)junk,(h(8,ii),ii=1,classcount)    !# interception capacity aug mm')	  
        elseif(line(j)(1:11).eq.':IntCap_Sep')then
          read(line(j),*,iostat=ios)junk,(h(9,ii),ii=1,classcount)    !# interception capacity sep mm')	  
        elseif(line(j)(1:11).eq.':IntCap_Oct')then
          read(line(j),*,iostat=ios)junk,(h(10,ii),ii=1,classcount)   !# interception capacity oct mm')	  
        elseif(line(j)(1:11).eq.':IntCap_Nov')then
          read(line(j),*,iostat=ios)junk,(h(11,ii),ii=1,classcount)   !# interception capacity nov mm')	  
        elseif(line(j)(1:11).eq.':IntCap_Dec')then
          read(line(j),*,iostat=ios)junk,(h(12,ii),ii=1,classcount)   !# interception capacity dec mm')	  
        endif
      end do
      if(ios.ne.0)then
        print*,'error in :Interception table'
        print*
        write(98,*)'Error: read error in Interception table'
        write(98,*)'Error: program aborted in read_par_parser @ 803'
        stop 'Program aborted in read_par_parser @ 803'
      endif

!     rev. 9.8.94  Nov.  20/13  - NK: Added check on interception capacity for water
!     This is to make sure no one inadvertently sets a high value on the
!     interception capacity for water.
      Do i=1,12
        ii=classcount-1
        h(i,ii)=amin1(h(i,ii),0.000001)
      end do


!     set the parameters for the coupled & parent wetland equal
      if(wetflg.eq.'y')then
      if(nclass(classcount-3).eq.nclass(classcount-2))then
          ds(classcount-2)           =ds(classcount-3)
          dsfs(classcount-2)         =dsfs(classcount-3)
          rec(classcount-2)          =rec(classcount-3)
          ak(classcount-2)           =ak(classcount-3)
          akfs(classcount-2)         =akfs(classcount-3)
          retn(classcount-2)         =retn(classcount-3)
          ak2(classcount-2)          =ak2(classcount-3)
          ak2fs(classcount-2)        =ak2fs(classcount-3)
          r3(classcount-2)           =r3(classcount-3)
          r3fs(classcount-2)         =r3fs(classcount-3)
          r4(classcount-2)           =r4(classcount-3)
          fpet(classcount-2)         =fpet(classcount-3)
          ftall(classcount-2)        =ftall(classcount-3)
          flint(classcount-2)        =flint(classcount-3)
          fcap(classcount-2)         =fcap(classcount-3)
          ffcap(classcount-2)        =ffcap(classcount-3)
          spore(classcount-2)        =spore(classcount-3)
!     rev. 9.8.84  Sep.  15/13  - NK: Added fratio to list of equal values for bog & fen
          fratio(classcount-2)        =fratio(classcount-3)
          do j=1,12
            if(h(j,classcount-2).ne.h(j,classcount-3))then
              write(98,98000)classcount-2,classcount-3,j
98000         format('WARNING: Int. values for classes',i5,'and',i5,'for month',i5,'are set equal')             
              h(j,classcount-2)=h(j,classcount-3)
            endif
          end do
          fm((classcount-2))          =fm(classcount-3)
          base((classcount-2))        =base(classcount-3)
          fmn((classcount-2))         =fmn(classcount-3)
          uadj((classcount-2))        =uadj(classcount-3)
          tipm((classcount-2))        =tipm(classcount-3)
          rho((classcount-2))         =rho(classcount-3)
          whcl((classcount-2))        =whcl(classcount-3)
          alb(classcount-2)          =alb(classcount-3)
          sublim_factor(classcount-2)=sublim_factor(classcount-3)
          sublim_rate(classcount-2)=sublim_rate(classcount-3)
          idump((classcount-2))       =idump(classcount-3)
          snocap((classcount-2))      =snocap(classcount-3)
          nsdc((classcount-2))        =nsdc(classcount-3)
          sdcsca(2,(classcount-2))    =sdcsca(2,classcount-3)
          sdcd(2,(classcount-2))      =sdcd(2,classcount-3)
      endif
    endif
!     rev. 9.8.08  Nov.  18/11  - NK: added fratio for interception hight optimization
!     correct the interception capacity by the multiplier fratio
!     This is done in sub so it does not affect writing the par file with the coupler
               if(debug_output)write(63,*)'EndInterceptionParLine - line ',EndInterceptionParLine

    if(ios.ne.0)then
        print*,'error in :InterceptionParameterTable'
      print*
      write(*,'(''paused, type [enter] to abort'')') 
      read (*,*) 
      stop 'Program aborted in read_par_parser @ 509'
      endif
      
!     rev. 10.2.30 Aug.  21/18  - NK: Added error for ftall(water) 
!     Check ftall has value for water class
      if(ftall(classcount-1).le.0.1.or.ftall(classcount-1).gt.1.5)then
          print*
          print*,'Warning:'
          print*,'ftall for the water class has an inappropriate value'
          print*,'It should be 0.1 < ftall(water) < 1.5'
          print*
          print*,'Your value =',ftall(classcount-1)
          write(98,*)'Warning: ftall for the water class has an inappropriate value'
          write(98,*)'Warning: It should be 0.1 < ftall(water) < 1.5'
          write(98,*)'Warning: Your value =',ftall(classcount-1)
          write(98,*)'Warning: A good initial value is fpet for water = 1.0'
!          write(98,*)'Error: program aborted in read_par_parser @ 889'
!          print*,'To continue with ftall(water) = 0.9, hit enter'
!          pause 'or Ctrl C  to quit'
!          ftall(classcount-1)=0.9
      endif

      do ii=1,classcount
!       adjust potential evaporation from vegetation by amount
        hmax=0.0  
!     rev. 9.2.22  Nov.  15/05  - NK: Fixed hmax bug in rdpar 
!             hmax=0 was erroniously moved down between 
!                    ver 9.2.03 & 9.2.2
        do i=1,12
!     rev  9.1.25  Sep.  11/02  - Added A11 as bare ground equiv. vegn height  
!          h(i,ii)=h(i,ii)+a11    ! taken out and replaced Jan. 25/06 nk
          h(i,ii)=amax1(a11,h(i,ii)) !  Jan. 25/06 nk
          hmax=amax1(hmax,h(i,ii))
        end do 

        do i=1,12
!         problems here nk Nov. 20/02
!         rev. 9.1.32  Nov.  20/02  - Fixed fpetmon() wrt. h()
!         h(month,class) is defined as the max interception storage. 
!         So for soil evaporation it has to become an index 
!         fpetmo() is that index - it modifies ftall(ii) to account for veg height
          fpetmo(i,ii)=ftall(ii)*h(i,ii)/hmax
!          fpetmo(i,ii)=ftall(ii)*h(i,ii)
        end do
!   if(debug_output)write(51,1003)'hmax ',hmax
!   if(debug_output)write(51,1003)'h(  )',(h(i,ii),i=1,12)
!   if(debug_output)write(51,1003)'fpetm',(fpetmo(i,ii),i=1,12)
      end do

!     SET POTENTIALS FOR EACH SOIL TYPE:
!     moved from soilinit   Mar. 14/07  nk
      do ii=1,classcount
!        NOTE FOR AK<0.0 WE HAVE WATER AREA & POT ISN'T USED
         if(ak(ii).gt.0.0)then
!           AK IS IN MM/HR
!           BUT FORMULA IS FOR MM/SEC
            xxx1=-alog10(ak(ii)/3600.)
            pot(ii)=250.*xxx1+100.
         endif
!        REV. 7.5 SEPERATE SNOW COVERED AND BARE GROUND
         if(akfs(ii).gt.0.0)then
!           ak is in mm/hr
!           but formula is for mm/sec
            xxx1=-alog10(akfs(ii)/3600.)
            potfs(ii)=250.*xxx1+100.
         endif
      end do

      if(iopt.eq.2)print*, 'param 10a'

!     Default monthly ET values - to be replaced by input
      do j=1,12
      do ii=1,classcount
        evap(ii,j)=0.0
      end do
    end do
    
      do j=MonthlyEvapParLine+1,EndMonthlyEvapParLine-1
        if(line(j)(1:14).eq.':Montly_ET_Jan')then
          read(line(j),*,iostat=ios)junk,(evap(ii,1),ii=1,classcount)     !# monthly evapotranspiration jan mm')	  
        elseif(line(j)(1:14).eq.':Montly_ET_Feb')then
          read(line(j),*,iostat=ios)junk,(evap(ii,2),ii=1,classcount)     !# monthly evapotranspiration feb mm')	  
        elseif(line(j)(1:14).eq.':Montly_ET_Mar')then
          read(line(j),*,iostat=ios)junk,(evap(ii,3),ii=1,classcount)     !# monthly evapotranspiration mar mm')	  
        elseif(line(j)(1:14).eq.':Montly_ET_Apr')then
          read(line(j),*,iostat=ios)junk,(evap(ii,4),ii=1,classcount)     !# monthly evapotranspiration apr mm')	  
        elseif(line(j)(1:14).eq.':Montly_ET_May')then
          read(line(j),*,iostat=ios)junk,(evap(ii,5),ii=1,classcount)     !# monthly evapotranspiration may mm')	  
        elseif(line(j)(1:14).eq.':Montly_ET_Jun')then
          read(line(j),*,iostat=ios)junk,(evap(ii,6),ii=1,classcount)     !# monthly evapotranspiration jun mm')	  
        elseif(line(j)(1:14).eq.':Montly_ET_Jul')then
          read(line(j),*,iostat=ios)junk,(evap(ii,7),ii=1,classcount)     !# monthly evapotranspiration jul mm')	  
        elseif(line(j)(1:14).eq.':Montly_ET_Aug')then
          read(line(j),*,iostat=ios)junk,(evap(ii,8),ii=1,classcount)     !# monthly evapotranspiration aug mm')	  
        elseif(line(j)(1:14).eq.':Montly_ET_Sep')then
          read(line(j),*,iostat=ios)junk,(evap(ii,9),ii=1,classcount)     !# monthly evapotranspiration sep mm')	  
        elseif(line(j)(1:14).eq.':Montly_ET_Oct')then
          read(line(j),*,iostat=ios)junk,(evap(ii,10),ii=1,classcount)    !# monthly evapotranspiration oct mm')	  
        elseif(line(j)(1:14).eq.':Montly_ET_Nov')then
          read(line(j),*,iostat=ios)junk,(evap(ii,11),ii=1,classcount)    !# monthly evapotranspiration nov mm')	  
        elseif(line(j)(1:14).eq.':Montly_ET_Dec')then
          read(line(j),*,iostat=ios)junk,(evap(ii,12),ii=1,classcount)    !# monthly evapotranspiration dec mm')	  
        endif
      end do

               if(debug_output)write(63,*)'EndMonthlyEvapParLine - line ',EndMonthlyEvapParLine

    if(ios.ne.0)then
        print*,'error in :MonthlyEvaporatioTabble'
      print*
      write(*,'(''paused, type [enter] to abort'')') 
      read (*,*) 
      stop 'Program aborted in read_par_parser @ 593'
    endif

      do j=OptimizationParLine+1,EndOptimizationParLine-1
        if(line(j)(1:5).eq.':numa')then
        read(line(j),*,iostat=ios)junk,numa      !# PS optimization 1=yes 0=no')
        elseif(line(j)(1:5).eq.':nper')then
          read(line(j),*,iostat=ios)junk,nper      !# opt 1=delta 0=absolute')
        elseif(line(j)(1:3).eq.':kc')then
          read(line(j),*,iostat=ios)junk,kc        !# no of times delta halved')
        elseif(line(j)(1:5).eq.':maxn')then
          read(line(j),*,iostat=ios)junk,maxn      !# max no of trials')
        elseif(line(j)(1:7).eq.':ddsflg')then
          read(line(j),*,iostat=ios)junk,dds_flag  !# 0=single run  1=DDS ')
        elseif(line(j)(1:7).eq.':errflg')then
          read(line(j),*,iostat=ios)junk,errflg    !# 1=wMSE 2=SSE 3=wSSE 4=VOL ')
        endif
      end do
      
      if(ios.ne.0)then
        print*,'error in :OptimizationSwitches'
        print*
        write(98,*)'Error: red error in :OptimizationSwitches'
        write(98,*)'Error: program aborted in read_par_parser @ 1007'
        stop 'Program aborted in read_par_parser @ 1008'
      endif
      if(debug_output)write(63,*)'EndOptimizationParLine - line ',EndOptimizationParLine

      do j=APILimitsLine+1,EndAPILimitsLine-1
        if(line(j)(1:6).eq.':a5dlt')then
          read(line(j),*,iostat=ios)junk,a5dlt
        elseif(line(j)(1:6).eq.':a5low')then
          read(line(j),*,iostat=ios)junk,a5low
        elseif(line(j)(1:6).eq.':a5hgh')then
          read(line(j),*,iostat=ios)junk,a5hgh
        endif
      end do
               if(debug_output)write(63,*)'EndAPILimitsLine - line ',EndAPILimitsLine

!     Hydrological parameters
      do j=HydrologicalParLimitsLine+1,EndHydrologicalParLimitsLine-1
        if(line(j)(1:6).eq.':akdlt')then
          read(line(j),*,iostat=ios)junk,(akdlt(ii),ii=1,classcount)
        elseif(line(j)(1:6).eq.':aklow')then
          read(line(j),*,iostat=ios)junk,(aklow(ii),ii=1,classcount)
        elseif(line(j)(1:6).eq.':akhgh')then
          read(line(j),*,iostat=ios)junk,(akhgh(ii),ii=1,classcount)
        elseif(line(j)(1:8).eq.':akfsdlt')then
          read(line(j),*,iostat=ios)junk,(akfsdlt(ii),ii=1,classcount)
        elseif(line(j)(1:8).eq.':akfslow')then
          read(line(j),*,iostat=ios)junk,(akfslow(ii),ii=1,classcount)
        elseif(line(j)(1:8).eq.':akfshgh')then
          read(line(j),*,iostat=ios)junk,(akfshgh(ii),ii=1,classcount)
        elseif(line(j)(1:7).eq.':recdlt')then
          read(line(j),*,iostat=ios)junk,(recdlt(ii),ii=1,classcount)
        elseif(line(j)(1:7).eq.':reclow')then
          read(line(j),*,iostat=ios)junk,(reclow(ii),ii=1,classcount)
        elseif(line(j)(1:7).eq.':rechgh')then
          read(line(j),*,iostat=ios)junk,(rechgh(ii),ii=1,classcount)
        elseif(line(j)(1:6).eq.':r3dlt')then
          read(line(j),*,iostat=ios)junk,(r3dlt(ii),ii=1,classcount)
        elseif(line(j)(1:6).eq.':r3low')then
          read(line(j),*,iostat=ios)junk,(r3low(ii),ii=1,classcount)
        elseif(line(j)(1:6).eq.':r3hgh')then
          read(line(j),*,iostat=ios)junk,(r3hgh(ii),ii=1,classcount)
        elseif(line(j)(1:8).eq.':fpetdlt')then
          read(line(j),*,iostat=ios)junk,(fpetdlt(ii),ii=1,classcount)
        elseif(line(j)(1:8).eq.':fpetlow')then
          read(line(j),*,iostat=ios)junk,(fpetlow(ii),ii=1,classcount)
        elseif(line(j)(1:8).eq.':fpethgh')then
          read(line(j),*,iostat=ios)junk,(fpethgh(ii),ii=1,classcount)
        elseif(line(j)(1:9).eq.':ftalldlt')then
          read(line(j),*,iostat=ios)junk,(ftalldlt(ii),ii=1,classcount)
        elseif(line(j)(1:9).eq.':ftalllow')then
          read(line(j),*,iostat=ios)junk,(ftalllow(ii),ii=1,classcount)
        elseif(line(j)(1:9).eq.':ftallhgh')then
          read(line(j),*,iostat=ios)junk,(ftallhgh(ii),ii=1,classcount)
        elseif(line(j)(1:10).eq.':fratiodlt')then
          read(line(j),*,iostat=ios)junk,(fratiodlt(ii),ii=1,classcount)
        elseif(line(j)(1:10).eq.':fratiolow')then
          read(line(j),*,iostat=ios)junk,(fratiolow(ii),ii=1,classcount)
        elseif(line(j)(1:10).eq.':fratiohgh')then
          read(line(j),*,iostat=ios)junk,(fratiohgh(ii),ii=1,classcount)
        elseif(line(j)(1:8).eq.':retndlt')then
          read(line(j),*,iostat=ios)junk,(retndlt(ii),ii=1,classcount)
        elseif(line(j)(1:8).eq.':retnlow')then
          read(line(j),*,iostat=ios)junk,(retnlow(ii),ii=1,classcount)
        elseif(line(j)(1:8).eq.':retnhgh')then
          read(line(j),*,iostat=ios)junk,(retnhgh(ii),ii=1,classcount)
        elseif(line(j)(1:7).eq.':ak2dlt')then
          read(line(j),*,iostat=ios)junk,(ak2dlt(ii),ii=1,classcount)
        elseif(line(j)(1:7).eq.':ak2low')then
          read(line(j),*,iostat=ios)junk,(ak2low(ii),ii=1,classcount)
        elseif(line(j)(1:7).eq.':ak2hgh')then
          read(line(j),*,iostat=ios)junk,(ak2hgh(ii),ii=1,classcount)
        elseif(line(j)(1:9).eq.':ak2fsdlt')then
          read(line(j),*,iostat=ios)junk,(ak2fsdlt(ii),ii=1,classcount)
        elseif(line(j)(1:9).eq.':ak2fslow')then
          read(line(j),*,iostat=ios)junk,(ak2fslow(ii),ii=1,classcount)
        elseif(line(j)(1:9).eq.':ak2fshgh')then
          read(line(j),*,iostat=ios)junk,(ak2fshgh(ii),ii=1,classcount)
!     this due to a typo in the write par file
        elseif(line(j)(1:9).eq.':ak2scdlt')then
          read(line(j),*,iostat=ios)junk,(ak2fsdlt(ii),ii=1,classcount)
        elseif(line(j)(1:9).eq.':ak2fslow')then
          read(line(j),*,iostat=ios)junk,(ak2fslow(ii),ii=1,classcount)
        elseif(line(j)(1:9).eq.':ak2schgh')then
          read(line(j),*,iostat=ios)junk,(ak2fshgh(ii),ii=1,classcount)
        endif
      end do
               if(debug_output)write(63,*)'EndHydrologicalParLimitsLine - line ',EndHydrologicalParLimitsLine

      if(ios.ne.0)then
        print*,'error in HydrologicalParLimits'
        print*
        write(98,*)'Error: read error in HydrologicalParLimits'
        write(98,*)'Error: program aborted in read_par_parser @ 1100'
        stop 'Program aborted in read_par_parser @ 1101'
      endif

!     globalsnow parameters
      do j=GlobalSnowParLimitsLine+1,EndGlobalSnowParLimitsLine-1
        if(line(j)(1:12).eq.':fmadjustdlt')then
          read(line(j),*,iostat=ios)junk,fmadjustdlt
        elseif(line(j)(1:12).eq.':fmadjustlow')then
          read(line(j),*,iostat=ios)junk,fmadjustlow
        elseif(line(j)(1:12).eq.':fmadjusthgh')then
          read(line(j),*,iostat=ios)junk,fmadjusthgh
        elseif(line(j)(1:10).eq.':fmalowdlt')then
          read(line(j),*,iostat=ios)junk,fmalowdlt
        elseif(line(j)(1:10).eq.':fmalowlow')then
          read(line(j),*,iostat=ios)junk,fmalowlow
        elseif(line(j)(1:10).eq.':fmalowhgh')then
          read(line(j),*,iostat=ios)junk,fmalowhgh
        elseif(line(j)(1:11).eq.':fmahighdlt')then
          read(line(j),*,iostat=ios)junk,fmahighdlt
        elseif(line(j)(1:11).eq.':fmahighlow')then
          read(line(j),*,iostat=ios)junk,fmahighlow
        elseif(line(j)(1:11).eq.':fmahighhgh')then
          read(line(j),*,iostat=ios)junk,fmahighhgh
        elseif(line(j)(1:12).eq.':gladjustdlt')then
          read(line(j),*,iostat=ios)junk,gladjustdlt
        elseif(line(j)(1:12).eq.':gladjustlow')then
          read(line(j),*,iostat=ios)junk,gladjustlow
        elseif(line(j)(1:12).eq.':gladjusthgh')then
          read(line(j),*,iostat=ios)junk,gladjusthgh
        endif
      end do
               if(debug_output)write(63,*)'EndGlobalSnowParLimitsLine - line ',EndGlobalSnowParLimitsLine

      if(ios.ne.0)then
        print*,'error in GlobalSnowParLinmits'
        write(98,*)'Error: read error in GlobalSnowParLinmits'
        write(98,*)'Error: program aborted in read_par_parser @ 1137'
        stop 'Program aborted in read_par_parser @ 1138'
      endif

      do j=SnowParLimitsLine+1,EndSnowParLimitsLine-1
        if(line(j)(1:6).eq.':fmdlt')then
          read(line(j),*,iostat=ios)junk,(fmdlt(ii),ii=1,classcount)
        elseif(line(j)(1:6).eq.':fmlow')then
          read(line(j),*,iostat=ios)junk,(fmlow(ii),ii=1,classcount)
        elseif(line(j)(1:6).eq.':fmhgh')then
          read(line(j),*,iostat=ios)junk,(fmhgh(ii),ii=1,classcount)
        elseif(line(j)(1:8).eq.':basedlt')then
          read(line(j),*,iostat=ios)junk,(basdlt(ii),ii=1,classcount)
        elseif(line(j)(1:8).eq.':baselow')then
          read(line(j),*,iostat=ios)junk,(baslow(ii),ii=1,classcount)
        elseif(line(j)(1:8).eq.':basehgh')then
          read(line(j),*,iostat=ios)junk,(bashgh(ii),ii=1,classcount)
        elseif(line(j)(1:7).eq.':subdlt')then
          read(line(j),*,iostat=ios)junk,(subdlt(ii),ii=1,classcount)
        elseif(line(j)(1:7).eq.':sublow')then
          read(line(j),*,iostat=ios)junk,(sublow(ii),ii=1,classcount)
        elseif(line(j)(1:7).eq.':subhgh')then
          read(line(j),*,iostat=ios)junk,(subhgh(ii),ii=1,classcount)
        elseif(line(j)(1:8).eq.':sdcddlt')then
          read(line(j),*,iostat=ios)junk,(sdcddlt(ii),ii=1,classcount)
        elseif(line(j)(1:8).eq.':sdcdlow')then
          read(line(j),*,iostat=ios)junk,(sdcdlow(ii),ii=1,classcount)
        elseif(line(j)(1:8).eq.':sdcddhgh')then
          read(line(j),*,iostat=ios)junk,(sdcdhgh(ii),ii=1,classcount)
        endif
      end do
      if(ios.ne.0)then
        print*,'error in error in SnowParLimits'
        write(98,*)'Error: read error in error in SnowParLimits'
        write(98,*)'Error: program aborted in read_par_parser @ 1173'
        stop 'Program aborted in read_par_parser @ 1173'
      endif
      if(debug_output)write(63,*)'EndSnowParLimitsLine - line ',EndSnowParLimitsLine

      if(sublimflg.eq.2)then          ! for sublim_factor only
        do ii=1,classcount
          sublow(ii)=sublow(ii)/24
          subhgh(ii)=subhgh(ii)/24
        end do
      endif

!     Routing parameters
      do j=RoutingParLimitsLine+1,EndRoutingParLimitsLine-1
        if(line(j)(1:7).eq.':flzdlt')then
          read(line(j),*,iostat=ios)junk,(flzdlt(i),i=1,nrvr)
        elseif(line(j)(1:7).eq.':flzlow')then
          read(line(j),*,iostat=ios)junk,(flzlow(i),i=1,nrvr)
        elseif(line(j)(1:7).eq.':flzhgh')then
          read(line(j),*,iostat=ios)junk,(flzhgh(i),i=1,nrvr)
        elseif(line(j)(1:7).eq.':pwrdlt')then
          read(line(j),*,iostat=ios)junk,(pwrdlt(i),i=1,nrvr)
        elseif(line(j)(1:7).eq.':pwrlow')then
          read(line(j),*,iostat=ios)junk,(pwrlow(i),i=1,nrvr)
        elseif(line(j)(1:7).eq.':pwrhgh')then
          read(line(j),*,iostat=ios)junk,(pwrhgh(i),i=1,nrvr)
        elseif(line(j)(1:7).eq.':r2ndlt')then
          read(line(j),*,iostat=ios)junk,(r2ndlt(i),i=1,nrvr)
        elseif(line(j)(1:7).eq.':r2nlow')then
          read(line(j),*,iostat=ios)junk,(r2nlow(i),i=1,nrvr)
        elseif(line(j)(1:7).eq.':r2nhgh')then
          read(line(j),*,iostat=ios)junk,(r2nhgh(i),i=1,nrvr)
        elseif(line(j)(1:9).eq.':thetadlt')then
          read(line(j),*,iostat=ios)junk,(thetadlt(i),i=1,nrvr)
        elseif(line(j)(1:9).eq.':thetalow')then
          read(line(j),*,iostat=ios)junk,(thetalow(i),i=1,nrvr)
        elseif(line(j)(1:9).eq.':thetahgh')then
          read(line(j),*,iostat=ios)junk,(thetahgh(i),i=1,nrvr)
        elseif(line(j)(1:9).eq.':kconddlt')then
          read(line(j),*,iostat=ios)junk,(kconddlt(i),i=1,nrvr)
        elseif(line(j)(1:9).eq.':kcondlow')then
          read(line(j),*,iostat=ios)junk,(kcondlow(i),i=1,nrvr)
        elseif(line(j)(1:9).eq.':kcondhgh')then
          read(line(j),*,iostat=ios)junk,(kcondhgh(i),i=1,nrvr)
        elseif(line(j)(1:9).eq.':rlakedlt')then
          read(line(j),*,iostat=ios)junk,(rlakedlt(i),i=1,nrvr)
        elseif(line(j)(1:9).eq.':rlakelow')then
          read(line(j),*,iostat=ios)junk,(rlakelow(i),i=1,nrvr)
        elseif(line(j)(1:9).eq.':rlakehgh')then
          read(line(j),*,iostat=ios)junk,(rlakehgh(i),i=1,nrvr)
          
        elseif(line(j)(1:8).eq.':pooldlt')then
          read(line(j),*,iostat=ios)junk,(pooldlt(i),i=1,nrvr)
        elseif(line(j)(1:8).eq.':poollow')then
          read(line(j),*,iostat=ios)junk,(poollow(i),i=1,nrvr)
        elseif(line(j)(1:8).eq.':poolhgh')then
          read(line(j),*,iostat=ios)junk,(poolhgh(i),i=1,nrvr)
          
        elseif(line(j)(1:7).eq.':r1ndlt')then
          read(line(j),*,iostat=ios)junk,(r1ndlt(i),i=1,nrvr)
        elseif(line(j)(1:7).eq.':r1nlow')then
          read(line(j),*,iostat=ios)junk,(r1nlow(i),i=1,nrvr)
        elseif(line(j)(1:7).eq.':r1nhgh')then
          read(line(j),*,iostat=ios)junk,(r1nhgh(i),i=1,nrvr)
          
        elseif(line(j)(1:7).eq.':aa2dlt')then
          read(line(j),*,iostat=ios)junk,(aa2dlt(i),i=1,nrvr)
        elseif(line(j)(1:7).eq.':aa2low')then
          read(line(j),*,iostat=ios)junk,(aa2low(i),i=1,nrvr)
        elseif(line(j)(1:7).eq.':aa2hgh')then
          read(line(j),*,iostat=ios)junk,(aa2hgh(i),i=1,nrvr)
          
        elseif(line(j)(1:7).eq.':aa3dlt')then
          read(line(j),*,iostat=ios)junk,(aa3dlt(i),i=1,nrvr)
        elseif(line(j)(1:7).eq.':aa3low')then
          read(line(j),*,iostat=ios)junk,(aa3low(i),i=1,nrvr)
        elseif(line(j)(1:7).eq.':aa3hgh')then
          read(line(j),*,iostat=ios)junk,(aa3hgh(i),i=1,nrvr)
          
        elseif(line(j)(1:9).eq.':widepdlt')then
          read(line(j),*,iostat=ios)junk,(widepdlt(i),i=1,nrvr)
        elseif(line(j)(1:9).eq.':wideplow')then
          read(line(j),*,iostat=ios)junk,(wideplow(i),i=1,nrvr)
        elseif(line(j)(1:9).eq.':widephgh')then
          read(line(j),*,iostat=ios)junk,(widephgh(i),i=1,nrvr)
          
        endif
      end do
      if(ios.ne.0)then
        print*,'error in error in RoutingParLimits'
        write(98,*)'Error: read error in error in RoutingParLimits'
        write(98,*)'Error: program aborted in read_par_parser @ 1228'
        stop 'Program aborted in read_par_parser @ 1229'
      endif
      if(debug_output)write(63,*)'EndRoutingParLimitsLine - line ',EndRoutingParLimitsLine

      !     rev. 9.8.01  Jul.  21/11  - NK: added ragmet optimization to dds setup
!     Global parameters
      if(GlobalParLimitsLine.gt.0)then
        do j=GlobalParLimitsLine+1,EndGlobalParLimitsLine-1
          if(line(j)(1:10).eq.':rlapsedlt')then
            read(line(j),*,iostat=ios)junk,rlapsedlt
          elseif(line(j)(1:10).eq.':rlapselow')then
            read(line(j),*,iostat=ios)junk,rlapselow
          elseif(line(j)(1:10).eq.':rlapsehgh')then
            read(line(j),*,iostat=ios)junk,rlapsehgh
            
          elseif(line(j)(1:10).eq.':tlapsedlt')then
            read(line(j),*,iostat=ios)junk,tlapsedlt
          elseif(line(j)(1:10).eq.':tlapselow')then
            read(line(j),*,iostat=ios)junk,tlapselow
          elseif(line(j)(1:10).eq.':tlapsehgh')then
            read(line(j),*,iostat=ios)junk,tlapsehgh

            !     rev. 10.4.31 Dec.  22/20  = NK fix climate end of year
          elseif(line(j)(1:16).eq.':rainsnowtempdlt')then
            read(line(j),*,iostat=ios)junk,rainsnowtempdlt
          elseif(line(j)(1:16).eq.':rainsnowtemplow')then
            read(line(j),*,iostat=ios)junk,rainsnowtemplow
          elseif(line(j)(1:16).eq.':rainsnowtemphgh')then
            read(line(j),*,iostat=ios)junk,rainsnowtemphgh
            
          elseif(line(j)(1:11).eq.':radinfldlt')then
            read(line(j),*,iostat=ios)junk,radinfldlt
          elseif(line(j)(1:11).eq.':radinfllow')then
            read(line(j),*,iostat=ios)junk,radinfllow
          elseif(line(j)(1:11).eq.':radinflhgh')then
            read(line(j),*,iostat=ios)junk,radinflhgh
            
          elseif(line(j)(1:13).eq.':smoothdisdlt')then
            read(line(j),*,iostat=ios)junk,smoothdistdlt
          elseif(line(j)(1:13).eq.':smoothdislow')then
            read(line(j),*,iostat=ios)junk,smoothdistlow
          elseif(line(j)(1:13).eq.':smoothdishgh')then
            read(line(j),*,iostat=ios)junk,smoothdisthgh
          endif
        end do
        if(ios.ne.0)then
          print*,'read error in GlobalParLimits'
          write(98,*)'Error: read error in error in GlobalParLimits'
          write(98,*)'Error: program aborted in read_par_parser @ 1258'
          stop 'Program aborted in read_par_parser @ 1259'
        endif
        if(debug_output)write(63,*)'EndGlobalParLimitsLine - line ',EndGlobalParLimitsLine
      
!     REV. 10.1.30 May   08/16  - NK: Added smoothdist warning in read_par_parser
!       check smoothing distance - if too large, program ragmet will not write data
        if(smoothdist.gt.radinfl/2.0)then
          print*
          print*,'WARNING: smoothing dist > radius of influence/2'
          print*,'Likely will lead to problems in RAGMET.exe'
          print*,'and TMP.exe'
          print*,'with no data written to file & program hanging'
          print*,'Program continues with ddsflg = 1'
          if(dds_flag.eq.0)then
            write(*,*)'Pls fix par file (Ctrl C) or hit enter to continue'	
            write(*,'(''paused, type [enter] to abort'')') 
            read (*,*) 
          endif
        endif  
      else
        rlapsedlt=-1.0
        tlapsedlt=-1.0
        rainsnowtempdlt=-1.0
        radinfldlt=-1.0
        smoothdistdlt=-1.0
        print*,'WARNING:'
        print*,'New section of GlobalParLimits not found'
        print*,'in the par file.'
        print*,'rlapsedlt,tlapsedlt,rainsnowtempdlt,radinfldlt & smoothdist'
        print*,'set - -1.0  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
        print*
      endif

! Assignments & checking
! Assignments & checking
! Assignments & checking

      if(unitNum.eq.51)then   ! write a few blank lines for readability
      write(line(i),*)  
      write(line(i),*)
        return
    endif
      
      close(unitNum,status='keep')

    if(iopt.gt.0)print*,'Closed unit=',unitNum
    if(iopt.gt.0)print*,fln(flnNum)(1:50),'read'
      
!  do n=1,noresv
!      if(rlake_o(n).gt.1.0E-20)then
!          print*,'WARNING:'
!          print*,'rlake value for river class ',n,' too large'
!          stop 'Program aborted in read_par_parser @ 1246'
!      endif
!  end do
          
!     rev. 10.5.01 Jan.  06/23  = NK Added check on rlake pars to avoid instability
!     rlake shoule be < 1.0E-08 for pond routing OR > 0.1 for r2n correction    
      do i=1,nrvr
          if(rlake_o(i).gt.1.0E-08.and.rlake_o(i).lt.0.1)then
              write(98,*)'Error: rlake for basin # ',i,'=',rlake_o(i),' >1.0E-08 '
              write(98,*)'Error: rlake lowered to 1.0E-08 to avoid instability'
              print*,'Error: rlake for basin # ',i,'=',rlake_o(i),' >1.0E-08 '
              print*,'rlake raised to 1.0E-08 to avoid instability'
              rlake_o(i)=1.0E-08
              rlake(i)=1.0E-08
              print*,'Please change the par file to avoid this message' 
!              pause 'In read_par_parser @ 1346'
          endif
      end do
      do i=1,nrvr
          if(rlakehgh(i).gt.1.0E-08)then
              print*,'Error:'
              print*,'rlakehgh for basin # ',i,' =',rlakehgh(i),' >1.0E-08 '
              print*,'rlakehgh lowered to 1.0E-08 to avoid instability'
              print*,'Please change the par file to '
              print*,'switch to pond routing and avoid this message' 
              rlakehgh(i)=1.0E-08
          endif
      end do
    
!     rev. 10.5.16 June  29/23  = NK re-assign the parameters if ibn has changed for any grid
!     grid the parameters for watroute & CHARM
      do n=1,naa
        flz(n)=flz_o(ibn(n))
        pwr(n)=pwr_o(ibn(n))
        r1n(n)=r1n_o(ibn(n))
        r2n(n)=r2n_o(ibn(n))
        mndr(n)=mndr_o(ibn(n))
        aa2(n)=aa2_o(ibn(n))
        aa3(n)=aa3_o(ibn(n))
        aa4(n)=aa4_o(ibn(n))
        theta(n)=theta_o(ibn(n))
        widep(n)=widep_o(ibn(n))
        kcond(n)=kcond_o(ibn(n))
        pool(n)=pool_o(ibn(n))
        rlake(n)=rlake_o(ibn(n))
        fpFactor(n)=fpFactor_o(ibn(n))
        
        
!       For NCA with channel size minimized to have swale flow only
        if(aa2(n).le.0.0000001)then
            theta(n)=-1.0*abs(theta(n))
            kcond(n)=-1.0*abs(kcond(n))
            rlake(n)=-1.0*abs(rlake(n))
        endif
        
      end do
      if(unitNum.eq.51)then    !  never gets here - see above!!!!!!!!!!!!!!!!!!
         write(line(i),*)
       write(line(i),*)line
         write(line(i),*)'diff  ',(diff(i),i=1,12)
         write(line(i),*)'humid ',(hu(i),i=1,12)
         write(line(i),*)'pres  ',(pres(i),i=1,12)
         write(line(i),*)
      endif

!  if(dds_flag.eq.1)iopt=0
      
      do ii=1,classcount-2
        if(ak(ii).le.0.0)then
            print*,'ak(',ii,')<= 0.0',ak(ii)
            write(98,*)'Error: ak(',ii,')<= 0.0',ak(ii)
            write(98,*)'Error: program aborted in read_par_parser @ 1347'
            stop 'program aborted in read_par_parser @ 1348'
        endif
      end do
      if(ak(classcount).le.0.0)then
            print*,'ak(',classcount,')<= 0.0',ak(classcount)
            write(98,*)'Error: ak(',classcount,')<= 0.0',ak(classcount)
            write(98,*)'Error: program aborted in read_par_parser @ 1354'
            stop 'program aborted in read_par_parser @ 1355'
      endif
      do ii=1,classcount-2
        if(akfs(ii).le.0.0)then
            print*,'akfs(',ii,')<= 0.0',akfs(ii)
            write(98,*)'Error: akfs(',ii,')<= 0.0',akfs(ii)
              write(98,*)'Error: program aborted in read_par_parser @ 1361'
            stop 'program aborted in read_par_parser @ 1361'
        endif
      end do
      if(akfs(classcount).le.0.0)then
            print*,'akfs(',classcount,')<= 0.0',akfs(classcount)
            write(98,*)'Error: akfs(',classcount,')<= 0.0',akfs(classcount)
            write(98,*)'Error: program aborted in read_par_parser @ 1368'
            stop 'program aborted in read_par_parser @ 1369'
      endif
      

      fratioflg1=.false.
      do ii=1,classcount
        if(fratio(ii).lt.0.001.or.fratio(ii).gt.10.0)then
          fratioflg1=.true.
          print*,'WARNING:'
          print*,'fratio for class',ii,'=',fratio(ii),' out of range 0.01 - 10.0'
        endif
      end do
      if(fratioflg1)then
        print*,'Please correct h & fratio values'
            write(98,*)'Error: Please correct h & fratio values'
            write(98,*)'Error: fratio for class',ii,'=',fratio(ii),' out of range 0.01 - 10.0'
            write(98,*)'Error: program aborted in read_par_parser @ 1385'
        stop 'Program aborted in read_par_parser @ 1386'
      endif
          
      fratioflg1=.false.
      do ii=1,classcount
        if(fratiohgh(ii).gt.10.0)then
          fratioflg1=.true.
          print*,'WARNING:'
          print*,'fratiohgh for class',ii,' above range 0.01 - 10.0'
        endif
      end do
      if(fratioflg1)then
        print*,'Please correct fratiohgh values'
            write(98,*)'Error: Please correct fratiohgh values'
            write(98,*)'Error: fratiohgh for class',ii,' above range 0.01 - 10.0'
            write(98,*)'Error: program aborted in read_par_parser @ 1401'
        stop 'Program aborted in read_par_parser @ 1402'
      endif
          
      fratioflg1=.false.
      do ii=1,classcount
        if(fratiolow(ii).lt.0.01)then
          fratioflg1=.true.
          print*,'WARNING:'
          print*,'fratiolow for class',ii,' below range 0.01 - 10.0','Value =',fratiolow(ii)
        endif
      end do
      if(fratioflg1)then
        print*,'Please correct fratiolow values'
            write(98,*)'Error: Please correct fratiolow values'
            write(98,*)'Error: fratiolow for class',ii,' below range 0.01 - 10.0','Value =',fratiolow(ii)
            write(98,*)'Error: program aborted in read_par_parser @ 1418'
        stop 'Program aborted in read_par_parser @ 1418'
      endif
          
!     rev. 10.4.50 Feb.  07/22  = NK Use of A7 in melt.f and remove constraint on A7
      if(a7.lt.0.50.or.a7.gt.0.99)then
        print*,' Value for weighting factor - old vs. new sca value'
        print*,' a7 not between 0.5-0.99in the par file'
!        print*,' a7=0.50 assumed - matches spl8 results'
!        write(63,9811)
        write(63,9812)
        write(63,9813)
        write(63,9814)
!        write(63,9815)
!        write(63,9811)
9811    format(' ')
9812    format('WARNING - in rdpar')
9813    format('WARNING -Value for weighting factor - old vs. new sca value')
9814    format('WARNING -a7 not between 0.5-0.99in the par file')
9815    format('WARNING -a7=0.50 assumed - matches spl8 results')
!        a7=0.5
!        tempflg='y'
      endif
!     rev. 9.6.01  Mar.  01/10  - NK: DDS capability added
      if(a3.lt.0.0.and.dds_flag.eq.1)then
        print*,'a3 < 0   Must be = or > 0 to use DDS'
        print*,'Please correct the par file'
        a3=0.05
        print*,'a3=0.05 assumed'
        print*
      endif
      if(a4.lt.0.0.and.dds_flag.eq.1)then
        print*,'a4 < 0   Must be = or > 0 to use DDS'
        print*,'Please correct the par file'
        a4=0.03
        print*,'a3=0.05 assumed'
        print*
      endif

!     rev. 9.03 jan 07/01 NK
      if(a12.le.0.0.and.smrflg.eq.'y')then
        print*,' Value for a12 not defined in the paramerter file'
        print*,' a12=0.001 assumed - matches previous results'
        print*,' Paused in rdpar after reading aNN values'
        print*
        tempflg='y'
      endif
      
      if(tempflg.eq.'y')then
        print*,' Hit enter to accept default & continue'
        write(*,'(''paused, type [enter] to continue'')') 
        read (*,*) 
      endif
       
      tempflg='n'
      if(a9.gt.5.0.or.a9.lt.0.3)then
!       change range from 0.33-2.0 to 0.33-5.0 Oct. 20/03 NK
        print*
        print*,' a9 (heat deficit / swe ratio outside range 0.333 - 5.0'
        print*,' a9 = .333 assumed'
        print*,' If this default value is unsuitable, please'
        print*,' edit the bsnm.par file and enter proper value'
        print*,' between 0.333 and 5.0'
        print*
        a9=0.333
        if(debug_output)write(51,5110)a9
        tempflg='y'
      endif

      if(iopt.eq.2)print*, 'in rdpar @ 637'

      if(a10.gt.3.0.or.a10.lt.0.5)then
        print*
        print*,' a10 uz discharge exponent outside range 0.5 - 3.0'
        print*,' a10 = 1.0 assumed'
        print*,' If this default value is unsuitable, please'
        print*,' edit the bsnm.par file and enter proper value'
        print*,' between 0.5 and 3.0'
        print*
        a10=1.0
    if(debug_output)write(51,5111)a10
        tempflg='y'
      endif
      if(tempflg.eq.'y')then
        print*
        print*,' To continue with the default value(s) hit enter'
        print*,' To prevent this message from reoccurring, please enter'
        print*,' values for a9 and/or a10 in the parameter .par file'
        print*
        write(*,'(''paused, type [enter] to continue'')') 
        read (*,*) 
      endif
5110  format(//,' Value for a9 heat deficit to swe outside range',f10.3,' assumed')
5111  format(//,' Value for a10 uz exponent outside range',f10.3,' assumed')

      if(iopt.eq.2)print*, 'in rdpar @ 660'
      
!       ERROR CHECK FOR REVISED LZS CLASSES
        do i=1,nrvr
!       rev. 9.1.37  Mar.  22/03  - Option to turn off leakage by setting LZF < 0.0
          if(flz_o(i).lt.0.0)flz_o(i)=0.1e-30
          if(dds_flag.eq.0.and.numa.eq.0)then
            if(pwr_o(i).lt.0.1.or.pwr_o(i).gt.4.0)then   !  Feb. 22/06  nk
              print*
              print*,'The value of pwr for river class ',nrvr
              print*,'is ',pwr(nrvr),' which' 
              print*,'is not in the range of 0.1 - 4.0'
              print*,'Please check all values'
              print*,'This msg not printed for DDS & PS runs'
              print*
              write(*,'(''paused, type [enter] to continue'')') 
              read (*,*) 
            endif
          endif
        end do

        if(iopt.eq.2)print*, 'in rdpar @ 712'
      
          do i=1,nrvr
            if(widep_o(i).le.0.0)then
              print*,'Width/depth ratio must be specified when using'
              print*,'Manning n as the roughness parameter.'
              print*,'Please fix the par file'
              print*
              write(98,*)'Error: Width/depth ratio must be specified when using'
              write(98,*)'Error: Manning n as the roughness parameter.'
              write(98,*)'Error: Please fix the par file'
              write(98,*)'Error: program aborted in read_par_parser @ 1543'
              stop 'Program aborted in rdpar @ 1543'
            endif
          end do
      
      if(wetflg.eq.'y')then
!       parameter checking
        do i=1,nrvr
          if(abs(theta_o(i)).lt.0.00001)then
            print*,' Value for theta for river class ',i
            print*,' too close to zero. Please provide reasonable'
            print*,' value  (0.01 to 0.5)'
            write(98,*)'Error:  Value for theta for river class ',i
            write(98,*)'Error:  too close to zero. Please provide reasonable'
            write(98,*)'Error:  value  (0.01 to 0.5)'
            tempflg='y'
          endif
          if(abs(widep_o(i)).lt.1.0)then
            print*,' Value for width/depth ratio for river class',i
            print*,' too close to zero. Please provide reasonable'
            print*,' value say (1 to 100)'
            write(98,*)'Error:  Value for width/depth ratio for river class',i
            write(98,*)'Error:  too close to zero. Please provide reasonable'
            write(98,*)'Error:  value say (1 to 100)'
            tempflg='y'
          endif
          if(abs(kcond_o(i)).lt.0.99e-09)then
            print*,' Value for kcond for river class ',i
            print*,' too close to zero. Please provide reasonable'
            print*,' value  (1.0e-09 to 1.00)'
            write(98,*)'Error:  Value for kcond for river class ',i
            write(98,*)'Error:  too close to zero. Please provide reasonable'
            write(98,*)'Error:  value  (1.0e-09 to 1.00)'
            tempflg='y'
          endif
        end do
        if(tempflg.eq.'y')then
            write(98,*)'Error: program aborted in read_par_parser @ 1580'
            stop 'Program aborted in rdpar @ 1581'
        endif
    endif       ! wetflg = y

    
    

! Conversion to f90 stopped here. The rest of the error messsages shoule only happen during program setup    
    
    
    
!     rev. 9.2.35  Mar.  22/06  - NK: Glacier flow bypasses wetlands
!     glacier_class_number will be used in tracer
      glacier_class_number=0  !default value

      do i=1,classcount      ! added 09/03/05  nk
        if(nclass(i)(1:7).eq.'GLACIER')nclass(i)='glacier   '
        if(nclass(i)(1:7).eq.'glacier')nclass(i)='glacier   '
        if(nclass(i)(1:7).eq.'Glacier')nclass(i)='glacier   '
        if(nclass(i)(1:8).eq.' glacier')nclass(i)='glacier   '
        if(nclass(i)(1:8).eq.' GLACIER')nclass(i)='glacier   '
        if(nclass(i)(1:8).eq.' Glacier')nclass(i)='glacier   '
        if(nclass(i)(1:7).eq.'glacier')glacier_class_number=i
        if(nclass(i)(1:7).eq.'WETLAND')nclass(i)='wetland   '
        if(nclass(i)(1:7).eq.'wetland')nclass(i)='wetland   '
        if(nclass(i)(1:7).eq.'Wetland')nclass(i)='wetland   '
        if(nclass(i)(1:8).eq.' wetland')nclass(i)='wetland   '
        if(nclass(i)(1:8).eq.' WETLAND')nclass(i)='wetland   '
        if(nclass(i)(1:8).eq.' Wetland')nclass(i)='wetland   '
        if(nclass(i)(1:5).eq.'WATER')nclass(i)='water     '
        if(nclass(i)(1:5).eq.'water')nclass(i)='water     '
        if(nclass(i)(1:5).eq.'Water')nclass(i)='water     '
        if(nclass(i)(1:6).eq.' water')nclass(i)='water     '
        if(nclass(i)(1:6).eq.' WATER')nclass(i)='water     '
        if(nclass(i)(1:6).eq.' Water')nclass(i)='water     '
      end do

!     check added Jul. 31/06 nk
      if(glacier_class_number.eq.0.and.itrace.eq.1.and.trcflg.eq.'y')then
          trcflg='n'
          print*,'Tracer is set to track glacier melt'
          print*,'and the tracer flag (trcflg) is set to "y"'
          print*,'but there is no glacier class'
          print*,'trcflg is set to "n" '
          print*
          write(*,'(''paused, type [enter] to continue'')') 
          read (*,*) 
      endif
                  
!     rev. 9.08.01 Apr.   3/01  - check wetland designation in rdpar
      if(wetflg.eq.'y')then
      if(nclass(classcount-2)(1:7).eq.'wetland')then
!       everything just fine
!         only class classcount-2 is acceptable as a coupled wetland
      else
        print*
        print*,'Class no. ',classcount-2,' is expected to be'
        print*,'the wetland class and be labelled as such'
        print*,'Name found ***',nclass(classcount-2),'***'
        print*,'It should be ***wetland   ***'
        print*,'Check the tpl file also'
        do i=1,classcount
            print*,i,nclass(i)
        end do
        print*,'Please be sure that the order of the parameters'
        print*,'and the land cover data are in the same order'
        print*,'Also: no  /  allowed in par file'
        print*,'Please fix the parameter file and try again'
        print*
        stop 'Program aborted in read_par_parser at line ~217'
      endif
      endif
         
      do i=1,classcount
!   if(rec(i).gt.1.00.or.rec(i).lt.0.0)then
        if(rec(i).lt.0.0)then
!         note: larger value will give NaN for x4 in runof6 
          print*,'the value for rec can not be less than 0.00'
          print*,'class=',i,' value=',rec(i)
          print*,'Please change the value in the par file'
          print*
          stop 'Program aborted in read_par_parser @ 1040'
        endif
      end do
      ijunko=1

      if(ak(classcount-1).ge.0.0.or.akfs(classcount-1).ge.0.0)then
        print*
        print*,'Error:'
        print*,'Water class not specified'
        print*,'ak and akfs NOT given a -ve value'
        print*,'ak(',classcount-1,') value found =',ak(classcount-1)
        print*,'akfs(',classcount-1,') value found =',akfs(classcount-1)
        print*,'Please correct the par file'
        print*,'No of classed expected =',classcount
        print*
!     REV. 10.1.23 Jan.  28/16  - NK: Added abort when water class not specified
        stop 'Program aborted in read_par_parser @ 1514'
      endif
      
!   rev. 9.5.65  Sep.  26/09  - NK: lapse rate changed from dC per 100 m to dC per m
!  if(rlapse.gt.0.1)then
!    print*,'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
!    print*,'WARNING:'
!    print*, 'Value for RLAPSE is too large'
!    print*, 'RLAPSE used to be mm per 100m'
!    print*, 'this has been changed to mm per km'
!    print*, 'Please correct the par file'
!    print*,'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
!    stop 'Program aborted in rdpar @ 1185'
!  endif

      if(iopt.eq.2)print*, 'param 9'
      if(fmadjust.gt.0.0)then
         if(fmalow.eq.0.0.or.fmahigh.eq.0.0)then
            print*,'fmadjust,fmalow,fmahigh/',fmadjust,fmalow,fmahigh
            print*,'fmadjust is set to 0.0'
            print*,' because there are no limits set in bsnm.par'
            fmadjust=0.0
            write(*,'(A)',advance='no')' hit enter to continue'
            read(*,*)
         endif
      endif
      
        if(fpet(classcount-1).le.0.0)then
          print*
          print*,'WARNING'
          print*,'Value for fpet for water class =',fpet(classcount-1)
          print*,'This value is multiplied by the pet for water'
          print*,'so for this value there will be no evaporation'
          print*,'from water bodies'
          print*,'Reasonable values are 0.5-1.0'
          print*,'Please change the par file accordingly'
          print*
!     stop 'Program aborted in rdpar @ 1356'
        endif
      
        do i=1,classcount
          if(ffcap(i).le.0.0)then
            print*,'pwp(',i,') Permanent Wilting Point set too low.'
            print*,'pwp(',i,')  Must be > 0.0'
            print*,'pwp(',i,') set to 0.01'
            print*
            ffcap(i)=0.01
          endif
        end do
      
        do ii=1,classcount
          if(spore(ii).le.0.0)then
            print*,'spore(',ii,')=',spore(ii)
            print*,' This value can not be .le. 0.0'
            print*,' Please correct the parameter file'
            print*
            stop 'Program aborted in rdpar at ~1254'
          endif
        end do
      
!     rev. 9.7.29  Jul.  07/11  - NK: Add sublim_rate to set sublimation rate/day to pat file
        if(sublimflg.eq.1)then          ! for sublim_factor only
          do i=1,classcount
            if(sublim_factor(i).gt.0.75)then
             if(numa.eq.0.and.dds_flag.eq.0)then    !added Apr.03/10 nk
                print*,'Warning: sublimation factor for class',i
                print*,'and is given as ',sublim_factor(i)
                print*,'is set too high. Sublim reduced to 0.75'
              endif
            endif
            sublim_factor(i)=amin1(0.75,sublim_factor(i))
          end do
        endif
      
         if(tempa3.lt.0.0001)then
          tempa3=0.0001
      if(debug_output)write(51,*)'WARNING:'
      if(debug_output)write(51,*)'Temp3 <0.0001  set temp3=0.001'
         if(debug_output)write(51,*)
        endif
        
      do ii=1,classcount
        if(spore(ii).le.0.0)then
          print*,'ii,spore',ii,spore(ii)
          print*,' line title =',title(39)
          print*,' value for spore can not be <= 0.0'
          STOP 'Program stopped in rdpar @ 9443'
        endif
      end do

      do ii=1,classcount

!       adjust potential evaporation from vegetation by amount
        hmax=0.0  
!     rev. 9.2.22  Nov.  15/05  - NK: Fixed hmax bug in rdpar 
!             hmax=0 was erroniously moved down between 
!                    ver 9.2.03 & 9.2.2
        do i=1,12
!     rev  9.1.25  Sep.  11/02  - Added A11 as bare ground equiv. vegn height  
!          h(i,ii)=h(i,ii)+a11    ! taken out and replaced Jan. 25/06 nk
         h(i,ii)=amax1(a11,h(i,ii)) !  Jan. 25/06 nk
         hmax=amax1(hmax,h(i,ii))
        end do 

        do i=1,12
!         problems here nk Nov. 20/02
!         rev. 9.1.32  Nov.  20/02  - Fixed fpetmon() wrt. h()
!         h(month,class) is defined as the max interception storage. 
!         So for soil evaporation it has to become an index 
!         fpetmo() is that index - it modifies ftall(ii) to account for veg height

          fpetmo(i,ii)=ftall(ii)*h(i,ii)/hmax
!          fpetmo(i,ii)=ftall(ii)*h(i,ii)

        end do
        if(iopt.gt.1)then
       if(debug_output)write(51,1003)'hmax ',hmax
       if(debug_output)write(51,1003)'h(  )',(h(i,ii),i=1,12)
       if(debug_output)write(51,1003)'fpetm',(fpetmo(i,ii),i=1,12)
 1003   format(a5,99f5.2)
        endif
      end do

      if(iopt.eq.2)print*, 'param 10a'

!     TO PREVENT DIVISION BY ZERO IN ETIN
      do ii=1,classcount
        do i=1,12
          h(i,ii)=amax1(0.01,h(i,ii))
        end do
      end do

       do ii=1,classcount
        if(ak2(ii).gt.1.0)then
          write(*,1024)ii
1024      format(' AK2(',i5,') =',e12.3)
          STOP ' Please change value for AK2- must be lt 1'
        endif
        if(ak2fs(ii).gt.1.0)then
          write(*,1027)ii
1027      format(' AK2fs(',i5,') =',e12.3)
          STOP ' Please change value for AK2fs - must be lt 1'
        endif
      end do
    
      do i=1,nbsn
        if(flz_o(i).le.0.0)then
          write(*,8429)flz_o(i)
 8429     format('In read_par,flz_o(',i2,')=',f12.7,' wrong parameter value - change flz(i) to +ve value')
        endif
      end do
       
!     SET POTENTIALS FOR EACH SOIL TYPE:
!     moved from soilinit   Mar. 14/07  nk
      do ii=1,classcount
!        NOTE FOR AK<0.0 WE HAVE WATER AREA & POT ISN'T USED
         if(ak(ii).gt.0.0)then
!           AK IS IN MM/HR
!           BUT FORMULA IS FOR MM/SEC
            xxx1=-alog10(ak(ii)/3600.)
            pot(ii)=250.*xxx1+100.
         endif
!        REV. 7.5 SEPERATE SNOW COVERED AND BARE GROUND
         if(akfs(ii).gt.0.0)then
!           ak is in mm/hr
!           but formula is for mm/sec
            xxx1=-alog10(akfs(ii)/3600.)
            potfs(ii)=250.*xxx1+100.
         endif
      end do
      
!     rev. 9.8.68  Jun   17/13  - NK: Added dds_override file
!     A file with only a zero in col 1 and line 1 will allow the 
!     copying of a par file from an ongoing DDS run and be run
!     in a normal production mode
      if(dds_flag.gt.0)then
        inquire(FILE='ddsflg_override.txt',EXIST=exists)
        if(exists)then
          open(unit=99,file='ddsflg_override.txt',iostat=ios)
            read(99,*)dds_flag
            print*,'---------------------------------------------'
            print*,'DDS_flag set to `0` for normal production run'
            print*
          close(unit=99,status='keep')
        endif
      endif

!      if(.not.debug_output)then
!          if(dds_flag.eq.1.or.netCDFflg)then
!          if(dds_flag.eq.1)then
!              iopt=0
!              iopt99=.false.
!          endif
!          if(iopt.ge.1)iopt99=.true.
!      endif
      
!     rev. 10.2.22 May   10/18  - NK: Set wetland classes uncoupled and coupled to same parameters
      if(dds_flag.eq.1)then
        if(nclass(ii-3)(1:7).eq.'wetland'.and.nclass(ii-4)(1:7).eq.'wetland')then    ! i.e. 2 wetland classes with same name   
          fratio(ii-3)=fratio(ii-4)
          fm(ii-3)=fm(ii-4)
          base(ii-3)=base(ii-4)
          sublim_rate(ii-3)=sublim_rate(ii-4)
          print*,'fratio, fm, base & subl set equal for wl classes'
        endif
      endif

!     This was fixed between version 9.1.11 and 9.1.13
!     It used to be all sin(...). End result: pet too high

!     rev. 9.9.03  Dec.  15/13  - NK: Change to gridded latitude for etharg
      if(.not.allocated(sinlat))then
        allocate(sinlat(na),coslat(na),tanlat(na),stat=iAllocate)
        if(iAllocate.ne.0) STOP 'Error with allocation of rivtype in rdpar'
      endif
      
      if(flgevp2.eq.2.0)then
        if(lat.eq.0)then
          print*
          print*,'You have lat = 0.0  and flgevp2 = 2.0'
          print*,'in the par file'
          print*,'flgevp2 = 2.0 requires a value for the lat'  
          print*,'to use original monthly_climate_normals file'
          print*,'the value of lat in the par file must be set as'
          print*,'the mid north-south latitude of the watershed'
          print*
          stop 'Program aborted in read_par_parser @ L1713'
        endif
        sinlat(1)=sin(3.1416*lat/180.0)
        coslat(1)=cos(3.1416*lat/180.0)
        tanlat(1)=tan(3.1416*lat/180.0)
!       this is the original way for using a fixed lat (centre of the grid usually)        
!d       print*,'   n            lat           sinlat         coslat   ',
!d    *              '       tanlat'
        do n=2,na
          sinlat(n)=sinlat(1)
          coslat(n)=coslat(1)
          tanlat(n)=tanlat(1) 
          if(debug_output)write(63,*)n,lat,sinlat(n),coslat(n),tanlat(n)
        end do
!d       print*,'   n            lat           sinlat         coslat   ',
!d    *              '       tanlat'
      elseif(flgevp2.eq.4.0)then
        print*,coordsys1
        if(coordsys1(1:3).eq.'UTM'.or.coordsys1(1:8).eq.'CARTESIAN')then
          if(lat.eq.0.0)then
            print*
            print*,'lat = 0 in the par file'
            print*,'lat must have the value of the latitude of the'
            print*,'midpoint of the watershed for UTM or Cartesian '
            print*,'coordinates'
            print*
            stop 'Program aborted in read_par_parser @ L1739'
          endif
        endif
        print*
        print*,'NEW:   flgevp2 = 4 option'
        print*,'the general mid-basin lat parameter in the par file'
        print*,'will be replaced by the lat for each grid.'
        print*
        if(coordsys1(1:7).eq.'LATLONG')then
!           print*,'   n            lat           sinlat         ',
!      *              'coslat          tanlat'
          do n=1,na
!           calculate the functions for each grid based on grid latitude        
            lat=yorigin+(float(yyy(n))-0.5)*ydelta
            sinlat(n)=sin(3.1416*lat/180.0)
            coslat(n)=cos(3.1416*lat/180.0)
            tanlat(n)=tan(3.1416*lat/180.0)
!             print*,n,lat,sinlat(n),coslat(n),tanlat(n)
          end do
!         print*,'   n            lat           sinlat         coslat   ',
!      *              '      tanlat'
        else       ! coordsys = UTN or CARTESIAN
!     rev. 9.9.03  Dec.  15/13  - NK: Change to gridded latitude for etharg
          print*
          print*,'NEW <<<<'
          print*,'variable lat is calculated based on the latitude'
          print*,'of the origin given in the par file as',lat
!           calculate the functions for each grid based on grid latitude   
          lat_origin=lat              ! given in the par file
          lat_delta=ydelta/111120.0   ! 1 degree lat = 1852*60 =111120 m     
          lat_max=-90.0
          do n=1,na
            lat=lat_origin+(float(yyy(n))-0.5)*lat_delta
            lat_max=amax1(lat_max,lat)
            sinlat(n)=sin(3.1416*lat/180.0)
            coslat(n)=cos(3.1416*lat/180.0)
            tanlat(n)=tan(3.1416*lat/180.0)
!             print*,n,lat,sinlat(n),coslat(n),tanlat(n)
          end do
          lat=lat_origin              ! give back its value in the par file 
          print*,'north-most latitude:',lat_max
          print*,'Please check'
          print*
        endif
!    print*,'NEW local latitude in read_par_parser @ 1722'
      endif   ! (flgevp2.eq.    )
            
      if(debug_output)write(98,*)'Info: read_par_parser before return @ 9989'

 9989 RETURN

      END SUBROUTINE read_par_parser

