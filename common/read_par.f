      SUBROUTINE read_par(unitNum,flnNum)

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
!  rev. 9.2.09  Sep.  11/05  - NK: removed write_par.for from rdpar.for
!  rev. 9.5.43  Oct.  27/08  - NK: changed bottom part of par file to be free format
!  rev. 9.5.45  Dec.  16/08  - NK: added various error calculations - user's choice with errflg
!  rev. 9.6.01  Mar.  01/10  - NK: DDS capability added
!  rev. 9.6.06  Apr.  18/10  - NK: added glacier adjust for optimization
!  rev. 9.8.01  Jul.  21/11  - NK: added ragmet optimization to dds setup
!
!  s/r created Sept. 11/05 to separate writing from reading 
!  in the rdpar.for subroutine. 
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
	INTEGER       :: nrvr1,nchr,iprtflg,linecount,i,unitNum,flnNum,
     *                 ios,iverflg,ix,iallocate,ijunko,ii,j,numb
      real*4        :: hmax,e1,ajunk,kcondflg,xxx1
	character(1)  :: firstpass
	character(256):: line
	logical       :: exists


c      DATA ntest/-10906/qstr/'param'/nchr/5/


      data firstpass/'y'/

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

c      i=max(nrvr,classcount)
c      allocate(ajunk(i),stat=iAllocate)
c      if(iAllocate.ne.0) STOP 'Error with allocation of ajunk in rdpar'

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
        if(iAllocate.ne.0) STOP
     *   'Error with allocation of rivtype in rdpar'
      endif

!     TS - ALLOCATIONS OF AREA4A ARRAYS
!     classcount for the number of land cover classes
!     nbsn for the number of channel or basin types
!     moved here from rdshed 27/07/06 because needed for bsn.for  nk
!     rev. 9.5.27  Apr.  15/08  - NK: fixed allocation for chnl in rdpar
      if(.NOT.allocated(ds))then
      allocate(
     *  ds(classcount),dsfs(classcount),chnl(5),  ! note: chnl = always 5
     *  r3(classcount),r4(classcount),r3fs(classcount),rec(classcount),
     *  ak(classcount),akfs(classcount),
     *  r3low(classcount),r3fslow(classcount),reclow(classcount),
     *  aklow(classcount),akfslow(classcount),ak2fslow(classcount),
     *  r3hgh(classcount),r3fshgh(classcount),rechgh(classcount),
     *  akhgh(classcount),
     *  akfshgh(classcount),ak2fshgh(classcount),r3dlt(classcount),
     *  r3fsdlt(classcount),recdlt(classcount),akdlt(classcount),
     *  akfsdlt(classcount),ak2fsdlt(classcount),retn(classcount),
     *  ak2(classcount),
     *  retnlow(classcount),ak2low(classcount),
     *  retnhgh(classcount),ak2hgh(classcount),
     *  retndlt(classcount),ak2dlt(classcount),
     *  retfs(classcount),ak2fs(classcount),fpet(classcount),
     *  fpetdlt(classcount),fpetlow(classcount),fpethgh(classcount),
     *  ftall(classcount),ftalldlt(classcount),ftalllow(classcount),
     *  ftallhgh(classcount),nclass(classcount),
     *  fratio(classcount),fratiodlt(classcount),
     *  fratiolow(classcount),fratiohgh(classcount),
     *  pot(classcount),potfs(classcount),   ! moved from read_shed_ef may 15/07 nk
     *  iiclass(classcount*2),h(12,classcount),fpetmo(12,classcount),
     *                        stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *   'Warning: error with allocation of area4a arrays @ 186'
	endif

!     TS - ALLOCATIONS OF AREAOPTSA ARRAYS
      if(.NOT.allocated(fmdlt))then
        allocate(
     *  fmdlt(classcount),fmlow(classcount),fmhgh(classcount),
     *  fmndlt(classcount),
     *  fmnlow(classcount),fmnhgh(classcount),uajdlt(classcount),
     *  uajlow(classcount),
     *  uajhgh(classcount),mbsdlt(classcount),mbslow(classcount),
     *  mbshgh(classcount),
     *  basdlt(classcount),baslow(classcount),bashgh(classcount),
     *  subdlt(classcount),sublow(classcount),subhgh(classcount),
     *  stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *   'Warning: error with allocation of areaoptsa arrays @ 196'
      endif

!     TS - ALLOCATIONS OF AREAMELTA ARRAYS (PARTIAL)
!     SNW,DSN,TTEMP,TMX,TMN,EL ALLOCATED IN SHEDA.FOR
!     SDCD,SDCSCA ALLOCATED IN RDSDCA.FOR
!     allocate(snowc(na,classcount),dsnow(na),tempv(na),tempvmin(na),
!     moved here from rdshed 27/07/06 because needed for bsn.for  nk
      if(.NOT.allocated(base))then
        allocate(base(classcount),fm(classcount),fmn(classcount),
     *  whcl(classcount),tipm(classcount),uadj(classcount),
     *  rho(classcount),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     * 'Warning: error with allocation of areamelta arrays in read_par'
	endif

!     TS - ALLOCATIONS OF AREAETA ARRAYS (PARTIAL)
!     RAD ALLOCATED IN SHEDA.FOR
!     TS - ADDED ALLOCATIONS FOR EVAP-SEPARATION PARAMS (22/03/06)
!     TS: CHANGED ALLOCATIONS OF alb,pet, evap TO classcount (27/03/06)
!     rev. 9.1.80  Mar.  31/05  - NK: added sublimation   (sublim)
!     moved here from rdshed 27/07/06 because needed for bsn.for  nk
      if(.NOT.allocated(evap))then
        allocate(evap(classcount,12),sublim_factor(classcount),
     *  sublim_rate(classcount),flint(classcount),
     *  diff(12),hu(12),ffcap(classcount),fcap(classcount),
     *  spore(classcount),alb(classcount),pres(12),stat=iAllocate)
 !      *rh(na),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *   'Warning: error with allocation of evap arrays in read_par'
	endif



c        allocate(nsdc(classcount),snocap(classcount),idump(classcount),
c     *        stat=iAllocate)
c        if(iAllocate.ne.0) STOP
c     *   'Error with allocation of areamelta arrays in spl9 @ 158'


      if(.not.allocated(r2nlow))then
        allocate(r2nlow(nrvr),r2nhgh(nrvr),r2ndlt(nrvr),        
     *         thetadlt(nrvr),thetalow(nrvr),thetahgh(nrvr),
     *         widepdlt(nrvr),wideplow(nrvr),widephgh(nrvr),
     *         kconddlt(nrvr),kcondlow(nrvr),kcondhgh(nrvr),
     *         flzlow(nrvr),flzhgh(nrvr),flzdlt(nrvr),
     *         pwrlow(nrvr),pwrhgh(nrvr),pwrdlt(nrvr),
     *         rlakelow(nrvr),rlakehgh(nrvr),rlakedlt(nrvr),
     *         stat=iAllocate)
        if(iAllocate.ne.0)then
          print*,'Error with allocation of area4 arrays in rdpar @ 524'
          print*,'Try setting numa = 0 in the par file'
          print*
          stop 'program aborted in rdpar @ 626'
        endif
	endif

c      endif

      if(.not.allocated(flz_o))then
        allocate(flz_o(nrvr),pwr_o(nrvr),r1n_o(nrvr),        
     *         r2n_o(nrvr),mndr_o(nrvr),aa2_o(nrvr),
     *         aa3_o(nrvr),aa4_o(nrvr),theta_o(nrvr),
     *         widep_o(nrvr),kcond_o(nrvr),pool_o(nrvr),
     *         rlake_o(nrvr),
     *  stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *        'Error with allocation of area4 arrays in rdpar -2'
	endif

      if(.not.allocated(sdcd))then
c       nnts=nsdc(classcount)
        nnts=2
        allocate(sdcd(nnts,classcount),sdcsca(nnts,classcount),
     *                  stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *   'Error with allocation of areamelta arrays in read_par' 
      endif

c      print*,'Allocations done in rdpar',classcount,nrvr

      endif  !firtspass

!     ADDED THIS FROM TODD'S ETPAR.FOR - FRANK S: NOV/97 
      alamb=2.478
      den=999.
      alpha=1.35
      dds_flag=0 ! default for old par files
      type1=float(itype)
c      itype=int(type1)

!     rev. 9.8.08  Nov.  18/11  - NK: added fratio for interception hight optimization
      do ii=1,classcount
        fratio(ii)=1.0
        fratiodlt(ii)=-1.0
        fratiolow(ii)=0.1
        fratiohgh(ii)=10.0
      end do

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

      inquire(FILE=fln(flnNum),EXIST=exists)
	if(.not.exists)then
	  print*,fln(flnNum)
	  print*,'not found in read_par @ 245'
	  pause 'program will abort with enter'
	  stop 'aborted in read_par'
	endif

      open(unitNum,file=fln(flnNum),
     *        	status='unknown',iostat=ios)
      if(ios.ne.0)then
        print*,' Error opening',fln(flnNum)
	  print*,'on unit=',unitNum
        print*
        stop ' Program aborted in read_par_10 @ ~104'
	endif

d	print*,'Opened unit',unitNum,'fln=',fln(flnNum)(1:30)	

      numb=0  ! so new.par can be copied to bsnm.par for batch job

	rewind (unit=unitNum)
256   format(a256)
      read(unitNum,256)line
d     print*,line(1:60)
      read(unitNum,256)line
d     print*,line(1:60)
      read(unitNum,256)line
d     print*,line(1:60)

      read(unitNum,*,iostat=ios)junk,iopt          !# debug level
	if(.not.debugflg)iopt=0                      !# needed for sensitivity
      read(unitNum,*,iostat=ios)junk,itype         !# channel type - floodplain/no 
      read(unitNum,*,iostat=ios)junk,itrace        !# Tracer choice
      read(unitNum,*,iostat=ios)junk,a1            !# ice cover weighting factor
      read(unitNum,*,iostat=ios)junk,a2            !# Manning`s correction for instream lake
      read(unitNum,*,iostat=ios)junk,a3            !# error penalty coefficient
      read(unitNum,*,iostat=ios)junk,a4            !# error penalty threshold
      read(unitNum,*,iostat=ios)junk,a5            !# API coefficien
      read(unitNum,*,iostat=ios)junk,a6            !# Minimum routing time step in seconds
      read(unitNum,*,iostat=ios)junk,a7            !# weighting - old vs. new sca value
      read(unitNum,*,iostat=ios)junk,a8            !# min temperature time offset
      read(unitNum,*,iostat=ios)junk,a9            !# max heat deficit /swe ratio')
      read(unitNum,*,iostat=ios)junk,a10           !# exponent on uz discharce function')
      read(unitNum,*,iostat=ios)junk,a11           !# bare ground equiv. veg height for ev')
      read(unitNum,*,iostat=ios)junk,a12           !# min precip rate for smearing')
      read(unitNum,*,iostat=ios)junk,fmadjust      !# snowmelt ripening rate')
      read(unitNum,*,iostat=ios)junk,fmalow        !# min melt factor multiplier')
      read(unitNum,*,iostat=ios)junk,fmahigh       !# max melt factor multiplier')
      read(unitNum,*,iostat=ios)junk,gladjust      !# glacier melt factor multiplier')
      read(unitNum,*,iostat=ios)junk,rlapse        !# precip lapse rate fraction/m')
      read(unitNum,*,iostat=ios)junk,tlapse        !# temperature lapse rate dC/100m')
      read(unitNum,*,iostat=ios)junk,elvref        !# reference elevation')
      read(unitNum,*,iostat=ios)junk,rainSnowTemp  !# rain/snow temperature')
      read(unitNum,*,iostat=ios)junk,radinfl       !# radius of influence km')
      read(unitNum,*,iostat=ios)junk,smoothDist    !# smoothing diatance km')
      read(unitNum,*,iostat=ios)junk,flgevp2       !# 1=pan;2=Hargreaves;3=Priestley-Taylor')
      read(unitNum,*,iostat=ios)junk,albe          !# albedo????')
      read(unitNum,*,iostat=ios)junk,tempa2        !# albedo????')
      read(unitNum,*,iostat=ios)junk,tempa3        !# albedo????')
      read(unitNum,*,iostat=ios)junk,tton          !# albedo????')
      read(unitNum,*,iostat=ios)junk,lat           !# albedo????')
      read(unitNum,*,iostat=ios)junk,chnl(1)       !# manning`s n multiplier')
      read(unitNum,*,iostat=ios)junk,chnl(2)       !# manning`s n multiplier')
      read(unitNum,*,iostat=ios)junk,chnl(3)       !# manning`s n multiplier')
      read(unitNum,*,iostat=ios)junk,chnl(4)       !# manning`s n multiplier')
      read(unitNum,*,iostat=ios)junk,chnl(5)       !# manning`s n multiplier')
      read(unitNum,*,iostat=ios)line
d     print*,line(1:60)
	if(ios.ne.0)then
        print*,'error in :GlobalParameters'
	  print*
	  pause 'hit enter to abort the program'
	  stop 'Program aborted in read_par @ 322'
	endif

      read(unitNum,*,iostat=ios)line
d     print*,line(1:60)
      read(unitNum,*,iostat=ios)line    !'# Routing Parameters'
d     print*,line(1:60)
      read(unitNum,*,iostat=ios)junk,nbsn  
d	print*,junk(1:15),nbsn
      if(nbsn.ne.nrvr)then
        print*
        print*,'No of river classes in the par file does not match'
        print*,'the number in the shd file'
        print*,'shd file number =',nrvr
        print*,'par file number =',nbsn
        print*
        stop 'Fatal error: program aborted in read_par @ 173'
	endif
      read(unitNum,*,iostat=ios)junk,(rivtype(i),i=1,nrvr) !#RiverClassName, 
      read(unitNum,*,iostat=ios)junk,(flz_o(i),i=1,nrvr)   !# lower zone oefficient')	  
      read(unitNum,*,iostat=ios)junk,(pwr_o(i),i=1,nrvr)   !# lower zone exponent')	  
      read(unitNum,*,iostat=ios)junk,(r1n_o(i),i=1,nrvr)   !# overbank Manning`s n')
      read(unitNum,*,iostat=ios)junk,(r2n_o(i),i=1,nrvr)   !# channel Manning`s n')	  
      read(unitNum,*,iostat=ios)junk,(mndr_o(i),i=1,nrvr)  !# meander channel length multiplier')	  
      read(unitNum,*,iostat=ios)junk,(aa2_o(i),i=1,nrvr)   !# channel area intercept = min channel xsect area')	  
      read(unitNum,*,iostat=ios)junk,(aa3_o(i),i=1,nrvr)   !# channel area coefficient')	  
      read(unitNum,*,iostat=ios)junk,(aa4_o(i),i=1,nrvr)   !# channel area exponent')	  
      read(unitNum,*,iostat=ios)junk,(theta_o(i),i=1,nrvr) !# wetland or bank porosity')	  
      read(unitNum,*,iostat=ios)junk,(widep_o(i),i=1,nrvr) !# channel width to depth ratio')	  
      read(unitNum,*,iostat=ios)junk,(kcond_o(i),i=1,nrvr) !# wetland/bank lateral conductivity')	  
      read(unitNum,*,iostat=ios)junk,(pool_o(i),i=1,nrvr)  !# average area of zero flow pools')	  
      read(unitNum,*,iostat=ios)junk,(rlake_o(i),i=1,nrvr) !# in channel lake retardation coefficient')	  
      read(unitNum,*,iostat=ios)line
d     print*,line(1:60)
	if(ios.ne.0)then
        print*,'error in :RoutingParameters'
	  print*
	  pause 'hit enter to abort the program'
	  stop 'Program aborted in read_par @ 356'
	endif

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,classcount_local  !# No land cover classes'
!     check added Mar. 08/11 nk
	if(classcount_local.ne.classcount)then
	  print*,'WARNING:'
	  print*,'no of land cover classes in this par file does not match'
	  print*,'the number in the shd file'
	  print*,'please check'
	  print*
d	  if(iopt.ge.1)then
d         pause 'Hit enter to abort'
d	    stop 'progran aborted in read_par_10 @ 348'
d	  endif
	endif

d	print*,'classcount=',classcount_local
      read(unitNum,*,iostat=ios)junk,(nclass(i),i=1,classcount)       !# class name')	  

      read(unitNum,*,iostat=ios)junk,(ds(i),i=1,classcount)           !# depression storage bare ground mm')	  
      read(unitNum,*,iostat=ios)junk,(dsfs(i),i=1,classcount)         !# depression storage snow covered area mm')	  
      read(unitNum,*,iostat=ios)junk,(rec(i),i=1,classcount)          !# interflow coefficient')	  
      read(unitNum,*,iostat=ios)junk,(ak(i),i=1,classcount)           !# infiltration coefficient bare ground')	  
      read(unitNum,*,iostat=ios)junk,(akfs(i),i=1,classcount)         !# infiltration coefficient snow covered ground')	  
      read(unitNum,*,iostat=ios)junk,(retn(i),i=1,classcount)         !# upper zone retention mm')	  
      read(unitNum,*,iostat=ios)junk,(ak2(i),i=1,classcount)          !# recharge coefficient bare ground')	  
      read(unitNum,*,iostat=ios)junk,(ak2fs(i),i=1,classcount)        !# recharge coefficient snow covered ground')	  
      read(unitNum,*,iostat=ios)junk,(r3(i),i=1,classcount)           !# overland flow roughness coefficient bare ground')	  
      read(unitNum,*,iostat=ios)junk,(r3fs(i),i=1,classcount)         !# overland flow roughness coefficient snow covered grnd')	  
      read(unitNum,*,iostat=ios)junk,(r4(i),i=1,classcount)           !# overland flow roughness coefficient impervious area')	  
      read(unitNum,*,iostat=ios)junk,(fpet(i),i=1,classcount)         !# interception evaporation factor * pet')	  
      read(unitNum,*,iostat=ios)junk,(ftall(i),i=1,classcount)        !# reduction in PET for tall vegetation')	  
      read(unitNum,*,iostat=ios)junk,(flint(i),i=1,classcount)        !# interception flag  1=on  <1=off')	  
      read(unitNum,*,iostat=ios)junk,(fcap(i),i=1,classcount)         !# not used - replaced by retn (retention)')	  
      read(unitNum,*,iostat=ios)junk,(ffcap(i),i=1,classcount)        !# wilting point - mm of water in uzs')	  
      read(unitNum,*,iostat=ios)junk,(spore(i),i=1,classcount)        !# soil porosity')	  
      read(unitNum,*,iostat=ios)line   !:EndHydrologicalParameters
d     print*,line(1:60)
	if(ios.ne.0)then
        print*,'error in :HydrologicalParameters'
	  print*
	  pause 'hit enter to abort the program'
	  stop 'Program aborted in read_par @ 401'
	endif

      read(unitNum,*,iostat=ios)   
      read(unitNum,*,iostat=ios)   !:SnowParameters'
      read(unitNum,*,iostat=ios)junk,(fm(ii),ii=1,classcount)         !# melt factor mm/dC/hour')	  
      read(unitNum,*,iostat=ios)junk,(base(ii),ii=1,classcount)       !# base temperature dC')	  
      read(unitNum,*,iostat=ios)junk,(fmn(ii),ii=1,classcount)        !# -ve melt factor')	  
      read(unitNum,*,iostat=ios)junk,(uadj(ii),ii=1,classcount)       !# not used')	  
      read(unitNum,*,iostat=ios)junk,(tipm(ii),ii=1,classcount)       !# coefficient for ati')	  
      read(unitNum,*,iostat=ios)junk,(rho(ii),ii=1,classcount)        !# snow density')	  
      read(unitNum,*,iostat=ios)junk,(whcl(ii),ii=1,classcount)       !# fraction of swe as water in ripe snow')	  
      read(unitNum,*,iostat=ios)junk,(alb(i),i=1,classcount)          !# albedo')	  

!     rev. 9.7.29  Jul.  07/11  - NK: Add sublim_rate to set sublimation rate/day to pat file
      read(unitNum,*,iostat=ios)junk,(sublim_factor(i),i=1,classcount)!# sublimation factor ratio')	  
      sublimflg=1  ! default for sublim_factor
      if(junk(1:12).eq.':sublim_rate')then
        sublimflg=2
	  do ii=1,classcount
	    sublim_rate(ii)=sublim_factor(ii)/24.0
	    sublim_factor(ii)=0.0
	  end do
	endif

      read(unitNum,*,iostat=ios)junk,(idump(ii),   ii=1,classcount)   !# receiving class for snow redistribution')	  
      read(unitNum,*,iostat=ios)junk,(snocap(ii),  ii=1,classcount)   !# max swe before redistribution')	  
      read(unitNum,*,iostat=ios)junk,(nsdc(ii),    ii=1,classcount)   !# no of points on scd curve - only 1 allowed')	  
      read(unitNum,*,iostat=ios)junk,(sdcsca(2,ii),ii=1,classcount)   !# snow covered area - ratio=1.0')	  
      read(unitNum,*,iostat=ios)junk,(sdcd(2,ii),  ii=1,classcount)   !# swe for 100% snow covered area')	  
      read(unitNum,*,iostat=ios)line      !:EndSnowParameters
d     print*,line(1:60)
	if(ios.ne.0)then
        print*,'error in :SnowParameters'
	  print*
	  pause 'hit enter to abort the program'
	  stop 'Program aborted in read_par @ 436'
	endif

      do ii=1,classcount
        sdcsca(1,ii)=0.0
        sdcd(1,ii)=0.0
      end do

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

 
c      read(unitNum,98229)title(45),(diff(i),i=1,12)
c98229 format(':diff,            ',12(f12.3,',')	  
c      read(unitNum,98230)title(46),(hu(i),i=1,12)
c98230 format(':hu,              ',12(f12.3,',')	  
c      read(unitNum,98231)title(47),(pres(i),i=1,12)
c98231 format(':pres,            ',12(f12.3,',')	  

!     *************************************************************** 
      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)!# Interception Capacity Table ')

      read(unitNum,*,iostat=ios)junk,(h(1,ii),ii=1,classcount)    !# interception capacity jan mm')	  
      read(unitNum,*,iostat=ios)junk,(h(2,ii),ii=1,classcount)    !# interception capacity feb mm')	  
      read(unitNum,*,iostat=ios)junk,(h(3,ii),ii=1,classcount)    !# interception capacity mar mm')	  
      read(unitNum,*,iostat=ios)junk,(h(4,ii),ii=1,classcount)    !# interception capacity apr mm')	  
      read(unitNum,*,iostat=ios)junk,(h(5,ii),ii=1,classcount)    !# interception capacity may mm')	  
      read(unitNum,*,iostat=ios)junk,(h(6,ii),ii=1,classcount)    !# interception capacity jun mm')	  
      read(unitNum,*,iostat=ios)junk,(h(7,ii),ii=1,classcount)    !# interception capacity jul mm')	  
      read(unitNum,*,iostat=ios)junk,(h(8,ii),ii=1,classcount)    !# interception capacity aug mm')	  
      read(unitNum,*,iostat=ios)junk,(h(9,ii),ii=1,classcount)    !# interception capacity sep mm')	  
      read(unitNum,*,iostat=ios)junk,(h(10,ii),ii=1,classcount)   !# interception capacity oct mm')	  
      read(unitNum,*,iostat=ios)junk,(h(11,ii),ii=1,classcount)   !# interception capacity nov mm')	  
      read(unitNum,*,iostat=ios)junk,(h(12,ii),ii=1,classcount)   !# interception capacity dec mm')	  
      read(unitNum,*,iostat=ios)line    !:EndInterceptionCapacityTable
d     print*,line(1:60)
	if(ios.ne.0)then
        print*,'error in :InterceptionParameterTable'
	  print*
	  pause 'hit enter to abort the program'
	  stop 'Program aborted in read_par @ 509'
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
c        write(51,1003)'hmax ',hmax
c        write(51,1003)'h(  )',(h(i,ii),i=1,12)
c        write(51,1003)'fpetm',(fpetmo(i,ii),i=1,12)
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

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)!:MonthlyEvapotranspirationTable ')
      read(unitNum,*,iostat=ios)junk,(evap(ii,1),ii=1,classcount)     !# monthly evapotranspiration jan mm')	  
      read(unitNum,*,iostat=ios)junk,(evap(ii,2),ii=1,classcount)     !# monthly evapotranspiration feb mm')	  
      read(unitNum,*,iostat=ios)junk,(evap(ii,3),ii=1,classcount)     !# monthly evapotranspiration mar mm')	  
      read(unitNum,*,iostat=ios)junk,(evap(ii,4),ii=1,classcount)     !# monthly evapotranspiration apr mm')	  
      read(unitNum,*,iostat=ios)junk,(evap(ii,5),ii=1,classcount)     !# monthly evapotranspiration may mm')	  
      read(unitNum,*,iostat=ios)junk,(evap(ii,6),ii=1,classcount)     !# monthly evapotranspiration jun mm')	  
      read(unitNum,*,iostat=ios)junk,(evap(ii,7),ii=1,classcount)     !# monthly evapotranspiration jul mm')	  
      read(unitNum,*,iostat=ios)junk,(evap(ii,8),ii=1,classcount)     !# monthly evapotranspiration aug mm')	  
      read(unitNum,*,iostat=ios)junk,(evap(ii,9),ii=1,classcount)     !# monthly evapotranspiration sep mm')	  
      read(unitNum,*,iostat=ios)junk,(evap(ii,10),ii=1,classcount)    !# monthly evapotranspiration oct mm')	  
      read(unitNum,*,iostat=ios)junk,(evap(ii,11),ii=1,classcount)    !# monthly evapotranspiration nov mm')	  
      read(unitNum,*,iostat=ios)junk,(evap(ii,12),ii=1,classcount)    !# monthly evapotranspiration dec mm')	  
      read(unitNum,*,iostat=ios)line     !:EndMonthlyEvapotranspirationTable
d     print*,line(1:60)
	if(ios.ne.0)then
        print*,'error in :MonthlyEvaporatioTabble'
	  print*
	  pause 'hit enter to abort the program'
	  stop 'Program aborted in read_par @ 593'
	endif

      read(unitNum,*,iostat=ios)   !:OptimizationSwitches
      read(unitNum,*,iostat=ios)               !# Optimization: switches,limits & initial values')
	read(unitNum,*,iostat=ios)junk,numa      !# PS optimization 1=yes 0=no')
      read(unitNum,*,iostat=ios)junk,nper      !# opt 1=delta 0=absolute')
      read(unitNum,*,iostat=ios)junk,kc        !# no of times delta halved')
      read(unitNum,*,iostat=ios)junk,maxn      !# max no of trials')
      read(unitNum,*,iostat=ios)junk,dds_flag  !# 0=single run  1=DDS ')
      read(unitNum,*,iostat=ios)junk,errflg    !# 1=wMSE 2=SSE 3=wSSE 4=VOL ')
      read(unitNum,*,iostat=ios)line    !:EndOptimizationSwitches
d     print*,line(1:60)
	if(ios.ne.0)then
        print*,'error in :OptimizationSwitches'
	  print*
	  pause 'hit enter to abort the program'
	  stop 'Program aborted in read_par @ 610'
	endif

!     Hydrological parameters
      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)    !:APILimits
      read(unitNum,*,iostat=ios)junk,a5dlt
      read(unitNum,*,iostat=ios)junk,a5low
      read(unitNum,*,iostat=ios)junk,a5hgh
c      read(unitNum,*,iostat=ios)junk,a5
      read(unitNum,*,iostat=ios)line    !:EndAPILimits
d     print*,line(1:60)

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)    !:HydrologicalParLimits
      read(unitNum,*,iostat=ios)junk,(nclass(i),i=1,classcount)
      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,(akdlt(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(aklow(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(akhgh(i),i=1,classcount)
c      read(unitNum,*,iostat=ios)junk,(ak(i),i=1,classcount)

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,(akfsdlt(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(akfslow(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(akfshgh(i),i=1,classcount)
c      read(unitNum,*,iostat=ios)junk,(akfs(i),i=1,classcount)

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,(recdlt(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(reclow(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(rechgh(i),i=1,classcount)
c      read(unitNum,*,iostat=ios)junk,(rec(i),i=1,classcount)

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,(r3dlt(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(r3low(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(r3hgh(i),i=1,classcount)
c      read(unitNum,*,iostat=ios)junk,(r3(i),i=1,classcount)

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,(fpetdlt(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(fpetlow(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(fpethgh(i),i=1,classcount)
c      read(unitNum,*,iostat=ios)junk,(fpet(i),i=1,classcount)

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,(ftalldlt(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(ftalllow(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(ftallhgh(i),i=1,classcount)
c      read(unitNum,*,iostat=ios)junk,(ftall(i),i=1,classcount)

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,(retndlt(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(retnlow(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(retnhgh(i),i=1,classcount)
c      read(unitNum,*,iostat=ios)junk,(retn(i),i=1,classcount)

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,(ak2dlt(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(ak2low(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(ak2hgh(i),i=1,classcount)
c      read(unitNum,*,iostat=ios)junk,(ak2(i),i=1,classcount)

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,(ak2fsdlt(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(ak2fslow(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(ak2fshgh(i),i=1,classcount)
c      read(unitNum,*,iostat=ios)junk,(ak2fs(i),i=1,classcount)
      read(unitNum,*,iostat=ios)line    !:EndHydrologicalParLimits
d     print*,line(1:60)
	if(ios.ne.0)then
        print*,'error in HydrologicalParLimits'
	  print*
	  pause 'hit enter to abort the program'
	  stop 'Program aborted in read_par @ 685'
	endif

!     snow parameters
      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)     !:GlobalSnowParLimits
      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,fmadjustdlt
      read(unitNum,*,iostat=ios)junk,fmadjustlow
      read(unitNum,*,iostat=ios)junk,fmadjusthgh
c      read(unitNum,*,iostat=ios)junk,fmadjust

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,fmalowdlt
      read(unitNum,*,iostat=ios)junk,fmalowlow
      read(unitNum,*,iostat=ios)junk,fmalowhgh
c      read(unitNum,*,iostat=ios)junk,fmalow

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,fmahighdlt
      read(unitNum,*,iostat=ios)junk,fmahighlow
      read(unitNum,*,iostat=ios)junk,fmahighhgh
c      read(unitNum,*,iostat=ios)junk,fmahigh

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,gladjustdlt
      read(unitNum,*,iostat=ios)junk,gladjustlow
      read(unitNum,*,iostat=ios)junk,gladjusthgh
c      read(unitNum,*,iostat=ios)junk,gladjust
      read(unitNum,*,iostat=ios)line    !:EndGlobalSnowParLimits
d     print*,line(1:60)
	if(ios.ne.0)then
        print*,'error in GlobalSnowParLinmits'
	  print*
	  pause 'hit enter to abort the program'
	  stop 'Program aborted in read_par @ 720'
	endif

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)    !:SnowParLimits
      read(unitNum,*,iostat=ios)junk,(nclass(i),i=1,classcount)
      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,(fmdlt(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(fmlow(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(fmhgh(i),i=1,classcount)
c      read(unitNum,*,iostat=ios)junk,(fm(i),i=1,classcount)

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,(basdlt(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(baslow(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(bashgh(i),i=1,classcount)
c      read(unitNum,*,iostat=ios)junk,(base(i),i=1,classcount)

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,(subdlt(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(sublow(i),i=1,classcount)
      read(unitNum,*,iostat=ios)junk,(subhgh(i),i=1,classcount)
c      read(unitNum,*,iostat=ios)junk,(sublim_factor(i),i=1,classcount)
      read(unitNum,*,iostat=ios)line   !:EndSnowParLimits
d     print*,line(1:60)
	if(ios.ne.0)then
        print*,'error in SnowParLimits'
	  print*
	  pause 'hit enter to abort the program'
	  stop 'Program aborted in read_par @ 749'
	endif

      if(sublimflg.eq.2)then          ! for sublim_factor only
        do ii=1,classcount
          sublow(ii)=sublow(ii)/24
          subhgh(ii)=subhgh(ii)/24
        end do
      endif

!     Routing parameters
      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)    !:RoutingParLimits
      read(unitNum,*,iostat=ios)junk,(rivtype(i),i=1,nrvr)
      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,(flzdlt(i),i=1,nrvr)
      read(unitNum,*,iostat=ios)junk,(flzlow(i),i=1,nrvr)
      read(unitNum,*,iostat=ios)junk,(flzhgh(i),i=1,nrvr)
c      read(unitNum,*,iostat=ios)junk,(flz_o(i),i=1,nrvr)

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,(pwrdlt(i),i=1,nrvr)
      read(unitNum,*,iostat=ios)junk,(pwrlow(i),i=1,nrvr)
      read(unitNum,*,iostat=ios)junk,(pwrhgh(i),i=1,nrvr)
c      read(unitNum,*,iostat=ios)junk,(pwr_o(i),i=1,nrvr)

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,(r2ndlt(i),i=1,nrvr)
      read(unitNum,*,iostat=ios)junk,(r2nlow(i),i=1,nrvr)
      read(unitNum,*,iostat=ios)junk,(r2nhgh(i),i=1,nrvr)
c      read(unitNum,*,iostat=ios)junk,(r2n_o(i),i=1,nrvr)

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,(thetadlt(i),i=1,nrvr)
      read(unitNum,*,iostat=ios)junk,(thetalow(i),i=1,nrvr)
      read(unitNum,*,iostat=ios)junk,(thetahgh(i),i=1,nrvr)
c      read(unitNum,*,iostat=ios)junk,(theta_o(i),i=1,nrvr)

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,(kconddlt(i),i=1,nrvr)
      read(unitNum,*,iostat=ios)junk,(kcondlow(i),i=1,nrvr)
      read(unitNum,*,iostat=ios)junk,(kcondhgh(i),i=1,nrvr)
c      read(unitNum,*,iostat=ios)junk,(kcond_o(i),i=1,nrvr)

      read(unitNum,*,iostat=ios)
      read(unitNum,*,iostat=ios)junk,(rlakedlt(i),i=1,nrvr)
      read(unitNum,*,iostat=ios)junk,(rlakelow(i),i=1,nrvr)
      read(unitNum,*,iostat=ios)junk,(rlakehgh(i),i=1,nrvr)
c      read(unitNum,*,iostat=ios)junk,(rlake_o(i),i=1,nrvr)
      read(unitNum,*,iostat=ios)line    !:EndRoutingParLimits
d     print*,line(1:60)
	  if(ios.ne.0)then
          print*,'error in RoutingParLimits '
	    print*
	    pause 'hit enter to abort the program'
	    stop 'Program aborted in read_par @ 797'
	  endif

!     rev. 9.8.01  Jul.  21/11  - NK: added ragmet optimization to dds setup
!     Global parameters
      read(unitNum,*,iostat=ios)junk
      read(unitNum,*,iostat=ios)junk    !:GlobalParLimits
      if(ios.ne.0)then
        print*,'end of file reached after routing parameters'
        print*,'Global parameter limits not found'
      endif
      if(junk(1:16).eq.':GlobalParLimits')then
        read(unitNum,*,iostat=ios)
        read(unitNum,*,iostat=ios)junk,rlapsedlt
        read(unitNum,*,iostat=ios)junk,rlapselow
        read(unitNum,*,iostat=ios)junk,rlapsehgh

        read(unitNum,*,iostat=ios)
        read(unitNum,*,iostat=ios)junk,tlapsedlt
        read(unitNum,*,iostat=ios)junk,tlapselow
        read(unitNum,*,iostat=ios)junk,tlapsehgh

        read(unitNum,*,iostat=ios)
        read(unitNum,*,iostat=ios)junk,radinfldlt
        read(unitNum,*,iostat=ios)junk,radinfllow
        read(unitNum,*,iostat=ios)junk,radinflhgh

        read(unitNum,*,iostat=ios)
        read(unitNum,*,iostat=ios)junk,smoothdistdlt
        read(unitNum,*,iostat=ios)junk,smoothdistlow
        read(unitNum,*,iostat=ios)junk,smoothdisthgh
        read(unitNum,*,iostat=ios)line    !:EndGlobalSnowParLimits
	  if(ios.ne.0)then
          print*,'error in GlobalParLimits '
	    print*
	    pause 'hit enter to abort the program'
	    stop 'Program aborted in read_par @ 833'
	  endif

d       print*,line(1:60)
      else
        rlapsedlt=-1.0
        tlapsedlt=-1.0
        radinfldlt=-1.0
        smoothdistdlt=-1.0
        print*,'WARNING:'
        print*,'rlapsedlt,tlapsedlt,radinfldlt & smoothdist'
        print*,'set - -1.0  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
        print*
      endif
       

! Assignments & checking
! Assignments & checking
! Assignments & checking

      if(unitNum.eq.51)then     ! write a few blank lines for readability
	  write(unitNum,*)  
	  write(unitNum,*)
        return
	endif
      
      close(unitNum,status='keep')

	if(iopt.gt.0)print*,'Closed unit=',unitNum
	if(iopt.gt.0)print*,fln(flnNum)(1:50),'read'
	
	
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
        end do

      if(unitNum.eq.51)then
         write(unitNum,*)
	   write(unitNum,*)line
         write(unitNum,*)'diff  ',(diff(i),i=1,12)
         write(unitNum,*)'humid ',(hu(i),i=1,12)
         write(unitNum,*)'pres  ',(pres(i),i=1,12)
         write(unitNum,*)
      endif

      if(dds_flag.eq.1)iopt=0
      
      if(a7.lt.0.50.or.a7.gt.0.99)then
        print*,' Value for weighting factor - old vs. new sca value'
        print*,' a7 not between 0.5-0.99in the par file'
        print*,' value found is',a7
        print*,' a7=0.50 assumed - matches spl8 results'
        write(51,9811)
        write(51,9812)
        write(51,9813)
        write(51,9814)
        write(51,9815)
        write(51,9811)
9811    format(' ')
9812    format(' WARNING - in rdpar')
9813    format('Value for weighting factor - old vs. new sca value')
9814    format('a7 not between 0.5-0.99in the par file')
9815    format('a7=0.50 assumed - matches spl8 results')
        a7=0.5
        tempflg='y'
      endif
!     rev. 9.6.01  Mar.  01/10  - NK: DDS capability added
      if(a3.lt.0.0.and.dds_flag.eq.1)then
        print*,'a3 < 0   Must be = or > 0 to use DDS'
        print*,' value found is',a3
        print*,'Please correct the par file'
        a3=0.05
        print*,'a3=0.05 assumed'
        print*
      endif
      if(a4.lt.0.0.and.dds_flag.eq.1)then
        print*,'a4 < 0   Must be = or > 0 to use DDS'
        print*,' value found is',a4
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
        pause ' Program paused in rdpar @ 480'
      endif
       
      tempflg='n'
      if(a9.gt.5.0.or.a9.lt.0.3)then
!       change range from 0.33-2.0 to 0.33-5.0 Oct. 20/03 NK
        print*
        print*,' a9 (heat deficit / swe ratio outside range 0.333 - 5.0'
        print*,' value found is',a9
        print*,' a9 = .333 assumed'
        print*,' If this default value is unsuitable, please'
        print*,' edit the bsnm.par file and enter proper value'
        print*,' between 0.333 and 5.0'
        print*
        a9=0.333
        write(51,5110)a9
        tempflg='y'
      endif

      if(iopt.eq.2)print*, 'in rdpar @ 637'

      if(a10.gt.3.0.or.a10.lt.0.5)then
        print*
        print*,' a10 uz discharge exponent outside range 0.5 - 3.0'
        print*,' value found is',a10
        print*,' a10 = 1.0 assumed'
        print*,' If this default value is unsuitable, please'
        print*,' edit the bsnm.par file and enter proper value'
        print*,' between 0.5 and 3.0'
        print*
        a10=1.0
        write(51,5111)a10
        tempflg='y'
      endif
      if(tempflg.eq.'y')then
        print*
        print*,' To continue with the default value(s) hit enter'
        print*,' To prevent this message from reoccurring, please enter'
        print*,' values for a9 and/or a10 in the parameter .par file'
        print*
        pause '  Program paused in rdpar @ 512'
      endif
5110  format(//,' Value for a9 heat deficit to swe outside range',/
     *     f10.3,' assumed')
5111  format(//,' Value for a10 uz exponent outside range',/
     *     f10.3,' assumed')

      if(iopt.eq.2)print*, 'in rdpar @ 660'
      
!       ERROR CHECK FOR REVISED LZS CLASSES
        do i=1,nrvr
!       rev. 9.1.37  Mar.  22/03  - Option to turn off leakage by setting LZF < 0.0
          if(flz_o(i).lt.0.0)flz_o(i)=0.1e-30
          if(pwr_o(i).lt.0.1.or.pwr_o(i).gt.4.0)then   !  Feb. 22/06  nk
            print*
            print*,'The value of pwr for land cover class ',nrvr
            print*,'is ',pwr(nrvr),' which' 
            print*,'is not in the range of 0.1 - 4.0'
            print*,'Please check all values'
            print*
            if(dds_flag.eq.0.and.numa.eq.0)then
              pause 'Hit enter to continue with this value (@ 607)'
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
              stop 'Program aborted in rdpar @ 715'
            endif
          end do
      
      if(wetflg.eq.'y')then
!       parameter checking
        do i=1,nrvr
          if(abs(theta_o(i)).lt.0.00001)then
            print*,' Value for theta for river class ',i
            print*,' too close to zero. Please provide reasonable'
            print*,' value  (0.01 to 0.5)'
            tempflg='y'
          endif
          if(abs(widep_o(i)).lt.1.0)then
            print*,' Value for width/depth ratio for river class',i
            print*,' too close to zero. Please provide reasonable'
            print*,' value say (1 to 100)'
            tempflg='y'
          endif
          if(abs(kcond_o(i)).lt.0.99e-09)then
            print*,' Value for kcond for river class ',i
            print*,' too close to zero. Please provide reasonable'
            print*,' value  (1.0e-09 to 1.00)'
            tempflg='y'
          endif
        end do
        if(tempflg.eq.'y')stop 'Program aborted in rdpar @ 595'
      endif       ! wetflg = y

!     rev. 9.2.35  Mar.  22/06  - NK: Glacier flow bypasses wetlands
!     glacier_class_number will be used in tracer
      glacier_class_number=0  !default value

      do i=1,classcount      ! added 09/03/05  nk
        if(nclass(i).eq.'GLACIER   ')nclass(i)='glacier   '
        if(nclass(i).eq.'Glacier   ')nclass(i)='glacier   '
        if(nclass(i).eq.' glacier  ')nclass(i)='glacier   '
        if(nclass(i).eq.' GLACIER  ')nclass(i)='glacier   '
        if(nclass(i).eq.' Glacier  ')nclass(i)='glacier   '
        if(nclass(i).eq.'glacier   ')glacier_class_number=i
        if(nclass(i).eq.'WETLAND   ')nclass(i)='wetland   '
        if(nclass(i).eq.'Wetland   ')nclass(i)='wetland   '
        if(nclass(i).eq.' wetland  ')nclass(i)='wetland   '
        if(nclass(i).eq.' WETLAND  ')nclass(i)='wetland   '
        if(nclass(i).eq.' Wetland  ')nclass(i)='wetland   '
        if(nclass(i).eq.'WATER     ')nclass(i)='water     '
        if(nclass(i).eq.'Water     ')nclass(i)='water     '
        if(nclass(i).eq.' water    ')nclass(i)='water     '
        if(nclass(i).eq.' WATER    ')nclass(i)='water     '
        if(nclass(i).eq.' Water    ')nclass(i)='water     '
      end do

!     check added Jul. 31/06 nk
      if(glacier_class_number.eq.0.and.itrace.eq.1.
     *                           and.trcflg.eq.'y')then
          trcflg='n'
          print*,'Tracer is set to track glacier melt'
          print*,'and the tracer flag (trcflg) is set to "y"'
          print*,'but there is no glacier class'
          print*,'trcflg is set to "n" '
          print*
          pause 'Continue by hitting enter'
      endif
                  
!     rev. 9.08.01 Apr.   3/01  - check wetland designation in rdpar
      if(wetflg.eq.'y')then
      if(nclass(classcount-2).eq.'wetland   ')then
!       everything just fine
!         only class classcount-2 is acceptable as a coupled wetland
      else
        print*,' Class no. ',classcount-2,' is expected to be'
        print*,' the wetland class and be labelled as such'
        print*,' Name found ***',nclass(classcount-2),'***'
        print*,' It should be ***wetland   ***'

        print*,' Please be sure that the order of the parameters'
        print*,' and the land cover data are in the same order'
        print*,' Please fix the parameter file and try again'
        print*
        stop ' Program aborted in rdpar at line ~217'
      endif
      endif
         
      do i=1,classcount
c       if(rec(i).gt.1.00.or.rec(i).lt.0.0)then
        if(rec(i).lt.0.0)then
!         note: larger value will give NaN for x4 in runof6 
          print*,'the value for rec can not be less than 0.00'
          print*,'class=',i,' value=',rec(i)
          print*,'Please change the value in the par file'
          print*
          stop 'Program aborted in rdpar @ 1040'
        endif
      end do
      ijunko=1

      if(ak(classcount-1).ge.0.0.or.akfs(classcount-1).ge.0.0)then
        print*,'Water class not specified'
        print*,'ak and akfs NOT given a -ve value'
        print*,'ak   value found =',ak(classcount-1)
        print*,'akfs value found =',akfs(classcount-1)
        print*,'Please correct the par file'
        print*,'No of classed expected =',classcount-1
        print*
      endif
      
c!   rev. 9.5.65  Sep.  26/09  - NK: lapse rate changed from dC per 100 m to dC per m
c      if(rlapse.gt.0.1)then
c        print*,'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
c        print*,'WARNING:'
c        print*, 'Value for RLAPSE is too large'
c        print*, 'RLAPSE used to be mm C per 100m'
c        print*, 'this has been changed to mm per km'
c        print*, 'Please correct the par file'
c        print*,'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
c        pause 'Program pause in read_par @ 1143'
c      endif

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
c         stop 'Program aborted in rdpar @ 1356'
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
          write(51,*)'WARNING:'
          write(51,*)'Temp3 <0.0001  set temp3=0.001'
          write(51,*)
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
        write(51,1003)'hmax ',hmax
        write(51,1003)'h(  )',(h(i,ii),i=1,12)
        write(51,1003)'fpetm',(fpetmo(i,ii),i=1,12)
 1003   format(a5,99f5.2)
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
        if(flz(i).le.0.0)then
          write(*,8429)flz(i)
 8429     format('In read_par,flz(',i2,')=',f12.7,
     * ' wrong parameter value - change flz(i) to +ve value')
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

!     This was fixed between version 9.1.11 and 9.1.13
!     It used to be all sin(...). End result: pet too high
      sinlat=sin(3.1416*lat/180.0)
      coslat=cos(3.1416*lat/180.0)
      tanlat=tan(3.1416*lat/180.0)
     
      
      
      
      
      
      
      
      
      if(iopt.eq.2)print*, 'read_par before return @ 9989'

 9989 RETURN

      END SUBROUTINE read_par

