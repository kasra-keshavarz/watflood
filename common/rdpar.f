      SUBROUTINE rdpar(iprtflg,iverflg,e1)

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
!  version numbers are added for version 7.0 and later. this allows
!  parameter files to backward compatible when more parameters are 
!  added to the file in the future. in version 7.0 some headers are 
!  added as well for readability.
!
!   REV. 7.5 seperate snow covered and bare ground - Jul/95
!   REV. 7.75 - May.  27/96 -  added ak2fs in rdpar & runof5
!
!   REV. 8.32 - June  13/97 -  bypassed non-flagged parameters in OPT
!   REV. 8.41 - July  21/97 -  added tipm to the optimization table
!   REV. 8.62 - Dec.  30/97 -  fixed param s/r comb'd et & par flgs
!   REV. 8.74 - Mar.  31/98 -  reinvented fs stuff in opt
!
!   REV. 9.00 - March 2000  -  TS: CONVERSION TO FORTRAN 90
!   REV. 9.03 - Nov   2000  -  TS: ADDED WATFLOOD SWAMP ROUTING 
!   rev. 9.07    Mar.  14/01  - fixed use of opt par's  for numa=0  
!   rev. 9.08.01 Apr.   3/01  - check wetland designation in rdpar
!   rev. 0.1.04  Oct.   4/01  - added A7 for weighting old/new sca in melt
!   rev. 0.1.05  Oct.   4/01  - new format parameter file
!   rev  9.1.20  Jun.  25/02  - Added A10 as the power on the UZ discharge function
!   rev  9.1.25  Sep.  11/02  - Added A11 as bare ground equiv. vegn height  
!   rev. 9.1.37  Mar.  22/03  - Option to turn off leakage by setting LZF < 0.0
!   rev. 9.1.40  Apr.  24/03  - Min time step A6 read in strfw over rides the A6 from the par file
!   rev. 9.1.76  Mar.  09/05  - NK: separated glacier parameters in par file
!   rev. 9.2.09  Sep.  11/05  - NK: removed write_par.for from rdpar.for
!   rev. 9.2.11  Sep.  11/05  - NK: added Manning's n  r1n & r2n
!   rev. 9.2.16  Oct.  10/05  - NK: Fixed bug for widep in rdpar
!   rev. 9.2.35  Mar.  22/06  - NK: Glacier flow bypasses wetlands
!   rev. 9.2.37  Mar.  31/06  - NK: Removed impervious area as special class
!   rev. 9.4.07  May.  15/07  - NK: converted opt to gridded routing parameters
!   rev. 9.5.06  Feb.  05/08  - NK: added pool and pool_o in rdpar & route
!   rev. 9.5.27  Apr.  15/08  - NK: fixed allocation for chnl in rdpar
!   rev. 9.5.43  Oct.  27/08  - NK: changed bottom part of par file to be free format
!   rev. 9.5.45  Dec.  16/08  - NK: added various error calculations - user's choice with errflg
!   rev. 9.5.72  Oct.  12/09  - NK: fixed bug in rdpar setting init values for fpet & ftal
!   rev. 9.6.01  Mar.  01/10  - NK: rlake parameter added for Manning n correction
!   rev. 9.6.02  Mar.  15/10  - NK: add sublimation to optimization
!   rev. 9.6.06  Apr.  18/10  - NK: added glacier adjust for optimization
!   rev. 9.7.19  Jan.  18/11  - NK: Added sensitivity analysis
!   rev. 9.7.20  Jan.  31/11  - NK: Moved open statement for rdpar to rdpar/f
!   rev.  check
!
!
!  iprtflg - if eq 2 write a new parameter file after optimization
!  r1      - the roughness of the flood plain
!  r2      - the roughness of the largest part of the channel
!  r1      - factor for raising r2 ie if r1=2 then f.p. roughness
!            is 2 times channel roughness
!  zz      - is an exponent in calculating the mannings  n
!  h       - crop height
!  dds_flag- flag to run the DDS optimization =0 no dds opt  =1 dds opt
!  itrace  - tracer flag 100=GW, 4=3-comp (SW, IF, GW)
!             or 5=6-comp (3-comp's + melt fractions).
!  e1       - void ratio
!  ak      - permeability
!  sm      -soil moisture(average for month)
!  nbsn    - no of basins with different r2 to be optimized
!          - must be smaller than nrvr
!  classcount_local  - no of classes to be optimized max=4
!  errflg  - melt debug function class selection number 
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

      real*4, dimension(:),   allocatable :: ajunk
      real*4, dimension(:),   allocatable :: par_temp


      CHARACTER(79) :: comment
      CHARACTER(10) :: time
      CHARACTER(8)  :: cday
      CHARACTER(128):: qstr
      CHARACTER(60) :: junk
      CHARACTER(30) :: filename1
      CHARACTER(1)  :: errorflg,answer,linetype,tempflg,firstpass
      LOGICAL       :: flzflg,pwrflg,R1nflg,R2nflg,mndrflg,
     *                 aa2flg,aa3flg,aa4flg,widepflg
      INTEGER(kind=2) :: result1
      INTEGER       :: nrvr1,ntest,nchr,iprtflg,linecount,i,
     *                 ios,iverflg,ix,iallocate,ijunko,ii,j,n
      real*4        :: hmax,e1,kcondflg,xxx1
      integer       :: classcount_local

      DATA ntest/-10906/qstr/'param'/nchr/5/
      data firstpass/'y'/

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      if(firstpass.eq.'y')then

!     rev. 9.7.19  Jan.  18/11  - NK: Added sensitivity analysis
      i=max(nrvr,classcount)
      if(.NOT.allocated(ajunk))then
        allocate(ajunk(i),stat=iAllocate)
        if(iAllocate.ne.0)STOP 'Error with allocation of ajunk in rdpar'
	endif

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
     *  ds(classcount),dsfs(classcount),chnl(5), ! note: chnl = always 5
     *  r3(classcount),r4(classcount),r3fs(classcount),rec(classcount),
     *  ak(classcount),akfs(classcount),
     *  r3low(classcount),r3fslow(classcount),reclow(classcount),
     *  aklow(classcount),akfslow(classcount),ak2fslow(classcount),
     *  r3hgh(classcount),r3fshgh(classcount),rechgh(classcount),
     *  akhgh(classcount),
     *  akfshgh(classcount),ak2fshgh(classcount),r3dlt(classcount),
     *  r3fsdlt(classcount),recdlt(classcount),akdlt(classcount),
     *  akfsdlt(classcount),ak2fsdlt(classcount),retn(classcount),
     *  ak2(classcount),retnlow(classcount),ak2low(classcount),
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
     *  fmndlt(classcount),fmnlow(classcount),fmnhgh(classcount),
     *  uajdlt(classcount),uajlow(classcount),uajhgh(classcount),
     *  mbsdlt(classcount),mbslow(classcount),mbshgh(classcount),
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
     *   'Warning: error with allocation of areamelta arrays in spl9'
	endif
!     TS - ALLOCATIONS OF AREAETA ARRAYS (PARTIAL)
!     RAD ALLOCATED IN SHEDA.FOR
!     TS - ADDED ALLOCATIONS FOR EVAP-SEPARATION PARAMS (22/03/06)
!     TS: CHANGED ALLOCATIONS OF alb,pet, evap TO classcount (27/03/06)
!     rev. 9.1.80  Mar.  31/05  - NK: added sublimation   (sublim)
!     moved here from rdshed 27/07/06 because needed for bsn.for  nk
      if(.NOT.allocated(evap))then
        allocate(evap(classcount,12),flint(classcount),
     *  sublim_factor(classcount),sublim_rate(classcount),
     *  diff(12),hu(12),ffcap(classcount),fcap(classcount),
     *  spore(classcount),alb(classcount),pres(12),stat=iAllocate)
 !      *rh(na),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *   'Warning: error with allocation of areaeta arrays in spl9'
	endif

c      print*,'Allocations done in rdpar',classcount,nrvr

      endif  !firtspass



!     rev. 9.7.20  Jan.  31/11  - NK: Moved open statement for rdpar to rdpar/f
!     basin/bsnm.par
      open(unit=32 ,file=fln(2) ,status='old',iostat=ios)
      if(ios.ne.0)then
        print*,'Problems on unit 32'
        write(*,99172)fln(2)
        write(98,99172)fln(2)
99172   format(' Warning: Error opening or reading fln:',a30/
     *  ' Probable cause: missing basin\bsnm.par input file'/
     *  ' OR: use backslash \ only in event files')
        print*,'iostat code =',ios
        STOP 'program aborted in spl.for @ 372'
      endif

!     PASS OVER THE COMMENT LINES ON THE PARAMETER FILE:
!     SEE SAME CODE ON PARAM.FOR
      linecount=0
      do i=1,20
        read(32,3010,iostat=ios)linetype,comment
        if(ios.ne.0)then
        print*,'Problems on unit 32'
          print*,' error reading the parameter file comment lines'
          print*
          stop 'Program aborted in SPL @ 372.1'
        endif
        if(linetype.eq.'#')linecount=linecount+1
      end do

d     if(iopt.eq.2)print*, 'passed 747 in rdpar'

      rewind 32

      do i=1,linecount
        read(32,3010)
      end do

d     if(iopt.eq.2)print*, 'passed 756 in rdpar'

!     ADDED THIS FROM TODD'S ETPAR.FOR - FRANK S: NOV/97 
      alamb=2.478
      den=999.
      alpha=1.35
      dds_flag=0 ! default for old par files

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


        title( 7) ='dsfs '
        title(10) ='akfs '
        title(18) ='r3fs '
        title(88) ='a10   '
        title(99) ='a9   '
        title(100)='a11  '

!     REWIND B/C SOME DATA WAS READ EARLIER
      rewind 32

!     PASS OVER THE COMMENT LINES ON THE PARAMETER FILE:
!     SEE SAME CODE IN SPL9.FOR
      linecount=0
      do i=1,20
        read(32,3010,iostat=ios)linetype,notes(i)
        if(ios.ne.0)then
          print*,' Problem reading comment lines'
        endif
        if(linetype.eq.'#')linecount=linecount+1
      end do

      rewind 32
      do i=1,linecount
        read(32,3010,iostat=ios)
        if(ios.ne.0)then
          print*,' Problem reading comment lines'
        endif
      end do

!     if we can read the next line, we have the old format  
      if(linecount.eq.0)rewind 32
      read(32,3000,iostat=ios)title(1),iopt,itype,numa,nper,kc,maxn,ver
      write(51,3000)title(1),iopt,itype,numa,nper,kc,maxn,ver

      if(iopt_start.eq.99)then
        iopt=99
        print*,'iopt temporarily changed to 99 via argument'
        print*
      endif

      if(ios.eq.0)then
        iverflg=0
        if(ver.lt.9.4)iverflg=1  ! a new par file will be written
        read(32,1000,iostat=ios)title(2),e1,ix,errflg
        if(ios.ne.0)then
          print*,' Old format par file found but'
          print*,' problems reading data line 2'
          print*
          stop ' Program aborted at 151 in rdpar'
        endif
        if(errflg.le.0.or.errflg.gt.2)then 
          write(51,5007)errflg
          errflg=0
        endif
        read(32,5000,iostat=ios)title(3),classcount_local,nbsn
        if(ios.ne.0)then
          print*,' Old format par file found but'
          print*,' problems reading data line 3'
          print*
          stop ' Program aborted at 152 in rdpar'
        endif
        if(nbsn.gt.nrvr)then
          print*,'The number of river par sets to be optimized is '
          print*,'larger than the number of rivers.'
          print*,'No of river classes from the .shd file = ',nrvr
          print*,'No of river classes from the par file  = ',nbsn
          print*,'Please change the par file'
          print*,'Program will proceed with nbsn set to nrvr'
          nbsn=nrvr
          print*
          pause 'Hit enter to continue  --  hit Ctrl C to abort'
        endif

        title(71)='ver'
        title(72)='iopt'
        title(73)='itype'
        title(74)='numa'
        title(75)='nper'
        title(76)='kc'
        title(77)='maxn'
        title(78)='ddsfl'
        title(79)='itrce'
        title(80)='errfl'
        title(81)='typeo'
        title(82)='nbsn'
        title(83)='mndr'
        title(84)='aa2'
        title(85)='aa3'
        title(86)='aa4'
        title(87)='theta'
        title(88)='widep'
        title(89)='kcond'
        title(90)='pool'
        title(70)='rlake'

        title(91)='a1'
        title(92)='a2'
        title(93)='a3'
        title(94)='a4'
        title(95)='a5'
        title(96)='a6'
        title(97)='a7'
        title(98)='a8'
        title(99)='a9'
        title(67)='a10'
        title(68)='a11'
        title(69)='a12'
        notes(71)='parameter file version number'
        notes(72)='debug level'
        notes(73)='flood plain type 0-yes  1=no'
        notes(74)='optimization 0=no 1=yes'
        notes(75)='opt delta 0=absolute 1=ratio'
        notes(76)='no of times delta halved'
        notes(77)='max no of trials'
        notes(78)='dds_flag 0=no  1=yes'
        notes(79)='trace flag - picks tracer'
        notes(80)='0=rms 1=r**2 2=dv 3=dv/Qmean'
        notes(81)='no of land classes optimized(part 2)'
        notes(82)='no of river classes optimized (part 2)'
        notes(91)='ice cover weighting factor'
        notes(92)='Manning`s n correction for instream lakes'
        notes(93)='N-S Efficiency penalty coefficient'
        notes(94)='N-S Efficiency penalty threshold'
        notes(95)='API coefficient'
        notes(96)='Minimum routing time step in seconds'
        notes(97)='weighting factor - old vs. new sca value'
        notes(98)='min temperature time offset'
        notes(99)='max heat deficit /swe ratio'
        notes(67)='exponent on uz discharce function'
        notes(68)='bare ground equivalent veg height for evap`n'
        notes(69)='min precip rate for smearing'
      else

!       new format file found

!!!!!!!!NOTE: void ratio e1 is dropped from the par file.
!             needs to be changed to use spore(ii) etv.
!             for now, use por = average(spore(ii))
!             done in runof6.for, see por=.......
        e1=-1.0     ! flag to use average in runof6.for

        backspace 32
        read(32,9807,iostat=ios)title(71),ver,notes(71)
        if(iopt.ge.1)print*,title(71),ver,notes(71)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading ver in the parameter file data'
          print*
          stop '@ ~227a in rdpar'
        endif

        if(ver.lt.9.4)iverflg=1

        read(32,9805,iostat=ios)title(72),iopt,notes(72)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading iopt in the parameter file data'
          print*
          stop '@ ~227b in rdpar'
        endif

        if(iopt_start.eq.99)then
          iopt=99
          print*,'iopt temporarily changed to 99 via argument'
          print*
        endif

        read(32,9805,iostat=ios)title(73),itype,notes(73)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading itype in the parameter file data'
          print*
          stop '@ ~227c in rdpar'
        endif
        read(32,9805,iostat=ios)title(74),numa,notes(74)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading numa in the parameter file data'
          print*
          stop '@ ~227d in rdpar'
        endif
        read(32,9805,iostat=ios)title(75),nper,notes(75)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading nper in  the parameter file data'
          print*
          stop '@ ~227e in rdpar'
        endif
        read(32,9805,iostat=ios)title(76),kc,notes(76)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading kc in  the parameter file data'
          print*
          stop '@ ~227f in rdpar'
        endif
        read(32,9805,iostat=ios)title(77),maxn,notes(77)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading maxn in  the parameter file data'
          print*
          stop '@ ~227g in rdpar'
        endif
        read(32,9805,iostat=ios)title(78),dds_flag,notes(78)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading errflg in the parameter file data'
          print*
          stop '@ ~227h in rdpar'
        endif
        read(32,9805,iostat=ios)title(79),itrace,notes(79)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading itrace in the parameter file data'
          print*
          stop '@ ~227i in rdpar'
        endif
        read(32,9805,iostat=ios)title(80),errflg,notes(80)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading errflg in the parameter file data'
          print*
          stop '@ ~227j in rdpar'
        endif
        if(dds_flag.eq.1)print*,'errfg=',errflg
        read(32,9805,iostat=ios)title(81),classcount_local,notes(81)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading classcount_local in the par file data'
          print*
          stop '@ ~227k in rdpar'
        endif
        read(32,9805,iostat=ios)title(82),nbsn,notes(82)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading nbsn in  the parameter file data'
          print*
          stop '@ ~227l in rdpar'
        endif
        if(nbsn.gt.nrvr)then
          print*,'The number of river par sets to be optimized is '
          print*,'larger than the number of rivers.'
          print*,'No of river classes from the .shd file = ',nrvr
          print*,'No of river classes from the par file  = ',nbsn
          print*,'Please change the par file'
          print*,'Program will proceed with nbsn set to nrvr'
          nbsn=nrvr
          print*
          pause 'Hit enter to continue  --  hit Ctrl C to abort'
        endif
      endif

      if(dds_flag.eq.1)iopt=0
     
      type1=float(itype)

!     Note: r2n,theta,widep,kcond,flz & pwr are allocated in read_shed_ef
!           needed to write a new par file if version is old.

c      if(numa.gt.0.or.par_file_version.lt.par_file_version_latest
c     *              .or.dds_flag.eq.1)then


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

        if(iopt.eq.2)print*, 'param 2'

c       if(ver.lt.9.10)then  !<<<<<<check version number in previous file
        if(ver.lt.9.1)then  !<<<<<<check version number in previous file

        read(32,1001,iostat=ios)title(4),a1,a2,a3,a4,a5  
        read(32,1001,iostat=ios)title(5),a6,a7,a8,a9,a10,a11,a12
        if(ios.ne.0)then
          print*,' Old format par file found but'
          print*,' problems reading data line 3 or 4'
          print*
          stop ' Program aborted at 152 in rdpar'
        endif
      else
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading mndr in  the parameter file data'
          print*
          stop '@ ~227l in rdpar'
        endif
        read(32,9806,iostat=ios)title(91),a1,notes(91) 
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading a1 in  the parameter file data'
          print*
          stop '@ ~227l in rdpar'
        endif
        read(32,9806,iostat=ios)title(92),a2,notes(92)  
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading a2 in  the parameter file data'
          print*
          stop '@ ~227l in rdpar'
        endif
        read(32,9806,iostat=ios)title(93),a3,notes(93)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading a3 in  the parameter file data'
          print*
          stop '@ ~227l in rdpar'
        endif
        read(32,9806,iostat=ios)title(94),a4,notes(94)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading a4 in  the parameter file data'
          print*
          stop '@ ~227l in rdpar'
        endif
        read(32,9806,iostat=ios)title(95),a5,notes(95)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading a5 in  the parameter file data'
          print*
          stop '@ ~227l in rdpar'
        endif
        read(32,9806,iostat=ios)title(96),a6,notes(96)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading a6 in  the parameter file data'
          print*
          stop '@ ~227l in rdpar'
        endif
        read(32,9806,iostat=ios)title(97),a7,notes(97)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading a7 in  the parameter file data'
          print*
          stop '@ ~227l in rdpar'
        endif
        read(32,9806,iostat=ios)title(98),a8,notes(98)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading a8 in  the parameter file data'
          print*
          stop '@ ~227l in rdpar'
        endif
        read(32,9806,iostat=ios)title(99),a9,notes(99)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading a9 in  the parameter file data'
          print*
          stop '@ ~227l in rdpar'
        endif
        read(32,9806,iostat=ios)title(67),a10,notes(67)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading a10 in  the parameter file data'
          print*
          stop '@ ~227l in rdpar'
        endif
        read(32,9806,iostat=ios)title(68),a11,notes(68)
        a11=amax1(0.01,a11)   ! bare ground equiv. veg height
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading a11 in  the parameter file data'
          print*
          stop '@ ~227l in rdpar'
        endif
        read(32,9806,iostat=ios)title(69),a12,notes(69)
        if(ios.ne.0)then
          print*,'Parameter file ver 9.1 or greater'
          print*,'Error reading a12 in  the parameter file data'
          print*
          stop '@ ~227l in rdpar'
        endif
      endif

      if(ios.ne.0)then
        print*,'Parameter file ver 9.1 or greater'
        print*,'Error in the a1 - a12 section of the parameter file'
        print*
        print*,title(91),a1,notes(91) 
        print*,title(92),a2,notes(92)  
        print*,title(93),a3,notes(93)
        print*,title(94),a4,notes(94)
        print*,title(95),a5,notes(95)
        print*,title(96),a6,notes(96)
        print*,title(97),a7,notes(97)
        print*,title(98),a8,notes(98)
        print*,title(99),a9,notes(99)
        print*,title(67),a10,notes(67)
        print*,title(68),a11,notes(68)
        print*,title(69),a12,notes(69)
        stop '@ ~256 in rdpar'
      endif

      if(iopt.eq.2)print*, 'in rdpar @ 581'

      tempflg='n'
!     a7 - weighting factor for old vs. new sca value default=.90
!     rev. 9.1.04  Oct. 12/01 nk
      if(a7.lt.0.50.or.a7.gt.0.99)then
        print*,' Value for weighting factor - old vs. new sca value'
        print*,' a7 not between 0.5-0.99in the par file'
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
!        a7=0.5
!        tempflg='y'
      endif
      STOP 88888
      if(iopt.eq.2)print*, 'in rdpar @ 605'

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
        pause ' Program paused in rdpar @ 480'
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
        write(51,5110)a9
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

      if(iopt.eq.2)print*, 'in rdpar @ 660'


!     REVISED TO CONSOLIDATE CLASSES
      if(ver.ge.8.99)then
!       read the column headers
        read(32,5009,iostat=ios)(rivtype(i),i=1,nrvr)
        if(ios.ne.0)then
          print*,'Error reading rivtype in the parameter file data'
          print*
          stop '@ ~571 in rdpar'
        endif

!       read the lowerzone function
        read(32,1001,iostat=ios)title(14),(flz_o(i),i=1,nrvr)
        if(ios.ne.0)then
          print*,'Error reading flz in the parameter file data'
          print*
          stop '@ ~578 in rdpar'
        endif
!       findout how many river classes there are:
        do i=1,nrvr
          if(ajunk(i).gt.0.0)then
            nrvr1=i
          endif
        end do

        if(iopt.eq.2)print*, 'in rdpar @ 686'

!       check no of rivers is ok
       read(32,1001,iostat=ios)title(15),(pwr_o(i),i=1,nrvr)
       
        if(ios.ne.0)then
          print*,'Error reading pwr in the parameter file data'
          print*

          stop '@ ~747 in rdpar'
        endif

!       ERROR CHECK FOR REVISED LZS CLASSES
        do i=1,nrvr
!       rev. 9.1.37  Mar.  22/03  - Option to turn off leakage by setting LZF < 0.0
          if(flz_o(i).lt.0.0)flz_o(i)=0.1e-30
          if(pwr_o(i).lt.0.1.or.pwr_o(i).gt.4.0)then   !  Feb. 22/06  nk
            print*
            print*,'The value of pwr for river class ',i
            print*,'is ',pwr_o(i),' which' 
            print*,'is not in the range of 0.1 - 4.0'
            print*,'Please check all values'
            print*
            if(dds_flag.eq.0.and.numa.eq.0)then
              pause 'Hit enter to continue with this value (@ 607)'
            endif
          endif
        end do

        if(iopt.eq.2)print*, 'in rdpar @ 712'

!      pause 33
!       read the over bank multiplier and river roughness parameter
        manningflg='n'
        read(32,1001,iostat=ios)title(16)  
        if(title(16).eq.'R1n  '.or.title(16).eq.'r1n  ')then
!         rev. 9.2.11  Sep.  11/05  - NK: added Manning's n  r1n & r2n
!         Manning's n is used instead of r1
          manningflg='y'

          backspace 32
          read(32,1001,iostat=ios)title(16),(r1n_o(i),i=1,nrvr)   
          if(ios.ne.0)then
            print*,'Error reading r1n in the parameter file data'
            print*
            stop '@ ~625 in rdpar'
          endif

          read(32,1001,iostat=ios)title(17),(r2n_o(i),i=1,nrvr)
          if(ios.ne.0)then
            print*,'Error reading r2n in the parameter file data'
            print*
            stop '@ ~631 in rdpar'
          endif

        else

          print*,'R1 and R2 values no longer accepted'
          print*,'Please change values to Manning`s n'
          print*,'and titles to R1n and R2n respectively'
          Print*,'R2n = R2/10 approximately'
          print*,'Values for width/depth ratio also required'
          print*
          stop 'Program aborted in rdpar @ 845'
        endif
      endif

      if(iopt.eq.2)print*, 'in rdpar @ 817'

      if(ver.ge.9.1)then
        read(32,1001,iostat=ios)title(83),(mndr_o(i),i=1,nrvr)
        if(ios.ne.0)then
          print*,'Error reading mndr in the parameter file data'
          print*
          stop '@ ~675 in rdpar'
        endif

!       read the three bankfull function coefficients
        read(32,1001,iostat=ios)title(84),(aa2_o(i),i=1,nrvr)
        if(ios.ne.0)then
          print*,'Error reading aa2 in the parameter file data'
          print*
          stop '@ ~682 in rdpar'
        endif

        read(32,1001,iostat=ios)title(85),(aa3_o(i),i=1,nrvr)
        if(ios.ne.0)then
          print*,'Error reading aa3 in the parameter file data'
          print*
          stop '@ ~688 in rdpar'
        endif

        read(32,1001,iostat=ios)title(86),(aa4_o(i),i=1,nrvr)
        if(ios.ne.0)then
          print*,'Error reading aa4 in  the parameter file data'
          print*
          stop '@ ~694 in rdpar'
        endif

        read(32,1001,iostat=ios)title(87),(theta_o(i),i=1,nrvr)
        if(ios.ne.0)then
          print*,'Error reading theta in the parameter file data'
          print*
          stop '@ ~700 in rdpar'
        endif

        read(32,1001,iostat=ios)title(88),(widep_o(i),i=1,nrvr)
        if(ios.ne.0)then
          print*,'Error reading widep in the parameter file data'
          print*
          stop '@ ~706 in rdpar'
        endif

        if(manningflg.eq.'y')then
          do i=1,nrvr
            if(widep_o(i).le.0.0)then
              print*,'Width/depth ratio must be specified when using'
              print*,'Manning n as the roughness parameter.'
              print*,'Please fix the par file'
              print*
              stop 'Program aborted in rdpar @ 715'
            endif
          end do
        endif

        read(32,1001,iostat=ios)title(89),(kcond_o(i),i=1,nrvr)
        if(ios.ne.0)then
          print*,'Problems reading kcond in the parameter'
          print*,' '
          stop 'Program aborted in rdpar @ 723'
        endif

!   rev. 9.5.06  Feb.  05/08  - NK: added pool and pool_o in rdpar & route
!   rev. 9.6.01  Mar.  01/10  - NK: rlake parameter added for Manning n correction
        if(ver.ge.9.4)then
          read(32,1001,iostat=ios)title(90),(pool_o(i),i=1,nrvr)
          if(ios.ne.0)then
            print*,'Problems reading pool in the parameter'
            print*,' '
            stop 'Program aborted in rdpar @ 880'
          endif
          read(32,1001,iostat=ios)title(90),(rlake_o(i),i=1,nrvr)
          if(ios.ne.0)then
            print*,'Problems reading rlake in the parameter'
            print*,' '
            stop 'Program aborted in rdpar @ 880'
          endif
        else
!         default value is 0.0 - no pools & riffles
          do i=1,nrvr
            pool_o(i)=0.0    !default
            rlake_o(i)=a2    !default from top of par file
          end do
        endif

!     rev. 9.4.07  May.  15/07  - NK: converted opt to gridded routing parameters
!     rev. 9.5.06  Feb.  05/08  - NK: added pool and pool_o in rdpar & route
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

      else    ! ver.lt.9.1
        print*,'parameter file versions before ver. 9.1'
        print*,'no longer accepted'
        print*,'Please update'
        print*
        stop 'Program aborted in rdpar @ 930'
      endif

      tempflg='n'
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

      if(iopt.eq.2)print*, 'param 5'

      read(32,5009,iostat=ios)(nclass(i),i=1,classcount)
      
      if(ios.ne.0)then
        print*,'Error reading nclass in the parameter file data'
        print*
        stop '@ ~694 in rdpar'
      endif

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
      if(iopt.eq.2)print*, 'param 6'
!      if(ver.ge.9.0)then
!        read(32,*)(iiclass(i),i=1,classcount)
!        print*,(iiclass(i),i=1,classcount)
!      endif
      read(32,1001,iostat=ios)title(6),(ds(i),i=1,classcount)   
      if(ios.ne.0)then
        print*,'Error reading depression storage (ds)'
        print*,' in the parameter file data'
        print*,ds
        print*
        stop '@ ~825 in rdpar'
      endif
      if(ver.ge.7.5)read(32,1001,iostat=ios) 
     *                  title(7),(dsfs(i),i=1,classcount)   
      if(ios.ne.0)then
        print*,'Error reading dsfs in the parameter file data'
        print*
        stop '@ ~832 in rdpar'
      endif

!     The true meaning of REC is that for a slope of 1.0, rec will
!     be the fraction of water taken out of the UZ. This amount is 
!     reduced by slope sl1   nk Jul. 21/04

      read(32,1001,iostat=ios)title(8),(rec(i),i=1,classcount)  
      if(ios.ne.0)then
        print*,'Error reading rec in the parameter file data'
        print*,rec
        print*
        stop '@ ~843 in rdpar'
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

      if(iopt.eq.2)print*, 'param 7'

      if(ijunko.eq.0)then
        STOP 'Program aborted in rdpar.for @ ~116'
      endif
      read(32,1001,iostat=ios)title(9),(ak(i),i=1,classcount)   
      if(ios.ne.0)then
        print*,'Error reading ak in the parameter file data'
        print*,ak
        print*
        stop '@ ~856 in rdpar'
      endif
      if(ver.ge.7.5)
     *    read(32,1001,iostat=ios)title(10),(akfs(i),i=1,classcount)
      if(ios.ne.0)then
        print*,'Error reading akfs in the parameter file data'
        print*,akfs
        print*
        stop '@ ~864 in rdpar'
      endif
      if(ak(classcount-1).ge.0.0.or.akfs(classcount-1).ge.0.0)then
        print*,'Water class not specified'
        print*,'ak and akfs NOT given a -ve value'
        print*,'ak   value found =',ak(classcount-1)
        print*,'akfs value found =',akfs(classcount-1)
        print*,'Please correct the par file'
        print*,'No of classed expected =',classcount
        print*
        stop 'Program aborted in rdpar @ 1104'
      endif

      read(32,1001,iostat=ios)title(11),(retn(i),i=1,classcount)   
      if(ios.ne.0)then
        print*,'Error reading retn in the parameter file data'
        print*
        stop '@ ~870 in rdpar'
      endif
      read(32,1001,iostat=ios)title(12),(ak2(i),i=1,classcount)   
      if(ios.ne.0)then
        print*,'Error reading ak2 in the parameter file data'
        print*
        stop '@ ~876 in rdpar'
      endif
      if(ver.ge.7.7)
     *    read(32,1001,iostat=ios)title(13),(ak2fs(i),i=1,classcount)
      if(ios.ne.0)then

        print*,'Error reading ak2fs in the parameter file data'
        print*
        stop '@ ~883 in rdpar'
      endif

!     REVISED PAR FILE TO GROUP CLASSES:
      if(ver.lt.8.99)then
        read(32,1001,iostat=ios)title(14),(flz(i),i=1,nrvr)
        if(ios.ne.0)then
          print*,'Error reading flz in the parameter file data'
          print*
          stop '@ ~896 in rdpar'
        endif
        read(32,1001,iostat=ios)title(15),(pwr(i),i=1,nrvr)
        if(ios.ne.0)then
          print*,'Error reading pwr in the parameter file data'
          print*
          stop '@ ~902 in rdpar'
        endif
!       ERROR CHECK FOR REVISED LZS CLASSES
        do i=1,nrvr
          if(flz(i).eq.0.0.or.pwr(i).eq.0.0)then
            print*,'Program revised to have nrvr classes for'
            print*,'flz & pwr - please fix parameter file'
            print*,' '
            STOP 'Program halted in rdpar.for line 123'
          endif
        end do
        read(32,1001,iostat=ios)title(16),(r1(i),i=1,nrvr)   
        if(ios.ne.0)then
          print*,'Error reading r1 in the parameter file data'
          print*
          stop '@ ~917 in rdpar'
        endif
        read(32,1001,iostat=ios)title(17),(r2(i),i=1,nrvr)
        if(ios.ne.0)then
          print*,'Error reading r2 in the parameter file data'
          print*
          stop '@ ~923 in rdpar'
        endif
      endif

!     rev. 9.2.37  Mar.  31/06  - NK: Removed impervious area as special class
!         now need the values below:
      if(ds(classcount).eq.0.0)ds(classcount)=1.0
      if(dsfs(classcount).eq.0.0)dsfs(classcount)=1.0
      if(ds(classcount).eq.0.0)ds(classcount)=1.0
      if(rec(classcount).eq.0.0)rec(classcount)=0.1
      if(ak(classcount).eq.0.0)ak(classcount)=1.0e-10
      if(akfs(classcount).eq.0.0)akfs(classcount)=1.0e-10
      if(ak2(classcount).eq.0.0)ak2(classcount)=1.0e-10
      if(ak2fs(classcount).eq.0.0)ak2fs(classcount)=1.0e-10
      if(r3(classcount).eq.0.0)r3(classcount)=0.04
      if(r3fs(classcount).eq.0.0)r3fs(classcount)=0.04
      if(r4(classcount).eq.0.0)r4(classcount)=10.0

      if(iopt.eq.2)print*, 'param 8'

      read(32,1001,iostat=ios)title(18),(r3(i),i=1,classcount) 
      if(ios.ne.0)print*,'last read',title(18),(r3(i),i=1,classcount)
      if(ver.ge.7.5)read(32,1001,iostat=ios) 
     *                  title(19),(r3fs(i),i=1,classcount) 
      if(ios.ne.0)print*,'last read',title(19),(r3fs(i),i=1,classcount)
      read(32,1001,iostat=ios)title(20),(r4(i),i=1,classcount)   
      if(ios.ne.0)print*,'last read',title(20),(r4(i),i=1,classcount)

      read(32,1001,iostat=ios)title(21),(chnl(i),i=1,5) 
      if(ios.ne.0)print*,'last read',title(21),(chnl(i),i=1,classcount)

      read(32,1001,iostat=ios)title(22),(fm(ii),ii=1,classcount)
      if(ios.ne.0)print*,'last read',title(22),(fm(i),i=1,classcount)
      read(32,1001,iostat=ios)title(23),(base(ii),ii=1,classcount)
      if(ios.ne.0)print*,'last read',title(23),(base(i),i=1,classcount)
      read(32,1001,iostat=ios)title(24),(fmn(ii),ii=1,classcount)
      if(ios.ne.0)print*,'last read',title(24),(fmn(i),i=1,classcount)
      read(32,1001,iostat=ios)title(25),(uadj(ii),ii=1,classcount)
      if(ios.ne.0)print*,'last read',title(25),(uadj(i),i=1,classcount)
      read(32,1001,iostat=ios)title(26),(tipm(ii),ii=1,classcount)
      if(ios.ne.0)print*,'last read',title(26),(tipm(i),i=1,classcount)
      read(32,1001,iostat=ios)title(27),(rho(ii),ii=1,classcount)
      if(ios.ne.0)print*,'last read',title(27),(rho(i),i=1,classcount)
      read(32,1001,iostat=ios)title(28),(whcl(ii),ii=1,classcount)
      if(ios.ne.0)print*,'last read',title(28),(whcl(i),i=1,classcount)
!     rev. 9.1.76  Mar.  09/05  - NK: separated glacier parameters in par file
      if(iopt.eq.2)print*, 'param 8a'
      if(ver.ge.9.3)then
        read(32,5004,iostat=ios)title(29),fmadjust
        if(ios.ne.0)print*,'last read =',title(29),fmadjust
        read(32,5004,iostat=ios)title(102),fmalow
        if(ios.ne.0)print*,'last read =',title(102),fmalow
        read(32,5004,iostat=ios)title(103),fmahigh
        if(ios.ne.0)print*,'last read =',title(103),fmahigh
        read(32,5004,iostat=ios)title(104),gladjust
        if(ios.ne.0)print*,'last read =',title(104),gladjust
        read(32,5004,iostat=ios)title(105),rlapse
!     rev. 9.5.65  Sep.  26/09  - NK: lapse rate changed from dC per 100 m to dC per m
        if(ios.ne.0)print*,'last read =',title(105),rlapse
        read(32,5004,iostat=ios)title(105),elvref
        if(ios.ne.0)print*,'last read =',title(106),elvref
        if(iopt.eq.2)print*, 'param 8b1'
      else
        read(32,5004,iostat=ios) 
     *    title(29),fmadjust,fmalow,fmahigh,gladjust,rlapse,elvref
        if(ios.ne.0)print*,'last read =',
     *    title(29),fmadjust,fmalow,fmahigh,gladjust,rlapse,elvref
!             ELEVATION IS IN HUNDREDS OF METERS
!             LAPSE RATE IS IN DEGREE C / 100 m - USUALLY 0.5 DEGREE C
        if(iopt.eq.2)print*, 'param 8b2'
        iverflg=1
      endif

!   rev. 9.5.65  Sep.  26/09  - NK: lapse rate changed from dC per 100 m to dC per m
      if(rlapse.gt.0.1)then
        print*,'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
        print*,'WARNING:'
        print*, 'Value for RLAPSE is too large'
        print*, 'RLAPSE used to be degree C per 100m'
        print*, 'this has been changed to degree C per metre'
        print*, 'Please correct the par file'
        print*,'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
        stop 'Program aborted in rdpar @ 1185'
      endif

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
      if(ios.ne.0)then
        print*,''
        print*,'Error in ds - whcl part of the parameter file'
        print*
        stop '@ ~469 in rdpar'
      endif

      if(ver.ge.8.5)then
        read(32,9804,iostat=ios)title(31),flgevp2,junk
!       REV. 8.71 - Feb.  24/98 -  ADDED EVPFLG2 TO INPEVTA.FOR
!       REV. 8.71 - deleted Nov. 6/98
        read(32,9802,iostat=ios)title(32),albe
        write(51,9802)title(32),albe
        write(51,*)'after reading albe'
        if(ios.ne.0)print*,'last read =',title(32),albe
!! TS: ADDED READ FOR IMPERV ALBEDO VALUE AND/OR DEFAULT VALUE: MAR 27/06
        read(32,9802,iostat=ios)title(33),(alb(i),i=1,classcount)
        if(ios.ne.0)print*,'last read =',title(33),alb
        if(alb(classcount).eq.0)then
          alb(classcount)=0.18
          write(51,*)
     *     '** Albedo for IMPERV class not found, assumed 0.18 **'
        endif     
        read(32,9802,iostat=ios)title(34),(fpet(i),i=1,classcount)
        if(ios.ne.0)print*,'last read =',title(34),fpet
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
        read(32,9802,iostat=ios)title(35),(ftall(i),i=1,classcount)
        if(ios.ne.0)print*,'last read =',title(35),ftall
        if(wetflg.eq.'y')then
          if(fpet(classcount-2).ne.fpet(classcount-3).or.
     *          ftall(classcount-2).ne.ftall(classcount-3))then
            print*
            print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print*,'WARNING'
            print*,'fpet & ftal need to be the same for the coupled'
            print*,'and non-coupled parent land cover class'
            print*,'i.e. the last two wetland classes'
            print*,'Program will set coupled class = source class'
            print*,'So.....DO NOT OPTIMIZE FPET OR FTALL FOR '
            PRINT*,'          THE COUPLED WETLAND CLASS'
            print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print*
            fpet(classcount-2)=fpet(classcount-3)
            ftall(classcount-2)=ftall(classcount-3)
          endif
        endif
        read(32,9802,iostat=ios)title(36),(flint(i),i=1,classcount)
        if(ios.ne.0)print*,'last read =',title(36),flint
        read(32,9802,iostat=ios)title(37),(fcap(i),i=1,classcount)
        if(ios.ne.0)print*,'last read =',title(37),fcap
        read(32,9802,iostat=ios)title(38),(ffcap(i),i=1,classcount)
        if(ios.ne.0)print*,'last read =',title(38),ffcap
        do i=1,classcount
          if(ffcap(i).le.0.0)then
            print*,'pwp(',i,') Permanent Wilting Point set too low.'
            print*,'pwp(',i,')  Must be > 0.0'
            print*,'pwp(',i,') set to 0.01'
            print*
            ffcap(i)=0.01
          endif
        end do
        read(32,9802,iostat=ios)title(39),(spore(i),i=1,classcount)
        if(ios.ne.0)print*,'last read =',title(39),spore
        do ii=1,classcount
          if(spore(ii).le.0.0)then
            print*,'spore(',ii,')=',spore(ii)
            print*,' This value can not be .le. 0.0'
            print*,' Please correct the parameter file'
            print*
            stop 'Program aborted in rdpar at ~1254'
          endif
        end do

!     rev. 9.2.37  Mar.  31/06  - NK: Removed impervious area as special class
!         now need the values below:
      if(fpet(classcount).eq.0.0)fpet(classcount)=1.0
      if(ftall(classcount).eq.0.0)ftall(classcount)=1.0
      if(flint(classcount).eq.0.0)flint(classcount)=1.0
      if(fcap(classcount).eq.0.0)fcap(classcount)=1.0
      if(ffcap(classcount).eq.0.0)ffcap(classcount)=1.0
      if(spore(classcount).eq.0.0)spore(classcount)=1.0
        
!     rev. 9.1.80  Mar.  31/05  - NK: added sublimation   (sublim)
!                  modified this Feb. 6/06
        read(32,9801,iostat=ios)title(40),tempa1
!       tempa1 no longer used. Line used to read sublimation factor
        if(title(40).eq.'sublm')then
          backspace 32
          read(32,9801,iostat=ios)title(40),
     *                       (sublim_factor(i),i=1,classcount)
          write(51,9801)title(40),(sublim_factor(i),i=1,classcount)
          if(ios.ne.0)print*,'last read =',title(40),sublim_factor
          do i=1,classcount
            if(sublim_factor(i).gt.0.5)then
             if(numa.eq.0.and.dds_flag.eq.0)then    !added Apr.03/10 nk
                print*,'Warning: sublimation factor for class',i
                print*,'and is given as ',sublim_factor(i)
                print*,'is set too high. Sublim reduced to 0.5'
              endif
            endif
            sublim_factor(i)=amin1(0.5,sublim_factor(i))
          end do
        else
          do i=1,classcount
            sublim_factor(i)=0.0
          end do
        endif
        sublimflg=1  !can only use sublim_factor with the old par file
        
        read(32,9801,iostat=ios)title(41),tempa2
        write(51,9801,iostat=ios)title(41),tempa2
        if(ios.ne.0)print*,'last read =',title(41),tempa2
        read(32,9801,iostat=ios)title(42),tempa3
        write(51,9801,iostat=ios)title(42),tempa3
        write(51,*)'after reading temp3'
        if(ios.ne.0)print*,'last read =',title(42),tempa3
        if(tempa3.lt.0.0001)then
          tempa3=0.0001
          write(51,*)'WARNING:'
          write(51,*)'Temp3 <0.0001  set temp3=0.001'
          write(51,*)
        endif
        read(32,9801,iostat=ios)title(43),tton
        write(51,9801,iostat=ios)title(43),tton
        if(ios.ne.0)print*,'last read =',title(43),tton
        if(tton.lt.0.000)tton=0.0
        read(32,9801,iostat=ios)title(44),lat
        if(ios.ne.0)print*,'last read =',title(44),lat
        read(32,9803,iostat=ios)title(45),(diff(i),i=1,12)
        if(ios.ne.0)print*,'last read =',title(45),diff
        read(32,9803,iostat=ios)title(46),(hu(i),i=1,12)
        if(ios.ne.0)print*,'last read =',title(46),hu
        read(32,9803,iostat=ios)title(47),(pres(i),i=1,12)
        if(ios.ne.0)print*,'last read =',title(47),pres
      endif
      if(iopt.eq.2)print*, 'param 10'
      if(ios.ne.0)then
        print*,''
        print*,'Error in the evaporation part of the parameter file'
        print*
        stop '@ ~497 in rdpar'
      endif

      do ii=1,classcount
        if(spore(ii).le.0.0)then
          print*,'ii,spore',ii,spore(ii)
          print*,' line title =',title(39)
          print*,' value for spore can not be <= 0.0'
          STOP 'Program stopped in rdpar @ 9443'
        endif
      end do

      read(32,5005,iostat=ios)heading(2)
!     read the interception capacities
      do ii=1,classcount
        read(32,*,iostat=ios)title(50+ii),(h(i,ii),i=1,12)
        write(51,1003)title(50+ii),(h(i,ii),i=1,12)
        if(ios.ne.0)then
          print*,'error reading the h(i,ii) values'
          print*,'Please ensure there is a line for each class'
          print*,'including the impervious class (last class)'
          print*,'There is now 1 more line with the new ensim'
          print*,'format shed file'
          print*
          stop 'Program aborted in rdpar @ 1312'
        endif
        if(ios.ne.0)then
          print*,''
          print*,'Error reading title(',50+ii,') in the parameter file'
          print*,'ii=',ii
          do j=1,ii
            print*,title(50+j),(h(i,j),i=1,12)
          end do
          print*,'last heading read =',heading(2)
          print*
          stop 'program aborted @ ~934 in rdpar'
        endif
      end do

!     rev. 9.2.37  Mar.  31/06  - NK: Removed impervious area as special class
!         now need the values below:
      do i=1,12
        h(i,classcount)=0.0
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

        if(ver.lt.7.5)then
          akfs(ii)=ak(ii)
          dsfs(ii)=ds(ii)
          r3fs(ii)=r3(ii)
        endif

        if(ver.lt.7.7)then
          ak2fs(ii)=ak2(ii)
        endif

        if(ak2fs(ii).gt.1.0)then
          write(*,1027)ii
1027      format(' AK2fs(',i5,') =',e12.3)
          STOP ' Please change value for AK2fs - must be lt 1'
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

      if(iopt.eq.2)print*, 'param 10b'
        read(32,5005)heading(3)

!     This was fixed between version 9.1.11 and 9.1.13
!     It used to be all sin(...). End result: pet too high
      sinlat=sin(3.1416*lat/180.0)
      coslat=cos(3.1416*lat/180.0)
      tanlat=tan(3.1416*lat/180.0)

!     PART II     OPTIMIZATION      SECTION
!     PART II     OPTIMIZATION      SECTION
!     PART II     OPTIMIZATION      SECTION
!     PART II     OPTIMIZATION      SECTION

      optflg=.false.
      if(numa.gt.0.or.dds_flag.eq.1)optflg=.true.

!     rev. 9.5.43  Oct.  27/08  - NK: changed bottom part of par file to be free format

ccccc      if(numa.ge.1.or.ver.lt.9.30.or.dds_flag.eq.1)then   !nk 28/10/07

!       this section needed so the par file can be updated with new.par

ccccc   if(ver.eq.par_file_version_latest)numa=1  !used as a flag if =1
!       dds_flag=1    !used as a flag if =1

!      if(numa.ge.1.or.iverflg.eq.1)then

!       IF NUMA >0 THEN WE ARE OPTIMIZING AND WE NEED THE LIMITS, ETC. 
!       READ THE OPTIMIZATION LIMITS AND STARTING VALUES
!       NOTE: FOR VER 6.2 AND LATER THE VALUES IN THE BOTTOM 
!       TABLE ARE USED AS THE INITIAL VALUES AND THE EARLIER VALUES ARE
!       IGNORED. THIS MAKES IT A LITTLE EASIER TO CHECK THE CONSTRAINTS.

!  REV. 8.32 - June  13/97 -  BYPASSED NON-FLAGGED PARAMETERS IN OPTA

!       WHEN A PARAMETER IS NOT BEING OPTIMIZED, THE VALUE FROM THE
!       TOP OF THE TABLE HAS TO BE USED. IN THE NEW.PAR FILE, THE TOP 
!       VALUE WILL REPLACE THE LOWER VALUE
        if(iopt.eq.2)print*, 'param 11'

!     rev. 9.7.02  Jun.  24/10  - NK: fixed bug in rdpar for classcount for imp area
!       added Apr. 07/10 nk
c        if(classcount_local.gt.classcount)then
        if(classcount_local.gt.classcount)then
          print*
          print*,'Fatal error:'
          print*,'No of classes in opt part of the par file is'
          print*,'greater than the no in the first part of the file'
          print*,'and the no of classes in the shd file.'
          print*,'Please fix the par file so the no of parameters to'
          print*,'optimized is <= the number in the shd & top part'
          print*,'of the par file. The value of typeo must be equal'
          print*,'to the no of classes in part 2 of the par file.'
          print*
          stop 'Program aborted in rdpar @ 1527'
        endif

!     rev. 9.07    Mar.  14/01  - fixed use of opt par's  for numa=0  
        if(classcount_local.ge.1)then
          do i=1,classcount_local
            read(32,*,iostat=ios)
     *        title(101),akdlt(i),aklow(i),akhgh(i),ajunk(i)
              if(iopt.eq.2)print*, 'param 11b'
          end do
          if(classcount_local.lt.classcount)then
!           this is added to give values tot he remainder for the new par file
            do i=classcount_local+1,classcount
              akdlt(i)=akdlt(classcount_local)
              aklow(i)=aklow(classcount_local)
              akhgh(i)=akhgh(classcount_local)
            end do
          endif
          do i=1,classcount_local
            if(akdlt(i).gt.0.0.and.optflg)ak(i)=ajunk(i)
          end do

          read(32,*,iostat=ios)
     *    (title(101),akfsdlt(i),akfslow(i),akfshgh(i),ajunk(i),
     *                         i=1,classcount_local)
          if(classcount_local.lt.classcount)then
!           this is added to give values tot he remainder for the new par file
            do i=classcount_local+1,classcount
              akfsdlt(i)=akfsdlt(classcount_local)
              akfslow(i)=akfslow(classcount_local)
              akfshgh(i)=akfshgh(classcount_local)
            end do
          endif
          do i=1,classcount_local
            if(akfsdlt(i).gt.0.0.and.optflg)akfs(i)=ajunk(i)
          end do

          read(32,*,iostat=ios)(title(101),
     *    recdlt(i),reclow(i),rechgh(i),ajunk(i),i=1,classcount_local)
          if(classcount_local.lt.classcount)then
!           this is added to give values tot he remainder for the new par file
            do i=classcount_local+1,classcount
              recdlt(i)=recdlt(classcount_local)
              reclow(i)=reclow(classcount_local)
              rechgh(i)=rechgh(classcount_local)
            end do
          endif
          do i=1,classcount_local
            if(recdlt(i).gt.0.0.and.optflg)rec(i)=ajunk(i)
          end do
          
          read(32,*,iostat=ios)(title(101),r3dlt(i),r3low(i),
     *               r3hgh(i),ajunk(i),i=1,classcount_local)
          if(classcount_local.lt.classcount)then
!           this is added to give values tot he remainder for the new par file
            do i=classcount_local+1,classcount
              r3dlt(i)=r3dlt(classcount_local)
              r3low(i)=r3low(classcount_local)
              r3hgh(i)=r3hgh(classcount_local)
            end do
          endif
          do i=1,classcount_local
            if(r3dlt(i).gt.0.0.and.optflg)r3(i)=ajunk(i)
          end do

!     rev. 9.5.72  Oct.  12/09  - NK: fixed bug in rdpar setting init values for fpet & ftal
          if(ver.ge.8.5)then
            read(32,*)
     *      (title(101),fpetdlt(i),fpetlow(i),fpethgh(i),ajunk(i),
     *                                           i=1, classcount_local)
          if(classcount_local.lt.classcount)then
!           this is added to give values tot he remainder for the new par file
            do i=classcount_local+1,classcount
              fpetdlt(i)=fpetdlt(classcount_local)
              fpetlow(i)=fpetlow(classcount_local)
              fpethgh(i)=fpethgh(classcount_local)
            end do
          endif
            do i=1,classcount_local
              if(fpetdlt(i).gt.0.0.and.optflg)fpet(i)=ajunk(i)
            end do
            read(32,*)
     *      (title(101),ftalldlt(i),ftalllow(i),ftallhgh(i),ajunk(i),
     *                                            i=1,classcount_local)
          if(classcount_local.lt.classcount)then
!           this is added to give values tot he remainder for the new par file
            do i=classcount_local+1,classcount
              ftalldlt(i)=ftalldlt(classcount_local)
              ftalllow(i)=ftalllow(classcount_local)
              ftallhgh(i)=ftallhgh(classcount_local)
            end do
          endif
            do i=1,classcount_local
              if(ftalldlt(i).gt.0.0.and.optflg)ftall(i)=ajunk(i)
            end do
          endif

          read(32,*,iostat=ios)(title(101),
     *    fmdlt(i),fmlow(i),fmhgh(i),ajunk(i),i=1,classcount_local)
          if(classcount_local.lt.classcount)then
!           this is added to give values tot he remainder for the new par file
            do i=classcount_local+1,classcount
              fmdlt(i)=fmdlt(classcount_local)
              fmlow(i)=fmlow(classcount_local)
              fmhgh(i)=fmhgh(classcount_local)
            end do
          endif
          do i=1,classcount_local
            if(fmdlt(i).gt.0.0.and.optflg)fm(i)=ajunk(i)
          end do

          read(32,*,iostat=ios)(title(101),
     *    basdlt(i),baslow(i),bashgh(i),ajunk(i),i=1,classcount_local)
          if(classcount_local.lt.classcount)then
!           this is added to give values tot he remainder for the new par file
            do i=classcount_local+1,classcount
              basdlt(i)=basdlt(classcount_local)
              baslow(i)=baslow(classcount_local)
              bashgh(i)=bashgh(classcount_local)
            end do
          endif
          do i=1,classcount_local
            if(basdlt(i).gt.0.0.and.optflg)base(i)=ajunk(i)
          end do

c          read(32,*,iostat=ios)(title(101),
c     *    fmndlt(i),fmnlow(i),fmnhgh(i),ajunk(i),i=1,classcount_local)
c          do i=1,classcount_local
c            if(fmndlt(i).gt.0.0.and.optflg)fmn(i)=ajunk(i)
c          end do

          read(32,*,iostat=ios)(title(101),
     *    subdlt(i),sublow(i),subhgh(i),ajunk(i),i=1,classcount_local)
          if(classcount_local.lt.classcount)then
!           this is added to give values tot he remainder for the new par file
            do i=classcount_local+1,classcount
              subdlt(i)=subdlt(classcount_local)
              sublow(i)=sublow(classcount_local)
              subhgh(i)=subhgh(classcount_local)
            end do
          endif
          do i=1,classcount_local
            if(subdlt(i).gt.0.0.and.optflg)sublim_factor(i)=ajunk(i)
          end do

          read(32,*,iostat=ios)
     *    (title(101),retndlt(i),retnlow(i),retnhgh(i),ajunk(i),
     *                                            i=1,classcount_local)
          if(classcount_local.lt.classcount)then
!           this is added to give values tot he remainder for the new par file
            do i=classcount_local+1,classcount
              retndlt(i)=retndlt(classcount_local)
              retnlow(i)=retnlow(classcount_local)
              retnhgh(i)=retnhgh(classcount_local)
            end do
          endif
          do i=1,classcount_local
            if(retndlt(i).gt.0.0.and.optflg)retn(i)=ajunk(i)
          end do

          read(32,*,iostat=ios)
     *    (title(101),ak2dlt(i),ak2low(i),ak2hgh(i),ajunk(i),
     *                                            i=1,classcount_local)
          if(classcount_local.lt.classcount)then
!           this is added to give values tot he remainder for the new par file
            do i=classcount_local+1,classcount
              ak2dlt(i)=ak2dlt(classcount_local)
              ak2low(i)=ak2low(classcount_local)
              ak2hgh(i)=ak2hgh(classcount_local)
            end do
          endif
          do i=1,classcount_local
            if(ak2dlt(i).gt.0.0.and.optflg)ak2(i)=ajunk(i)
          end do

          read(32,*,iostat=ios)
     *    (title(101),ak2fsdlt(i),ak2fslow(i),ak2fshgh(i),ajunk(i),
     *                                           i=1,classcount_local)
          if(classcount_local.lt.classcount)then
!           this is added to give values tot he remainder for the new par file
            do i=classcount_local+1,classcount
              ak2fsdlt(i)=ak2fsdlt(classcount_local)
              ak2fslow(i)=ak2fslow(classcount_local)
              ak2fshgh(i)=ak2fshgh(classcount_local)
            end do
          endif
          do i=1,classcount_local
            if(ak2fsdlt(i).gt.0.0.and.optflg)ak2fs(i)=ajunk(i)
          end do
          if(iopt.eq.2)print*, 'param 10c'

        endif
        if(iopt.eq.2)print*, 'param 12'

        if(ver.lt.9.0)then
!         old version without wetlands
          if(nbsn.ge.1)then
            read(32,*,iostat=ios)(title(101),
     *       flzdlt(i),flzlow(i),flzhgh(i),ajunk(i),i=1,nbsn)
            if(nbsn.lt.nrvr)then
!             this is added to give values tot he remainder for the new par file
              do i=nbsn+1,nrvr
                flzdlt(i)=flzdlt(nbsn)
                flzlow(i)=flzlow(nbsn)
                flzhgh(i)=flzhgh(nbsn)
              end do
            endif
            do i=1,nbsn
              if(flzdlt(i).gt.0.0.and.optflg)flz_o(i)=ajunk(i)
            end do
            read(32,*,iostat=ios)(title(101),
     *       pwrdlt(i),pwrlow(i),pwrhgh(i),ajunk(i),i=1,nbsn)  
            if(nbsn.lt.nrvr)then
!             this is added to give values tot he remainder for the new par file
              do i=nbsn+1,nrvr
                pwrdlt(i)=pwrdlt(nbsn)
                pwrlow(i)=pwrlow(nbsn)
                pwrhgh(i)=pwrhgh(nbsn)
              end do
            endif
            do i=1,nbsn
              if(pwrdlt(i).gt.0.0.and.optflg)pwr_o(i)=ajunk(i)
            end do
            read(32,*,iostat=ios)
     *        title(101),a5dlt,a5low,a5hgh,ajunk(1)
            if(a5dlt.gt.0.0.and.optflg)a5=ajunk(1)
            read(32,*,iostat=ios)
     *       (title(101),r2dlt(i),r2low(i),r2hgh(i),ajunk(i),i=1,nbsn)
            if(title(101).ne.'pwr  ')then
              print*,' Please note that bottom of par file needs to be '
              print*,' changed. Order should be pwr(nbsn),A5,A9,A10,A11'
              print*,' Hit enter to continue if not optimizing these'
        !      pausee ' else hit ^break'
            endif
            do i=1,nbsn
              if(r2dlt(i).gt.0.0.and.optflg)r2(i)=ajunk(i)
            end do
          endif
          if(iopt.eq.2)print*, 'param 10d'
        else    !    ver.ge.9.0
!         new version for wetlands
          if(nbsn.ge.1)then
            read(32,*,iostat=ios)(title(101),
     *       flzdlt(i),flzlow(i),flzhgh(i),ajunk(i),i=1,nbsn)
            if(nbsn.lt.nrvr)then
!             this is added to give values tot he remainder for the new par file
              do i=nbsn+1,nrvr
                flzdlt(i)=flzdlt(nbsn)
                flzlow(i)=flzlow(nbsn)
                flzhgh(i)=flzhgh(nbsn)
              end do
            endif
            do i=1,nbsn
              if(flzdlt(i).gt.0.0.and.optflg)flz_o(i)=ajunk(i)
            end do
            read(32,*,iostat=ios)(title(101),
     *       pwrdlt(i),pwrlow(i),pwrhgh(i),ajunk(i),i=1,nbsn)  
            if(nbsn.lt.nrvr)then
!             this is added to give values tot he remainder for the new par file
              do i=nbsn+1,nrvr
                pwrdlt(i)=pwrdlt(nbsn)
                pwrlow(i)=pwrlow(nbsn)
                pwrhgh(i)=pwrhgh(nbsn)
              end do
            endif
            do i=1,nbsn
              if(pwrdlt(i).gt.0.0.and.optflg)pwr_o(i)=ajunk(i)
            end do
!           rev. 9.2.11  Sep.  15/05  - NK: added Manning's n  r1n & r2n
            read(32,*,iostat=ios)(title(101),
     *       r2ndlt(i),r2nlow(i),r2nhgh(i),ajunk(i),i=1,nbsn)
            if(nbsn.lt.nrvr)then
!             this is added to give values tot he remainder for the new par file
              do i=nbsn+1,nrvr
                r2ndlt(i)=r2ndlt(nbsn)
                r2nlow(i)=r2nlow(nbsn)
                r2nhgh(i)=r2nhgh(nbsn)
              end do
            endif
            do i=1,nbsn
              if(r2ndlt(i).gt.0.0.and.optflg)r2n_o(i)=ajunk(i)
            end do

          endif

          if(ver.lt.9.2)then
            read(32,*,iostat=ios)
     *        title(101),a5dlt,a5low,a5hgh,ajunk(1)
            if(a5dlt.gt.0.0.and.optflg)a5=ajunk(1)
            read(32,*,iostat=ios)
     *        title(101),a9dlt,a9low,a9hgh,ajunk(1)
            if(a9dlt.gt.0.0.and.optflg)a9=ajunk(1)
            read(32,*,iostat=ios)
     *        title(101),a10dlt,a10low,a10hgh,ajunk(1)
            if(a10dlt.gt.0.0.and.optflg)a10=ajunk(1)
            read(32,*,iostat=ios)
     *        title(101),a11dlt,a11low,a11hgh,ajunk(1)
            if(a11dlt.gt.0.0.and.optflg)a11=ajunk(1)
          else      !   ver.ge.9.2
            read(32,*,iostat=ios)(title(101),
     *       thetadlt(i),thetalow(i),thetahgh(i),ajunk(i),i=1,nbsn)  
             if(nbsn.lt.nrvr)then
!             this is added to give values tot he remainder for the new par file
              do i=nbsn+1,nrvr
                thetadlt(i)=thetadlt(nbsn)
                thetalow(i)=thetalow(nbsn)
                thetahgh(i)=thetahgh(nbsn)
              end do
            endif
           do i=1,nbsn
              if(thetadlt(i).gt.0.0.and.optflg)theta_o(i)=ajunk(i)
            end do
            read(32,*,iostat=ios)(title(101),
     *       kconddlt(i),kcondlow(i),kcondhgh(i),ajunk(i),i=1,nbsn)  
!     rev. 9.2.16  Oct.  10/05  - NK: Fixed bug for widep in rdpar
            if(nbsn.lt.nrvr)then
!             this is added to give values tot he remainder for the new par file
              do i=nbsn+1,nrvr
                kconddlt(i)=kconddlt(nbsn)
                kcondlow(i)=kcondlow(nbsn)
                kcondhgh(i)=kcondhgh(nbsn)
              end do
            endif
            do i=1,nbsn
              if(kconddlt(i).gt.0.0.and.optflg)kcond_o(i)=ajunk(i)
            end do

!           pool() isnot optimizable

!     rev. 9.6.01  Mar.  01/10  - NK: rlake parameter added for Manning n correction
            if(ver.ge.9.4)then
              read(32,*,iostat=ios)(title(101),
     *          rlakedlt(i),rlakelow(i),rlakehgh(i),ajunk(i),i=1,nbsn)  
            if(nbsn.lt.nrvr)then
!             this is added to give values tot he remainder for the new par file
              do i=nbsn+1,nrvr
                rlakedlt(i)=rlakedlt(nbsn)
                rlakelow(i)=rlakelow(nbsn)
                rlakehgh(i)=rlakehgh(nbsn)
              end do
            endif
              do i=1,nbsn
                if(rlakedlt(i).gt.0.0.and.optflg)rlake_o(i)=ajunk(i)
              end do
            else
              do i=1,nbsn
                rlakedlt(i)=-0.1
                rlakelow(i)=0.0
                rlakehgh(i)=3.0
                rlake_o(i)=0.0
              end do
              iverflg=1
            endif
            read(32,*,iostat=ios)
     *        title(101),a5dlt,a5low,a5hgh,ajunk(1)
            if(a5dlt.gt.0.0.and.optflg)a5=ajunk(1)
          endif
!     rev. 9.6.06  Apr.  18/10  - NK: added glacier adjust for optimization
          if(ver.ge.9.5)then
            read(32,*,iostat=ios)
     *        title(101),fmadjustdlt,fmadjustlow,fmadjusthgh,ajunk(1)
              if(fmadjustdlt.gt.0.0.and.optflg)fmadjust=ajunk(1)
            read(32,*,iostat=ios)
     *        title(101),fmalowdlt,fmalowlow,fmalowhgh,ajunk(1)
              if(fmalowdlt.gt.0.0.and.optflg)gladjust=ajunk(1)
            read(32,*,iostat=ios)
     *        title(101),fmahighdlt,fmahighlow,fmahighhgh,ajunk(1)
              if(fmahighdlt.gt.0.0.and.optflg)gladjust=ajunk(1)
            read(32,*,iostat=ios)
     *        title(101),gladjustdlt,gladjustlow,gladjusthgh,ajunk(1)
              if(gladjustdlt.gt.0.0.and.optflg)gladjust=ajunk(1)
          else
            fmadjustdlt=-1.0
            fmadjustlow=0.5
            fmadjusthgh=1.0
            fmalowdlt=-0.1
            fmalowlow=0.1
            fmalowhgh=0.75
            fmahighdlt=-0.1
            fmahighlow=0.75
            fmahighhgh=1.5
            gladjustdlt=-0.1
            gladjustlow=0.5
            gladjusthgh=1.5
            iverflg=1
          endif

          if(iopt.eq.2)print*, 'param 10e'

        endif
        if(iopt.eq.2)print*, 'param 20'

c!     rev. 9.4.07  May.  15/07  - NK: converted opt to gridded routing parameters
c        do n=1,naa
c          flz(n)=flz_o(ibn(n))
c          pwr(n)=pwr_o(ibn(n))
c          r2n(n)=r2n_o(ibn(n))
c          theta(n)=theta_o(ibn(n))
c          kcond(n)=kcond_o(ibn(n))
c          if(ver.ge.9.4)then
c            pool(n)=pool_o(ibn(n))
c            rlake(n)=rlake_o(ibn(n))
c          else
c            pool(n)=0.0
c            rlake(n)=a2
c          endif
c        end do



!       NOTE:
!       TABLE WILL BE BLENDED DEPENDING ON WHETHER A PARAMETER WAS 
!       FLAGGED TO BE OPTIMIZED
!       VALUES FROM TOP PART OF TABLE IF NOT OPTIMIZED
!       VALUES FROM LOWER PART IF FLAGGED TO BE OPTIMIZED
!       THIS HAPPENS ONLY IN AN OPT RUN


!       DDELTA IS CHANGED DURING EXECUTION - NEED TEMP VARIABLE TO 
!       SAVE ORIGINAL VALUES:

ccccc      endif

!     close the parameter file
      close(unit=32) 

d      print*,'Closed unit 32  fln=',fln(2)

      if(iopt.eq.2)print*, 'param 11'

      if(numa.eq.0.and.abs(dds_flag).eq.0)then                     
!       Echos par file TO SIMOUT/SPL.TXT FOR THE RECORD
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c        call write_par(51,iverflg)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      endif

      if(iopt.eq.2)print*, 'param 12'

      if(ios.ne.0.and.dds_flag.eq.0)then
        print*,''
        print*,' An error has been encountered in the parameter file'
        print*,' Please see the simout/spl.txt file to find out how'
        print*,' the program got before encountering the error'
        print*
c        stop ' Program aborted in rdpar.for @ ~716'
      endif

      if(flz(1).le.0.0)then
        write(*,8429)flz
        STOP
      endif

!     POTENTIAL EVAPORATION:
!     DEFAULT VALUES:
!     CHANGED SO IT WOULD COMPILE IN UNIX - FRANK S OCT/97
!     TS: ADDED IMPERV=classcount EVAP VALUE (0 CAUSE NO GROWTH)
      do i=1,classcount       
        do j=1,12
          evap(i,j)=0.0
          if(i.eq.classcount) evap(i,j)=0.0
        end do
      end do

!     READ THE MONTHLY EVAP:
      filename1='basin/evap.dat'
      open(unit=99,file=filename1,status='unknown',err=99901)
      do n=1,classcount
        read(99,3101,end=3201)(evap(n,j),j=1,12)
        write(51,3101)(evap(n,j),j=1,12)
      end do
	close(unit=99,status='keep')

      if(iopt.eq.2)print*, 'param 13'

      GO TO 3202

99901 write(*,5100)
      write(51,5100)

      GO TO 3202

 3201 write(51,3203)

 3202 CONTINUE

c      close(unit=99,status='keep')

      if(iopt.eq.2)print*, 'param 13'

      firstpass='n'
      answer='n'

      if(ver.ge.9.5)RETURN

      print*
      print*
      print*,' Version (ver) =',ver
      print*,' An old parameter file format was found and a new'
      print*,' format file called basin/new.par will be created '
      
      print*
      print*,'New parameter file format now recommended for use'
      print*,'To avoid this message & annoying pause, rename '
      print*,'basin\new_par.csv   basin\bsnm_par.csv'
      print*,'where `bsnm` is the basin name'
      print*,'and change the par file name in the event file'
      print*
      print*,' However, the new.par file will NOT be compatible'
      print*,' with splx versions older than  9.6.04 '
      print*,' or dds connectors not written by NK'
      print*
      if(dds_flag.eq.0)then
        pause 'hit enter to continue with the old par file format'
      endif
 
      
      
      answer='y'

!     REPLACE AN OLD VERSION PAR FILE WITH THE NEW FORMAT

      if(snwflg.eq.'y') call rdsdc()


!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      filename(99)='basin\new_par.csv'
      call write_par_10(99,99)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if(iopt.eq.2)print*, 'param before return @ 9989'


! 9989 RETURN

      firstpass='n'
      RETURN

 9999 CONTINUE

      print*,' stopped in rdpar due to data problems - see spl.txt'
      print*,' spl.txt will show how much data was properly read'
      STOP 'program stopped in rdpar.for'

! FORMATS

 1000 format(a5,g10.3,2i5,a65)
 1001 format(a5,99e10.3)
 1002 format(a5,4e10.3)
 1003 format(a5,99f5.2)
 1004 format(a5,99e10.3)
 1005 format(a5,99e10.3)
 1006 format(a5,99e10.3)
 1007 format(a5,99e10.3)
 1008 format(a5,99e10.3)
 1009 format('     e       ix')
 1010 format(' maximun interception storage in mm')
 1011 format(' depression storage ds(i)')
 1012 format(' permeabilities ak(i)')
 1013 format(' roughness r3(i)  -  pervious area')
 1014 format(' velocity factor for channels/square   chnl(i)')
 1015 format(' roughness r4(i)  -  impervious area')
 1019 format(' a1 ... a12')
 1022 format(' interflow recession constant rec(i)')
 1023 format(' ','ti            delta        low         high   param') 
 1025 format(/' ','flgevp2 set to',f5.1,' in rdpar.for'/)

 1031 format(' lower zone discharge function flz(i)')
 1032 format(' lower zone power pwr(i)')
 1033 format(' flood plain roughness multiplier r1(i)')
 1034 format(' river roughness r2(i)')
 1035 format(' meander length multiplier mndr(i)')
 1036 format(' bankfull constant aa2(i)')
 1037 format(' bankfull constant aa3(i)')
 1038 format(' bankfull constant aa4(i)')
 1039 format(' wetland porosity theta(i)')
 1040 format(' wetland width/depth ratio widep(i)')
 1041 format(' wetland conductivity kcond(i)')

 3000 format(a5,6i5,25x,f10.0) 
 3010 format(a1,a79)
 3101 format(12f5.0)
 3203 format(' Warning: basin/evap.dat table incomplete'/
     *         '          zero values are inserted for evap.dat'/)
 5000 format(a5,2i5,a75)
c 5001 format(a5,17f10.3)
c 5002 format(17f10.3)
 5003 format(a5,4e10.3)
 5004 format(a5,6f10.3)
 5005 format(a80)
 5006 format(' warning: default values for retn, ak2, flz & pwr'/
     *'          are used. par file is out of date'/)
 5007 format(' warning: errflg  =',i5,' reset to 0 in rdpar')
 5009 format(5x,99a10)
 5010 format(a1)
        write(51,5111)a10

 5100 format(' in rdpar - problem opening basin/evap.dat file'/
     *         '          zero values are inserted for evap.dat'/)
5110  format(//,' Value for a9 heat deficit to swe outside range',/
     *     f10.3,' assumed')
5111  format(//,' Value for a10 uz exponent outside range',/
     *     f10.3,' assumed')
 5999 format(' mflag=',i5,'   errflg=',i5)
c 6003 format(a5,17e10.3)

 6011 format('# runtime    ',2(a2,':'),a2)
 6012 format('# rundate  ',a4,'-',a2,'-',a2)


 8429 format(//' in runof5,flz/',5f10.6/
     * ' wrong parameter value - change flz(i) to +ve value'/)
 9801 format(a5,99f10.0)
 9802 format(a5,99f10.2)
 9803 format(a5,12f5.1)
 9804 format(a5,f10.2,a60)
 9805 format(a5,i10,5x,a40)
 9806 format(a5,f10.0,5x,a40)
 9807 format(a5,f10.3,5x,a40)
 9808 format(//' Parameters used for this run:'/)
 9809 format('A new par file called basin\new.par has been created')

      firstpass='n'

      RETURN

      END SUBROUTINE rdpar

