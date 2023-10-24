      SUBROUTINE read_evt(date,conv,scale,smc5,nhg,nhf)

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

! THIS SUBROUTINE INPUTS EVENT PARTICULARS.


!     REV. 8.71 - Feb.  24/98 -    added flgevp2 to rdevt.for
!     REV. 8.78 - July   7/98 -    added scalesnw and scaletem to rdevt
!     REV. 8.82 - July  10/98 -    added runoff output option: routeflg
!     REV. 8.94 - Feb.  01/99 -    crseflg to read resume and snow course
!     REV. 8.96.1 May   12/99 -    added ireport for reporting interval
!     REV. 9.00   March  2000 -    TS: CONvert version TO FORTRAN 90      
!     rev. 9.1.34  Dec.  23/02  - Added ensim1flg - if ensimflg='a' for 1st id then 'y' for all events
!     rev. 9.1.54  aPR.  12/04  - NK: SEDFLG set for multiple events at event No. 1
!     rev. 9.1.77  Mar.  07/05  - NK: added .psm .gsm and .glz  files
!     rev. 9.1.78  Mar.  15/05  - NK: added WQD file to event file
!     rev. 9.3.11  Feb.  28/07  - NK: ch_par added / event file ver = 9.5
!     rev. 9.3.11  ???.  ??/07  - NK: initflg added / event file ver = 9.6
!     rev. 9.3.11  ???.  ??/07  - NK: deltat_report added / event file ver = 9.7
!     rev. 9.5.08  Feb.  08/08  - NK: new event parser
!     rev. 9.5.22  Mar.  12/08  - NK: added grdflg to print gridded flow, swe & evap
!     rev. 9.5.26  Apr.  04/08  - NK: added Julian day calc. to read_evt
!     rev. 9.5.48  Dec.  26/08  - NK: added event_fln() to allow unlimited events
!     rev. 9.5.57  Apr.  13/09  - NK: added ntrlflg for natural lake flows
!     rev. 9.5.58  Apr.  16/09  - NK: added nudgeflg for forcing gauge flows
!     rev. 9.5.60  Sep.  01/09  - NK: added deltat_report for lake_sd.csv file
!     rev. 9.5.79  Nov.  04/09  - NK: added resumflg='s' for read_soilinit ONLY
!     rev. 9.6.04  Apr.  05/10  - NK: fixed filename carry over in read_evt
!     rev. 9.8.22  Jul.  17/12  - NK: Added resetflg to reset cumm. precip Sept.1
!     rev. 9.8.23  Aug.  03/12  - NK: Added resinid1flg to use resinflg for id=1
!     rev. 9.8.24  Aug.  07/12  - NK: Added reading yyyymmdd_lvl.tb0 for lake levels
!     rev. 9.8.25  Sep.  26/12  - NK: Added warning for resumflg=y and ID > 1
!     rev. 9.8.47  Feb.  04/13  - NK: Headers added for spl & resin csv files
!     rev. 9.8.53  Mar.  20/13  - NK: Add Lake St. Joseph diversion algorithm to REROUT.f
!     rev. 9.8.81  Sep.  03/13  - NK: Add pafflg and update precip adjustment factors PAF!***********************************************************************
!     rev. 9.9.00  Dec.  08/13  - NK: Added Lake Evaporation model
!     rev. 9.9.28  Sep.  18/14  - NK: Added 'a' as option for ntrlflg & smrflg
!     rev. 9.9.47  Dec.  23/14  - NK: Added lakeflg for lake evaporation option
!     REV. 10.1.41 Oct   11/16  - NK: Added tb0flg to write lake_*.tb0 files
!     rev. 10.1.64 Jan.  26/17  - NK: Added XML output file 
!     rev. 10.2.28 Jul.  08/18  - NK: Revised for TSW = NBS without lake precip/evap'n

!     rev mar05/07

! - LIST OF ARGUMENTS:
 
!   o - date    char*14  event identifier
!   o - conv    real*4   conv version factor for precipitation inputs
!   o - smc     real*4   initial soil moisture content
!   o - fln( )  char*12  file names

! - LIST OF INPUT FILES:

! unit=31  fln(1) - basin file (bsnm_shd.r2c)
! unit=32  fln(2) - parameter file (bsnm.par)
! unit=33  fln(3) - point data location file (bsnm.pdl)
! unit=34  fln(4) - not used in spl
! unit=35  fln(5) - point precipitation file  yyyymmdd.rag
! unit=36  fln(6) - streamflow data (yyyymmdd.str)
! unit=37  fln(7) - reservoir release data (yyyymmdd.rel)
! unit=38  fln(8) - reservoir inflow data  yyyymmdd.rin
! unit=39  fln(9) - unadjusted radar file  yyyymmdd.rad
! unit=40  fln(10)- precipitation data (yymmdd.met)
! unit=41  fln(11)- radar scan (cappi) file
! unit=42  fln(12)- radar clutter file    yyyymmdd.clt
! unit=43  fln(13)- snow cover depletion curve   bsnmsdc
! unit=44  fln(14)- point temperatures yymmdd.tag
! unit=45  fln(15)- gridded temperatures yymmdd.tem
! unit=46  fln(16)- max temperatures yymmdd.tmx
! unit=47  fln(17)- min temperatures yymmdd.tmn
! unit=48  fln(18)- saily sbow files yymmdd.dsn
! unit=49  fln(19)- radiation data gridded yymmdd.flx
! unit=40  fln(20)- radiation data point yymmdd.prn
! units 51-99 reserved for output files
! filenames 100-199 reserved for event file names
! unit=251  fln(21)- humidity yymmdd.grh
! unit=252  fln(22)- 
! unit=253  fln(23)- longwave radiation yymmdd.glw
! unit=254  fln(24)- shortwave radiatin yymmdd.gsw
! unit=255  fln(25)- atmospheric pressure yymmdd.gpr
! unit=256  fln(26)- point relative humidity yyyymmdd.prh
! unit=257  fln(27)- 
! unit=258  fln(28)- point longwave radiation yyyymmdd.plw
! unit=259  fln(29)- point shortwave radiation yyyymmdd.psw
! unit=260  fln(30)- point atmospheric pressure yyyymmdd.ppr
! unit=261  fln(31) - gridded runoff files  yyyymmdd.rff
! unit=262  fln(32) - gridded recharge    yyyymmdd.rch
! unit=263  fln(33) - gridded leakage     yyyymmdd.lkg
! unit=264  fln(34)- gridded lake evaporation
! unit=265  fln(35)- snow course data file  yyyymmdd.crs 
! unit=266  fln(36)- gridded snow water equivalant  yyyymmdd.swe
! unit=267  fln(37)- gridded soil moisture         yyyymmdd.gsm
! unit=268  fln(38)- gridded lower zone storage    yyyymmdd.lzs
! unit=269  fln(39)- point soil moisture  .psm 
! unit=270  fln(40)- water quality data file  .wqd
! unit=271  fln(41)- routing parameter file  bsnm_ch_par.r2c
! unit=272  fln(42)- gridded wetland surface flux    - r2c file
! unit=273  fln(43)- gridded water surface flux file - 2rc file
! unit=274  fln(44)- point snow precip              yyyymmdd_snw.tbO
! unit=275  fln(45)- point drain (O18)              yyyymmdd_drn.tbO
! unit=276  fln(46)- point dsnow (O18)              yyyymmdd_dsn.tbO
! unit=277  fln(47)- gridded snow precip            yyyymmdd_snw.r2c
! unit=278  fln(48)- gridded drain (O18)            yyyymmdd_drn.r2c
! unit=279  fln(49)- gridded dsnow (O18)            yyyymmdd_dsn.r2c
! unit=280  fln(50)- point initial lake conditions  yyyymmdd_ill.pt2
! unit=281  fln(51)- gridded evaporation            yyyymmdd_evp.r2c
! unit=282  fln(52)- point diversion flow file      yyyymmdd_div.tb0
! unit=283  fln(53)- recorded lake level file       yyyymmdd_lvl.tb0
! unit=284  fln(54)- swe time series pillows & crs  yyyymmdd_swe.tb0
!     rev. 9.9.00  Dec.  08/13  - NK: Added Lake Evaporation model
! unit=285  fln(55)- point wind speed               yyyymmdd_spd.tb0
! unit=286  fln(56)- point wind direction           yyyymmdd_dir.tb0
! unit=287  fln(57)- gridded wind speed             yyyymmdd_spd.r2c
! unit=288  fln(58)- gridded wind direction         yyyymmdd_dir.r2c
! filenames 100-199 reserved for event file names



!   o - nhg     int      number of hours of rain gauge data
!   o - nhf     int      number of hours of flow data


! THE FOLLOWING FLAGS IF .EQ. 'y' WILL:
!  1 snwflg   - Use MELT ROUTINES
!  2 sedflg   - Use SEDIMENT ROUTINES
!  3 vapflg   - Use EVAPORATION ROUTINES
!  4 smrflg   - SMEAR PRECIP DATA
!  5 resinflg - READ THE yyyymmdd.rin FILE AND OUTPUT RESIN.TXT FILE
!  6 tbcflg   - WRITE A \spl\bsname\resume.txt FILE AT THE END OF RUN
!  7 resumflg - READ THE \spl\bsname\resume.txt FILE AT START OF RUN
!  8 contflg  - CONT THE STATISTICS FROM PREVIOUS RUN VIA RESUME.TXT
!  9 routeflg - WRITE THE \spl\bsnm\runof\yymmdd.txt FOR WATROUTE.EXE
!               write \spl\bsnm\rchrg\yymmdd.txt for MODFLOW
!               and write \spl\bsnm\lkage\yymmdd.lkg for watroute
! 10 crseflg  - READ SNOW COURSE DATA TO REPLACE RESUME FILE DATA
! 11 ensimflg - write the wfo file for ENSIM
! .  ensim1flg- if = 'a' wfo file will be written for 'all' events
! 12 picflg   - write the simout/pic.txt file for mapper
! 13 wetid1flg- run the wetland routing module - read in first event only
! 14 modelflg - if ='i' watroute with surface flow only 
!               if ='l' watroute for surface and groundwater routing
!               if ='r' watroote for surface to channel and drainage thru lz
! 15 shdflg   - replace the watershed file basin\bsnm.shd for next event
! 16 trcflg   - use the tracer module
! 17 frcflg   - use isotope fractionation
! 18 initflg
! 19 grdflg   - write a gridded outflow
! 20 ntrlflg
! 21 nudgeflg
! 22 resetflg
! 23 hdrflg
! 24 divertflg - y = use diversion data (default) g = generate diversion flow
! 25 pafflg    - generate & use precip adjustment factors (PAF)

!
!***********************************************************************

      USE area_watflood
      use area_debug
      use EF_Module
      implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      type(EventParam) :: event
      type(FrameRecord) :: frameRec
      character*4096 line, subString
      character*128 keyword, value ,keyword_no_evt
      integer lineLen, keyLen, wordCount
      logical rStat, foundEndevent,foundFirstFrame,foundFrame
      logical       :: file58,file59,file60   !  used to check the existance of wind files
      integer*4 unitNum, flnNum, iStat
      CHARACTER(14) :: date
      character(1)  ::   wetid1flg,ensim1flg,sedid1flg,
     *                   routeid1flg,modelid1flg,trcid1flg,
     *                   snwid1flg,frcid1flg,resinid1flg,
     *                   divertid1flg,pafid1flg,grd1flg,
     *                   ntrlid1flg,lakeid1flg,smrid1flg,
     *                   iceid1flg,icerivid1flg,icelakeid1flg,
     *                   tb0id1flg,xmlid1flg,netCDFoutID1flg
      character(30) :: junk
      real*4   :: smc5(16),conv,scale
c      real(4)  :: scaleallsnw,readscalesnw,readscale
      real*4   :: readscaletemp
      integer  :: nhg,nhf,i,ios,n,nkeep,iAllocate,deltat_report_id1
      logical  :: exists,firstpass

!     rev. 9.5.26  Apr.  04/08  - NK: added Julian day calc. to read_evt
      INTEGER        :: mohours(24),ju_mon(24),dayrad(12),daysum(12)
      DATA mohours/744,672,744,720,744,720,744,744,720,744,720,744,
     *             744,672,744,720,744,720,744,744,720,744,720,744/
      DATA ju_mon/1, 32, 60, 91,121,152,182,213,244,274,305,335,
     *          366,397,425,456,486,517,547,578,609,639,670,700/
      data dayrad/31,28,31,30,31,30,31,31,30,31,30,31/
      data daysum/31,59,90,120,151,181,212,243,273,304,334,365/

      data firstpass/.true./

!     set defaults defaults defaults defaults defaults defaults defaults defaults
!     set defaults defaults defaults defaults defaults defaults defaults defaults
!     set defaults defaults defaults defaults defaults defaults defaults defaults
!     if these flags are not in subsequent event files
!     rev. 9.8.22  Jul.  17/12  - NK: Added resetflg to reset cumm. precip Sept.1
!     rev. 10.4.42 Sep.  29/21  = NK Added flowresetflg to the event file
      if(firstpass)then     ! set defaults
	  grdflg='n'
	  ntrlflg='n'
	  nudgeflg='n'
	  resetflg=.false.
	  hdrflg0='n'
	  idskip=0      ! default: zero events skipped for spinup
c	  divertflg='y' ! default - use diversion files if present
	  divertflg='n' ! default - must be set to 'y' in the event file
	  pafflg='n'    ! default - do not use or generate PAF's
	  iceflg='n'
        event%iceflg='n'
	  icerivflg='n'
        event%icerivflg='n'
	  icelakeflg='n'
        event%icelakeflg='n'
!       REV. 10.1.41 Oct   11/16  - NK: Added tb0flg to write lake_*.tb0 files
!       This is the default.  Add tb0flg = 'y' in the event file to activate this option
        tb0flg='n'
        event%tb0flg='n'
        xmlflg='n'
        event%xmlflg='n'
        nbsflg='n'
        event%nbsflg='n'
        fewsflg='n'
        event%fewsflg='n'
        netCDFoutflg='n'
        event%netCDFoutflg='n'
	endif

	fliflg='n'    ! defaults
      event%fliflg='n'
      lakeflg='n'
      event%lakeflg='n'
      flowresetflg='n'   ! must be explicitly set in each event file
      
!     TW  MAR. 3/98 - GIVES ERROR LINKING ON SUN, SO COMMENT OUT
!      DATA flgevp2/-1.0/

d       if(iopt.eq.2)print*,' In rdevt, passed location  300'
c      if(debug_output)write(51,1300)fln(99)
d       if(iopt.eq.2)print*,' opening file name ',fln(99)
!     OPEN THE EVENT FILE:

      unitNum=99

!     Initialize parameters
      CALL InitEventParam(event)

      INQUIRE(FILE=fln(99),EXIST=exists)

      IF(exists)THEN
        open(unit=99,file=fln(99),iostat=ios)
d        if(iopt.ge.1)print*,'Opened event file ',fln(99)(1:40)
d        if(iopt.eq.2)print*,' In rdevt, passed location  30001'
        if(ios.ne.0)then
          write(63,99921)fln(99)
          if(numa.le.0)write(51,99921)fln(99)(1:40)
99921     format(' error opening unit 99, file name= ',a30)
          print*, 'iostat= ',ios
          print*  
          if(dds_flag)then
            print*,'Continuing will not calculate the error function'  
            pause ' program paused in rdevt'
          else
            STOP ' program terminated in rdevt'
          endif
        endif
      else
        print*,'Attempting to open the file ',fln(99)(1:40)
        print*,'but it is not found (in this location)'
	  print*,'in read_evt @ 188  unit=',unitNum
	  print*,'Possibel cause: Not in proper working directory'
        print*
        if(dds_flag)then
          print*,'Continuing will not calculate the error function'  
          pause ' program paused in rdevt'
        else
          STOP ' program terminated in rdevt'
        endif
      endif

d      if(iopt.eq.2)print*,' In rdevt, passed location  301'

!     rev. 9.6.04  Apr.  05/10  - NK: fixed filename carry over in read_evt
!     ensure filenames do not get carried over into subsequent files
!     if there are blanks in the event file 
      do i=1,60
        event%infln(i)='dummy'
      end do
      event%writestr='n'
      writeflg(2)   =.false.
      
!     INPUT EVENT PARTICULARS:

      do while(.not.ios.lt.0)

        read(unit=unitNum, FMT='((A))', iostat=ios) line      ! read a line

d        write(63,99991)line
99991	   format(a60)

c                if(ios .eq. -1)then
c        print*, 'ERROR: Premature EndOfFile encountered'
c                      STOP ' Stopped in read_rain_ef'
c                end if


        rStat = Detab(line)                       ! replace tabs with spaces
        line = adjustl(line)          ! Get rid of leading white space
        lineLen = len_trim(line)            ! Find the length excluding trailing spaces
        if(line(1:1) .eq. ':')then
          wordCount = SplitLine(line, keyword, subString) ! find the keyword
          rStat = ToLowerCase(keyword)
          KeyLen = len_trim(keyword)


!         get numerical values
          if(keyword(1:KeyLen) .eq. ':year')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':month')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':day')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':hour')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':initsoilmoisture')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':rainconvfactor')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':eventprecipscalefactor')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':precipscalefactor')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
c          elseif(keyword(1:KeyLen) .eq. ':eventscalesnowfactor')then
          elseif(keyword(1:KeyLen) .eq. ':eventsnowscalefactor')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':snowscalefactor')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':eventtempscalefactor')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':tempscalefactor')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':disaggregate')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':hoursraindata')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':hoursflowdata')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':deltat_report')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':spinupevents')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)


!         get the flags
          elseif(keyword(1:KeyLen) .eq. ':snwflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':sedflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':vapflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':smrflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':resinflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':tbcflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':resumflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':contflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':routeflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':crseflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':ensimflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':picflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':wetflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':modelflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':shdflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':trcflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':frcflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':initflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':grdflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':ntrlflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':nudgeflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':flowresetflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':resetflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
!     rev. 9.8.47  Feb.  04/13  - NK: Headers added for spl & resin csv files
          elseif(keyword(1:KeyLen) .eq. ':hdrflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
!     rev. 9.8.53  Mar.  20/13  - NK: Add Lake St. Joseph diversion algorithm to REROUT.f
          elseif(keyword(1:KeyLen) .eq. ':divertflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
!     rev. 9.8.81  Sep.  03/13  - NK: Add pafflg and update precip adjustment factors PAF!***********************************************************************
          elseif(keyword(1:KeyLen) .eq. ':pafflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':fliflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
!     rev. 9.9.47  Dec.  23/14  - NK: Added lakeflg for lake evaporation option
          elseif(keyword(1:KeyLen) .eq. ':lakeflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':iceflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':icerivflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':icelakeflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':tb0flg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':xmlflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':nbsflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':fewsflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':netCDFoutflg')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':writestr')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)


!         get the file names
          elseif(keyword(1:KeyLen) .eq. ':basinfilename')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':parfilename')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':pointdatalocations')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':pointprecip')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':streamflowdatafile')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':reservoirreleasefile')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':reservoirinflowfile')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':radarfile')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddedrainfile')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':rawradarfile')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':clutterfile')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen).eq.':snowcoverdepletioncurve')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':pointtemps')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddedtemperaturefile')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. '::griddeddailydifference')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddednetradiation')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':pointnetradiation')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddedhumidity')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddedwind   ')then    ! = old key word
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddedwindspd')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddedwinddir')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddedlongwave')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddedshortwave')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddedatmpressure')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':pointhumidity')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':pointwind   ')then    ! = old key word
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':pointwindspd')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':pointwinddir')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':pointlongwave')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':pointshortwave')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':pointatmpressure')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddedrunoff')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddedrecharge')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddedleakage')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddedsnowfile')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':snowcoursefile')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddedinitsnowweq')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen).eq.':griddedinitsoilmoisture')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddedinitlzs')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':pointsoilmoisture')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':waterqualitydatafile')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':climateNormalsDiff')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':channelparfile')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':pointsnow')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':pointdrain')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':pointdsnow')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddedsnow')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddeddrain')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddeddsnow')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddedevaporation')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':initlakelevel')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':diversionflowfile')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':observedlakelevel')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':observedswe')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
              
!     rev. 10.2.11 Dec.  18/17  - NK: 4 files added for BLEND.exe    
          elseif(keyword(1:KeyLen) .eq. ':pointHoulyPrecip')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':pointDailyPrecip')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddedHoulyPrecip')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':griddedDailyPrecip')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
              
              
              
          elseif(keyword(1:KeyLen) .eq. ':flowinit')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
          elseif(keyword(1:KeyLen) .eq. ':filetype')then
              continue !do nothing - no problem
          elseif(keyword(1:KeyLen) .eq. ':fileversionno')then
              continue !do nothing - no problem
          elseif(keyword(1:KeyLen) .eq. ':noeventstofollow')then
              iStat = ParseEventParam(event,keyword,keyLen,subString)
!             Once we read this line, we're done reading the filenames.
!             Go ahead and read the chained events.
              read(99,*)line    ! read the #  line
d             write(63,99991)line
              nch=event%nch
d              print*,'nch=',nch,event%nch
!               ni = no of events
!               default is 1 but if nch > 0 then the total no of events
!               will be nch+1

!             but wait a minute
!     rev. 9.8.26  Sep.  26/12  - NK: Added error check on # chained files for id>1'
              if(id.gt.1.and.nch.gt.0)then
                error_msg=1
                nch=0
              endif

              if(nch.ge.1)then
	          ni=nch+1    ! number of events in this run
!     rev. 9.5.48  Dec.  26/08  - NK: added event_fln() to allow unlimited events
	          if(.NOT.allocated(event_fln))then 
                  allocate(event_fln(ni),stat=iAllocate)
                  if (iAllocate.ne.0) STOP 
     *    'Warning: error with allocation of event_fln in read_evt' 
	          endif


!               READ SUBSEQUENT EVENT NAMES TO MODELLED:
!               READ SUBSEQUENT EVENT NAMES TO MODELLED:
                do n=2,ni
                  nkeep=nch
                  read(99,1300,iostat=ios)event_fln(n)
                  if(iopt.gt.1)write(51,1300)event_fln(n)
d                  if(iopt.ge.1)write(63,1300,iostat=ios)event_fln(n)
                  if(numa.le.0)write(51,1300)event_fln(n)
                  if(ios.ne.0)then
                    print*,'Error reading file no ',nkeep,' in rdevt'
                    print*,'possible cause: number of files listed '
                    print*,'in the event file is shorter'
                    print*,'than the specified ',ni-1
                    print*
                    stop 'Program abortedin rdevt @ 460'
                  endif
                  if(iopt.eq.2)print*,n,event_fln(n)
                end do
              endif
         
d              if(iopt.ge.1)print*,'Reached end of ',fln(99)(1:40)
          else
              iStat = ParseEventParam(event,keyword,keyLen,subString)
              if(iStat .lt. 0) then
                   write(98,'(2(A))') 'ERROR parsing ', fln(99)
                   write(98,'(2(A))') '   in line: ',line                             
                   write(98,*)'Error: Stopped in read_evt @ 577'
                   write(*,'(2(A))') 'ERROR parsing ', fln(99)
                   write(*,'(2(A))') '   in line: ',line                             
                   STOP ' Stopped in read_evt @ 577'
                   return
              else if(iStat .eq. 0) then
                  write(98,'((A), (A))')
     *           'Info: Unrecognized keyword line in evt file: '
                  write(98,*)'Info: line'
              endif
          end if
        end if
	end do  ! (.not.ios.lt.0)

!     rev. 9.8.47  Feb.  04/13  - NK: Headers added for spl & resin csv files
      write(line,10000)event%year1,'/',event%mo1,'/',event%day1
      if(iopt.ge.1)then    !  NK Mar. 20/13
        write(63,10000)event%year1,'/',event%mo1,'/',event%day1
10000   format(i4,a1,i2,a1,i2)
      endif
      read(line,10001)start_date
10001 format(a10)            

      if(iopt.ge.1)then
        if(debug_output)write(51,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        if(debug_output)write(51,*)'data as read from the event file:'
        if(debug_output)write(51,*)event%year1
        if(debug_output)write(51,*)event%mo1
        if(debug_output)write(51,*)event%day1
        if(debug_output)write(51,*)event%hour1
        if(debug_output)write(51,*)event%smc_init
        if(debug_output)write(51,*)event%conv
        if(debug_output)write(51,*)event%scale
        if(debug_output)write(51,*)event%readscale
        if(debug_output)write(51,*)event%scalesnw
        if(debug_output)write(51,*)event%readscalesnw
        if(debug_output)write(51,*)event%scaletem
        if(debug_output)write(51,*)event%readscaletem
        if(debug_output)write(51,*)event%nhg
        if(debug_output)write(51,*)event%nhf
        if(debug_output)write(51,*)event%deltat_report
        if(debug_output)write(51,*)event%nch
!     flags as read
        if(debug_output)write(51,*)'flags as read from the event file:'
        if(debug_output)write(51,*)'event%snwflg   ',event%snwflg  
        if(debug_output)write(51,*)'event%sedflg   ',event%sedflg  
        if(debug_output)write(51,*)'event%vapflg   ',event%vapflg  
        if(debug_output)write(51,*)'event%smrflg   ',event%smrflg  
        if(debug_output)write(51,*)'event%resinflg ',event%resinflg
        if(debug_output)write(51,*)'event%tbcflg   ',event%tbcflg  
        if(debug_output)write(51,*)'event%resumflg ',event%resumflg
        if(debug_output)write(51,*)'event%contflg  ',event%contflg 
        if(debug_output)write(51,*)'event%routeflg ',event%routeflg
        if(debug_output)write(51,*)'event%crseflg  ',event%crseflg 
        if(debug_output)write(51,*)'event%ensimflg ',event%ensimflg
        if(debug_output)write(51,*)'event%picflg   ',event%picflg  
        if(debug_output)write(51,*)'event%wetflg   ',event%wetflg  
        if(debug_output)write(51,*)'event%modelflg ',event%modelflg
        if(debug_output)write(51,*)'event%shdflg   ',event%shdflg  
        if(debug_output)write(51,*)'event%trcflg   ',event%trcflg  
        if(debug_output)write(51,*)'event%frcflg   ',event%frcflg  
        if(debug_output)write(51,*)'event%initflg  ',event%initflg 
        if(debug_output)write(51,*)'event%grdflg   ',event%grdflg  
        if(debug_output)write(51,*)'event%ntrlflg  ',event%ntrlflg 
        if(debug_output)write(51,*)'event%nudgeflg ',event%nudgeflg
        if(debug_output)
     *   write(51,*)'event%flowresetflg ',event%flowresetflg
        if(debug_output)write(51,*)'event%resetflg ',event%resetflg 
        if(debug_output)write(51,*)'event%hdrflg   ',event%hdrflg0 
        if(debug_output)write(51,*)'event%divertflg',event%divertflg 
        if(debug_output)write(51,*)'event%pafflg   ',event%pafflg 
        if(debug_output)write(51,*)'event%fliflg   ',event%fliflg 
        if(debug_output)write(51,*)'event%lakeflg   ',event%lakeflg 
        if(debug_output)write(51,*)'event%iceflg   ',event%iceflg 
        if(debug_output)write(51,*)'event%icerivflg   ',event%icerivflg 
        if(debug_output)
     *   write(51,*)'event%icelakeflg   ',event%icelakeflg 
        if(debug_output)write(51,*)'event%tb0flg   ',event%tb0flg 
        if(debug_output)write(51,*)'event%xmlflg   ',event%xmlflg 
        if(debug_output)write(51,*)'event%nbsflg   ',event%nbsflg 
        if(debug_output)write(51,*)'event%fewsflg   ',event%fewsflg 
!       netCDFoutflg   not written
      endif

!     NOTE:
!     rev. 10.1.93 Aug   17/17  - NK: allow year1 etc. to be passed for each event
!     For CHARM, for it's first pass year1 = set to year now and then timer keeps track of the years
!     For RAGMET & TMP, year1, mo1, day1 & hour1 are the values for the first day of each event
!                       and are neede to start the time stamping for the output files.
c      if(firstpass)then
          year1         =event%year1
          mo1           =event%mo1
          day1          =event%day1
          hour1         =event%hour1
	    jul_day1      =ju_mon(mo1)+day1-1
	    if(mod(year1,4).eq.0.and.mo1.ge.3)then
	    jul_day1=jul_day1+1
          endif
c      endif
          
          
      do i=1,16
        smc5(i)      =event%smc_init
	end do

      conv          =event%conv         !rain conversion factor
      scale         =event%scale        !event precip conversion factor
      readscale     =event%readscale    !precip scale factor - all events
      scalesnw      =event%scalesnw     !event snow scale factor - seems to be used i snow init
      readscalesnw  =event%readscalesnw !snow scale factor - all events
      scaletem      =event%scaletem     !event temp scale factor
      readscaletemp  =event%readscaletem!temp scale factor - all events
      smearfactor   =event%disaggregate !minimum disaggreation amount
      nhg           =event%nhg
      nhf           =event%nhf
      deltat_report =event%deltat_report
      nch           =event%nch
      snwflg        =event%snwflg
      sedflg        =event%sedflg
      vapflg        =event%vapflg
      smrflg        =event%smrflg
      resinflg      =event%resinflg
      tbcflg        =event%tbcflg
      resumflg      =event%resumflg
      contflg       =event%contflg
      routeflg      =event%routeflg
      crseflg       =event%crseflg
      ensimflg      =event%ensimflg
      picflg        =event%picflg
      wetflg        =event%wetflg
      modelflg      =event%modelflg
      shdflg        =event%shdflg
      trcflg        =event%trcflg
      frcflg        =event%frcflg
      initflg       =event%initflg
      grdflg        =event%grdflg
      ntrlflg       =event%ntrlflg
      nudgeflg      =event%nudgeflg
      flowresetflg  =event%flowresetflg
      divertflg     =event%divertflg
      pafflg        =event%pafflg 
      fliflg        =event%fliflg 
      lakeflg       =event%lakeflg 
      iceflg        =event%iceflg
      icerivflg     =event%icerivflg
      icelakeflg    =event%icelakeflg
      tb0flg        =event%tb0flg
      xmlflg        =event%xmlflg
      nbsflg        =event%nbsflg
      fewsflg       =event%fewsflg
      netCDFoutflg  =event%netCDFoutflg
      if(event%writestr.eq.'y')writeflg(2)   =.true.    ! set = .false in CHARM
!     special cases:
      if(event%divertflg.eq.'g')divertflg=event%divertflg
!     rev. 9.8.22  Jul.  17/12  - NK: Added resetflg to reset cumm. precip Sept.1
      if(firstpass)then
       if(event%resetflg.eq.'y'.or.event%resetflg.eq.'Y')resetflg=.true.
       if(event%hdrflg0.eq.'y'.or.event%hdrflg0.eq.'Y')hdrflg0='y'
      endif
                     
      if(iopt.ge.1)then   
        if(debug_output)write(51,*)'    file #    file name (fln(#)'
        do i=1,70
	      if(debug_output)write(51,*)i,event%infln(i)
c	      write(63,*)i,event%infln(i)
	  end do 
      endif      

!     INPUT FILES
!     INPUT FILES
!     INPUT FILES
!     INPUT FILES
!     assign the names
!     rev. 10.2.46 Jan.  21/19  - NK: Added netCDF capability

      if(id.eq.1)then      
	  fln(1)=event%infln(1)
        call find_filetype(1)
        FLtype(1)=filetype
      else
        if(shdflg.eq.'y')then
	    fln(1)=event%infln(1)
          call find_filetype(1)
          FLtype(1)=filetype
c	    print*,i,event%infln(1)
        endif
      endif
      
      do i=2,60
	  fln(i)=event%infln(i)
        call find_filetype(i)
        FLtype(i)=filetype
c	  print*,i,event%infln(i)
      end do       
!     Note: fln(61) is in use for swe                  
      do i=62,65
	  fln(i)=event%infln(i)
        call find_filetype(i)
        FLtype(i)=filetype
c	  print*,i,event%infln(i)
	end do       
	fln(68)=event%infln(68)
        call find_filetype(68)
        FLtype(68)=filetype
c	  print*,68,event%infln(68)

!     rev. 9.9.00  Dec.  08/13  - NK: Added Lake Evaporation model
!       Basic idea is that lakeflg is either on or off throughout the run
!       check to see of wind files are present: 
c      if(lakeflg.eq.'y')then
c!       check to see if wind files exist:                
c        if(program_name(1:3).eq.'spl')then
c          file58=.false.
c          file59=.false.
c          INQUIRE(FILE=fln(58),EXIST=exists)
c          if(exists)file58=.true.
c          INQUIRE(FILE=fln(59),EXIST=exists)
c          if(exists)file59=.true.
c          if(file58.and.file59)then
c            lakeflg='y' ! if true: no wind files
c          else
c            print*,'Fatal error:'
c            print*,'lakeflg = `y` but no wind files found'
c            print*,'Looking for  ',fln(58)(1:40)
c            print*,'Looking for  ',fln(59)(1:40)
c            print*
cc            stop 'Program aborted in read_evt @ 709'
c          endif
c        endif
c      endif

      if(firstpass)then
!       rev. 9.9.04  Dec.  17/13  - NK: Change over to gridded clamate normals fo diff
!       this can only be done in the first event
!       Basic idea is that diffflg is either on or off throughout the run
!     rev. 9.9.06  Jan.  08/14  - NK: Add daily differences to Harfreaves ETHarg.f
        INQUIRE(FILE=fln(60),EXIST=exists)
        if(exists)file60=.true.
        if(file60)diffflg=.true.
      elseif(id.ge.2)then
        INQUIRE(FILE=fln(60),EXIST=exists)
        if(.not.exists.and.file60)then
          print*
          print*,'daily differences found for previous event(s) but'
          print*,'no this one. Not a good idea so program terminated'
          stop 'Program aborted in read_evt @ L727'
        endif
      else
        continue   ! do nothing      
      endif
               
      if(id.eq.1)then
!       REMEMBER WHICH MONTH WE START IN:
        month1=mo1
        shd_fln=fln(1)
	  par_fln=fln(2)
        newevtflg='y'
	  idskip=event%spinupevents  ! use value only from first event
      endif  
              
!     turn these off for opt
      if(numa.ne.0)then
	  trcflg='n'
	  frcflg='n'
      endif

!       conv is used as a conversion factor for the precip data.
!       scale is used separately in each event to scale all the precip in this event
!       scalesnw is used to scale the snow course data in this event
!       scaletem is used to add or subtract the value to/from each temperature
!       readscale is used to scale precip for all events to follow
!       readscalesnw is used to scale all snow data to follow
!       readscaletemp is used to scale all temperature data to follow
!       >>>>> the readscale... OVERRIDES the scale... variables !!!!!!!!

!       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!         read the file names for this event:      
!       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!       Note: the file names can be numbered from 1 to nn
!       but the unit numbers can not. 
!       The unit numbers for the input files go from 31 to 50
!       and from 251 to nn
!       This leaves room for more output files from unit=100-199
!       and leaves room for more input files above unit=251
!       Also, in the future, all input unit numbers should be changed to 
!                                     unit=251-300   nk Mar. 11/04
      if(numa.le.0)write(51,1031)
!       In the new revised file format, the order of the files has been
!       changed to group the files by file type:
!       e.g. watershed data files, point data files, gridded files, etc. 

!       Note: discontinuity in the unit numbers.31->50 and 251->infinity
!       Assign unit numbers starting at 250
                          
d      if(iopt.eq.2)print*,' In rdevt, passed location  504'

!       DETERMINE IF IT IS A LEAP YEAR
        leapflg='0'
        if(mod(year1,4).eq.0)leapflg='1'

!       Note: discontinuity in the unit numbers.31->50 and 251->infinity
!       Assign unit numbers starting at 250

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

d      if(iopt.eq.2)print*,' In rdevt, passed location  305'

!     CLOSE INPUT FILE - IT IS USED AS A SCRATCH UNIT # LATER
 999  close(unit=99,status='keep')

!     these flags can only be set in the first event and carry on
!     for all subsequent events - they can not be changed on the fly
!     revised Apr. 27/07  nk
!     rev. 9.5.60  Sep.  01/09  - NK: added deltat_report for lake_sd.csv file
!     rev. 9.8.78  Jul   16/13  - NK: Fixed divertflg to have the first event file value
!     rev. 10.1.07 Dec.  02/15  - NK: Added ice_fctr(n) to route 
!     REV. 10.1.41 Oct   11/16  - NK: Added tb0flg to write lake_*.tb0 files
!     rev. 10.1.44 Dec.  02/15  - NK: Reworked icerivflg & icelakeflg
      if(id.eq.1)then
        snwid1flg=snwflg
        wetid1flg=wetflg    
	  trcid1flg=trcflg
        sedid1flg=sedflg
	  frcid1flg=frcflg
        routeid1flg=routeflg
        modelid1flg=modelflg
        if(.not.nbsflg.eq.'n'.and.deltat_report.ne.24)then
!           e.g.  = 0    or  =1         
            write(98,*)'Warning: deltat_report set to 24'
            write(98,*)'Warning: because nbsflg = y'
            write(*,*)'Warning: deltat_report set to 24'
            write(*,*)'Warning: because nbsflg = y'
            deltat_report=24
        endif
	  deltat_report_id1=deltat_report
	  resinid1flg=resinflg
	  divertid1flg=divertflg
	  pafid1flg=pafflg
	  ntrlid1flg=ntrlflg
	  lakeid1flg=lakeflg
	  iceid1flg=iceflg
	  icerivid1flg=icerivflg
	  icelakeid1flg=icelakeflg
	  tb0id1flg=tb0flg
	  xmlid1flg=xmlflg
        netCDFoutID1flg=netCDFoutflg
      else
	  snwflg=snwid1flg
        wetflg=wetid1flg    
	  trcflg=trcid1flg
        sedflg=sedid1flg
	  frcflg=frcid1flg
        routeflg=routeid1flg
        modelflg=modelid1flg
	  resinflg=resinid1flg
	  deltat_report=deltat_report_id1
	  divertflg=divertid1flg
	  pafflg=pafid1flg
	  lakeflg=lakeid1flg
	  iceflg=iceid1flg
	  icerivflg=icerivid1flg
	  icelakeflg=icelakeid1flg
	  tb0flg=tb0id1flg
	  xmlflg=xmlid1flg
        netCDFoutflg=netCDFoutID1flg
!     rev. 9.9.40  Nov.  19/14  - NK: Modified the 'a' option for ntrlflg
	  if(ntrlid1flg.eq.'a')then
	    ntrlflg='a'   ! overrides value in later events
!         otherwise, ntrlflg can switch between y & n for subsequent events	 
!         by what ever is in the event file - i.e left untouched   
	  endif
      endif

!     write the wfo file for all events if ensimflg='a' for first event
      if(id.eq.1)then
!       For convenience, these flag can be set in the 1st event for the whole run
!       normally, these flag can be changed from one event to another      
!     rev. 9.1.34  Dec.  23/02  - Added ensim1flg - if ensimflg='a' for 1st id then 'y' for all events
!       wfo file written for all events
        if(ensimflg.eq.'a'.or.ensimflg.eq.'A')then
          ensim1flg='a'
c        else                !took out Feb. 09/15
c          ensim1flg=' '
        endif
        if(grdflg.eq.'a'.or.grdflg.eq.'A')then
          grd1flg='a'
c        else                !took out Feb. 09/15
c          grd1flg=' '
        endif
!     rev. 9.9.28  Sep.  18/14  - NK: Added 'a' as option for ntrlflg & smrflg
        if(ntrlflg.eq.'a'.or.ntrlflg.eq.'A')then
          ntrlid1flg='a'
c        else                !took out Feb. 09/15
c          ntrlflg=' '
        endif
!     rev. 9.9.28  Sep.  18/14  - NK: Added 'a' as option for ntrlflg & smrflg
        if(smrflg.eq.'a'.or.smrflg.eq.'A')then
          smrid1flg='a'
c        else                !took out Feb. 09/15              
c          smrflg=' '
        endif
      endif

!     special cases for event #1 flags  = 'a'    
      if(ensim1flg.eq.'a')ensimflg='y'         ! <<<<<  
      if(grd1flg.eq.'a')grdflg='y'           ! <<<<<
      if(ntrlid1flg.eq.'a')ntrlflg='y'          ! <<<<< 
      if(smrid1flg.eq.'a')smrflg='y'          ! <<<<< 
      
      if(snwflg.eq.'Y')snwflg='y'
      if(snwid1flg.eq.'Y')snwid1flg='y'         ! <<<<<
      if(sedid1flg.eq.'Y')sedid1flg='y'         ! <<<<<
      if(vapflg.eq.'Y')vapflg='y'
      if(smrflg.eq.'Y')smrflg='y'
!     rev. 10.3.04 Jan.  29/20  = NK added smrflg = "t" for FEWS only      
      if(smrflg.eq.'T')smrflg='t'
!     rev. 10.4.57 Jan.  29/20  = NK added smrflg = "d" for Disaggredation on CHARM for Ostrich 
!     Only for daily timesteps - 24 hours      
      if(smrflg.eq.'D')smrflg='d'
      if(resinflg.eq.'Y')resinflg='y'
      if(resumflg.eq.'Y')resumflg='y'
!     rev. 9.5.79  Nov.  04/09  - NK: added resumflg='s' for read_soilinit ONLY
      if(resumflg.eq.'S')resumflg='s'  !reads soilinit only
      if(tbcflg.eq.'Y')tbcflg='y'
      if(contflg.eq.'Y')contflg='y'
      if(routeid1flg.eq.'Y')routeid1flg='y'     ! <<<<<
      if(crseflg.eq.'Y')crseflg='y'  
	if(ensimflg.eq.'Y')ensimflg='y'
	if(picflg.eq.'Y')picflg='y'
      if(wetid1flg.eq.'Y')wetid1flg='y'         ! <<<<<
      if(modelid1flg.eq.'R')modelid1flg='r'         ! <<<<<
      if(modelid1flg.eq.'I')modelid1flg='i'         ! <<<<<
      if(modelid1flg.eq.'D')modelid1flg='d'         ! <<<<<
      if(shdflg.eq.'Y')shdflg='y'
	if(trcid1flg.eq.'Y')trcid1flg='y'         ! <<<<<
	if(frcid1flg.eq.'Y')frcid1flg='y'         ! <<<<<
      if(grdflg.eq.'Y')grdflg='y'
      if(ntrlflg.eq.'Y')ntrlflg='y'
      if(nudgeflg.eq.'Y')nudgeflg='y'
      if(flowresetflg.eq.'Y')flowresetflg='y'
      if(pafflg.eq.'Y')pafflg='y'
      if(fliflg.eq.'Y')fliflg='y'
      if(lakeflg.eq.'Y')lakeflg='y'
      if(iceflg.eq.'Y')iceflg='y'
      if(icerivflg.eq.'Y')icerivflg='y'
      if(icelakeflg.eq.'Y')icelakeflg='y'
      if(tb0flg.eq.'Y')tb0flg='y'
      if(xmlflg.eq.'Y')xmlflg='y'
      if(nbsflg.eq.'Y')nbsflg='y'
      if(fewsflg.eq.'Y')fewsflg='y'
      if(netCDFoutID1flg.eq.'Y')netCDFoutflg='y'

        
      if(snwflg.eq.'N')snwflg='n'
      if(snwid1flg.eq.'N')snwid1flg='n'         ! <<<<<
      if(sedid1flg.eq.'N')sedid1flg='n'         ! <<<<<
      if(vapflg.eq.'N')vapflg='n'
      if(smrflg.eq.'N')smrflg='n'
      if(resinflg.eq.'N')resinflg='n'
      if(resumflg.eq.'N')resumflg='n'
      if(tbcflg.eq.'N')tbcflg='n'
      if(contflg.eq.'N')contflg='n'
      if(routeid1flg.eq.'N')routeid1flg='n'     ! <<<<<
      if(crseflg.eq.'N')crseflg='n'  
	if(ensimflg.eq.'N')ensimflg='n'
	if(picflg.eq.'N')picflg='n'
      if(wetid1flg.eq.'N')wetid1flg='n'         ! <<<<<
      if(modelid1flg.eq.'N')modelid1flg='n'     ! <<<<<
      if(shdflg.eq.'N')shdflg='n'
	if(trcid1flg.eq.'N')trcid1flg='n'         ! <<<<<
	if(frcid1flg.eq.'N')frcid1flg='n'         ! <<<<<
	if(grdflg.eq.'N')grdflg='n'
	if(ntrlflg.eq.'N')ntrlflg='n'
	if(nudgeflg.eq.'N')nudgeflg='n'
	if(flowresetflg.eq.'N')flowresetflg='n'
	if(pafflg.eq.'N')pafflg='n'
	if(fliflg.eq.'N')fliflg='n'
	if(lakeflg.eq.'N')lakeflg='n'
	if(iceflg.eq.'N')iceflg='n'
	if(icerivflg.eq.'N')icerivflg='n'
	if(icelakeflg.eq.'N')icelakeflg='n'
	if(tb0flg.eq.'N')tb0flg='n'
	if(xmlflg.eq.'N')xmlflg='n'
	if(nbsflg.eq.'N')nbsflg='n'
	if(fewsflg.eq.'N')fewsflg='n'
      if(netCDFoutflg.eq.'N')netCDFoutflg='n'
        
      if(snwflg.eq.' ')snwflg='n'
      if(snwid1flg.eq.' ')snwid1flg='n'         ! <<<<<
      if(sedid1flg.eq.' ')sedid1flg='n'         ! <<<<<
      if(vapflg.eq.' ')vapflg='n'
      if(smrflg.eq.' ')smrflg='n'
      if(resinflg.eq.' ')resinflg='n'
      if(resumflg.eq.' ')resumflg='n'
      if(tbcflg.eq.' ')tbcflg='n'
      if(contflg.eq.' ')contflg='n'
      if(routeid1flg.eq.' ')routeid1flg='n'     ! <<<<<
      if(crseflg.eq.' ')crseflg='n'  
	if(ensimflg.eq.' ')ensimflg='n'
	if(picflg.eq.' ')picflg='n'
      if(wetid1flg.eq.' ')wetid1flg='n'         ! <<<<<
      if(modelid1flg.eq.' ')modelid1flg='n'     ! <<<<<
      if(shdflg.eq.' ')shdflg='n'
	if(trcid1flg.eq.' ')trcid1flg='n'         ! <<<<<
	if(frcid1flg.eq.' ')frcid1flg='n'         ! <<<<<
	if(grdflg.eq.' ')grdflg='n'
	if(ntrlflg.eq.' ')ntrlflg='n'
	if(nudgeflg.eq.' ')nudgeflg='n'
	if(flowresetflg.eq.' ')flowresetflg='n'
	if(pafflg.eq.' ')pafflg='n'
	if(fliflg.eq.' ')fliflg='n'
	if(lakeflg.eq.' ')lakeflg='n'
	if(iceflg.eq.' ')iceflg='n'
	if(icerivflg.eq.' ')icerivflg='n'
	if(icelakeflg.eq.' ')icelakeflg='n'
	if(tb0flg.eq.' ')tb0flg='n'
	if(xmlflg.eq.' ')xmlflg='n'
	if(nbsflg.eq.' ')nbsflg='n'
	if(fewsflg.eq.' ')fewsflg='n'
      if(netCDFoutflg.eq.' ')netCDFoutflg='n'
	
d      if(iopt.eq.2)print*,'flags assigned in read_evt'

      if(modelid1flg.ne.'n')then
        if(initflg.eq.'y')then
          print*,'In the event file: initflg can not be `y`'
          print*,'if the modelflg = not `n`'
          print*,'i.e. if the modelflg = `i`, `r` or `l`'
          print*
          stop 'Program aborted in rdevt @ 257'
        endif
      endif

!     Added March 10/07  NK
!     rev. 10.1.67 Feb.  18/17  - NK: Ignore start year in subsequent event files 
      if(firstpass)then
	  year_start=year1
	  mo_start=mo1
	  day_start=day1
	  hour_start=hour1
!     rev. 10.2.46 Jan.  21/19  - NK: Added netCDF capability
        if(fewsflg.eq.'y')then
            netCDFflg=.true.
        else
            netCDFflg=.false.
        endif
      endif
!     REV. 10.1.42 Oct   20/16  - NK: Reinstated read_ice_factor.f as default if present
!     WATROUTE:
!     for watroute, the ice correction can not be calculated (no temepratures)
!     so the correction can only be made with the ice_factor.tb0 file      
      if(modelflg.ne.'n')iceflg='n'        ! for watroute, ice_fctr can NOT be calculated
      if(modelflg.ne.'n')icerivflg='n'     !i.e. when = i, l or r
      if(icerivflg.eq.'y')iceflg='y'       !i.e ice_fctr(n) will be calculated
      if(iceflg.eq.'y')icerivflg='y'       !i.e ice_fctr(n) will be used on rivers
!     If the ice factors can be read from the lake_ice_factor.tb0 file, then
!     they can be used with watroute AND for CHARM it will override the calculated values

      if(id.eq.1)then
!       rev. 9.1.45  Jun.  11/03  - runoff, recharge and leakage files added 
        if(routeflg.eq.'y'.and.modelflg.ne.'n')then
          print*
          print*,' You can not have the routeflg = y AND'
          print*,' the modelflg = r, i or l '
          print*
          print*,' Please pick one or the other as "n"'
          print*,'   and try again. Sorry.'
          print*
          stop ' Program aborted in rdevt @ 198'
        endif
      endif

!     rev  9.1.24  Sep.  11/02  - Added scaleallsnw to set snw scale in event #1
!     scaleall is the master multiplier
!     if larger than 0.0, for the first event then use it throughout the run
!     scaleall overrides everything!!!!
      if(id.eq.1)then
c        if(event%readscale.gt.0.0)then
        if(int(readscale).ne.1)then
          scaleall=readscale   ! scaleall set for all events to come
          scale=readscale
        else
          scaleall=0.00
          if(scale.lt.0.00001)scale=1.0
        endif
      else
!       override event scaling
        if(scaleall.gt.0.0)scale=scaleall
      endif
      if(scale.lt.0.00001)scale=1.0
      
!     scaleallsnw is the master multiplier for snow
!     if larger than 0.0, for the first event then use it throughout the run
!     It will override the eventsnowscalefactor
      if(id.eq.1)then
        if(readscalesnw.gt.0.0)then
          scaleallsnw=readscalesnw   ! scaleallsnw set for all events 
        else
          scaleallsnw=0.00
        endif
      endif
!     This will scale snw with evt #1 if readscalesnw > 0.0   
      if(scalesnw.ne.0.0.and.scaleallsnw.ne.0)then
          print*,'Warning'
          print*,':eventsnowscalefactor will be ignored'
          print*,'To use this factor, :snowscalefactor must be 0.0'
          print*,':snowscalefactor will be used if you proceed'
          pause 'hit enter to continue - @ 1192; Ctrl C to abort'
      endif
      if(scaleallsnw.gt.0.0)then
        scalesnw=scaleallsnw
        write(98,*)'Info: snow precip scaled by ',scalesnw
        write(98,*)'Info: as set in evt #1 :snowscalefactor'
        write(98,*)'Info: this will over write :eventsnowscalefacto'
      endif

!     rev. 10.3.06 Feb.  27/20  = NK Fixed temperature scaling factors for single events
!     scaletem is used in read_temp to scale the temperature
!     It's value is set here 
!     Keyword in evt file   = variable name here
!     :eventtempscalefactor = scaletem=event%scaletem          !event temp scale factor
!     :tempscalefactor      = readscaletemp=event%readscaletem !temp scale factor - all events
      if(id.eq.1)then
        if(readscaletemp.lt.-0.01.or.readscaletemp.gt.0.01)then
          scalealltem=readscaletemp   ! readscaletemp set for all events 
          write(98,*)'Info: all temps modified by ',scalealltem
        else
          scalealltem=0.00
        endif
      endif
      if(scaletem.lt.-0.01.or.scaletem.gt.0.01)then
!         Note:  this overrides the scalealltem set in first event bur will be used in 
!                all events where scaletem = 0.0          
	    scaletem=scaletem
          write(98,*)'Info: temps for event #',id,'changed by ',scaletem
      else
	    scaletem=scalealltem          ! used in read_temp to modify temperatures
	endif

d      if(iopt.eq.2)print*,' In rdevt, passed location 870'

!     rev. 9.8.25  Sep.  26/12  - NK: Added warning for resumflg=y and ID > 1
      if(id.gt.1.and.resumflg.ne.'n')then
        print*,'WARNING:'
        print*,'resumflg =',resumflg
        print*,'Carrying on with the model values'
        print*,'resumflg set to `n` & resume files ignored '
        print*
        write(98,*)'WARNING:'
        write(98,*)'resumflg =',resumflg
        write(98,*)'Carrying on with the model values'
        write(98,*)'resumflg set to `n` & resume files ignored '
        write(98,*)
        resumflg='n'
      endif

      if(iopt.ge.1)then
!       write this stuffto spl.txt for the record
        if(debug_output)write(51,*)'flags and data as used in the model'
        if(debug_output)write(51,*)'Some data/flags from evt 1 are used'
        if(debug_output)write(51,*)'#                       '        
        if(debug_output)write(51,*)':snwflg                 ',snwflg         
        if(debug_output)write(51,*)':sedflg                 ',sedid1flg         
        if(debug_output)write(51,*)':vapflg                 ',vapflg         
        if(debug_output)write(51,*)':smrflg                 ',smrflg         
        if(debug_output)write(51,*)':esimflg                ',resinflg   
        if(debug_output)write(51,*)':tbcflg                 ',tbcflg         
        if(debug_output)write(51,*)':resumflg               ',resumflg         
        if(debug_output)write(51,*)':contflg                ',contflg         
        if(debug_output)write(51,*)':routeflg              ',routeid1flg         
        if(debug_output)write(51,*)':crseflg                ',crseflg         
        if(debug_output)write(51,*)':ensimflg               ',ensimflg         
        if(debug_output)write(51,*)':picflg                 ',picflg         
        if(debug_output)write(51,*)':wetflg                 ',wetid1flg         
        if(debug_output)write(51,*)':modelflg              ',modelid1flg         
        if(debug_output)write(51,*)':shdflg                 ',shdflg         
        if(debug_output)write(51,*)':trcflg                 ',trcid1flg  
        if(debug_output)write(51,*)':frcflg                 ',frcid1flg  
        if(debug_output)write(51,*)':initflg                ',initflg    
        if(debug_output)write(51,*)':grdflg                 ',grdflg    
        if(debug_output)write(51,*)':ntrlflg                ',ntrlflg    
        if(debug_output)write(51,*)':nudgeflg               ',nudgeflg    
        if(debug_output)write(51,*)':flowresetflg         ',flowresetflg    
        if(debug_output)write(51,*)':divertflg              ',divertflg    
        if(debug_output)write(51,*)':pafflg                 ',pafflg
        if(debug_output)write(51,*)'#                       '        
        if(debug_output)write(51,*)':intSoilMoisture   ',(smc5(i),i=1,5)         
        if(debug_output)write(51,*)':rainConvFactor         ',conv         
        if(debug_output)write(51,*)':eventPrecipScaleFactor ',scale         
        if(debug_output)write(51,*)':precipScaleFactor      ',readscale         
        if(debug_output)write(51,*)':eventSnowScaleFactor   ',scalesnw        
        if(debug_output)write(51,*)':snowScaleFactor      ',readscalesnw        
        if(debug_output)write(51,*)':eventTempScaleFactor   ',scaletem         
        if(debug_output)write(51,*)':tempScaleFactor     ',readscaletemp          
        if(debug_output)write(51,*)'#                       '        
        if(debug_output)write(51,*)':hoursRainData          ',nhg         
        if(debug_output)write(51,*)':hoursFlowData          ',nhf    
        if(debug_output)write(51,*)':deltat_report       ',deltat_report     
        if(debug_output)write(51,*)'#                       '
      endif
      
      firstpass=.false.
        
      RETURN

! FORMATS:

 1000 format(26x,a31)
 1010 format(26x,a14,20a1)
 1020 format(' ', /'you are trying to link too many events'/
     *'      the maximum is 500, you entered,',i5,/
     *'      program aborted. fix the event file and try again')
 1030 format(' ','Unit no. =',i3,' file no',i3,' = ',a35)
 1031 format(' Input files from event.evt')
 1040 format(' Event no. ',i5)
 1050 format(' ')
 1100 format(26x,6f5.2)
 1110 format(26x,6f5.2)
 1200 format(26x,3i5)
 1300 format(a30)
 6080 format(a14)
 6081 format(i4)

99001 format(a30)
99002 format(a30,a1)
99003 format(a30,f12.2)
99004 format(a30,i10)
99005 format(a30,a31)

      END SUBROUTINE read_evt

