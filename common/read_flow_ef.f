      SUBROUTINE read_flow_ef(hdrflg,date)

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
!*****************************************************************************
!  READ_FLOW_EF - written Aug/06 by Dave Watson, CHC
!	- Derived from rdflow written by Nick Kouwen
!	This subroutine reads streamflow hydrograph (STR) file 
!	(tb0 format)
!*****************************************************************************

!     rev. 9.5.58  Apr.  16/09  - NK: added nudgeflg for forcing gauge flows
!     rev. 9.5.62  Sep.  04/09  - NK: new tb0 file for DW routing
!     rev. 9.5.64  Sep.  16/09  - NK: corrected nudging wrt first event
!     rev. 9.7.07  Sep.  05/10  - NK: increased allowed # flow stations from 128 to 512
!     rev. 9.8.37  Oct.  27/12  - NK: added section to read_flow_ef to check # columns = no
      use area_watflood

! TB0 data module
      USE EF_Module

      implicit none
      TYPE(FlowParam) :: header
      TYPE(FlowColumnMetaData) :: colHeader

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      CHARACTER(1) :: hdrflg,firstpass,usefirstpass

      INTEGER :: i,j,k,l,n,ios,no_old,iAllocate,ideallocate,
     *            linecount,dl,lineno,idlast,WFruntime
      LOGICAL exists,nudge,msgFlag
      character(14) :: date,junk
!     rev. 10.2.58 Jul.  17/19  - NK Remap nopt for changed flow stationlocations in FEWS
      real(4),  dimension(:),  allocatable :: xNudge,yNudge
      integer,  dimension(:),  allocatable :: nNudge,jNudge,iNudge
      integer,  dimension(:),  allocatable :: nopt_temp
      character(7), dimension(:), allocatable :: gage_temp

      data firstpass/'y'/

! Parameter type definitions
      integer(4) unitNum, flnNum, iStat

! Local variables
      character*128 keyword, value
      character*4096 line, subString, tmpString
      integer lineLen, keyLen, wordCount
      logical rStat, lineType, foundEndHeader, insideColMetaData

! Set unit and fln number
      unitNum = 36
      flnNum = 6
!     rev. 10.2.17 Feb.  26/18  - NK: Added WFruntime to read_flow_ef
      WFruntime=0

!     rev. 9.9.35  Oct.  20/14  - NK: Added keyword & file checks
      if(FLtype(6).eq.'tb0')then
      
      inquire(FILE=fln(6),EXIST=exists)
	if(exists) then
	  open(unit=unitNum,file=fln(flnNum),status='old',iostat=ios)
        if(IsFileTypeTB0(fln(6))) then
d         print*,'reading flow header'      
!         Open the file
          open(unit=unitNum,file=fln(flnNum),status='old',iostat=ios)
          if(ios.ne.0)then
            print*,'Problems opening ',fln(flnNum)(1:40)
            print*
            STOP ' Stopped in read_flow_ef @ 67'
          endif
        else
!         call rdflow('0',date)   
          print*,'Old format .str files not accepted'
          print*,'Please create yyyymmdd_str.tb0 files & rerun'
          stop 'Program aborted in sub @ 192'
        endif
	  if(iopt.ge.1)print*,'Opened unit'
     *                      ,unitNum,' filename  ',fln(flnNum)(1:40)
	else
	  print*, 'Attempting to open file name ',fln(flnNum)(1:40)
	  print*, 'Unit number ',unitNum
	  print*, 'but it is not found (in this location)'
	  STOP 'Program aborted in read_flow @64'
      endif
      
!     Initialize default values
      CALL InitFlowParam(header)	

!     Search for and read tb0 file header
      linecount=0
      line(1:1) = '#'
      foundEndHeader = .false.
      insideColMetaData = .false.

      do WHILE((.NOT.foundEndHeader) .AND.
     &	    ((line(1:1) .eq. '#') .OR.
     &		(line(1:1) .eq. ':') .OR.
     &		(LEN_TRIM(line) .eq. 0))) 	
                  linecount=linecount+1
!	if(iopt.eq.2)print*,'reading line ',linecount,' in read_flow_ef'
        read(UNIT=unitNum, FMT='((A))', iostat=ios) line	! read a line
        if(ios .eq. -1)then
            write(6,'((A))') 'ERROR: Premature EndOfFile encountered'
            STOP ' Stopped in read_flow_ef'
        end if

        rStat = Detab(line)				! replace tabs with spaces
        line = ADJUSTL(line)		! Get rid of leading white space
        lineLen = LEN_TRIM(line)		! Find the length excluding trailing spaces

        if(line(1:1) .eq. ':')then
            wordCount = SplitLine(line, keyword, subString)	! find the keyword
            rStat = ToLowerCase(keyword)
            KeyLen = LEN_TRIM(keyword)

            if(keyword(1:KeyLen) .eq. ':endheader')then
                foundEndHeader = .TRUE.

            else if(keyword(1:KeyLen) .eq. ':columnmetadata')then
                insideColMetaData = .TRUE.
            else if(keyword(1:KeyLen) .eq. ':endcolumnmetadata')then
                insideColMetaData = .FALSE.
            else if(insideColMetaData) then
                iStat = ParseFlowColumnMetaData(colHeader,keyword,
     &							keyLen,subString)
                if(iStat .lt. 0) then
                    write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
                    write(*,'(2(A))') '   in line: ',line					
                    STOP ' Stopped in read_flow_ef'
                    return
                endif
            else
                iStat = ParseFlowParam(header,keyword,keyLen,subString)
                if(iStat .lt. 0) then
                    write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
                    write(*,'(2(A))') '   in line: ',line					
                    STOP ' Stopped in read_flow_ef'
                    return
                else if(iStat .eq. 0) then
!					write(*,'((A), (A))')  'Unrecognized keyword line: ',
!     &										line
                endif
            end if
        end if
        
!     rev. 10.2.17 Feb.  26/18  - NK: Added WFruntime to read_flow_ef
        if(line(1:10).eq.'#WFruntime')then
            read(line,*)junk,WFruntime
        endif
      end do	
!     ***************************************
!	Finished reading the header
!     ***************************************
        
!     Assign the variables from the types
      kt =  header%tb0p%deltaT    !data time step in hours
      irdt =  header%routingDeltaT    ! routing time step in hours
      flowUnitConv = header%tb0p%unitConv

!     FILL IN FLOW DATA flowfillflg='y'
!     this is when you have say 24 hr flow data but you want to plot hourly
      flowfillflg='n'
      if(GetBoolean(header%fillFlag) .eq. .true.) flowfillflg = 'y' 

!     no    =     no of streamflow stations
      no = colHeader%tb0cmd%colCount

d     print*,'no of flow stations found =',no

!     rev. 9.7.07  Sep.  05/10  - NK: increased allowed # flow stations from 128 to 512
      if(no.gt.1024)then
        print*,'no of flow stations found =',no
        print*,'maximun nuber allowed is 1024'
        print*,'Please contact developer(s) to have allowance'
        print*,'to have allowanceincreased'
        print*
        stop 'Program aborted in read_flow_ef @ 149'
      endif

!     nl    =     no of hours of streamflow data
!     Scan lines of data
      rewind unitNum
!	nl = CountDataLinesAfterHeader(unitNum)

      nl = CountDataLinesAfterHeader(unitNum)*kt

d     print*,'no columns in the str file', no
d     print*,'no of data lines in the str file ',nl
d     print*
      if(nl.gt.8784)then
        print*,'WARNING'
        print*,'No of hours in the str file = ',nl
        print*,'This is > no of hours in one year'
        print*,'Check deltaT and # lines in the file'
        print*
      endif
      
!     rev. 10.2.17 Feb.  26/18  - NK: Added WFruntime to read_flow_ef
!     Undocumented option to reduce run time to < nl
!     this can be used for florecasting prctice - 
!     so long str file can be left but shorter run made      
      if(WFruntime.gt.0)nl=min(WFruntime,nl)
      

!******************************************************************
!******************************************************************
!******************************************************************
!     added (or replaced Nov. 10/10 nk
      mhrd=nl
!******************************************************************
!******************************************************************
!******************************************************************

!     fixed  nk  Nov. 22/06

      rewind unitNum
      CALL GoToStartOfData(unitNum)
      
        if(firstpass.eq.'y')then
          allocate(qhyd(no,nl),stat=iAllocate)
          if(iAllocate.ne.0) STOP
     *		'Error with allocation of qhyd arrays in in read_flow_ef @ 234'
        endif
      
      else          ! FLtype(6).eq.'tb0'
!     rev. 10.2.46 Jan.  21/19  - NK: Added netCDF capability
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call read_ts_flow_nc(6)                                        ! netCDF
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!           FILL IN FLOW DATA 
            flowfillflg='y'
      endif        ! FLtype(6).eq.'tb0'
      
!       rev. 9.1.68  Dec.  19/04  - NK: rewrote rdflow c/w memory allocation 
!       allocate stuff      
        if(firstpass.eq.'y')then
            no_old=no  ! # flow stations
            nl_max=nl  ! # # hours of recorded data
!             ystr & xstr allocation is here or in read_ts_nc:            
              if(.not.allocated(ystr))allocate(ystr(no),xstr(no))
              allocate(nnsum(no),ppsum(no),
     *	    nxtbasin(no),statnerr(no),suberr(no),subarea(no),
     *		iflowgrid(no),qbar(no),iy(no),jx(no),            
     *	    frc(no,4),hydfctr(no),synfctr(no),area(no),
     *		gage(no),nlow(no),nopt(no),nopt1(no),
     *	    qsyn(no,nl),statnash(no),
     *        flowflag(no),mean_flow_error(no),
     *        us_sta(no),ds_sta(no),
     *		qhyd_dly(no,nl/24),qsyn_dly(no,nl/24),
     *		qhyd_dly_sum(no),qsyn_dly_sum(no),
!     rev. 9.9.30  Sep.  30/14  - NK: fixed allocation for qhyd_mly & qsyn_mly
c     *		qhyd_mly(no,max(1,nl/672)),qsyn_mly(no,max(1,nl/672)),
!     rev. 10.1.02 Oct.  09/15  - NK: Fixed allocation for qhyd_mly & qsyn_mly
     *		qhyd_mly(no,12),qsyn_mly(no,12),
c     *		areasum(no),
     *		qhyd_mly_sum(no),qsyn_mly_sum(no),stat=iAllocate)
            if(iAllocate.ne.0) STOP
     *		'Error with allocation of area5a arrays in in read_flow_ef'
!     moved to flow init
c!         TS - ALLOCATION OF AREANASHA ARRAYS
c			allocate(aa(no),bb(no),cc(no),ashnum(no),ashden(no),
c     *		rsquare(no),nq(no),nqc(no),qbarobs(no),stat=iAllocate)
c			if(iAllocate.ne.0) STOP
c     *		'Error with allocation of areanasha arrays in sub'

!     REV. 10.1.34 Jul   05/16  - NK: Added Obs. & Model mean flows to wfo file
            allocate(wfo_qhyd(na,nl),wfo_qsyn(na),stat=iAllocate)
            if(iAllocate.ne.0) STOP
     *		'Error with allocation of wfo_qhyd arrays in read_flow_ef'

!     rev. 10.2.16 Feb.  14/18  - NK: Added bankfull flow calculation
            allocate(bankfull(ycount,xcount),stat=iAllocate)
            if(iAllocate.ne.0) STOP
     *		'Error with allocation of bankfull array in read_flow_ef'

!     rev. 
            allocate(xNudge(no),yNudge(no),nNudge(no),gage_temp(no),
     *               jnudge(no),iNudge(no),nopt_temp(no),stat=iAllocate)
            if(iAllocate.ne.0) STOP
     *		'Error with allocation of xNudge & yNudge in read_flow_ef'


        else                    !  firstpass
!         check to see memory allocated is adequate      

            if(no.ne.no_old)then
                print*,'No of streamflow stations has been changed in'
                print*,'in file ',fln(6)
                print*,'This is not permitted'
                print*
                stop 'Program aborted in rdstr @ 172'
            endif
            if(nl.gt.nl_max)then
                nl_max=nl

!           the file length is longer than any of the previous events so 
!           more memory has to be allocated


!           DEALLOCATION OF ARRAYS FROM AREA10A:
                deallocate(qhyd,qsyn,qhyd_dly,qsyn_dly,
     *			qhyd_mly,qsyn_mly,wfo_qhyd,stat=iDeallocate)
                if(iDeallocate.ne.0)then
                    print*,'Error with deallocation of qhyd & qsyn'
                    print*
                    stop 'Program aborted in read_flow_ef @ 180'
                endif

                  print*,'reallocation for more flows -',nl

                allocate(qhyd(no,nl),qsyn(no,nl),
     *			qhyd_dly(no,nl/24),qsyn_dly(no,nl/24),
c     *			qhyd_mly(no,max(1,nl/672)),qsyn_mly(no,max(1,nl/672)),
     *			qhyd_mly(no,max(12,nl/672)),qsyn_mly(no,max(12,nl/672)),
     *                  stat=iAllocate)
                if(iAllocate.ne.0) STOP
     *			'Allocation Error: area5b arrays in read_flow_ef@230'

!     REV. 10.1.34 Jul   05/16  - NK: Added Obs. & Model mean flows to wfo file
                allocate(wfo_qhyd(na,nl),stat=iAllocate)
                if(iAllocate.ne.0) STOP
     *		    'Error with allocation of wfo_qhyd arrays in read_flow_ef'


!     rev. 9.5.62  Sep.  04/09  - NK: new tb0 file for DW routing
!                 these can only be deallocated if allocated in the first place
                  if(allocated(lake_elv))then
                deallocate(lake_elv,lake_stor,
     *            lake_inflow,net_lake_inflow,lake_outflow,
     *            net_lake_outflow,del_stor,
     *			qdwpr,stat=iDeallocate)
                if(iDeallocate.ne.0)then
                    print*,'Error with deallocation of resv stuff'
                    print*
                    stop 'Program aborted in read_flow_ef @ 190'
                endif

!     rev. 10.1.65 Jan.  28/17  - NK: Fixed allocate lake_elv from read_flow 
!                     moved to read_resv_ef.f
!                allocate(lake_elv(noresv,nl),lake_stor(noresv,nl),
                allocate(lake_stor(noresv,nl),
     *            lake_inflow(noresv,nl),net_lake_inflow(noresv,nl),
     *            lake_outflow(noresv,nl),
     *            net_lake_outflow(noresv,nl),del_stor(noresv,nl),
     *                  qdwpr(noresv,nl),
     *                  stat=iAllocate)
                if(iAllocate.ne.0) STOP
     *			'Allocation Error:  arrays in read_flow_ef@239'
                  endif
            endif
      endif                   !  firstpass
!       rev. 9.1.68  Dec.  19/04  - NK: rewrote rdflow c/w memory allocation 
!     rev. 10.2.46 Jan.  21/19  - NK: Added netCDF capability
      
      if(FLtype(6).eq.'tb0')then
!       ASSIGN THE GAUGE LOCATIONS NAMES AND FUNCTION COEFFICIENTS
        do l=1,no
            gage(l) = colHeader%tb0cmd%colName(l) ! streamflow station name
            xstr(l) = colHeader%tb0cmd%colLocX(l) ! x coordinate
            ystr(l) = colHeader%tb0cmd%colLocY(l) ! y coordinate
            frc(l,1) = colHeader%colCoeff1(l) ! coefficient 1
            frc(l,2) = colHeader%colCoeff2(l) ! coefficient 2
            frc(l,3) = colHeader%colCoeff3(l) ! coefficient 3
            frc(l,4) = colHeader%colCoeff4(l) ! coefficient 4
            nopt(l) = colHeader%colValue1(l) ! nopt
        end do
        deallocate(colHeader%tb0cmd%colName)
        deallocate(colHeader%tb0cmd%colLocX)
        deallocate(colHeader%tb0cmd%colLocY)
        deallocate(colHeader%colCoeff1)
        deallocate(colHeader%colCoeff2)
        deallocate(colHeader%colCoeff3)
        deallocate(colHeader%colCoeff4)
        deallocate(colHeader%colValue1)
!       endif moved to here  from below so nudging will work for FEWS        
      endif    ! FLtype(6).eq.'tb0'

!       rev. 9.2.21  Nov.  11/05  - NK: 
!       Set nopt in first event .str file 
!       This allows a pick of the stations for opt 
!       in the first event only
!       All nopt will be the same for all events.
!       Does not allow opt stations to change during the run. For that, each
!       .str file has to be edited.
        
      if(firstpass.eq.'y')then 
        nudgeflg1=nudgeflg
        usefirstpass='n'
!     REV. 10.1.25 Feb.  21/16  - NK: Added nudge_flags.txt 
!       it's easiest to make this file by copying the header stuff from a str file
!       and transposing & rearranging the columns
        if(nudgeflg.ne.'n')then   ! note: something other than 'n'
          INQUIRE(FILE='strfw\nudge_flags.xyz',EXIST=exists)
          if(exists)then
!           if this file exists, it replaces the values in the firts event only
!           For FEWS, thise file is the only way to nudge              
            open(unit=99,file='strfw\nudge_flags.xyz',status='old')
            do l=1,no
              read(99,*)xNudge(l),yNudge(l),nopt_temp(l),gage_temp(l)
            end do
!     rev. 10.2.58 Jul.  17/19  - NK Remap nopt for changed flow stationlocations in FEWS
!           check the nudge locations match the flow locations
            msgFlag=.true.
            if(iopt99)write(98,*)
     *   '      flow station locations      nudge locations '      
            if(iopt99)write(98,*)
     *   'rch#  long     lat     iy   jx    long     lat   iN   jN   fl'
            do l=1,no
!             find the grid #                
              iy(l)=int((ystr(l)-yorigin)/ydelta)+1
              jx(l)=int((xstr(l)-xorigin)/xdelta)+1
              iNudge(l)=int((yNudge(l)-yorigin)/ydelta)+1
              jNudge(l)=int((xNudge(l)-xorigin)/xdelta)+1
              if(iopt99)write(98,98000)l,xstr(l),ystr(l),jx(l),iy(l),
     *                     xNudge(l),yNudge(l),
     *                     jNudge(l),iNudge(l),nopt_temp(l)
98000          format(i5,2f8.3,2i5,2f8.3,3i5)              
            end do
            do l=1,no
                if(jNudge(l).ne.jx(l).or.iNudge(l).ne.iy(l))then
                  if(msgflag)then  
                    write(98,98001)'WARNING:',
     *               'Nudge location for some or all stations ',
     *               ' do not match the flow station',
     *               ' location in the strfw location.'
                    if(iopt99)write(98,*)
     *                      '          location        x           y'
98001               format(a8,a40,a30,a32)                    
                    msgFlag=.false.
                  endif
                  if(iopt99)then
                  write(98,*)'flow ',l,jx(l),iy(l)
                  write(98,*)'nudge',l,jNudge(l),iNudge(l),nopt_temp(l)
                  endif
                endif
            end do
            close(unit=99,status='keep')
            write(98,*)'Finsished reading strfw\nudge_flags.xyz'
          endif
            
          if(.not.msgFlag)then
!         Reset the flags using the lat-long in the nudge_flags.xyz file
              write(98,*)
     *        'Reset the flags using the lat-long in nudge_flags.xyz' 
              do l=1,no
                  nopt(l)=1
              end do
              if(iopt99)write(98,*)
     *         'flow station #  nudge station #  nudgeflg'   
              do i=1,no        ! nudge file
                  do l=1,no    ! flow file
                      if(jNudge(i).eq.jx(l).and.iNudge(i).eq.iy(l))then
                          nopt(l)=nopt_temp(i)
                          if(iopt99)write(98,*)l,i,nopt(l)
                      endif
                  end do
              end do
          endif

            write(98,*)'Nudge flags as used in CHARM:'
            write(98,*)'   Flow station #  Nudgeflg'
            do l=1,no
                write(98,*)l,flow_sta_name(l),nopt(l)
            end do

!     rev. 10.2.58   end

            open(unit=99,file='debug\nudged_stations.xyz',
     *                              status='unknown')
            do l=1,no
              if(nopt_temp(l).eq.2)then
                write(99,*)xstr(l),ystr(l),l,gage_temp(l)
              endif
            end do
            close(unit=99,status='keep')
          endif

        if(nudgeflg.eq.'1')usefirstpass='y' 
        do l=1,no
          nopt1(l)=nopt(l)
!     rev. 9.5.64  Sep.  16/09  - NK: corrected nudging wrt first event
          if(nopt(l).eq.-1)usefirstpass='y' 
!         use first event nopt if even one value = -1
          nopt(l)=abs(nopt1(l))
        end do

!       flag the stations that will be nudged:
        nudge=.false.
        do l=1,no
          if(nopt(l).eq.2)then
            print*,'value1(',l,')=2 - nudging in 1st event is risky'
            if(nopt(l).eq.2)nudge=.true.
          endif
        end do
        if(nudge)then
          print*,'Flows nudged at designated flow stations'
          print*,'as designated in event #1'
        endif
        if(dds_flag.eq.0)print*,'usefirstpass =',usefirstpass
      else
!       for subsequent events use the first event opt flags
c        print*,'usefirstpass$ =',usefirstpass
        if(usefirstpass.eq.'y')then
!         used only if some stations want to be nudged
          nudge=.false.
          do l=1,no
            nopt(l)=abs(nopt1(l))
            if(nopt(l).eq.2)nudge=.true.
          end do
          if(nudge)then
            print*,'Flows nudged at designated flow stations'
            print*,'as designated in event #1'
          endif
          
        else
!         do nothing - nopt values from current event files are used
!         use nopt as read in this event's flow file
          nudge=.false.
          do l=1,no
            if(nopt(l).eq.2)nudge=.true.
          end do
          if(nudge)then
            print*,'Flows nudged at designated flow stations'
            print*,'as designated in the current event'
          endif
        endif
        endif
        
!          moved up        
!      endif    ! FLtype(6).eq.'tb0'

      if(nudgeflg1.eq.'a')then
!       nudge all the flows - this overrides everything
        do l=1,no
        nopt(l)=2
        end do
        write(98,*)'Warning: Flows nudged at ALL flow stations'
      endif

d	write(51,*)nopt


!     IF STEP < 1 km THEN A LOCAL COORDINATE GRID HAS TO BE USED
!     NOT THE UTM GRID
!     THIS WAS TAKEN OUT BY FRANK S FEB/99 TO ALLOW FOR GRIDS
!     DEFINED ON A LAT/LONG BASIS

!     FLOW STATION COORDINATES ARE NOW CONVERTED TO COMPUTATION
!     GRID COORDINATES

      if(id.eq.1)then
!     rev. 9.2.33  Feb.  14/06  - NK: /str stations from first event ONLY!!
        do l=1,no
!			convert iy and jx to local coordinates
!			Jul. 21/04  nk
            iy(l)=int((ystr(l)-yorigin)/ydelta)+1
            jx(l)=int((xstr(l)-xorigin)/xdelta)+1
        end do
          if(iopt.ge.1)then
          write(55,*)
          write(55,*)'Just after converting streamflow stations to '
          write(55,*)'grid coordinates'
          do l=1,no
            write(55,1778)id,l,iy(l),jx(l)
          end do
          write(55,*)
          endif
      endif

      lineno=0
      if(FLtype(6).eq.'tb0')then
!       WATFLOOD COLUMN FORMAT
        do j=kt,nl,kt   ! changed nk  sept. 29/06

!     rev. 9.8.37  Oct.  27/12  - NK: added section to read_flow_ef to check # columns = no
!           Added this section to check # cols = # stations in the headet
!           Oct. 27/2012  NK
            read(UNIT=unitNum, FMT='((A))', iostat=ios) line	! read a line
            if(ios .eq. -1)then
               write(6,'((A))') 'ERROR: Premature EndOfFile encountered'
               STOP ' Stopped in read_flow_ef @ 266'
            end if
            lineno=lineno+1
            wordCount = countwords(line)
            if(no.ne.wordcount)then
              print*,'Error:'
              print*,'No of columns of data ',wordcount
              print*,'not equal to the no columns in the header',no
              print*,'in data line ',lineno
              no=wordcount
c			  stop 'Program aborted in read_flow @ 276'
            endif
            
            read(line,*,iostat=ios)(qhyd(l,j),l=1,no)
            
c      print*,j,qhyd(33,j),nopt(33)            

c			read(unitNum,*,iostat=ios)(qhyd(l,j),l=1,no)
            if(ios.ne.0)then
                print*, 'In strfw'
                print*,' NEW format flow file'
                print*,' problems reading the hydrograph at hour '
     *				,j/kt
                print*,' and column',l,'  no=',no
                print*,'Weird hidden characters in the file will
     *			 do this'
                print*,'last values read:'
                if(j.gt.2)then
                    write(*,206)(qhyd(l,j-2),l=1,no)
                endif
                if(j.gt.1)then
                    write(*,206)(qhyd(l,j-1),l=1,no)
                endif
                write(*,206)(qhyd(l,j),l=1,no)
                print*
                print*,'Look for this line:'
                print*,line(1:72)
                stop 'Program aborted in strfw @ 368'
                  endif

 9993             format(' error reading flows time step =',i3,/
     *	   ' possible cause: station list does not match # specified'/
     *		' or there may be tabs in the data')

      end do
      
!     rev. 9.1.10  Jan.  29/02  - flow nudging added for nopt(l)=2
!         fill in the blanks - this is needed for nudging
!         added nopt=2 to the conditional Feb. 20/08 -nk-
      if(kt.gt.1)then
        do l=1,no
          do j=kt,nl,kt   ! changed nk  sept. 29/06
            do i=j-kt+1,j-1
              qhyd(l,i)=qhyd(l,j)
c      if(l.eq.5)print*,i,j,qhyd(l,j)                  
            end do
          end do
        end do
      endif

!     apply the flow conversion factor e.g 0.02843 for cfs to cms
      if(flowunitconv.ne.0.0)then
        do j=kt,nl,kt
            do l=1,no
                qhyd(l,j)=qhyd(l,j)*flowunitconv
            end do
        end do
      endif

!	Convert stage hydrographs to flow hydrographs:
      do l=1,no
        if(frc(l,1).ne.0)then
            do j=kt,nl,kt
!			BUG: changed nl from mhtot  NK  Nov. 13/03
                if(qhyd(l,j).gt.0) then
                    if(qhyd(l,j).ge.frc(l,4))then
                        qhyd(l,j)=frc(l,3)
     *		      		+frc(l,1)*(qhyd(l,j)-frc(l,4))**frc(l,2)
                    else
                        qhyd(l,j)=max(frc(l,3),0.0001)
                    endif
                endif
            end do
        endif
      end do

!     FILL IN STREAMFLOW DATA UNTIL END OF RAINGFAL DATA:
      do l=1,no 
        do j=mhrd+kt,nl,kt
            if(qhyd(l,j).lt.0.0)qhyd(l,j)=-1.0  
        end do  
      end do
      endif               ! FLtype(6).eq.'tb0'
      
!     REV. 10.1.34 Jul   05/16  - NK: Added Obs. & Model mean flows to wfo file
      do n=1,na
        do j=1,nl
          wfo_qhyd(n,j)=-999.0
        end do
      end do
!     map the observed flow to the grid (rank) number      
      do l=1,no
        i=iy(l)
        j=jx(l)
!     REV. 10.1.39 Sep   15/16  - NK: Fixed stations outside the watershed for tb0
!       these exclusions are needed when the gauges are outside the watershed         
        if(i.ge.1.and.j.ge.1.and.i.le.ycount.and.j.le.xcount)then
          n=s(i,j)
          if(n.ge.1)then
c          do k=kt,nl,kt
            do k=1,nl
!     REV. 10.1.38 Jul   28/16  - NK: Added noDataValue to WFO & tb0 files
              if(qhyd(l,k).lt.0.0)qhyd(l,k)=-999.0   
              wfo_qhyd(n,k)=qhyd(l,k)
            end do
          endif
        endif
      end do

!     FILL IN FLOW DATA
!     this is when you have say 24 hr flow data but you want to plot hourly
!     This can not be done if you want to extend the routing time step

      if(irdt.le.kt)then   
        if(flowfillflg.eq.'y'.or.flowfillflg.eq.'Y')then
            do l=1,no
                do j=nl,kt,-kt                 ! 744,24,-24
                    do i=j,j-kt+1,-1            ! 744,721,-1
                        qhyd(l,i)=qhyd(l,j)
                    end do 
                end do
            end do
            kt=1
        endif
      else

!           do we need anything here???????????????????????????????????????????

      endif

!     write the converted hydrographs (stage to flow):
      if(iopt.gt.0.and.frc(1,1).ne.0)then
        write(55,5001)
        do j=kt,nl,kt
            write(55,202)(qhyd(l,j),l=1,no)
        end do
      endif

      do l=1,no
        qbar(l)=0.
!		QBAR(I,D,L)=MEAN FLOW DURING STORM PERIOD AT STA.=L,STORM=I
        do j=kt,mhtot,kt
            qbar(l)=qbar(l)+qhyd(l,j)
        end do
        qbar(l)=qbar(l)/float(mhtot/kt)
      end do

!     rev. 9.1.10  Jan.  29/02  - flow nudging added for nopt(l)=2
!     identify grid numbers with flow stations
      do l=1,no
        if(iy(l).le.0.or.jx(l).le.0)then
!         this can happen if the stream gauge is outside the waterhsed
!         as when subwaterhseds are modelled as separate watersheds
!         added Mar 14/04 nk.      
            iflowgrid(l)=0
        elseif(iy(l).gt.ycount.or.jx(l).gt.xcount)then
!         this can happen if the stream gauge is outside the waterhsed
!         as when subwaterhseds are modelled as separate watersheds
!         added Jul. 06/05 nk.          
            iflowgrid(l)=0
        else
            iflowgrid(l)=s(iy(l),jx(l))
        endif
      end do

!     rev. 9.1.68  Dec.  19/04  - NK: rewrote rdflow c/w memory allocation 
!     moved from sub
      if(FLtype(6).eq.'tb0')then
      close(unit=unitNum,status='keep',iostat=ios)
      if(ios.ne.0)then
        print*,'Problem closing unit 36 fln=',fln(6)
        print*
        stop ' program aborted in rdflow @ 623'
      endif
      endif


!     rev. 9.8.51  Mar.  11/13  - NK: Link skiplines in s/r stats to value1 in the str file
      if(id.eq.1)then
        skipflg=.true.
        skiplines=0
!       if even one value1 .ne. 0 then don't skip stats calculation        
        do l=1,no
          if(nopt(l).ne.0)skipflg=.false.  
        end do
        if(skipflg)skiplines=nl/kt
!       skiplines is the number of data lines to be passed over for stats.
        idlast=1
      endif     
      if(idlast.eq.id-1.and.skipflg)then            
!       keep adding to the skiplines as long as consecutive events have
!       all value1's = 0        
!       as sooon as we find one value1 .ne. 0 no more skipping
        do l=1,no
          if(nopt(l).ne.0)skipflg=.false.  
        end do
        if(skipflg)skiplines=skiplines+nl/kt
        idlast=id
      endif
      if(iopt.ge.1)then
        print*,'ID,# data lines skipped for stats:',id,skiplines
      endif
      firstpass='n'
      RETURN

! FORMATS

c  201 format(1024f10.0)
  202 format(1024f10.3)
  203 format(' id=',i5,' l= ',i5,a12,2f12.3,4e10.3)
  205 format(2f12.3,a12,f12.3,4e10.3)
  206 format(1024f10.3)
 1778 format(' id,l,iy(l),jx(l)',6i5)
 1801 format(14x,10(f5.0,a1),i6)
 1802 format(14x,11(f5.0,a1))
 5000 format(' echo recorded hydrographs:')
 5001 format(/' echo converted hydrographs:')
 5004 format(a20,a10)
 5005 format(a20,i5)
 6226 format(' error encountered opening unit=37 fln=',a31/)

      END SUBROUTINE read_flow_ef



