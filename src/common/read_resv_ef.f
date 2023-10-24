      SUBROUTINE read_resv_ef()

!***********************************************************************
!    Copyright (C) 1987-2018 by Nicholas Kouwen and Dave Watson (NRC)  
        
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
!*****************************************************************************
!  READ_RESV_EF - written Sep/06 by Dave Watson, CHC
!     - Derived from rdresv written by Nick Kouwen
!     This subroutine reads the reservoir release (REL) file 
!     (tb0 format)
!*****************************************************************************

!     rev. 9.5.14  Feb.  26/08  - NK: padded rel file for missing data
!     rev. 9.5.23  Mar.  12/08  - NK: fixed allocation error in read_resv_ef
!     rev. 9.5.24  Mar.  18/08  - NK: fixed missing data in read_resl_ef.f
!     rev. 9.5.29  May.  26/08  - NK: fixed initialization in read_resv_ef
!     rev. 9.5.36  Oct.  01/08  - NK: fixed ires bug for unevent dx & dy in read_resv
!     rev. 9.5.41  Oct.  22/08  - NK: read in reservoir coefficients each event
!     rev. 9.5.49  Dec.  31/08  - NK: changed conditional to read releases in rerout
!     rev. 9.5.57  Apr.  13/09  - NK: added ntrlflg for natural lake flows
!     rev. 9.5.59  Jul.  26/09  - NK: added fpet_lake for each lake in ill file
!     rev. 9.5.60  Sep.  01/09  - NK: added deltat_report for lake_sd.csv file
!     rev. 9.5.62  Sep.  04/09  - NK: new tb0 file for DW routing
!     rev. 9.5.81  Jan.  16/09  - NK: allow reservoirs outside watershed in resv file
!     rev. 9.8.57  Apr.  12/13  - NK: Added lakeEflg to stop lake evaporation whan levels very low
!     rev. 9.9.20  Jul.  24/14  - NK: Added dead storage for lakes "stroe_dead"
!     rev. 10.2.27 Jul.  08/18  - NK: replaced lake_inflow_sum with temp_flow_sum
!******************
      use area_watflood

! TB0 data module
      USE EF_Module
      implicit none
      TYPE(ResvParam) :: header
      TYPE(ResvColumnMetaData) :: colHeader

      Integer  :: ios,j,k,i,n,l,jz,jj
      integer  :: noresv_firstpass,miss_count,
     *            iAllocate,iDeallocate
      integer  :: anno(100)
      real*4   :: factor
      real*4,dimension(:,:),allocatable :: ff

!     rev. 9.1.55  Jun.  12/04  - NK: write new str files to strfw\newfmt folder.
      LOGICAL exists,read_res,fudge
      CHARACTER(1)  :: firstpass

      data firstpass/'y'/
      data fudge/.false./

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

! Parameter type definitions
      integer*4 unitNum, flnNum, iStat

! Local variables
      character*4096 line, subString, tmpString
      character*128 keyword, value 
c      character*12  outlet_type
      integer lineLen, keyLen, wordCount
      logical rStat, lineType, foundEndHeader, insideColMetaData

      data read_res/.true./
      if(.not.read_res)return

!     rev. 9.5.62  Sep.  04/09  - NK: new tb0 file for DW routing
!     no need to read rel file when dumping out flows for flow1D
!     all we need is the reach locations and to initialize the variables.

      if(routeflg.eq.'q'.and.firstpass.ne.'y')return

! Set unit and fln number
      unitNum = 37
      flnNum = 7

c      if(id.eq.1.and.ntrlflg.eq.'y'.or.ntrlflg.ne.'y')then
!       when using the natural flow option, 
!       we can bypass reading subsequent rel files
!       Just need to reallocate if events get longer

!     rev. 9.9.35  Oct.  20/14  - NK: Added keyword & file checks
      inquire(FILE=fln(flnNum),EXIST=exists)
      if(.not.exists)then
!       maxr is the # of reached in the shd file - see flowinit
        if(maxr.eq.0)then
!         no problem - rel file not needed        
          if(numa.eq.0.and.dds_flag.eq.0)then
            print*,'WARNING:'
            print*,'Reservoir release (rel) file '
            print*,fln(flnNum)(1:60)
            print*,'is NOT found'
            print*,'Program continues with no lakes or resevoirs'
            print*
          endif
          noresv=0
          noresvi=0
          nrel=0
          read_res=.false.
          return
        else
          print*
          print*,'Fatal error:   file not found'
          print*,fln(flnNum)(1:72)
          print*,'but there are',maxr,' reaches in the shd file'
          print*,'Please check the event\event.evt file'
          print*,'for correct keyword and/or file name'
          print*,'for the rel file'
          print*
          stop 'Program aborted in read_resv_ef @91'
        endif
      endif  
      if(IsFileTypeTB0(fln(flnNum))) then
d       print*,'reading flow header'      
!       Open the file
        open(unit=unitNum,file=fln(flnNum),status='old',iostat=ios)
        if(ios.ne.0)then
           print*,'Problems opening ',fln(flnNum)(1:40)
           print*
           STOP ' Stopped in read_res_ef @ 100'
        endif
        noresv=1   ! assume there is at least one reservoir if there is a file
      else
         print*,'Old format .rel files not accepted'
         print*,'Please create yyyymmdd_str.tb0 files & rerun'
         stop 'Program aborted in read_resrl_ef @ 124'
      endif

      if(noresv.ge.1)then

! Initialize default values
        CALL InitResvParam(header)    

d        if(iopt.eq.2)print*,'in read_resv_ef passed 76'

! Search for and read tb0 file header
        line(1:1) = '#'
        foundEndHeader = .false.
        insideColMetaData = .false.

        do WHILE((.NOT.foundEndHeader) .AND.
     &        ((line(1:1) .eq. '#') .OR.
     &        (line(1:1) .eq. ':') .OR.
     &        (LEN_TRIM(line) .eq. 0)))   

          read(UNIT=unitNum, FMT='((A))', iostat=ios) line    ! read a line
          if(ios .eq. -1)then
              write(6,'((A))') 'ERROR: Premature EndOfFile encountered'
              STOP ' Stopped in read_resv_ef'
          end if

          rStat = Detab(line)             ! replace tabs with spaces
          line = ADJUSTL(line)        ! Get rid of leading white space
          lineLen = LEN_TRIM(line)        ! Find the length excluding trailing spaces

          if(line(1:1) .eq. ':')then
              wordCount = SplitLine(line, keyword, subString) ! find the keyword
              rStat = ToLowerCase(keyword)
              KeyLen = LEN_TRIM(keyword)

              if(keyword(1:KeyLen) .eq. ':endheader')then
                  foundEndHeader = .TRUE.

              else if(keyword(1:KeyLen) .eq. ':columnmetadata')then
                  insideColMetaData = .TRUE.
              else if(keyword(1:KeyLen) .eq. ':endcolumnmetadata')then
                  insideColMetaData = .FALSE.
              else if(insideColMetaData) then
                  iStat = ParseResvColumnMetaData(colHeader,keyword,
     &                                                keyLen,subString)
                  if(iStat .lt. 0) then
                      write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
                      write(*,'(2(A))') '   in line: ',line                   
                      STOP ' Stopped in read_resv_ef'
                      return
                  endif
              else
                  iStat= ParseResvParam(header,keyword,keyLen,subString)
                  if(iStat .lt. 0) then
                      write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
                      write(*,'(2(A))') '   in line: ',line                   
                      STOP ' Stopped in read_resv_ef'
                      return
                  else if(iStat .eq. 0) then
!                     write(*,'((A), (A))')  'Unrecognized keyword line: ',
!     &                                       line
                  endif
              end if
          end if
        end do    
!***************************************
!       Finished reading the header
!***************************************

d        if(iopt.eq.2)print*,'in read_resv_ef passed 132'


! Assign the variables from the types
        ktr =  header%tb0p%deltaT    !data time step in hours
        factor = header%tb0p%unitConv ! conversion to cms
        noresv = colHeader%tb0cmd%colCount !no of reservoirs

!     nrel    =     no of hours of data
! Scan lines of data
        rewind unitNum

!        added *ktr     Nov. 27/06  nk
        nrel = CountDataLinesAfterHeader(unitNum)*ktr


        rewind unitNum
        CALL GoToStartOfData(unitNum)

      endif

c     print*,'in read_resv_ef',nl

      if(iopt.ge.1)then
        do n=1,naa
          if(ireach(n).gt.0)then
!           rev. 9.5.57  Apr.  13/09  - NK: added ntrlflg for natural lake flows
            if(ireach(n).gt.noresv.and.ntrlflg.ne.'y')then
!           this check is bypasses when running natural flows
c           if(ireach(n).gt.noresv)then
              print*,' There are more reservoirs in the map & shed file'
              print*,' than there are outlet locations in the .rel file'
              print*,' Please fix either file so the number will match'
              print*,' No of reaches found = at least',ireach(n)
              print*,' No of reservours in rel file =',noresv
	        print*,' grid #',n,' ireach(n)=',ireach(n)
              print*
c            stop ' Aborted in read_resv_ef @207'
            endif
          endif
        end do
      endif

c      endif  ! (id.eq.1.and.ntrlflg.eq.'y'.or.ntrlflg.ne.'y')

!       allocate stuff      
      if(firstpass.eq.'y')then
          firstpass='n'
          noresv_firstpass=noresv
!         nl comes from the .str file and is the # hour of the event

          nrel_max=max0(nrel,nl)   !?????????????????????????????????????????
          
          if(nrel.gt.nl)then
            print*,'rel file longer than the str file'
            print*,'data past ',nl,' hours ignored'
          endif
!         but we need to provide enough memory to simulate a whole event
!         sometimes users specify the duration in the rel to be just 1 hr.
!         when a rule is given. However, we need memory of all the variables
!     rev. 9.5.59  Jul.  26/09  - NK: added fpet_lake for each lake in ill file
          if(noresv.gt.0)then
            allocate(b1(noresv),b2(noresv),b3(noresv),
     *      b4(noresv),b5(noresv),b6(noresv),b7(noresv),
     *      ires(noresv),jres(noresv),yres(noresv),xres(noresv),      !nk 18/06/04
     *      res_next(noresv),lake_tmp(noresv),
     *      resname(noresv),fpet_lake(noresv),lakeEflg(noresv),
     *                          store_dead(noresv),stat=iAllocate)
            if(iAllocate.ne.0) STOP
     *       'Error with allocation in read_resv_ef @205'

!     rev. 10.1.65 Jan.  28/17  - NK: Fixed allocate lake_elv from read_flow 
!     In lst we may need an extra reservoir location to look after 
!     reservoirs not in the watershed - lake_elv(noresv+1,nl)
      
            allocate(qrel(noresv,nrel_max),qdwprmd(noresv),
     *      qdwpr(noresv,nl),q24(noresv),lake_elv(noresv+1,nl),
     *      lake_stor(noresv,nl),lake_outflow(noresv,nl), 
     *      del_stor(noresv,nl),lake_inflow(noresv,nl),
     *      net_lake_inflow(noresv,nl),
     *      net_lake_outflow(noresv,nl),
     *      qstream_sum(noresv,nl),strloss_sum(noresv,nl),
     *      qmin(noresv),safe_max(noresv),DecayT(noresv),
     *                     stat=iAllocate)
            if(iAllocate.ne.0) STOP
     *       'Error with allocation in read_resv_ef @217'

!     rev. 9.5.60  Sep.  01/09  - NK: added deltat_report for lake_sd.csv file
            allocate(temp_flow_sum(noresv),stat=iAllocate)
            if(iAllocate.ne.0) STOP
     *       'Error with allocation - temp_flow_sum read_resv_ef @233'
     
!     rev. 9.9.38  Nov.  12/14  - NK: Added LKdepth to ill file
!     rev. 10.1.96 Sep   11/17  - NK: Added variable lake depth calculation lake_elv()-LKinvert()
c                  if(.not.allocated(LKdepth))then
            allocate(LKdepth(noresv),LKinvert(noresv),stat=iAllocate)
	      if(iAllocate.ne.0)STOP 
     *      	      'Error allocating LKdepth in read_resv.for'  
!     rev. 9.9.65  Apr.  03/15  - NK: Added rule s/r; resrl\rules.txt & ruleflg
            allocate(resRuleFlg(noresv),stat=iAllocate)
	      if(iAllocate.ne.0)STOP 
     *      	      'Error allocating resRuleFlg in read_resv.for'  
c                  endif

!     rev. 
            allocate(fpetLakeOverride(noresv),stat=iAllocate)
	      if(iAllocate.ne.0)STOP 
     *      	      'Error allocating fpetLakeOverride in read_resv.for'  
c                  endif

c          else
c            allocate(b1(1),b2(1),b3(1),b4(1),b5(1),b4(1),b5(1),
c     *      b6(1),b7(1),ires(1),jres(1),
c     *      resname(1),qrel(1,1),qdwpr(1,1),stat=iAllocate)
c            if(iAllocate.ne.0) STOP
c     *       'Error with allocation in read_resv_ef @178'
          endif

!         initialize lake inflows so nbs.tb0 looks ok
          do l=1,noresv
            do j=1,nl
              lake_inflow(l,j)=-1.0
!             assume there are rules for the reservoirs
!             but if no rules are found in rules.f these will be changed to .false.              
              resRuleFlg=.true.
            end do
            lakeEflg=.false.   ! no lake evap in first time step to be safe
          end do

          if(iopt.eq.2)print*,' after area5 allocate in read_resv_ef'

!     TS - ALLOCATION OF AREA6A ARRAYS (REMAINDER)
          allocate(inbsnflg(no+noresv),stat=iAllocate)
          if(iAllocate.ne.0) STOP              ! mod 03/05/02 nk
     *    'Error with allocation of area6a arrays in read_resv_ef @ 302'

!     TS: FOR ARRAY DELTA
          nnch=max0(1,nch+1)   ! added the +1 Feb. 25/08 -nk-

d          if(iopt.eq.2)print*,'In read_resv_ef @ 276'

!     TS - ALLOCATION OF AREA10A ARRAYS (REMAINDER)



!     rev. 9.5.68  Oct.  07/09  - NK: debugged read_resvin_ef.f
!         if additional allocation is needed for reservoirs, this is done in 
!         read_resvin.f, which is called for the first time after this s/r

          allocate(
c     *    delta(nnch,no+noresvi),     ! fixed bug 24/07/00
c     *    qpeakh(no+noresvi),qpeaks(no+noresvi),
c     *    dpeakh(no+noresvi),dpeaks(no+noresvi),
c     *    volsyn(no+noresvi),volhyd(no+noresvi),
     *    delta(nnch,no),     ! fixed bug 24/07/00
     *    qpeakh(no),qpeaks(no),
     *    dpeakh(no),dpeaks(no),
     *    volsyn(no),volhyd(no),
     *    stat=iAllocate)
          if(iAllocate.ne.0) STOP
     *       'Error with allocation of area10a arrays in read_resv_ef'

d          if(iopt.eq.2)print*,'In read_resv_ef @ 299'

!         initialize ERROR FUNCTIONS = 0
          do i=1,ni
            do l=1,no+noresvi
              delta(i,l)=0.0
            end do
          end do
d          if(iopt.eq.2)print*,'In read_resv_ef @ 307'

!         ASSIGN THE GAUGE LOCATIONS NAMES AND FUNCTION COEFFICIENTS
!         this is done only during the first pass if coefficient values 
!         are set to -1 for subsequent events. This makes tweaking easy
!         as only the values in the first event need to be adjusted.
          do i=1,noresv
            resname(i) = colHeader%tb0cmd%colName(i) ! reservoir name
            xres(i) = colHeader%tb0cmd%colLocX(i) ! x coordinate
            yres(i) = colHeader%tb0cmd%colLocY(i) ! y coordinate
            b1(i) = colHeader%colCoeff1(i) ! coefficient 1
            b2(i) = colHeader%colCoeff2(i) ! coefficient 2
            b3(i) = colHeader%colCoeff3(i) ! coefficient 3
            b4(i) = colHeader%colCoeff4(i) ! coefficient 4
            b5(i) = colHeader%colCoeff5(i) ! coefficient 5

          end do

!     rev. 10.1.82 May   09/17  - NK: Added reservoir_fudge_factors.csv
          INQUIRE(FILE='resrl\reservoir_fudge_factors.csv',EXIST=exists)
          if(exists)then
            open(unit=99,file='resrl\reservoir_fudge_factors.csv',
     *                   status='old',iostat=ios)
            if(ios.ne.0)then
                print*,'Error:'
                print*,'file resrl\reservoir_fudge_factors.csv exists'
                print*,'but can not be opened.'
                print*,'Open elsewhere ??'
                stop 'Program abortd in tread_resv @ 384'
            endif      
            read(99,*)
            read(99,*)
            read(99,*)
            j=0
            allocate(ff(100,noresv),stat=iallocate)
              if(iallocate.ne.0)then
                  print*,'problem allocating ff in read_resv'
                  stop 'Program aborted in readresv @ 398'
              endif
            do while(.not.eof(99))
              j=j+1  
              if(j.gt.100)then
                  print*,'Hey: you`ve gone over the limit '
                  print*,'with the # of lines in '
                  print*,'reservoir_fudge_factors.csv'
                  print*,'Max = 100'
                  stop 'Program aborted in read_resv @ 407'
              endif
              read(99,*,iostat=ios)anno(j),(ff(j,l),l=1,noresv)
              if(ios.ne.0)then
                print*,'Problems reading reservoir_fudge_factors.csv'
c                stop 'Program aborted in readresv @ 405'
              endif
              fudge=.true.
              write(*,99999)anno(j),(ff(j,l),l=1,noresv)
99999         format(i5,<noresv>f8.3)              
            end do
            print*,year1,anno(1)
            if(year1.ne.anno(1))then
                print*,'Error:'
                print*,'The first year in '
                print*,'resrl\reservoir_fudge_factors.csv'
                print*,'does not match the first year of the simulation'
                stop 'Progam aborted in read_resv @ 423'
            endif
      
          endif  ! exist:    reservoir_fudge_factors.csv
      
      else      !   firstpass
!         subsequent passes
!         check to see memory allocated is adequate     

!     rev. 9.5.57  Apr.  13/09  - NK: added ntrlflg for natural lake flows
        if(ntrlflg.ne.'y')then
 
d          if(iopt.eq.2)print*,'In read_resv_ef @ 364'
c          if(noresv.ne.noresv_firstpass.and.ntrlflg.ne.'y')then
          if(noresv.ne.noresv_firstpass)then
            print*,'No of reservoirs has been changed in'
            print*,'in file ',fln(7)
            print*,'This is not permitted'
            print*
            stop 'Program aborted in read_resvo_ef @ 264'
          endif
d          if(iopt.eq.2)print*,'In read_resv_ef @ 266'

!     rev. 9.5.57  Apr.  13/09  - NK: added ntrlflg for natural lake flows
c          if(ntrlflg.eq.'n')then

!     rev. 9.5.41  Oct.  22/08  - NK: read in reservoir coefficients each event
!     this is because their mode of operation might have changed such as
!     natural to regulated or vice versa. 
!     but not if ntrlflg='y' when we use only the coefficients in the first event
            do i=1,noresv
              b1(i) = colHeader%colCoeff1(i) ! coefficient 1
              b2(i) = colHeader%colCoeff2(i) ! coefficient 2
              b3(i) = colHeader%colCoeff3(i) ! coefficient 3
              b4(i) = colHeader%colCoeff4(i) ! coefficient 4
              b5(i) = colHeader%colCoeff5(i) ! coefficient 5
            end do
c         endif

      endif   ! firstpass

!     rev. 9.5.57  Apr.  13/09  - NK: added ntrlflg for natural lake flows
!     these reallocations are needed 
      if(nl.gt.nrel_max)then
          nrel_max=max0(nrel,nl,nl_max)
!         the event is longer than any of the previous events so 
!         more memory has to be allocated
!         DEALLOCATION OF ARRAYS FROM AREA10A:
d          if(iopt.eq.2)print*,'in read_resv_ef @ 213'

!         DEALLOCATION OF ARRAYS FROM AREA5A:
          deallocate(lake_elv,stat=iDeallocate)     
          if (iDeallocate.ne.0) print*,    
     *        'Error with deallocation of lake_elv arrays'

!         DEALLOCATION OF ARRAYS FROM AREA5A:
          deallocate(qrel,qdwpr,  !lake_elv,
     *    lake_stor,lake_inflow,lake_outflow,del_stor,net_lake_inflow,
     *    net_lake_outflow,qstream_sum,strloss_sum,stat=iDeallocate)     
          if (iDeallocate.ne.0) print*,    
     *        'Error with deallocation of area5a arrays'

!         re-allocate for larger arrays
          allocate(qrel(noresv,nrel_max), 
     *    qdwpr(noresv,nrel_max),lake_elv(noresv+1,nrel_max),
     *    lake_stor(noresv,nrel_max),lake_outflow(noresv,nrel_max),
     *    del_stor(noresv,nrel_max),net_lake_inflow(noresv,nrel_max),
     *    net_lake_outflow(noresv,nrel_max),
     *    lake_inflow(noresv,nrel_max),
     *    qstream_sum(noresv,nrel_max),strloss_sum(noresv,nrel_max),
     *                           stat=iAllocate)
          if(iAllocate.ne.0) STOP
     *       'Error with allocation of area10a arrays in read_resv'

      if(iopt.ge.1)
     *    print*,'qdwpr has been reallocated with ',noresv,nrel_max

!         re-initialize lake inflows so nbs.tb0 looks ok
          do i=1,noresv
            do j=1,nl
              lake_inflow(i,j)=-1.0
            end do
          end do

!     rev. 9.5.23  Mar.  12/08  - NK: fixed allocation error in read_resv_ef
      endif

!     rev. 9.5.57  Apr.  13/09  - NK: added ntrlflg for natural lake flows
        if(ntrlflg.eq.'y')return

!       REASSIGN THE GAUGE LOCATIONS NAMES AND FUNCTION COEFFICIENTS
!       This part is used only if coefficient values area >0
!       Used if values change over time & need to be reassigned.
        do i=1,noresv
          resname(i) = colHeader%tb0cmd%colName(i) ! reservoir name
          xres(i) = colHeader%tb0cmd%colLocX(i) ! x coordinate
          yres(i) = colHeader%tb0cmd%colLocY(i) ! y coordinate
!         if -ve values are entered for all entries for one lake
!         only values in the first event are used
!         This makes it a lot easier to tweak as only one file
!         needs to be edited.
          if(colHeader%colCoeff1(i).gt.0.0.or.
     *      colHeader%colCoeff2(i).gt.0.0.or.
     *      colHeader%colCoeff3(i).gt.0.0.or.
     *      colHeader%colCoeff4(i).gt.0.0.or.
     *      colHeader%colCoeff5(i).gt.0.0)then
            b1(i) = colHeader%colCoeff1(i) ! coefficient 1
            b2(i) = colHeader%colCoeff2(i) ! coefficient 2
            b3(i) = colHeader%colCoeff3(i) ! coefficient 3
            b4(i) = colHeader%colCoeff4(i) ! coefficient 4
            b5(i) = colHeader%colCoeff5(i) ! coefficient 5
c            if(allocated(colHeader%colCoeff6))then
c               b6(i) = colHeader%colCoeff6(i) ! coefficient 6
c            else
c              b6(i)=0.0
c            endif
c            if(allocated(colHeader%colCoeff7))then
c             b7(i) = colHeader%colCoeff7(i) ! coefficient 7
c            else
c              b7(i)=0.0
c            endif

          endif
        end do

      endif      ! firstpass
!       rev. 9.1.69  Dec.  19/04  - NK: rewrote rdresv c/w memory allocation 

      if(iopt.eq.2)print*,'In read_resv_ef @ 295'

      if(noresv.eq.0)return   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<-----

      deallocate(colHeader%tb0cmd%colName)
      deallocate(colHeader%tb0cmd%colLocX)
      deallocate(colHeader%tb0cmd%colLocY)
      deallocate(colHeader%colCoeff1)
      deallocate(colHeader%colCoeff2)
      deallocate(colHeader%colCoeff3)
      deallocate(colHeader%colCoeff4)
      deallocate(colHeader%colCoeff5)
      if(allocated(colHeader%colCoeff6))then
        deallocate(colHeader%colCoeff6)
      endif
      if(allocated(colHeader%colCoeff7))then
        deallocate(colHeader%colCoeff7)
      endif

d      if(iopt.eq.2)print*,'in read_resv_ef passed 274'


      if(noresv.gt.0)then
!         FIND THE LOCAL COORDINATES FOR THE RESERVOIRS
!         THE VALUES FOR IRES AND JRES ARE CHANGED

          do i=1,noresv
!           convert to local coordinate unit system for new .res file
!     rev. 9.2.26  Dec.  23/05  - NK: Fixed reservoir outlet location bug 
            jres(i)=int((xres(i)-xorigin)/xdelta)+1
!     rev. 9.5.36  Oct.  01/08  - NK: fixed ires bug for unevent dx & dy in read_resv
            ires(i)=int((yres(i)-yorigin)/ydelta)+1
          end do
          if(iopt.ge.1)then
            write(53,*)
            write(53,*)'In read_resv_ef:'
            write(53,*)'Note: set iopt = 0 to not write this.'
            write(53,1011)
            write(53,1013)(i,ires(i),jres(i),b1(i),b2(i),b3(i),
     *                    b4(i),b5(i),resname(i),i=1,noresv)
          endif

!         THE ORDER OF READING THE COORDINATES OF THE RESERVOIRS
!         MUST BE THE SAME AS READING THE CORRESPONDING FLOWS
!         IN S/R REROUT.
!         READ RELEASES
!         THE RESERVOIR OUTFLOWS ARE ALL READ IN THE FIRST TIME
!         REROUT IS CALLED. THEY ARE THEN STORED AND USED EACH TIME
!         REROUT IS CALLED.
!         IF NATURAL CONTROL, FLOWS ARE SET TO -1.0 ABOVE

!         initialize releases
          do k=1,noresv
            do j=1,nrel
!     REV. 10.1.38 Jul   28/16  - NK: Added noDataValue to WFO & tb0 files
              qrel(k,j)=-999.0
            end do
          end do

d         print*,'# of reservoirs in the rel file =',noresv
d         print*,'# releases read in: nrel =',nrel

c          if(b1(1).eq.0.0)then

!         Read the whole record and then deal with missing data

          if(nrel.gt.0)then
!           case with reservoir releases
!     rev. 9.1.13  Mar.  23/02  - fixed resv. timing, moved to beginning of dt
!           so the releases are constant for the day if ktr= 24 say
!           do j=ktr,nrel,ktr
            do j=1,nrel,ktr
              read(unitNum,*,iostat=ios)(qrel(k,j),k=1,noresv)
c              write(*,*)(qrel(k,j),k=1,noresv)
d              if(iopt.eq.2)print*,j,(qrel(k,j),k=1,noresv)
                if(ios.ne.0)then
                  write(98,*)' Error on unit=37,fln=',fln(7)
                  write(98,*)' Trying to read releases hour =',j
                  print*,' Error on unit=37,fln=',fln(7)
                  print*,' Trying to read releases'
                  print*,' ios= ',ios
                  if(ios.eq.-1)then
                    write(98,*)'End of file in fln= ',fln(7)
                    write(98,*)
     *                   'Possibly last line does not have a return'
                    print*,'End of file in fln= ',fln(7)
                    print*,'Possibly last line does not have a return'
                    print*
                  else
                    print*
                    STOP ' program aborted in read_resv_ef.for'
                  endif
                endif
            end do

!     rev. 10.1.82 May   09/17  - NK: Added reservoir_fudge_factors.csv
            if(fudge)then
            
!           Adjust faulty reservoir outflows 
              do l=1,noresv
                do j=1,nrel,ktr
                  if(qrel(l,j).gt.0.0)qrel(l,j)=qrel(l,j)*ff(id,l)
                end do
              end do      
            endif  

!     rev. 9.5.29  May.  26/08  - NK: fixed initialization in read_resv_ef
!     rev. 9.9.56  Feb.  04/15  - NK: Fixed missing initial rel data in read_resv
!           split this do loop into 2: first read all the data & then fix mising data
!           read it:
            if(dds_flag.eq.0)then
            do k=1,noresv
              if(qrel(k,1).lt.0.0.and.b1(k).le.0.0)then
                if(routeflg.eq.'q')then
!                 do nothing                  
                else
!                 Don't bother with this msg for DDS or when writing '
!                 a file for the flow1D            
                  print*
                  print*,'WARNING <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
                  print*,'b1(1) = 0.0 and first rel < 0.0 for '
                  print*,'reservoir #',k
                  print*,'first release should not be < 0.0'
                  print*,'releases set to first +ve release found'
                endif
c                qrel(k,ktr)=0.0
                j=1
                do while(qrel(k,j).lt.0.0.and.j.le.nrel-ktr)
d                 print*,'resv#,hour,q_rel',k,j,qrel(k,j)
                  j=j+ktr
                  if(j.gt.nrel)then
                    print*,'WARNING <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
                    print*,'no release flows found for reservoir #',k
!     rev. 9.9.63  Apr.  06/15  - NK: Changed reas_resv to carry on with lask known release(s)
c                    stop 'Program aborted in read_resv @ 621'
                  endif
                end do
                miss_count=j
                
d               print*,'Back-filled releases:',miss_count/ktr                
                do j=1,miss_count,ktr
                  qrel(k,j)=qrel(k,miss_count)
d                 print*,'resv#,hour,q_rel',k,j,qrel(k,j)
                end do
                if(routeflg.ne.'q')then
                  print*,'WARNING <<<<<<<<<<<<<<<<<<<<<'
                  print*,'Initial ',(miss_count-1)/ktr,'-ve flows in' 
                  print*,'the rel file for reservoir #',k 
                  print*,'replaced by first +ve value found'
                  print*,'or the last know value in a previous event',
     *                                       qrel(k,miss_count)
                  print*
                endif
              endif
            end do
            endif   ! ddsflg
          endif
      
!         fix it:
          if(nrel.gt.0)then
!           case with reservoir releases
!     rev. 9.1.13  Mar.  23/02  - fixed resv. timing, moved to beginning of dt
!           so the releases are constant for the day if ktr= 24 say
!           do j=ktr,nrel,ktr
            do j=1,nrel,ktr

!     rev. 9.5.24  Mar.  18/08  - NK: fixed missing data in read_resl_ef.f
!     rev. 9.9.56  Feb.  04/15  - NK: Fixed missing initial rel data in read_resv
!               fill in the gaps to create hourly data
                if(ktr.gt.1)then
                  do jj=j+1,j+ktr-1
                    do k=1,noresv
                      qrel(k,jj)=qrel(k,jj-1)
                    end do
                  end do
                endif
!               fill in missing data (-ve data)
                do k=1,noresv
!     rev. 9.5.29  May.  26/08  - NK: fixed initialization in read_resv_ef
                  if(qrel(k,j).lt.0.0.and.b1(k).le.0.0)then
                    do jj=j,j+ktr-1
                      if(jj.gt.1)then                  
                        qrel(k,jj)=qrel(k,jj-1)
                      endif  
                    end do
                  endif
                end do
              end do
            endif     !  if(nrel.gt.0)


d          if(iopt.ge.1)then
d            write(53,*)'In read_resv_ef'
d            write(53,*)'Reservoir releases this event:'
cd            j=ktr
d            write(53,53001)(resname(k)(1:9),k=1,noresv)
53001        format(<noresv>a10)
d            do j=1,nrel,ktr
d              write(53,53000)(qrel(k,j),k=1,noresv)
53000          format(<noresv>f10.3)
d            end do
d          endif
      endif !         if(noresv.gt.0)

d      if(iopt.eq.2)print*,'in read_resv_ef passed 187'

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 
      istep=al/1000
d      if(iopt.eq.2)print*,'in read_resv_ef passed 272'

      if(noresv.gt.0)then                             !999999999999999
!       FIND THE LOCAL COORDINATES FOR THE RESERVOIRS
!       THE VALUES FOR IRES AND JRES ARE CHANGED

!       THE ORDER OF READING THE COORDINATES OF THE RESERVOIRS
!       MUST BE THE SAME AS READING THE CORRESPONDING FLOWS
!       IN S/R rdresv.

!       READ RELEASES
!       THE RESERVOIR OUTFLOWS ARE ALL READ IN THE FIRST TIME
!       rdresv IS CALLED. THEY ARE THEN STORED AND USED EACH TIME
!       rdresv IS CALLED.

!       IF NATURAL CONTROL, FLOWS ARE SET TO -1.0 ABOVE
d        if(iopt.eq.2)print*,'in read_resv_ef passed 304'

      endif           !   if(noresv.gt.0)

d      if(iopt.eq.2)print*,'in read_resv_ef passed 551'

!     WE CAN'T HAVE -VE RELEASES WHEN WE START
!     rev. 10.2.18 Mar.  12/18  - NK: Fixed array fault in read_resv_ef
      if(nrel.gt.0)then
        do k=1,noresv
          if(qrel(k,1).lt.0.0)qrel(k,ktr)=0.000001
        end do
      endif

d      if(iopt.eq.2)print*,'in read_resv_ef passed 555'
d      print*,'no releases read in  nrel =',nrel,' nl =',nl

!     SET FUTURE RELEASES = LAST NON-NEGATIVE RELEASE
!     REALLY, WE'RE WORKING IN HOURLY INTERVALS ALTHOUGH THE      
!     RELEASES MAY BE READ IN ONLY WHEN THE RES OUTFLOW IS CHANGED.

!     rev. 9.5.14  Feb.  26/08  - NK: padded rel file for missing data
!     fill in missing data if rel file is shorter than the str file
!     nl is the length in nrs of the str file
!     nrel is the length of the rel file
!     added Feb. 26/08  -nk-


      if(nrel.gt.0)then
!       fill data at end of file only if there are releases
        if(nrel.lt.nl)then
          do j=nrel+1,nl
            do k=1,noresv
              qrel(k,j)=qrel(k,j-1)
            end do
          end do
          if(numa.eq.0.and.dds_flag.eq.0)then
            print*
            print*,'WARNING:'
            print*,'reservoir release file is shorter than the str file'
            print*,'missing data set = to last recorder release'
            print*
          endif
          do j=2,nl
            do k=1,noresv
              if(qrel(k,j).lt.0.0)qrel(k,j)=qrel(k,j-1)
            end do
          end do
        endif
      endif                        

d      if(iopt.eq.2)print*,'in read_resv_ef passed 683'


d      if(iopt.ge.2.and.iopt.le.10)then
d        do j=1,nl
d          write(53,6801)k,nrel
d          write(53,*)j,(qrel(k,j),k=1,noresv)
d        end do   
d      endif

d      if(iopt.eq.2)print*,'in read_resv_ef passed 678'

!check to see that this does not have to go to flowint !!!!!!!!!!1
!     REVISED JAN 17/96

      if(id.eq.1)then
        do k=1,noresv
!         rev. 9.1.72  Dec.  28/04  - NK: fix bug in rdresv setting reach # 
!     rev. 9.5.81  Jan.  16/09  - NK: allow reservoirs outside watershed in resv file
            
c          if(inbsnflg(no+k).eq.1)then    !oops - inbsnflg not defined yet NK
            i=ires(k)
            j=jres(k)
            if(i.ge.1.and.i.le.ycount.and.j.ge.1.and.j.le.xcount)then
              n=s(i,j)
            else
              n=0
            endif
            
c          if(inbsnflg(no+k).eq.1)then      !inbsnflg not defined yet!!!!
!     v. 8.99n  Dec. 31/2001-     fixed nat. res initial flow (JW)
!           this was fixed in spl8 on dec. 31/01
            if(b1(k).le.0.0)then
              if(n.eq.0)then
                print*,'reservoir no',k,' at row',i,' column',j 
                print*,'is not in the basin or'
                print*,'is in a grid with zero area: please check'
                print*
c                STOP 'program aborted in read_resv_ef @ 718'
              else
                qo1(n)=qrel(k,1)
                qo2(n)=qrel(k,1)
              endif
            endif
!           FOR RESERVOIRS WITH RULES, IREACH IS READ IN FOR THE 
!           DWOPER FORMAT BUT WHEN NO DWOPER REACHES ARE SPECIFIED
!           IREACH HAS TO BE DEFINED HERE.              
!           read_resv_ef WON'T BE CALLED AGAIN IF IREACH = 0 
c          endif   ! inbsnflg
!         rev. 9.1.72  Dec.  28/04  - NK: fix bug in rdresv setting reach # 
          IF(N.EQ.0.and.dds_flag.eq.0)THEN
            print*
            print*,'WARNING'
            print*,'For reservoir/lake no ',k
c            print*,'the inbsnflg = ',inbsnflg(no+k)
            print*,'for row = ',i,' & col = ',j
            print*,'for lat = ',yres(k),' & long = ',xres(k)
            print*,'but the rank seems to be ',n
            print*,'i.e. the grid is not part of the watershed'
            print*,'and the reach # can not be assigned'
            print*,'OK if running a sub-watershed'
c            stop 'Program aborted in read_resv @ l 736'
          endif
          if(n.gt.0)then
            ireach(n)=k
          endif
        end do
      endif

d      if(iopt.eq.2)print*,'in read_resv_ef passed 612'

      close(unit=unitNum,status='keep')


!     REV. 10.1.17 Jan.  11/16  - NK: Added fpetLakeOverride factor
      
      INQUIRE(FILE='basin\fpet_lake_override.csv',EXIST=exists) 
      if(exists)then
        open(unit=99,file='basin\fpet_lake_override.csv',
     *              status='old',iostat=ios)
        if(ios.ne.0)then
          print*,'Problems opening file basin\fpet_lake_override.csv'
          stop 'Program aborted in read_resv_ef @ 853'
        endif
        read(99,*)
        do while(.not.eof(99))
          read(99,*,iostat=ios)i,fpetLakeOverride(i)
          if(ios.ne.0)then
            print*,'Problems reading basin\fpet_lake_override.csv'
            stop 'Program aborted in read_resv_ef @ 860'
          endif
        end do
        if(i.ne.noresv)then
          print*
          print*,'Error: found basin\fpet_lake_override.csv'
          print*,'but the number of entrees does not match the no of'
          print*,'lakes and/or reservoirs i.e. # reaches'
          stop 'Program aborted in read_resv_ef.f @ 858'
        endif
        close(unit=99)
      endif
          




  999 RETURN                 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<-----

! FORMATS

  500 format(256f10.3)
 1011 format(' ',3x,'  i  ires(i) jres(i)    b1(i)     b2(i)',
     * '      b3(i)       b4(i)      b5(i)      resvname')
 1013 format(' ',3x,i3,2i8,5g12.3,a12/)
 6801 format('   read_resv_ef: reservoir no =',i3,' mhtot =',i5)
 6802 format('   ',i5,<noresv>g12.3)

      END SUBROUTINE read_resv_ef

