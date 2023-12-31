      SUBROUTINE read_divert(unitNum,flnNum)

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

!    You should have recieved a copy of the GNU Lesser General Public License
!    along with WATFLOOD.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
     
     
!*****************************************************************************
!  read_divert_ef - written Jan/09 by Nick Kouwen, UW
!     - Derived from rdresv written by Dave Watson CHC
!     This subroutine reads the diversion (div) file 
!     (tb0 format)
!*****************************************************************************

!     rev. 9.5.52  Jan.  20/09  - NK: added reading yyyymmdd_div.pt2 for diversions
!     rev. 9.5.81  Jan.  05/11  - NK: allow reservoirs outside watershed in resv file
!     rev. 9.7.24  Apr.  20/11  - NK: Added diverflg to indicate if a diversion is in grid
!     rev. 10.2.59 Aug.  17/19  - NK Convert read_flow_ef & read_divert to F90
!******************
      use area_watflood
      use area_debug

! TB0 data module
      USE EF_Module

      implicit none
      TYPE(DivParam) :: header
      TYPE(DivColumnMetaData) :: colHeader

      Integer  :: ios,j,k,i,n,l,jz,jj,ng,nt,skip
      integer  :: nodivert_firstpass,ndiv_max,iAllocate,iDeallocate,ntake,ngive
      real*4   :: factor

!     rev. 9.1.55  Jun.  12/04  - NK: write new str files to strfw\newfmt folder.
      LOGICAL exists,diversion_local,errchk
      logical       :: climateFlg
      CHARACTER(1)  :: firstpass
      integer       ::  today,startHour,julian_day
      integer       :: cdayYear,cdayMonth,cdayDay

      CHARACTER(10) :: ctime
      CHARACTER(8)  :: cday

      data firstpass/'y'/
      data factor/1.0/      ! unit conversion factor
      data climateFlg/.false./
      data errchk/.false./

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

! Parameter type definitions
      integer*4 unitNum, flnNum, iStat

! Local variables
      character*4096 line, subString, tmpString
      character*128 keyword, value 
      character*12  outlet_type
      integer lineLen, keyLen, wordCount
      logical rStat, lineType, foundEndHeader, insideColMetaData

      if(FLtype(flnNum).eq.'tb0')then
          ! Open the file
          INQUIRE(FILE=fln(flnNum),EXIST=exists)
          if(exists)then
                open(unit=unitNum,file=fln(flnNum),status='old',iostat=ios)
	          if(iopt.ge.1)print*,'opened unit=',unitNum,fln(flnNum)(1:40)
                if(ios.ne.0)then
                      print*,'Problems opening ',fln(flnNum)
                      print*
                      STOP ' Stopped in read_divert_ef'
                endif
            nodivert=1   ! assume there is at least one diversion if there is a file
          else
            if(numa.eq.0)then
              print*,'WARNING:'
              print*,'Diversion (div) file ',fln(flnNum)(1:40)
              print*,'is NOT found'
              print*,'Program continues with no diversions'
              print*
            endif
            nodivert=0
            ndiv=0
          endif

    !d      if(iopt.eq.2)print*,'in read_divert_ef passed 70'


          if(nodivert.ge.1)then


    ! Initialize default values
            CALL InitDivParam(header)   

    !d       if(iopt.ge.2)print*,'in read_divert_ef passed 76'

    ! Search for and read tb0 file header
            line(1:1) = '#'
            foundEndHeader = .false.
            insideColMetaData = .false.

            do WHILE((.NOT.foundEndHeader) .AND.((line(1:1) .eq. '#') .OR.(line(1:1) .eq. ':') .OR.(LEN_TRIM(line) .eq. 0)))     

              read(UNIT=unitNum, FMT='((A))', iostat=ios) line      ! read a line
    !c      print*,line(1:72)
              if(ios .eq. -1)then
                 write(6,'((A))') 'ERROR: Premature EndOfFile encountered'
                 STOP ' Stopped in read_divert_ef'
              end if

              rStat = Detab(line)                       ! replace tabs with spaces
              line = ADJUSTL(line)          ! Get rid of leading white space
              lineLen = LEN_TRIM(line)            ! Find the length excluding trailing spaces

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

    !             this is not needed I think as there are no coefficients
                 else if(insideColMetaData) then
                    iStat = ParseDivColumnMetaData(colHeader,keyword,keyLen,subString)
                    if(iStat .lt. 0) then
                       write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
                       write(*,'(2(A))') '   in line: ',line                             
                       STOP ' Stopped in read_divert_ef'
                       return
                    endif
                 else
                    iStat = ParseDivParam(header,keyword,keyLen,subString)
                    if(iStat .lt. 0) then
                          write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
                          write(*,'(2(A))') '   in line: ',line                             
                          STOP ' Stopped in read_divert_ef'
                          return
                    else if(iStat .eq. 0) then
    !                     write(*,'((A), (A))')  'Unrecognized keyword line: ',
    !     &                                                   line
                    endif
                 end if
              end if
            end do    
    !***************************************
    !       Finished reading the header
    !***************************************

    !d      if(iopt.ge.2)print*,'in read_divert_ef passed 132'


    ! Assign the variables from the types
            ktr =  header%tb0p%deltaT    !data time step in hours
    !     rev. 9.8.54  Apr.  02/13  - NK: deltat conversion seconds to hours
    !       In the past, deltat has bee in hours but GK wants them in seconds
    !       This converts deltat to hours as needed in spl if a large deltat is found        
            if(ktr.gt.1000)then
              ktr=ktr/3600
            endif
            factor = header%tb0p%unitConv ! conversion to cms
            if(factor.lt.abs(1.0E-24))factor=1.0  ! for missing value in the div file
            nodivert = colHeader%tb0cmd%colCount !no of reservoirs

    !     ndiv    =     no of hours of data
    ! Scan lines of data
            rewind unitNum

    !        added *ktr     Nov. 27/06  nk
            ndiv = CountDataLinesAfterHeader(unitNum)*ktr
            rewind unitNum
            CALL GoToStartOfData(unitNum)

          endif

      else    ! FLtype(flnNum).eq.'tb0'
      
!     rev. 10.2.46 Jan.  21/19  - NK: Added netCDF capability
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!            call read_ts_diversion_nc(flnNum)                                    ! netCDF
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!           FILL IN FLOW DATA 
            flowfillflg='y'
      
      endif   ! FLtype(flnNum).eq.'tb0'
      
      
!       allocate stuff      
      if(firstpass.eq.'y')then
        nodivert_firstpass=nodivert
        diversion_local=.true.
!       nl comes from the .str file and is the # hour of the event

        ndiv_max=max0(ndiv,nl)
        if(ndiv.gt.nl)then
          print*,'div file longer than the str file'
          print*,'data past ',nl,' hours ignored'
        endif

!       but we need to provide enough memory to simulate a whole event
!       sometimes users specify the duration in the rel to be just 1 hr.
!       when a rule is given. However, we need memory of all the variables
!     rev. 10.2.59 Aug.  17/19  - NK Convert read_flow_ef & read_divert to F90
        if(FLtype(flnNum).eq.'tb0')then
            if(nodivert.gt.0)then
              allocate(divertname(nodivert),&
              xsource(nodivert),ysource(nodivert),&
              xrecieve(nodivert),yrecieve(nodivert),&
              jtake(nodivert),itake(nodivert),&
              jgive(nodivert),igive(nodivert),&
              gridsource(nodivert),gridrecieve(nodivert),&
              qdivert(nodivert,ndiv),&
              divert_inbasin_flg(nodivert),&
              upstrda(nodivert),stat=iAllocate)

              if(iAllocate.ne.0)then
                print*,'Error with allocation in read_divert_ef @172'
                print*,'firstpass =',firstpass
                STOP 'Program aborted in read_divert @ 198'
              endif
    !         initialize upstrda
              do l=1,nodivert
                upstrda(l)=0.0
              end do
                    
            endif
    !     rev. 9.9.43  Nov.  26/14  - NK: Allocation for divertflg = 'g' 
            if(divertflg.eq.'g')return   ! ' only here to allocate

    !       ASSIGN THE GAUGE LOCATIONS NAMES AND FUNCTION COEFFICIENTS
    !       this is done only during the first pass if coefficient values 
    !       are set to -1 for subsequent events. This makes tweaking easy
    !       as only the values in the first event need to be adjusted.
            do n=1,nodivert
              divertname(n) = colHeader%tb0cmd%colName(n) ! diversion name
              xsource(n) = colHeader%tb0cmd%colLocX(n) ! x coordinate
              ysource(n) = colHeader%tb0cmd%colLocY(n) ! y coordinate
              xrecieve(n) = colHeader%tb0cmd%colLocX1(n) ! x coordinate
              yrecieve(n) = colHeader%tb0cmd%colLocY1(n) ! y coordinate
    !     rev. 9.9.55  Jan.  22/15  - NK: Added diversion upstream drainage area in div file
              if(allocated(colHeader%colValue1))then
    !           means upstream drainage area is present          
                upstrda(n) = colHeader%colValue1(n)     ! upstream drainage area
              else
                upstrda(n)=0.0
                print*,'Warning: upstream drainage area at diversion not found'
                print*,'Insert as value1'
              endif
            end do
        endif   ! FLtype(flnNum).eq.'tb0'
    
!     rev. 9.5.81  Jan.  16/09  - NK: allow reservoirs outside watershed in resv file
!     check to see if diversion locations are in the watershed
!     first see if it is in the grid
        do n=1,nodivert
          if(xsource(n).lt.xorigin.or.&
            xsource(n).gt.xorigin+xcount*xdelta.or.&
            ysource(n).lt.yorigin.or.&
            ysource(n).gt.yorigin+ycount*ydelta)then
            write(98,*)'WARNING: Diversion origin and/or destination not in watershed Diversion ',n,' is ignored'
          endif

!     rev. 9.9.41  Nov.  20/14  - NK: Added check if diversion = in-basin
          write(51,*) 
          write(51,*)'Check if diversion locations are in the grid'
          write(51,*)'There are ',na,' grids' 
          write(51,*)
	      divert_inbasin_flg(n)=.false.
          i=int((ysource(n)-yorigin)/ydelta)+1
	      j=int((xsource(n)-xorigin)/xdelta)+1
          
          if(FLtype(flnNum).eq.'tb0')then
	        write(51,*)'i,ycount,j,xcount',i,ycount,j,xcount
		    if(i.ge.1.and.i.le.ycount.and.j.ge.1.and.j.le.xcount)then
      	      if(s(i,j).ge.1.and.s(i,j).le.na)then
                gridsource(n)=s(i,j)
                write(51,*)'gridsource,n,i,j',gridsource(n),n,i,j,s(i,j)
              else
                gridsource(n)=-1
              endif
            else
              gridsource(n)=-1  
            endif
          else
!            for FEWS - there is no source grid - just take the outlet              
             gridsource(n)=na  
          endif   ! FLtype(flnNum).eq.'tb0'
          
!c          gridsource(n)=s(i,j)
            i=int((yrecieve(n)-yorigin)/ydelta)+1
            j=int((xrecieve(n)-xorigin)/xdelta)+1
	        if(i.ge.1.and.i.le.ycount.and.j.ge.1.and.j.le.xcount)then
!     rev. 10.5.13 May   26/23  = NK Fixed diversion outlet check:  naa -> na
      	      if(s(i,j).ge.1.and.s(i,j).le.na)then
                gridrecieve(n)=s(i,j)
              else
                gridrecieve(n)=-1
              endif
            else
              gridrecieve(n)=-1  
            endif
!          endif          ! FLtype(flnNum).eq.'tb0'
          
          write(51,*)'n,gridsource(n),gridrecieve(n)',n,gridsource(n),gridrecieve(n)
          
          if(gridsource(n).ge.1.and.gridrecieve(n).ge.1)then
            divert_inbasin_flg(n)=.true.
          else
              if(iopt99)Then
                if(gridsource(n).eq.0)Then
                  print*,'Diversion source # ',n,'is not in the domain'
                  write(98,*)'Diversion source # ',n,'is not in the domain'
                  errchk=.true.
                endif
                if(gridsource(n).eq.0)Then
                  print*,'Diversion target # ',n,'is not in the domain'
                  write(98,*)'Diversion target # ',n,'is not in the domain'
                  errchk=.true.
                  pause 'please fix'
                endif
                if(errchk)then
                    print*,'Please fix the diversion location(s) so '
                    print*,'are in the domain - i.e a model grid'
                    print*,'This diversion data will be ignored if you continue'
                    pause 'In read-divert @ 342'
                endif
              endif
          endif
         write(51,*)'n,divert_inbasin_flg(n)',n,divert_inbasin_flg(n)
         write(51,*)'gridsource ,x,y,n',n,xsource(n),ysource(n)
         write(51,*)'gridrecieve,x,y,n',n,xrecieve(n),yrecieve(n)
         write(51,*)
         write(51,*)
        end do
        
!     rev. 9.9.55  Jan.  22/15  - NK: Added diversion upstream drainage area in div file
!     REV. 10.1.40 Oct   11/16  - NK: Fixed bug in read_divert for missing u/s DA
        if(FLtype(flnNum).eq.'tb0')then
            do l=1,nodivert
              n=gridrecieve(l)
    !d         print*,'source grid =',l,' in grid # ',n          
    !         add the u/s area to the grid receiving the diversion
              if(n.ge.1)then
    !d            print*,n,da(n),upstrda(l),da(n)+upstrda(l)
                da(n)=da(n)+upstrda(l)
              endif
    !d          print*,'naa=',naa,'  in read_divert @ 309'
              i=0
    !     rev. 10.2.19 Mar.  13/18  - NK: Fixed array fault read_divert.f
              if(n.lt.naa.and.n.gt.1)then
                if(next(n).ge.1)then
                  da(next(n))=da(next(n))+upstrda(l)
    !d              print*,n,next(n),da(next(n))
                  n=next(n)
                endif
              endif
            end do
        endif   ! FLtype(flnNum).eq.'tb0'
      
        firstpass='n'

      endif      !   firstpass
      
!d      if(iopt.ge.2)print*,'End first pass in read_divert'
      
!     rev. 10.1.85 May   17/17  - NK: Level_station_location.xyx for iopt > 0 only
!     rev. 10.1.86 May   17/17  - NK: Diversion_location.xyx for iopt > 0 only
      if(iopt99)then
        open(unit=99,file='diversion_location.xyz',status='unknown')
        do n=1,nodivert
          write(99,*)xsource(n),ysource(n),-n
          write(99,*)xrecieve(n),yrecieve(n),n
        end do
        close(unit=99,status='keep')
      endif

	if(.not.diversion_local)then
!       rule: if no diversion in the first event, then not later either
	  diversion=.false.
	  print*,'WARNING  <<<<<<<<<<<<<<<<<<<<'
	  print*,'A diversion file was found but the data can not be used'
	  print*,'for some reason'
	  print*,'Possible problem: locations not in watershed grids'
	  print*,'Program continues without diversion flows'
      pause 'In read_divert @ 383'
	  print*
	  return
	endif

    if(FLtype(flnNum).eq.'tb0')then
!     subsequent passes
!     check to see memory allocated is adequate      
!d      if(iopt.eq.2)print*,'In read_divert_ef @ 258'
      if(nodivert.ne.nodivert_firstpass)then
	  print*
        print*,'No of diversions has been changed in'
        print*,'in file ',fln(7)
        print*,'This is not permitted'
        print*
        stop 'Program aborted in read_divert @ 264'
      endif
!d      if(iopt.eq.2)print*,'In read_divert_ef @ 266'
      
      
      
      
      if(nl.gt.ndiv_max.or.ndiv.gt.ndiv_max)then
        ndiv_max=max0(ndiv,nl,ndiv)
!       the event is longer than any of the previous events so 
!       more memory has to be allocated
!       DEALLOCATION OF ARRAYS FROM AREA10A:

!d       if(iopt.eq.2)print*,'in read_divert_ef @ 213'

!       DEALLOCATION OF ARRAYS FROM AREA5A:
        deallocate(qdivert,stat=iDeallocate)     
        if (iDeallocate.ne.0) print*,'Error with deallocation of area5a arrays'

!       re-allocate for larger arrays
        allocate(qdivert(nodivert,ndiv_max),stat=iAllocate)
        if(iAllocate.ne.0) STOP'Error with allocation of area10a arrays in read_divert'
!     rev. 9.5.23  Mar.  12/08  - NK: fixed allocation error in read_divert_ef
        endif

!       REASSIGN THE GAUGE LOCATIONS NAMES AND FUNCTION COEFFICIENTS
!       This part is used only if coefficient values area >0
!       Used if values change over time & need to be reassigned.
!     rev. 10.2.59 Aug.  17/19  - NK Convert read_flow_ef & read_divert to F90
        if(FLtype(flnNum).eq.'tb0')then
            do i=1,nodivert
              divertname(i) = colHeader%tb0cmd%colName(i) ! reservoir name
              xsource(i) = colHeader%tb0cmd%colLocX(i) ! x coordinate
              ysource(i) = colHeader%tb0cmd%colLocY(i) ! y coordinate
              xrecieve(i) = colHeader%tb0cmd%colLocX1(i) ! x coordinate
              yrecieve(i) = colHeader%tb0cmd%colLocY1(i) ! y coordinate
            end do
        endif

!       rev. 9.1.69  Dec.  19/04  - NK: rewrote rdresv c/w memory allocation 

!d      if(iopt.eq.2)print*,'In read_divert_ef @ 295'
      if(nodivert.eq.0)return   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<-----
!     rev. 10.2.59 Aug.  17/19  - NK Convert read_flow_ef & read_divert to F90
        if(FLtype(flnNum).eq.'tb0')then
          deallocate(colHeader%tb0cmd%colName)
          deallocate(colHeader%tb0cmd%colLocX)
          deallocate(colHeader%tb0cmd%colLocY)
          deallocate(colHeader%tb0cmd%colLocX1)
          deallocate(colHeader%tb0cmd%colLocY1)
        endif
!d      if(iopt.eq.2)print*,'in read_divert_ef passed 274'


      if(nodivert.gt.0)then
!       FIND THE LOCAL COORDINATES FOR THE RESERVOIRS
!       THE VALUES FOR idiv AND jdiv ARE CHANGED

        do i=1,nodivert
!         convert to local coordinate unit system for new .res file
          jtake(i)=int((xsource(i)-xorigin)/xdelta)+1
          itake(i)=int((ysource(i)-yorigin)/ydelta)+1
          jgive(i)=int((xrecieve(i)-xorigin)/xdelta)+1
          igive(i)=int((yrecieve(i)-yorigin)/ydelta)+1
        end do
        if(iopt.ge.1)then
          write(53,*)
          write(53,*)'In read_divert_ef:'
          write(53,*)'Note: set iopt = 0 to not write this.'
          write(53,1011)
          write(53,1013)(i,itake(i),jtake(i),igive(i),jgive(i),divertname(i),i=1,nodivert)
        endif
      endif    
    endif               ! FLtype(flnNum).eq.'tb0'

!       THE ORDER OF READING THE COORDINATES OF THE RESERVOIRS
!       MUST BE THE SAME AS READING THE CORRESPONDING FLOWS
!       IN S/R REROUT.
!       READ RELEASES
!       THE RESERVOIR OUTFLOWS ARE ALL READ IN THE FIRST TIME
!       REROUT IS CALLED. THEY ARE THEN STORED AND USED EACH TIME
!       REROUT IS CALLED.
!       IF NATURAL CONTROL, FLOWS ARE SET TO -1.0 ABOVE

!       initialize releases
        if(FLtype(flnNum).eq.'tb0')then
            do k=1,nodivert
              do j=1,ndiv
                qdivert(k,j)=-1.0
              end do
            end do
        endif
        
        if(ndiv.gt.0)then
!         case with reservoir releases
!         do j=ktr,ndiv,ktr

!        Assumes diversion if for a whole year
        if(fln(flnNum)(1:13).eq.'diver\climate')then
            skip=(jul_day_now)*24/ktr            
        endif
        print*,'lines skipped in div file: ',skip,skip*ktr
            
!     rev. 10.2.59 Aug.  17/19  - NK Convert read_flow_ef & read_divert to F90
!     rev. 9.1.13  Mar.  23/02  - fixed resv. timing, moved to beginning of dt
              if(FLtype(flnNum).eq.'tb0')then
                if(skip.gt.0)then
                  do j=1,skip
                    read(unitNum,*,iostat=ios)
                  end do
                endif
                do j=ktr,ndiv-skip*ktr,ktr
                  read(unitNum,*,iostat=ios)(qdivert(k,j),k=1,nodivert)
!                  write(*,*,iostat=ios)j,(qdivert(k,j),k=1,nodivert)
                  If(factor.ne.1.0)then
                      do k=1,nodivert
                          qdivert(k,j)=factor*qdivert(k,j)
                      end do
!                      print*,'diversion flows converted'
                  endif
                     
!    write(866,*)year_now,j            
!    write(866,*)year_now,j,qdivert(1,j)            
    !d            if(iopt.eq.2)print*,j,(qdivert(k,j),k=1,nodivert)
                  if(ios.ne.0)then
                    write(98,*)' Error on unit=',unitNum,' fln=',fln(flnNum)
                    write(98,*)' Trying to read diversions hour =',j
                    print*,' Error on unit=',unitNum,' fln=',fln(flnNum)
                    print*,' Trying to read diversions'
                    print*,' ios= ',ios
                    if(ios.eq.-1)then
                      write(98,*)'End of file in fln= ',fln(flnNum)
                      write(98,*)'Possibly last line does not have a return'
                      print*,'End of file in fln= ',fln(flnNum)
                      print*,'Possibly last line does not have a return'
                      print*
                    else
                      print*
                      STOP ' program aborted in read_divert_ef.for'
                    endif
                  endif
    !     rev. 9.5.24  Mar.  18/08  - NK: fixed missing data in read_resl_ef.f
    !             fill in the gaps to hourly data
                  if(ktr.gt.1)then
                    do k=1,nodivert
	                  do jj=ktr-1,1,-1
                          qdivert(k,j-jj)=qdivert(k,j)
!                          write(886,*)year_now,month_now,day_now,jj,qdivert(k,j)
	                  end do
                    end do
                  endif

    !             fill in missing data (-ve data)
                  do k=1,nodivert
    !     rev. 9.5.29  May.  26/08  - NK: fixed initialization in read_divert_ef
                    if(qdivert(k,ktr).lt.0.0)then
                      print*
                      print*,'WARNING <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
                      print*,'first diversion flow should not be < 0.0'
                      print*,'diversions set to 0.0 until +ve value found'
                      qdivert(k,ktr)=0.0
                    endif
    !     rev. 9.5.29  May.  26/08  - NK: fixed initialization in read_divert_ef
                    if(qdivert(k,j).lt.0.0)then
                      do jj=j,j+ktr-1
                        qdivert(k,jj)=qdivert(k,jj-1)
                      end do
                    endif
                  end do
                end do
              else              ! FLtype(flnNum).eq.'tb0'
              
!     rev. 10.2.59 Aug.  17/19  - NK Convert read_flow_ef & read_divert to F90
!     rev. 9.5.24  Mar.  18/08  - NK: fixed missing data in read_resl_ef.f
!             fill in the gaps to hourly data

                if(ktr.gt.1)then
                  do k=1,nodivert
                    do j=ktr,ndiv,ktr
	                do jj=ktr-1,1,-1
                        qdivert(k,j-jj)=qdivert(k,j)
                      end do
                    end do
                  end do
                endif
            endif             ! FLtype(flnNum).eq.'tb0'
        endif     !  if(ndiv.gt.0)

!      write(855,*)'year=',year_now
!      write(855,*)qdivert
        
!d     if(iopt.eq.2)print*,'in read_divert_ef passed 187'

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 
      istep=al/1000

      if(nodivert.gt.0)then                             !999999999999999
!       FIND THE LOCAL COORDINATES FOR THE RESERVOIRS
!       THE VALUES FOR idiv AND jdiv ARE CHANGED

!       THE ORDER OF READING THE COORDINATES OF THE RESERVOIRS
!       MUST BE THE SAME AS READING THE CORRESPONDING FLOWS
!       IN S/R rdresv.

!       READ RELEASES
!       THE RESERVOIR OUTFLOWS ARE ALL READ IN THE FIRST TIME
!       rdresv IS CALLED. THEY ARE THEN STORED AND USED EACH TIME
!       rdresv IS CALLED.

!       IF NATURAL CONTROL, FLOWS ARE SET TO -1.0 ABOVE
!d       if(iopt.eq.2)print*,'in read_divert_ef passed 304'

      endif           !   if(nodivert.gt.0)

!d     if(iopt.eq.2)print*,'in read_divert_ef passed 551'

!     WE CAN'T HAVE -VE flows WHEN WE START
      do k=1,nodivert
        if(qdivert(k,1).lt.0.0)qdivert(k,ktr)=0.001
      end do

!     SET FUTURE RELEASES = LAST NON-NEGATIVE RELEASE
!     REALLY, WE'RE WORKING IN HOURLY INTERVALS ALTHOUGH THE      
!     RELEASES MAY BE READ IN ONLY WHEN THE RES OUTFLOW IS CHANGED.

!     rev. 9.5.14  Feb.  26/08  - NK: padded rel file for missing data
!     fill in missing data if rel file is shorter than the str file
!     nl is the length in nrs of the str file
!     ndiv is the length of the rel file
!     added Feb. 26/08  -nk-

!     rev. 10.4.34 Aug.  01/21  = NK Added check on qdivert dimension in read-divert
!      if(ndiv.gt.0)then
      if(ndiv.gt.0.and.nl.le.mhtot)then
!       fill data at end of file only if there are values
        if(ndiv.lt.nl)then
          do j=ndiv+1,nl
            do k=1,nodivert
              qdivert(k,j)=qdivert(k,j-1)
            end do
          end do
          
          print*
          print*,'WARNING:'
          print*,'diversion file is shorter than the str file'
          print*,'missing data has been set = to last recorder release'
          print*
          do j=2,nl
            do k=1,nodivert
              if(qdivert(k,j).lt.0.0)qdivert(k,j)=qdivert(k,j-1)
            end do
          end do
        endif
      endif                        

      close(unit=unitNum,status='keep')

      firstpass='n'
  999 RETURN                 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<-----

! FORMATS

  500 format(256f10.3)
 1011 format(' ',3x,'  i  idiv(i) jdiv(i)    b1(i)     b2(i)    b3(i)     b4(i)')
 1013 format(' ',3x,i3,4i8,a12/)
 6801 format('   read_divert_ef: reservoir no =',i3,' mhtot =',i5)
 6802 format('   ',i5,256f8.2)

      END SUBROUTINE read_divert

