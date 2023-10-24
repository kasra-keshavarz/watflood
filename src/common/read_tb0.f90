!      SUBROUTINE read_tb0(unitNum,flnNum,hdrflg,date)
      SUBROUTINE read_tb0(unitNum,flnNum)

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
!  Converted to read_tb0 Sept. 05/21 by NK
!	This subroutine reads tb0 files 
!	(tb0 format)
!*****************************************************************************

      use area_watflood

! TB0 data module
	USE EF_Module

	implicit none
	TYPE(FlowParam) :: header
	TYPE(FlowColumnMetaData) :: colHeader

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      CHARACTER(1) :: hdrflg,firstpass,usefirstpass

      INTEGER :: l,i,j,ios,no_old,iAllocate,ideallocate,linecount,dl
      LOGICAL exists
      character(14) :: date

      data firstpass/'y'/

! Parameter type definitions
	integer*4 unitNum, flnNum, iStat

! Local variables
	character*128 keyword, value
	character*4096 line, subString, tmpString
	integer lineLen, keyLen, wordCount
	logical rStat, lineType, foundEndHeader, insideColMetaData

! Set unit and fln number
!	unitNum = 36
!	flnNum = 6

! Open the file
	INQUIRE(FILE=fln(flnNum),EXIST=exists)
	if(exists)then
		open(unit=unitNum,file=fln(flnNum),status='old',iostat=ios)
		if(ios.ne.0)then
			print*,'Problems opening ',fln(flnNum)
			print*
			STOP ' Stopped in read_flow_ef'
		endif
	else
		print*,'ERROR: the recorded stream hydrographs (STR) file ',&
              fln(flnNum)
		print*,'is NOT found'
		STOP ' Program STOPPED in read_flow_ef'
	endif

! Initialize default values
	CALL InitFlowParam(header)	

! Search for and read tb0 file header
      linecount=0
	line(1:1) = '#'
	foundEndHeader = .false.
	insideColMetaData = .false.

	do WHILE((.NOT.foundEndHeader) .AND.&
     	    ((line(1:1) .eq. '#') .OR.&
     		(line(1:1) .eq. ':') .OR.&
     		(LEN_TRIM(line) .eq. 0))) 	
                  linecount=linecount+1
!	if(iopt.eq.2)print*,'reading line ',linecount,' in read_flow_ef'
		read(UNIT=unitNum, FMT='((A))', iostat=ios) line	! read a line
		if(ios .eq. -1)then
			write(6,'((A))') 'ERROR: Premature EndOfFile encountered'
            print*,'in file ',fln(flnNum)(1:50)
            print*,'Delete the file & try again'
			STOP ' Stopped in read_flow_ef @ 99'
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
				iStat = ParseFlowColumnMetaData(colHeader,keyword,&
     												keyLen,subString)
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
!					write(*,'((A), (A))')  'Unrecognized keyword line: ',&
!     										line
				endif
			end if
		end if
	end do	
!***************************************
!	Finished reading the header
!***************************************


! Assign the variables from the types
!	kt =  header%tb0p%deltaT    !data time step in hours
	deltaT_temp =  header%tb0p%deltaT    !data time step in hours
	irdt =  header%routingDeltaT    ! routing time step in hours
	UnitConv_temp = header%tb0p%unitConv

!     FILL IN FLOW DATA flowfillflg='y'
!     this is when you have say 24 hr flow data but you want to plot hourly
      flowfillflg='n'
	if(GetBoolean(header%fillFlag) .eq. .true.) flowfillflg = 'y' 


!     no    =     no of streamflow stations = # columns
!	no = colHeader%tb0cmd%colCount
	xcount_temp = colHeader%tb0cmd%colCount

      print*,'no of flow stations found =',xcount_temp


!     nl    =     no of hours of streamflow data
! Scan lines of data
      rewind unitNum
!	nl = CountDataLinesAfterHeader(unitNum)

!      if(kt.eq.24)then
!        nl = CountDataLinesAfterHeader(unitNum)*kt
        ycount_temp = CountDataLinesAfterHeader(unitNum)*deltaT_temp     
!	else
!        nl = CountDataLinesAfterHeader(unitNum)*kt
!      endif

!     fixed  nk  Nov. 22/06

        print*
        print*,'In read_tb0 fln = ',fln(56)(1:40)
        print*,'no of data lines found =',ycount_temp


	rewind unitNum
	CALL GoToStartOfData(unitNum)

!       rev. 9.1.68  Dec.  19/04  - NK: rewrote rdflow c/w memory allocation 
!       allocate stuff      
			no_old=xcount_temp
			nl_max=ycount_temp
              
      if(.not.allocated(inarray))then        
	        allocate(inarray(xcount_temp,ycount_temp),stat=iAllocate)
			if(iAllocate.ne.0) STOP &
                  'Error with allocation of inarray'
			if(iAllocate.ne.0) STOP&
     		'Error with allocation of inarray in read_tb0 @ 61'
      endif    
          
!       rev. 9.1.68  Dec.  19/04  - NK: rewrote rdflow c/w memory allocation 
  
!       ASSIGN THE GAUGE LOCATIONS NAMES AND FUNCTION COEFFICIENTS
		do l=1,xcount_temp
			gage(l) = colHeader%tb0cmd%colName(l) ! streamflow station name
			xstr(l) = colHeader%tb0cmd%colLocX(l) ! x coordinate
			ystr(l) = colHeader%tb0cmd%colLocY(l) ! y coordinate
		end do
		deallocate(colHeader%tb0cmd%colName)
		deallocate(colHeader%tb0cmd%colLocX)
		deallocate(colHeader%tb0cmd%colLocY)

      print*,'xcount_temp or # cols =',xcount_temp
          
!       WATFLOOD COLUMN FORMAT
!		do j=kt,nl,kt   ! changed nk  sept. 29/06
		do j=deltaT_temp,ycount_temp,deltaT_temp   ! changed nk  sept. 29/06
			read(unitNum,*,iostat=ios)(inarray(j,l),l=1,xcount_temp)
!			write(*,*,iostat=ios)(inarray(j,l),l=1,xcount_temp)
			if(ios.ne.0)then
				print*, 'In strfw'
				print*,' NEW format flow file'
				print*,starttime
				print*,' problems reading the hydrograph at hour '&
     				,j/kt
				print*,' and column',l,'  no=',no
				print*,'Weird hidden characters in the file will do this'
				print*,'last values read:'
				if(j.gt.2)then
					write(*,206)(inarray(j-2,l),l=1,xcount_temp)
				endif
				if(j.gt.1)then
					write(*,206)(inarray(j-1,l),l=1,xcount_temp)
				endif
				write(*,206)(inarray(j,l),l=1,no)
				print*
		   		stop 'Program aborted in strfw @ 368'
	        endif

206           format(256f10.3)
9993			format(' error reading flows time step =',i3,/&
     	   ' possible cause: station list does not match # specified'/&
     		' or there may be tabs in the data')
		end do

!     rev. 9.1.68  Dec.  19/04  - NK: rewrote rdflow c/w memory allocation 
!     moved from sub
      close(unit=unitNum,status='keep',iostat=ios)
      if(ios.ne.0)then
		print*,'Problem closing unit 36 fln=',fln(6)
		print*
		stop ' program aborted in rdflow @ 623'
      endif


      firstpass='n'


      RETURN



      END SUBROUTINE read_tb0



