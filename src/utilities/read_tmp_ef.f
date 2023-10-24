      SUBROUTINE read_tmp_ef(unitNum,flnNum,nhg,conv,ng)

!***********************************************************************
!    Copyright (C) 1987-2018 by Nicholas Kouwen  and Dave Watson (NRC)
        
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
!  READ_TMP_EF - written Jul/08 by Nick Kouwen
!	- Derived from read_flow_ef written by Dave Watson, CHC
!	This subroutine reads point temperature yyyymmdd_tag.tb0 file 
!	(tb0 format)
!*****************************************************************************

      use area_watflood

! TB0 data module
	USE EF_Module

	implicit none
	TYPE(RagParam) :: header
	TYPE(RagColumnMetaData) :: colHeader

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      CHARACTER(1)    :: hdrflg,firstpass,usefirstpass
      character(4096) :: inline
      character(30)   :: name_local

      INTEGER :: l,i,j,n,ios,no_old,iAllocate,ideallocate,
     *            linecount,dl,nhg,ng,old_nhg,old_ng,old_ng_1,
     *            line_no
	real*4  ::  conv
      LOGICAL exists,dataflg

      integer, dimension(:),   allocatable :: throwout
      integer,   dimension(:),  allocatable :: station_n,time_count
      character, dimension(:),  allocatable :: sta_data_flg

      data firstpass/'y'/
      data dataflg/.false./

! Parameter type definitions
	integer*4 unitNum, flnNum, iStat

! Local variables
	character*128 keyword, value

c	character*4096 line, subString, tmpString
      character*32768 line   ! changed Feb. 29/12 nk
      character*32768 subString, tmpString

	integer lineLen, keyLen, wordCount, colcount
	logical rStat, lineType, foundEndHeader, insideColMetaData

      if(firstpass.eq.'y')then
        write(52,211)(gname(l)(1:4),l=1,min(200,ng))
        write(53,211)(gname(l)(1:4),l=1,min(200,ng))
	endif

! Open the file
	INQUIRE(FILE=fln(flnNum),EXIST=exists)
	if(exists)then
		open(unit=unitNum,file=fln(flnNum),status='old',iostat=ios)
		if(ios.ne.0)then
			print*,'Problems opening ',fln(flnNum)
			print*
			STOP ' Stopped in read_tmp_ef @ 50'
		endif
	else
		print*,'ERROR: the data file ',fln(flnNum)
		print*,'is NOT found'
		STOP ' Program STOPPED in read_tmp_ef @ 69'
	endif

! Initialize default values
	CALL InitRagParam(header)	

! Search for and read tb0 file header
      linecount=0
      colcount=0
	line(1:1) = '#'
	foundEndHeader = .false.
	insideColMetaData = .false.

	do WHILE((.NOT.foundEndHeader) .AND.
     &	    ((line(1:1) .eq. '#') .OR.
     &		(line(1:1) .eq. ':') .OR.
     &		(LEN_TRIM(line) .eq. 0))) 	
                  linecount=linecount+1
!	if(iopt.eq.2)print*,'reading line ',linecount,' in read_tmp_ef'
		read(UNIT=unitNum, FMT='((A))', iostat=ios) line	! read a line
		if(ios .eq. -1)then
			write(6,'((A))') 'ERROR: Premature EndOfFile encountered'
			STOP ' Stopped in read_tmp_ef @ 90'
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
				iStat = ParseRagColumnMetaData(colHeader,
     &						keyword,keyLen,subString)
				if(iStat .lt. 0) then
					write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
					write(*,'(2(A))') '   in line: ',line					
					STOP ' Stopped in read_flow_ef'
					return
				endif
			else
				iStat = ParseRagParam(header,keyword,keyLen,subString)
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
	end do	
!***************************************
!	Finished reading the header
!***************************************

! Assign the variables from the types
      deltat=header%tb0p%DeltaT  !data time step in hours
	conv=header%tb0p%unitConv  !conversion to get mm of precip
	ng = colHeader%tb0cmd%colCount  !no of precip stations
      name_local=header%tb0p%sourceKey
      
      if(name_local(1:7).eq.'CMC_glb')then
!       the CMC_glb model provides 10 km gridded data
!       so with this radius of influence we ensure that
!       each grid in the watershed gets data only from 
!       the nearest grids
!       It is assumed that this is the last event
        radinfl=10.0
        smoothdist=1.0          ! no smoothing really
        print*
        print*,'NEW <<<<<'
        print*,'radius of influence set to ',radinfl
        print*,'for CMC_glb data'
        print*,'This value is kept for the rest of this run'
        print*,'as it is assumed this is a forecast'
        print*
      endif


d     print*,'ng of flow stations found =',ng

! Scan lines of data
      rewind unitNum

      nhg = CountDataLinesAfterHeader(unitNum)*deltat

!     fixed  nk  Nov. 22/06

	rewind unitNum
	CALL GoToStartOfData(unitNum)

      print*,'# stations=',ng
	print*,'# hours of data=',nhg
	print*,'DeltaT=',deltat
      print*,'UnitConversion=',conv
	print*






!       rev. 9.1.68  Dec.  19/04  - NK: rewrote read_tbo c/w memory allocation 
!       allocate stuff      
c      if(firstpass.eq.'y')then
c	        allocate(rrain(ng,nhg),throwout(ng),stat=iAllocate)
c     *		
c			if(iAllocate.ne.0) STOP
c     *		'Error with allocation of area5a arrays in sub'
c!         TS - ALLOCATION OF AREANASHA ARRAYS
c			allocate(xsta(ng),ysta(ng),gname(ng),sta_elv(ng),
c     *                   ewg(ng),sng(ng),stat=iAllocate)
c			if(iAllocate.ne.0) STOP
c     *		'Error with allocation of areanasha arrays in sub'
c	endif

!     revised for multiple ID's  nk Apr.07/10
      if(allocated(rrain))then
!       check to see if there is more data and 
!       reallocate if necessary
        if(nhg.gt.old_nhg.or.ng.gt.old_ng)then
          deallocate(rrain,throwout,station_n,sta_data_flg)  
          allocate(rrain(ng,nhg),throwout(ng),
     *    station_n(ng),sta_data_flg(ng),
     *        stat=iAllocate)
          if(iAllocate.ne.0)then
            print*,'No of stations found =',ng
            print*,'No of hrs data found =',nhg
            print*,'Error with allocation of arrays in read_tmp_ef'
            STOP 'Program stopped in read_tmp_ef @ 165'
          endif
          print*,'rrain re-allocated for ',ng,' location(s)'
          old_nhg=nhg		
          old_ng=ng
        endif
      else
        allocate(rrain(ng,nhg),throwout(ng),
     *    station_n(ng),sta_data_flg(ng),stat=iAllocate)
        if(iAllocate.ne.0)then
          print*,'No of stations found =',ng
          print*,'No of hrs data found =',nhg
          print*,'Error with allocation of arrays in read_tmp_ef'
          STOP 'Program stopped in read_tmp_ef @ 176'
        endif
        print*,'rrain allocated for ',ng,' location(s)'
        old_nhg=nhg
        old_ng=ng
      endif

      if(allocated(xsta))then
!       check to see if there are more stations and 
!       reallocate if necessary
        if(ng.gt.old_ng_1)then
          deallocate(xsta,ysta,gname,sta_elv,ewg,sng)
          allocate(xsta(ng),ysta(ng),gname(ng),sta_elv(ng),
     *                   ewg(ng),sng(ng),stat=iAllocate)
        
          if(iAllocate.ne.0)then
            print*,'Error with allocation of station descripters'
            print*,' in read_tmp_ef'
            STOP 'Program aborted in read_tmp_ef @ 195'
          endif
          old_ng_1=ng
        endif
      else
        allocate(xsta(ng),ysta(ng),gname(ng),sta_elv(ng),
     *                   ewg(ng),sng(ng),stat=iAllocate)
        if(iAllocate.ne.0)then
          print*,'Error with allocation of station descripters'
          print*,' in read_tmp_ef'
          STOP 'Program aborted in read_tmp_ef @ 204'
        endif
        old_ng_1=ng
      endif
   

!       rev. 9.1.68  Dec.  19/04  - NK: rewrote read_tbo c/w memory allocation 
  
!       ASSIGN THE GAUGE LOCATIONS NAMES AND FUNCTION COEFFICIENTS
	do l=1,ng
	  gname(l) = colHeader%tb0cmd%colName(l) ! temperature station name
	  xsta(l) = colHeader%tb0cmd%colLocX(l) ! x coordinate
	  ysta(l) = colHeader%tb0cmd%colLocY(l) ! y coordinate
c	      sta_elv(l)=colHeader%tb0cmd%colValue1(l)     ! station elevation
	  sta_elv(l)=colHeader%tb0cmd%elevation(l)     ! station elevation
	end do

	deallocate(colHeader%tb0cmd%colName)
	deallocate(colHeader%tb0cmd%colLocX)
	deallocate(colHeader%tb0cmd%colLocY)
	deallocate(colHeader%tb0cmd%elevation)

!     This does the lenght calculations in terms of no of grids.
!     That way it works for utm & latlong - although for lat long the 
!     n-s and e-w distances are not properly accounted for if the grids
!     are not set up as nearly square. But that's the user's problem.
      do n=1,ng
        ysta(n)=(ysta(n)-yorigin)/ydelta        
        xsta(n)=(xsta(n)-xorigin)/xdelta
      end do
            
!     Read the data:            
            
        If(deltat.eq.0)then
          print*
          print*,'ERROR:'
          print*,'deltaT in ',fln(flnNum)(1:40)
          print*,'is 0'
          print*,'please fix the file'            
          stop 'Program aborted in read_tmp_ef @ 267'
        endif
        line_no=0    
!       WATFLOOD COLUMN FORMAT
		do j=deltat,nhg,deltat   ! changed nk  sept. 29/06

!           Added this section to check # cols = # stations in the headet
!           Oct. 27/2012  NK
	 	    read(UNIT=unitNum, FMT='((A))', iostat=ios) line	! read a line
		    if(ios .eq. -1)then
			    write(6,'((A))') 'ERROR: Premature EndOfFile encountered'
			    STOP ' Stopped in read_tmp_ef @ 266'
		    end if
			line_no=line_no+1
			wordCount = countwords(line)
			if(ng.ne.wordcount)then
			  print*,'Error:'
			  print*,'No of columns of data ',wordcount
			  print*,'not equal to the no columns in the header',ng
			  print*,'in data line ',line_no
			  ng=wordcount
			  pause 'Program pause in read_tmp @ 276'
!			  stop 'Program aborted in read_tmp @ 276'
			endif
		    
c			read(unitNum,*,iostat=ios)(rrain(l,j),l=1,ng)
			read(line,*,iostat=ios)(rrain(l,j),l=1,ng)
			
			if(ios.ne.0)then
				print*, 'In read_tmp_ef'
				print*,' NEW format precip file'
				print*,starttime
				print*,' problems reading the temperature at line '
     *				    ,line_no
c     *				,j/deltat
				print*,' and column',l,'  ng=',ng
				print*,'Weird hidden characters in the file will
     *			 do this'
				print*,'last values read:'
				if(line_no.gt.2)then
					write(*,206)(rrain(l,j-2*deltat),l=1,ng)
				endif
				if(line_no.gt.1)then
					write(*,206)(rrain(l,j-deltat),l=1,ng)
			endif
				write(*,206)(rrain(l,j),l=1,ng)
				print*
		   		stop 'Program aborted in read_tmp_ef @ 368'
	        endif
	    end do

!     The unit conversion factor is processed in the calling program

c	do l=1,ng
c	  throwout(l)=1
c	end do
c	throwout(2)=-1

c	do l=1,ng
c	  if(throwout(l).lt.0)then
c	    do j=deltat,nhg,deltat
c		  rrain(l,j)=-999.0
c		end do
c	  endif
c	end do

!     write the temp to tmp_info.txt:
c      if(iopt.gt.0)then
c		write(55,5001)
c		do j=deltat,nhg,deltat
c			write(55,202)(rrain(l,j),l=1,ng)
c			write(*,202)(rrain(l,j),l=1,ng)
c		end do
c      endif

	if(firstpass.eq.'y')then
 	    do l=1,ng
	      if(rrain(l,deltat).gt.-90.0)dataflg=.true.  !data in first timestep
		end do
	    if(.not.dataflg)then
	        print*,'No temperature values found in the first'
	        print*,'time step. Not allowed.'
			print*,'A temperature is required for the first time step'
	        Print*,'Please correct the yyyymmdd_tag.tb0 file'
	        Stop 'Program aborted in read_tmp_ef @ 222'
	    endif
	endif




!     tabulate data availability
!     March 31, 2011  nk
      do l=1,ng
        sta_data_flg(l)='n'
        station_n(l)=0
        do j=deltat,nhg,deltat
          if(rrain(l,j).ge.-90.0)then
            station_n(l)=station_n(l)+1
          endif
        end do
        if(station_n(l)*deltat.eq.nhg)then
          sta_data_flg(l)='y'
        elseif(station_n(l).gt.0)then
          sta_data_flg(l)='p'
        else
          sta_data_flg(l)='n'
        endif
	  station_n(l)=station_n(l)*deltat/24
      end do

      write(52,207)id,(station_n(l),l=1,min(200,ng))
      write(53,208)id,(sta_data_flg(l),l=1,min(200,ng))



!     rev. 9.1.68  Dec.  19/04  - NK: rewrote read_tbo c/w memory allocation 
!     moved from sub

      close(unit=unitNum,status='keep',iostat=ios)
      if(ios.ne.0)then
		print*,'Problem closing unit',unitNum,' fln=',flnNum
		print*
		stop ' program aborted in read_tbo @ 623'
      endif

      firstpass='n'
      RETURN

! FORMATS

  201 format(9999f10.0)
  202 format(9999f10.3)
  203 format(' id=',i5,' l= ',i5,a12,2f12.3,4e10.3)
  205 format(2f12.3,a12,f12.3,4e10.3)
  206 format(9999f10.3)
  207 format(i5,9999i5)
  208 format(i5,9999a5)
  209 format(i5,99996f10.3)
  210 format(5x,9999a10)
  211 format(5x,9999a5)
 1778 format(' id,l,iy(l),jx(l)',6i5)
 1801 format(14x,10(f5.0,a1),i6)
 1802 format(14x,11(f5.0,a1))
 5000 format(' echo recorded data:')
 5001 format(/' echo converted data:')
 5004 format(a20,a10)
 5005 format(a20,i5)
 6226 format(' error encountered opening unit=37 fln=',a31/)

      END SUBROUTINE read_tmp_ef



