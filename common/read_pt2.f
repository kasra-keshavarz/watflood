      SUBROUTINE read_pt2(unitNum,flnNum,nrows,ncols)

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
!     s/r created Jan. 22/08 NK
!     copyright (c) nick kouwen 2003-

!     rev. 9.7.09  Oct.  02/10  - NK: ensure fpet_lake is not assigned unintended values
!     rev. 9.9.45  Dec.  03/14  - NK: Revamped read_pt2 for general use

      use area_watflood

!///////////////////////// 
!// Added by Dave
	USE EF_module
!// End Dave addition
!/////////////////////////
      implicit none
	save

      real, dimension(:),     allocatable :: ewg_fl,sng_fl

      CHARACTER(12) :: date
	INTEGER       :: iAllocateStatus,iDeallocateStatus,
     *                 iAllocate,ios,ng
      INTEGER(2)       status1
      CHARACTER(1)     firstpass  
      character(20) :: junk,junk2
	character(80) :: junk60
      character(60) :: junk1  
      logical       :: exists
	integer*4     :: unitNum,flnNum,nrows,ncols,i,j

	DATA firstpass/'y'/

      if(id.ne.1)return

!     test if file already exists. 
      INQUIRE(FILE=fln(flnNum),EXIST=exists)
      IF(exists)THEN
        open(unit=unitNum,file=fln(flnNum),
     *                 status='unknown',iostat=ios)
        if(ios.ne.0)then
          print*,'fln=',fln(flnNum)(1:40)
          print*,'ios=',ios
          STOP 'Program read_pt2 aborted @ 46'
        endif 
	  if(iopt.ge.1)print*,'Opened file ',fln(flnNum)(1:40)
      else
        print*,'Attempting to open the file ',fln(flnNum)(1:40)
        print*,'but it is not found (in this location)'
        print*
        stop 'Program aborted in read_pt2 @ 51'
      endif

!     programmd for pt2 file type 14/07/06 nk

      call find_filetype(flnNum)

!     Read the yyyymmdd_crs.pt2 file

d	write(51,*)
d	write(51,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
d     write(51,*)'PT2 debug information:'
d	write(51,*)'Reading fln=',fln(flnNum)(1:40)
d	write(51,*)'Reading the header - data found:'
d	write(51,*)

      if(filetype.ne.'pt2')then
	  print*,'pt2 file type expected but '
	  print*,filetype,' found for lake initialization'
	  print*
	  stop 'Program aborted in read_pt2 @ 71'
	else
        do while(junk1.ne.':Name     ')
          read(unitNum,*)junk1
d          write(51,*)junk1
        end do
        read(unitNum,*)junk1
d        write(51,*)junk1
        read(unitNum,*,iostat=ios)junk,coordsys2
	  write(51,*)'Projection=',coordsys2
        if(coordsys2.eq.'UTM       ')then
          read(unitNum,*,iostat=ios)junk,zone2
d	    write(51,*)zone2
          read(unitNum,*,iostat=ios)junk,datum2
d	    write(51,*)'Ellipsoid=',datum2
        endif
        if(coordsys2.eq.'LATLONG   '.or.coordsys2.eq.'latlong   '
     *           .or.coordsys2.eq.'LatLong   '   )then
          read(unitNum,*,iostat=ios)junk,datum2
d	    write(51,*)'Ellipsoid=',datum2
        endif

d       write(63,*)'In read_pt2 reading the attribute names:'
        ncols=0
        junk1='                                     '
        do while(junk1(1:14).ne.':endHeader')
          read(unitNum,60000)junk1                !,junk2,junk
60000     format(a60)          
          if(junk1(1:14).eq.':AttributeType')ncols=ncols+1
d          write(63,*)ncols,junk1(1:60)            !,junk2(1:15),junk(1:15)
        end do
d        write(63,*)'# attributes =',ncols
        
!       we're now at the first line after the header - i.e. first data line        
        if(allocated(attributeName))deallocate(attributeName)
		allocate(attributeName(ncols),stat=iAllocate)
		if(iAllocate.ne.0)then
		  write(98,*)'Error: allocating attributeName in read_pt2 @ 104'
		  print*,'Problem allocating attributeName in read_pt2 @ 104'
		  STOP 
        endif
        if(allocated(attributeType))deallocate(attributeType)
		allocate(attributeType(ncols),stat=iAllocate)
		if(iAllocate.ne.0)then
		  write(98,*)'Error: allocating attributeType in read_pt2 @ 110'
		  print*,'Problem allocating attributeType in read_pt2 @ 110'
		  STOP 
        endif

!       how many locations are there? Need to read 
!       to the end of the file and rewind.  
        nrows=0
        do while(.not.eof(unitNum))
          nrows=nrows+1
          read(unitNum,*,iostat=ios)junk1
d          write(63,*)nrows,junk1
        end do
c        nrows=nrows-1
d       write(63,*)'No locations found =',nrows  
 
        if(allocated(sng_fl))then    
          deallocate(sng_fl,ewg_fl,gname,
     *              stat=iAllocateStatus)
          if (iAllocateStatus .ne. 0) STOP
     *      '**DeAllocation failed in read_pt2 @ 127**'
        endif	 
        allocate(sng_fl(nrows),ewg_fl(nrows),gname(nrows),
     *              stat=iAllocateStatus)
        if (iAllocateStatus .ne. 0) STOP
     *      '**Allocation failed in read_pt2 @ 131**'

        if(allocated(inarray))then    
          deallocate(inarray,stat=iAllocate)
		  if(iAllocate.ne.0) STOP 
     *	   'Error with deallocation in read_pt2 @ 141'
        endif
		allocate(inarray(nrows,nCols),stat=iAllocate)
		if(iAllocate.ne.0) STOP 
     *	   'Error with allocation in read_pt2 @ 144'
!               the +1 is needed for fpet_lake
!             Initialize inarray so garbage won't corrupt AET s/r
!             for lake avaporation if data is not in the init lake level file
!             nk. Sept. 3/09
        
        do i=1,nrows
	    do j=1,ncols
	      inarray(1,j)=-1.0
	    end do
	  end do
	  
        rewind(unitNum)
        i=0
        do while(junk1.ne.':endHeader'.and.junk1.ne.':EndHeader')
          read(unitNum,*)junk1  !reads "attribute" and counts how many
          
          if(junk1(1:14).eq.':AttributeName')then
            backspace unitNum
            read(unitNum,*)junk1,junk1,junk1  !reads "attribute" and counts how many
            i=i+1
            AttributeName(i)=junk1(1:25)  ! may have to get the string length here
d          write(51,*)i,AttributeName(i)
          endif
        end do

d        write(51,*)
d        write(51,*)'Echo  ',fln(flnNum)(1:60)
d        write(51,*)
        do i=1,nrows
          read(unitNum,*,iostat=ios)
     *       ewg_fl(i),sng_fl(i),gname(i),(inarray(i,j),j=1,nCols-1)
d	    write(51,51000)
d    *       ewg_fl(i),sng_fl(i),gname(i),(inarray(i,j),j=1,nCols-1)
51000     format(2f15.3,5x,a12,<ncols-1>f10.3)     
        end do

d        write(51,*)'ios= ',ios,' No of locations found= ',nrows
        close(unit=unitNum,status='keep')
d        write(51,*)
d	  write(51,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

	endif

      return

      END SUBROUTINE read_pt2
