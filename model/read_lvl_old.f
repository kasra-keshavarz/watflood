      SUBROUTINE read_lvl(unitNum,flnNum)

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

!    You should have received a copy of the GNU Lesser General Public License
!    along with WATFLOOD.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
     
!***************************************************************************
! PROGRAM BY: NK  June 2012
!
      use area_watflood
      USE EF_Module

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      implicit none
      save

      TYPE(LvlParam) :: header
      TYPE(LvlColumnMetaData) :: colHeader

      CHARACTER(1) :: new
      INTEGER      :: jan,ju,iallcnt2,iallocate,n,ii,i,j,k,l,nold
      integer      :: n_max,ideallocate,no_hdr_lines
      integer      :: ios,linecount,unitnum,flnNum,iStat
      integer      :: nlines,nlines_old
      REAL*4       :: time,ttime,taold,temmmp  
      character*16 :: swe_file              
      logical      :: exists,firstpass,errflg1

! Local variables
      character*128 keyword, value
      character*4096 line, subString, tmpString
      CHARACTER(1)  :: col0(13),col1(13),col2(26)
      integer lineLen, keyLen, wordCount
      logical rStat, lineType, foundEndHeader, insideColMetaData

      save
      
      DATA new/'f'/firstpass/.true./errflg1/.false./
      DATA iallcnt2/0/

      DATA col1/'b','d','f','h','j','l','n','p','r','t','v','x','z'/
      DATA col0/'c','e','g','i','k','m','o','q','s','u','w','y','a'/
      DATA col2/' ','a','b','c','d','e','f','g','h','i','j','k','l'
     *         ,'m','n','o','p','q','r','s','t','u','v','w','x','y'/

d      if(iopt.eq.2)print*,'in melt @ 0'

      print*,'~~~~~~~~~~~~~~~~~~~~~~~~'
      print*,'NEW  <<<<<'

! Open the file
      INQUIRE(FILE=fln(flnNum),EXIST=exists)
      if(exists)then
        open(unit=unitNum,file=fln(flnNum),status='old',iostat=ios)
        if(ios.ne.0)then
          print*,'Problems opening:',fln(flnNum)(1:40)
          print*
          STOP ' Stopped in read_lvl @ 39'
        endif
        print*,'Opened ',fln(flnNum)(1:40)
      else
        print*,fln(flnNum)(1:40),'not found'
        print*,'program continues without lake level analysis'
        lvlflg=.false.
        return
      endif

c      pause 2

! Initialize default values
      CALL InitlvlParam(header)	
      
c	pause 3
      
! Search for and read tb0 file header
      linecount=0
      line(1:1) = '#'
      foundEndHeader = .false.
      insideColMetaData = .false.
      no_hdr_lines=0

      write(51,*)
      write(51,*)'Level - observed level/elv. file'
      write(51,*)

      do WHILE((.NOT.foundEndHeader) .AND.
     &        ((line(1:1) .eq. '#') .OR.
     &        (line(1:1) .eq. ':') .OR.
     &        (LEN_TRIM(line) .eq. 0)))   
                  linecount=linecount+1
!     if(iopt.eq.2)print*,'reading line ',linecount,' in read_lvl'
          read(UNIT=unitNum, FMT='((A))', iostat=ios) line    ! read a line
          write(51,*)line(1:60)
          if(ios .eq. -1)then
              write(6,'((A))') 'ERROR: Premature EndOfFile encountered'
              STOP ' Stopped in read_lvl @ 88'
          end if
          no_hdr_lines=no_hdr_lines+1
          
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
                  iStat = ParseLvlColumnMetaData(colHeader,
     &                        keyword,keyLen,subString)
                  if(iStat .lt. 0) then
                      write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
                      write(*,'(2(A))') '   in line: ',line                   
                      STOP ' Stopped in read_lvl @113'
                      return
                  endif
              else
                  iStat = ParseLvlParam
     *                (header,keyword,keyLen,subString)
                  if(iStat .lt. 0) then
                      write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
                      write(*,'(2(A))') '   in line: ',line                   
                      STOP ' Stopped in read_lvl @122'
                      return
                  else if(iStat .eq. 0) then
!                     write(*,'((A), (A))')  'Unrecognized keyword line: ',
!     &                                       line
                  endif
              end if
          endif
          
      end do  
!***************************************
!     Finished reading the header
!***************************************

!     nolvl    = no of level data data
      ktlvl =  header%tb0p%deltaT    !data time step in hours
      nolvl    = colHeader%tb0cmd%colCount
      print*,'no of level columns found =',nolvl
      print*,'no header lines =',no_hdr_lines

d     print*
d     print*,'In read_lvl'
d     print*,colHeader%tb0cmd%colCount
d     print*,'No of level columns found =',nolvl
d     print*

c      pause 111

      
! Scan lines of data
      rewind unitNum
      nlines = CountDataLinesAfterHeader(unitNum)*ktlvl
      print*,'no of data lines found =',nlines
      
      rewind unitnum
      do n=1,no_hdr_lines
        read(unitnum,*)line
c        print*,line(1:60)
      end do
      

!     rev. 9.1.68  Dec.  19/04  - NK: rewrote rdflow c/w memory allocation 
!     rev. 9.8.41  Jan.  28/13  - NK: fixed bug in lst for level print statement
!     allocate stuff      
      if(firstpass)then
        n_max=n
        nold=nolvl
        nlines_old=nlines
        allocate(lvl_obs(nolvl,nlines),stat=iAllocate)
        allocate(lvl_calc(nolvl,nlines),stat=iAllocate)
        write(51,*)'lvl_obs dimensioned for:',nlines
        if(iAllocate.ne.0) STOP
     *		'Error with allocation of arrays in read_lvl_ef @ 192'
      else                    !  firstpass
!       check to see memory allocated is adequate      
        if(nolvl.ne.nold)then
          print*,'No of level stations has been changed in'
!     rev. 10.1.92 Aug   29/17  - NK: Fixed col check bug in read_lvl
          print*,'nold=',nold,' nolvl=',nolvl
          print*,'in file ',fln(53)(1:60)
          print*,'This is not permitted'
          print*
          stop 'Program aborted in read_lvl_ef @ 202'
        endif
        if(nlines.gt.nlines_old)then
          print*,'No of data lines have changed in'
          print*,'nlines_old=',nlines_old,' nlines=',nlines
          print*,'in file ',fln(53)(1:60)
          print*
c	  	if(n.gt.n_max)then
          n_max=n

!         the file length is longer than any of the previous events so 
!         more memory has to be allocated

!         DEALLOCATION OF ARRAYS FROM AREA10A:
          deallocate(lvl_obs,stat=iDeallocate)
          if(iDeallocate.ne.0)then
            print*,'Error with deallocation of lvl_obs'
            print*
            stop 'Program aborted in read_lvl_ef @ 220'
          endif
          allocate(lvl_obs(nolvl,nlines),
     *         				stat=iAllocate)
          if(iAllocate.ne.0) STOP
     *	      'Allocation Error: arrays in read_lvl_ef @228'
                           nlines_old=nlines
          write(51,*)'lvl_obs re-dimensioned for:',nlines
        endif
      endif                   !  firstpass

      
      if(.not.allocated(xlvl))then
!     rev. 10.1.69 Mar.  03/17  - NK: Changed allocation for lvl_reach in read_lvl.f 
          allocate(xlvl(nolvl),ylvl(nolvl),gname_lvl(nolvl),
     *                   ewg_lvl(nolvl),sng_lvl(nolvl),
     *                   lvl_reach(nolvl),stat=iAllocate)
        if(iAllocate.ne.0)then
          print*,'Error with allocation of station descripters'
          print*,' in read_lvl_ef'
          STOP 'Program aborted in read_lvl_ef @ 241'
        endif
      endif
  
!     ASSIGN THE GAUGE LOCATIONS NAMES AND FUNCTION COEFFICIENTS
      do n=1,nolvl
        gname_lvl(n) = colHeader%tb0cmd%colName(n) ! station name
        xlvl(n) = colHeader%tb0cmd%colLocX(n) ! x coordinate
        ylvl(n) = colHeader%tb0cmd%colLocY(n) ! y coordinate
      end do
      deallocate(colHeader%tb0cmd%colName)
      deallocate(colHeader%tb0cmd%colLocX)
      deallocate(colHeader%tb0cmd%colLocY)

c      pause '261'

      if(firstpass)then
        write(51,*)(gname_lvl(n),gname_lvl(n),n=1,nolvl)
91101   format(<2*nolvl>(1x,a12))      
      endif

c      pause '269'

!     write lake level location file
!     rev. 9.8.38  Nov.  13/12  - NK: changed name level_plotting.xyz > level_station_location.xyz
!     rev. 10.1.85 May   17/17  - NK: Level_station_location.xyx for iopt > 0 only
      if(iopt99)then
        open(unit=99,file='level_station_location.xyz',
     *            status='unknown',iostat=ios)
        do l=1,nolvl
          i=l
          j=1+i/13
          k=j
          if(l.gt.13)i=l-(i/13)*13
          if(i.eq.0)i=13
          if(mod(i,13).eq.0)j=j-1
          write(99,*)xlvl(l),ylvl(l),l,gname_lvl(l),
     *                  col2(j),col1(i),col2(k),col0(i)
c 6012   format(2f15.3,i5,1x,a12,3x,2a1,3x,2a1)
        end do
        close(unit=99,status='keep')
      endif

!     turn into local coordinates

      do n=1,nolvl
        ylvl(n)=(ylvl(n)-yorigin)/ydelta+1.0  
        xlvl(n)=(xlvl(n)-xorigin)/xdelta+1.0
c        print*,n,'xlvl',xlvl(n),'ylvl',ylvl(n)
        ewg_lvl(n)=int(xlvl(n))
        sng_lvl(n)=int(ylvl(n))
      end do

c      pause 297

      write(51,*)'Echo levels',ktlvl
      do j=ktlvl,nlines,ktlvl   
        read(unitNum,*,iostat=ios)(lvl_obs(n,j),n=1,nolvl)
        if(ios.ne.0)then
          print*, 'In read_lvl_ef @ 298'
c          write(*,*)j,(lvl_obs(n,j),n=1,4)
          print*,'Got as far as data line',j
        endif
      end do
!     rev. 10.2.31 Aug.  23/18  - NK: Echo recorded levels for iopt>=99 only 
      if(iopt99)then
          do j=ktlvl,nlines,ktlvl   
              write(51,*,iostat=ios)(lvl_obs(n,j),n=1,nolvl)
          end do
      endif
    
      close(unit=99,status='keep')

      close(unit=unitNum,status='keep',iostat=ios)
      if(ios.ne.0)then
		print*,'Problem closing unit 36 fln=',fln(6)
		print*
		stop ' program aborted in rdflow @ 623'
      endif

      if(firstpass)then
!       find out what lake (reach) each point is in
        write(51,*)
        write(51,*)'   level_station  grid#    reach#'
        write(*,*)'   level_station  grid#       reach#'
        do l=1,nolvl
          if(sng_lvl(l).le.ycount.and.ewg_lvl(l).le.xcount
     *             .and.sng_lvl(l).gt.1.and.ewg_lvl(l).gt.1)then   ! added Jan. 23/14 NK
            n=s(sng_lvl(l),ewg_lvl(l))
            
c            print*,l,n
c            print*,ireach(n)
c            print*,lvl_reach(l)
            lvl_reach(l)=ireach(n)
            if(lvl_reach(l).eq.0)then
              errflg1=.true.
c               write(*,*)l,n,lvl_reach(l)
            endif
            write(51,*)l,n,lvl_reach(l)
          endif
        end do
        
      if(errflg1)then
        print*
        print*,'WARNING'
        print*,'1 or more Level locations are not in lakes'
        print*,'designated as a reach - i.e. in a lake'
        print*,'Please correct input data:'
        print*,'Either add grid to a lake or move the level station'
        print*,'into a grid designated as a reach'
        print*,'Put the level_station_location.xyz file on the '
        print*,'map file with the reach entries showing' 
        print*,'If you don`t elect to make this change, the'
        print*,'date will be set to -999 for a dummy reach'
!       don't delete this stop!!  will crash in lst        
!        stop 'Program terminated in read_lvl_ef @ 312'
      endif
      endif

      firstpass=.false.  ! reset for next event to read header only
      
      return

      END SUBROUTINE read_lvl
