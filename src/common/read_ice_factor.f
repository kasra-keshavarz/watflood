      SUBROUTINE read_ice_factor()

!***********************************************************************
!    Copyright (C) 1987-2018 by Nicholas Kouwen and Dave Watson  
        
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
!  Adapted from
!  READ_RESV_EF - written Sep/06 by Dave Watson, CHC
!     - Derived from rdresv written by Nick Kouwen
!     This subroutine reads the ice_factor.tb0 file 
!     (tb0 format)
!*****************************************************************************

!******************
      use area_watflood

!     rev. 9.7.27  May.  26/11  - NK: Add lake_ice_factor
!     REV. 10.1.42 Oct   20/16  - NK: Reinstated read_ice_factor.f as default if present


! TB0 data module
      USE EF_Module
      implicit none
      TYPE(ResvParam) :: header
      TYPE(ResvColumnMetaData) :: colHeader

      Integer  :: ios,j,k,i,n,l,jz,jj
      integer  :: noresv_firstpass,
     *            iAllocate,iDeallocate
      real*4   :: factor

!     rev. 9.1.55  Jun.  12/04  - NK: write new str files to strfw\newfmt folder.
      LOGICAL exists
      CHARACTER(1)  :: firstpass

      data firstpass/'y'/

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

      if(firstpass.ne.'y')return

c! Set unit and fln number
      unitNum=99
        
      if(.not.allocated(lake_ice_factor))then
        allocate(lake_ice_factor(noresv,12),stat=iAllocate)
        if(iAllocate.ne.0)then
          print*,'Error with allocation - read_ice_factor  @259'
          print*,'noresvo=',noresv
          STOP
        endif
      else
          print*,'Lake_ice_factors allocated',noresv,12
          pause
      endif
      
!     assign default values        
	do l=1,noresv
	  do j=1,12
	    lake_ice_factor(l,j)=1.00
	  end do
	end do
	  
!     If the icelakeflg in the event file = n or not present lake ice factors 
!     will not be used. If the icelakeflg = y and there is no file, the ice_fctr
!     computed in s/r ice_factor will be used

      if(icelakeflg.ne.'y')return
      
! Open the file
      INQUIRE(FILE='resrl\ice_factor.tb0',EXIST=exists)
      if(exists)then
        open(unit=unitNum,file='resrl\ice_factor.tb0',
     *     		status='old',iostat=ios)
        if(ios.ne.0)then
          print*,'Problems opening resrl\ice_factor.tb0'
          print*
          STOP ' Stopped in read_ice_factor'
        endif
      else
        if(numa.eq.0.and.dds_flag.eq.0)then
          print*,'WARNING: Ice factor file resrl\ice_factor.tb0'
          print*,'WARNING: is NOT found. Program continues '
          print*,'WARNING: with no lake ice correction'
          write(98,*)'WARNING: Ice factor file resrl\ice_factor.tb0'
          write(98,*)'WARNING: is NOT found. Program continues '
          write(98,*)'WARNING: with no lake ice correction'
          print*
        endif
	  return
      endif

d      if(iopt.eq.2)print*,'in _ice_factor passed 70'

      if(noresv.ge.1)then

! Initialize default values
        CALL InitResvParam(header)    

d        if(iopt.eq.2)print*,'in read_ice_factor passed 95'

!       NOTE:
!       column count is done on the line :ColumnLocationX <<<<<!!!!!

! Search for and read tb0 file header
        line(1:1) = '#'
        foundEndHeader = .false.
        insideColMetaData = .false.

        do WHILE((.NOT.foundEndHeader) .AND.
     &        ((line(1:1) .eq. '#') .OR.
     &        (line(1:1) .eq. ':') .OR.
     &        (LEN_TRIM(line) .eq. 0)))   

          read(UNIT=unitNum, FMT='((A))', iostat=ios) line    ! read a line
	    write(53,*)line(1:72)
          if(ios .eq. -1)then
              write(6,'((A))') 'ERROR: Premature EndOfFile encountered'
              STOP ' Stopped in read_ice_factor @ 103'
          end if

          rStat = Detab(line)         ! replace tabs with spaces
          line = ADJUSTL(line)        ! Get rid of leading white space
          lineLen = LEN_TRIM(line)    ! Find the length excluding trailing spaces

          if(line(1:1).eq.':')then
              wordCount = SplitLine(line,keyword,subString) ! find the keyword
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
                      write(*,'(2(A))') 'ERROR parsing ice_factor.tb0'
                      write(*,'(2(A))') '   in line: ',line                   
                      STOP ' Stopped in read_ice_factor @ 132'
                      return
                  endif
              else
                  iStat= ParseResvParam(header,keyword,keyLen,subString)
                  if(iStat .lt. 0) then
                      write(*,'(2(A))') 'ERROR parsing ice_factor.tb0'
                      write(*,'(2(A))') '   in line: ',line                   
                      STOP ' Stopped in read_ice_factor @ 140'
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

d       if(iopt.eq.2)print*,'in read_ice_factor passed 157'
c       print*,noresv,colHeader%tb0cmd%colCount
! Assign the variables from the types
c       ktr =  header%tb0p%deltaT    !data time step in hours
	  if(noresv.ne.colHeader%tb0cmd%colCount)then
		  print*,'No of columns in the resrl\ice_factor.tb0 file'
	    print*,'does not match the number in the' 
	    print*,'yyyymmdd_rel.tb0 file'
	    print*,'colHeaderCount=',colHeader%tb0cmd%colCount
	    print*,'No reservoirs in rel file =',noresv
	    print*
	    stop 'Program aborted in read_ice_factor @ 166'
	  endif
!     nrel    =     no of hours of data
! Scan lines of data
        rewind unitNum
        k = CountDataLinesAfterHeader(unitNum)
	  if(k.ne.12)then
	    print*,'No of entries in the resrl\ice_factor.tb0 file'
	    print*,'must be 12 - one for each month' 
	    print*,'No found =',k
	    print*
	    stop 'Program aborted in read_ice_factor @ 178'
	  endif
        rewind unitNum
        CALL GoToStartOfData(unitNum)
      endif
d      if(iopt.eq.2)print*,'in read_resv_ef passed 274'
 
       if(k.gt.0)then
!        case with reservoir releases
!        do j=ktr,nrel,ktr
!     rev. 9.1.13  Mar.  23/02  - fixed resv. timing, moved to beginning of dt
         write(53,*)'month - ice factors each lake:'
         do j=1,12
           read(unitNum,*,iostat=ios)(lake_ice_factor(l,j),l=1,noresv)
           write(53,53001)j,(lake_ice_factor(l,j),l=1,noresv)
53001      format(i5,<noresv>f5.2)              
d          if(iopt.eq.2)print*,j,(lake_ice_factor(l,j),l=1,noresv)
           if(ios.ne.0)then
             write(98,*)' Error on resrl\ice_factor.tb0'
             write(98,*)' Trying to read ice_factor =',j
             print*,' Error on resrl\ice_factor.tb0'
             print*,' Trying to read read ice_factor =',j
             print*,' ios= ',ios
             if(ios.eq.-1)then
               write(98,*)'End of file in read ice_factor '
               write(98,*)'Possibly last line does not have a return'
               print*,'End of file in read ice_factor '
               print*,'Possibly last line does not have a return'
               print*
             else
               print*
               STOP ' program aborted in read_ice_factor.f @345'
             endif
           endif
         end do
       endif     !  if(nrel.gt.0)
      close(unit=unitNum,status='keep')
      
      print*
      print*,'NEW'

      icefactorfile=.true.

  999 RETURN                 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<-----

! FORMATS

  500 format(256f10.3)
 1011 format(' ',3x,'  i  ires(i) jres(i)    b1(i)     b2(i)',
     * '      b3(i)       b4(i)      b5(i)      resvname')
 1013 format(' ',3x,i3,2i8,5g12.3,a12/)
 6801 format('   read_resv_ef: reservoir no =',i3,' mhtot =',i5)
 6802 format('   ',i5,256f8.2)

      END SUBROUTINE read_ice_factor

