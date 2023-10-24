      SUBROUTINE rdsdc()

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

!  READ SNOW COVER DEPLETION CURVE DATA
!  - CALLED ONLY ON FIRST TIME STEP

!     REV. 9.00    Mar.  2000   - TS: CONVERTED TO FORTRAN 90 

!***********************************************************************

      use area_watflood
	implicit none

	INTEGER :: iallcnt6,linenumber,ios,ii,i,nnts,iallocate,j
      logical :: exists

      DATA iallcnt6/0/

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE


!     OPEN THE SDC DATA FILE:
      INQUIRE(FILE=fln(13),EXIST=exists)
	IF(.NOT.exists)THEN
	  print*
	  print*,'The sdc file is not found'
	  print*,'Looking for filename:',fln(13)
	  print*
	  stop 'Program aborted in rdsdc @ 59'
	endif
      open(unit=43,file=fln(13),status='unknown',iostat=ios)
      if(ios.ne.0)then    ! added Nov. 10/14  nk
        print*
        print*,'Unable to open file',fln(13)(1:40)
        print*,'Possible cause(s):'
        print*,'file in use by another application'
        print*,'or target directory does not exist'
        stop 'Program aborted in rdsdc.f @ 40'
      endif
      
!     READ THE SDC'S:  FIRST THE CAPACITIES
      linenumber=1
      read(43,*,iostat=ios)
      if(ios.ne.0)then
        print*,' Error reading line #',linenumber
        print*
        stop ' Program aborted in rdsdc @ 40'
      endif

!     NEW FORMAT - SNOWCAPS NOT IN TOP LINE
      do ii=1,classcount
         linenumber=linenumber+1
         read(43,5001,iostat=ios)nsdc(ii),idump(ii),snocap(ii)
         if(ios.ne.0)then
           print*,' Error reading line #',linenumber
           print*
           stop ' Program aborted in rdsdc @ 45'
         endif
         do i=1,nsdc(ii)
           linenumber=linenumber+1
           read(43,*,iostat=ios)
           if(ios.ne.0)then
             print*,' Error reading line #',linenumber
             print* 
             stop ' Program aborted in rdsdc @ 53'
           endif
         end do
      end do

      if(iallcnt6.eq.0)then
!     TS - ADDED ALLOCATION OF AREAMELTA ARRAYS (REMAINDER)
      nnts=nsdc(classcount)
      allocate(sdcd(nnts,classcount),sdcsca(nnts,classcount),
     *                  stat=iAllocate)
      if(iAllocate.ne.0) STOP
     *   'Error with allocation of areamelta arrays in rdsdca' 
      iallcnt6=1
	endif

      rewind(43)
!     READ THE SDC'S:  FIRST THE CAPACITIES
      linenumber=1
      read(43,*,iostat=ios)
      if(ios.ne.0)then
        print*,' Error reading line #',linenumber
        print*
        stop ' Program aborted in rdsdc @ 81'
      endif
      do ii=1,classcount
         linenumber=linenumber+1
         read(43,5001,iostat=ios)
         if(ios.ne.0)then
           print*,' Error reading line #',linenumber
           print*
           stop ' Program aborted in rdsdc @ 89'
         endif
         if(snocap(ii).le.0.0)then
            write(51,6001)
         endif
         if(snocap(ii).gt.8850)then
            write(*,6002)ii
            snocap(ii)=8850.0
            write(*,'(A)',advance='no')
     *                 'snowcap set to <8850 - hit enter to continue'
            linenumber=linenumber+1
            read(*,*,iostat=ios)
            if(ios.ne.0)then
              print*,' Error reading line #',linenumber
              print*
              stop ' Program aborted in rdsdc @ 104'
            endif
         endif
         if(nsdc(ii).gt.20) write(6,5010)
         do j=1,nsdc(ii)
           linenumber=linenumber+1
           read(43,5002,iostat=ios)sdcsca(j,ii),sdcd(j,ii)
           if(ios.ne.0)then
             print*,' Error reading line #',linenumber
             print*
             stop ' Program aborted in rdsdc@45'
           endif
           if(sdcd(j,ii).lt.0)write(6,5013)
           if(sdcsca(j,ii).lt.0) write(6,5014)
         end do
      end do
      close (unit=43)

!     CHECK SDC'S - MONOTONICALLY INCREASING? IN ORDER?
      do ii=1,classcount
!         do j=2,j2
         do j=2,nsdc(ii)
            if(sdcd(j,ii).le.sdcd(j-1,ii)) write(6,5012)
            if(sdcsca(j,ii).le.sdcsca(j-1,ii)) write(6,5011)        
         end do
      end do

! FORMATS:

 5000 format(16f10.3)
 5001 format(2i5,f10.0)
 5002 format(2f10.3)
 5010 format(1x,'maximum of 20 points allowed on sdc!')
 5011 format(1x,'sca must be input in order \(i.e. 0 thru 1.0\)')
 5012 format(1x,'sdc depths must be monotonically increasing')
 5013 format(1x,'sdc depths must be positive')
 5014 format(1x,'sdc sca values must be positive') 
 6001 format(1x,'snow redistribution disabled in basin\bsn.sdc'/) 
 6002 format(' ','snocap set too high - weird things happen'/
     *           ' max value set to 8850 for class=',i5/
     *           ' hit enter to continue but fix basin/xxx.sdc file'/)


      RETURN

      END SUBROUTINE rdsdc

