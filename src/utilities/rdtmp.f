      SUBROUTINE rdtmp(ng,nhg,conv,ixr,iyr)

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
!     include 'debug.for'

!********************************************************************
! RDTMP - THIS SUBROUTINE INPUTS temperature DATA
!
! Modified by Tricia Stadnyk - September 2000
! Converted common blocks to modules and added dynamically allocated
! run-time arrays as part of Fortran 90 conversion.
!
!
! - list of arguments:
!
!   o - rain( , )real*4    rain gauge precipitation data
!   i - ng       int       number of temp. stations
!   i - nhg      int       number of hours of data
!   i - conv     real*4    conversion factor
!   i - fln( )   char*12   file names
!   o - sta( )   cmplx    gauge locations (grid square numbers)
!   o - ewg( )   real*4   east-west utm coordinates of gauges
!   o - sng( )   real*4   north-south utm coordinates of gauges
!   o - ixr      int      number of east-west grid squares
!   o - iyr      int      number of north-south grid squares
!   o - radw     int      utm coordinate of west side of grid
!   o - rads     int      utm coordinate of south side of grid
!   o - gname( ) char*12  rain gauge names
!   i - fln( )   char*30  file names
!
!********************************************************************

c      USE area2
c     USE area12
c      USE area16
c     use areawfo

      use area_watflood

      CHARACTER(20) :: junk
      CHARACTER(25) :: junk1
	character(1)  :: junk2
	character(4096) :: line
	logical :: exists

      DATA firstpass/'y'/

c      if(firstpass.eq.'y')then
c        allocate(d(ng),stat=iAllocate)
c        if (iAllocate.ne.0) STOP 
c     *    'Warning: error with allocation of d() in weight' 
c        firstpass='n'
c     endif

! OPEN point temperature INPUT FILE  tag
      INQUIRE(FILE=fln(14),EXIST=exists)
      IF(exists)THEN
        open(unit=44,file=fln(14),status='old',iostat=ios)
        if(ios.ne.0)then
          print*,' Problems opening unit 44 fln=',fln(14)
          print*
          stop 'program aborted in rdtmp @ 48'
	  endif
	else
	  print*
	  print*,'Fatal error'
	  print*,'file ',fln(14),' not found'
	  print*,'file name should be something like: tempg\yyyymmdd.tag'
	  print*
	  print*,'This program requires a temp gauge data file (.tag)'
	  print*,'If you do not have gauge data, please create a dummy'
	  print*,'with all -99.9 for data for at least one station'
	  print*,'and make reference to it in the event.evt file'
	  print*
	  stop 'program aborted in rdtmp @ 57'
      endif


      write(51,*)'Opened  ',fln(14),'  in rdtmp.for'
	write(51,*)
      rewind 44
      read(44,2003)junk1
      write(51,*)'First line in ',fln(14)(1:50)
      write(51,*)junk1
      write(51,*)
      write(*,*)'First line in ',fln(14)(1:50)
      write(*,*)junk1
      write(*,*)


!     To run the weather model stuff:
!     Made up some new subroutines to read in IOP-2 data from other models
!     Did this by rewriting the rdrag s/r to read the files directly
      if(junk1.eq.'MOLOCH                   ')then
        print*,' Going to special case inprag_mo'
c        call inprag_mo(nhg,conv,ng)
        print*,' Back from inprag_mo  .... '
        return
      elseif(junk1.eq.'MESO-NH                  ')then
c        call inprag_nh(nhg,conv,ng)
        print*,' Back from inprag_nh  .... '
        return
      elseif(junk1.eq.'MC2                      ')then
!       cummulative
c        call inprag_mc(nhg,conv,ng)
        print*,' Back from inprag_mc  .... '
        return
      elseif(junk1.eq.'MM5RE                    ')then
!       cummulative
c        call inprag_mc(nhg,conv,ng)
        print*,' Back from inprag_mc  .... '
        return
      elseif(junk1.eq.'MM5E1                    ')then
!       cummulative
c        call inprag_mc(nhg,conv,ng)
        print*,' Back from inprag_mc  .... '
        return
      elseif(junk1.eq.'BOLAM3                   ')then
c        call inprag_nh(nhg,conv,ng)
        print*,' Back from inprag_mc  .... '
        return
      elseif(junk1.eq.'GAUGE                    ')then
c        call inprag_nh(nhg,conv,ng)
        print*,' Back from inprag_mc  .... '
        return
      elseif(junk1.eq.' -103.81  50.04 remo.422.')then
c        call inp_data_remo(nhg,44,ng)
        print*,' Back from inp_data_remo  .... '
	  conv=1.0
	  deltat=6

      print*,'nhg=',nhg

        return
      endif

      junk2=junk1  ! this just makes junk2 = first character of junk1

c      if(junk2.eq.'#                        ')then
      if(junk2(1:1).eq.'#')then
!     programmd for new file type 08/06/04 nk
        print*,'yyyymmdd.tag format file found'
!       new file type
        if(iopt.eq.2)print*,' In tmp @ checkpoint a1'
        
        read(44,44001)line
44001   format(a4096)
        print*,line(1:72)
        read(line,*)junk,fileformat        ! :FileTypece 19961001.evtragmet


        write(51,*)junk,fileformat        ! :FileType
!       want a different variable name for coordsys as in shed
        read(44,*)junk,coordsys1
        write(51,*)junk,coordsys1
        newformat=1

        read(44,*)junk,datum1
        write(51,*)junk,datum1
        read(44,*)junk,zone1
        write(51,*)junk,zone1
        read(44,*,iostat=ios)junk        !#
        write(51,*)junk
        read(44,*,iostat=ios)junk,startdate
        write(51,*)junk,startdate
        read(44,*,iostat=ios)junk,starttime
        write(51,*)junk,starttime
        read(44,*,iostat=ios)junk        !#
        read(44,*,iostat=ios)junk,ng     
        write(51,*)junk,ng
        read(44,*,iostat=ios)junk,nhg    ! no hours of data
        write(51,*)junk,nhg
        read(44,*,iostat=ios)junk,deltat     !data time step in hours
        write(51,*)junk,deltat
        read(44,*,iostat=ios)junk,conv   ! conversion to cms
        write(51,*)junk,conv
        read(44,*,iostat=ios)junk        !#
        write(51,*)junk

! TS - ALLOCATION OF AREA16A ARRAYS.
        write(51,*)'Allocation for',ng,' stations'

        print*,'iall=',iall
c        if(iall.eq.0)then
        if(.not.allocated(xsta))then
          allocate(smcrag(ng),xsta(ng),ysta(ng),
     *       sng(ng),sngmin(ng),ewg(ng),
     *       ewgmin(ng),gname(ng),stat=iAllocateStatus)
          iall=iall+1
          if(iAllocateStatus.ne.0) STOP
     *     '***Allocation of AREA16A arrays failed in INPRAGA @ 189***'

         else

!         allocations were previously made using # of gauges in
!         fln(3) (.pdl) but there are a different # of gauges in the
!         .rag file. So re-allocation has to be done:  nk 03/12/04         
! nk - DEALLOCATION OF AREA16A ARRAYS
          deallocate(smcrag,xsta,ysta,
     *       sng,sngmin,ewg,
     *       ewgmin,gname,
     *       stat=iDeallocateStatus)
          if (iDeallocateStatus.ne.0) STOP   
     *    '**Deallocation for AREA16A arrays in tmpdsta failed @ 202**' 

          allocate(smcrag(ng),xsta(ng),ysta(ng),
     *       sng(ng),sngmin(ng),ewg(ng),
     *       ewgmin(ng),gname(ng),stat=iAllocateStatus)
          iall=iall+1
          if(iAllocateStatus.ne.0) STOP
     *     '***Allocation of AREA16A arrays failed in INPRAGA***'
        endif

! TS - ALLOCATION STATEMENT TO DYNAMICALLY ALLOCATE ARRAY FOR AREA16A
! ALLOCATIONS FOR XSTA,YSTA,EWG,SNG,GNAME OCCUR IN INPGRDA

        if(.not.allocated(rrain))then
          allocate(rrain(ng,nhg),temp(ng),smc(ng),
     *         stat=iAllocateStatus)
          if (iAllocateStatus .ne. 0) STOP 
     *    '**Allocation failed for area16 in tmpdsta 220**'
        endif
     

c      if(.not.allocated(d))then
c        allocate(d(ng),stat=iAllocateStatus)
c        if (iAllocateStatus .ne. 0) STOP 
c     *    '**Allocation failed for d in tmpdsta l75**'
c     endif


      if(allocated(sta_elv))then
!       check to see if there are more stations and 
!       reallocate if necessary
        if(ng.gt.old_ng_1)then
          deallocate(sta_elv)
          allocate(sta_elv(ng),stat=iAllocate)
        
          if(iAllocate.ne.0)then
            print*,'Error with allocation of sta_elv'
            print*,' in rdtmp'
            STOP 'Program aborted in rdtmp_ef @ 195'
          endif
          old_ng_1=ng
        endif
      else
        allocate(sta_elv(ng),stat=iAllocate)
        if(iAllocate.ne.0)then
          print*,'Error with allocation of sta_elv'
          print*,' in rdtmp'
          STOP 'Program aborted in rdtmp_ef @ 204'
        endif
        old_ng_1=ng
      endif
   






!       READ THE GAUGE LOCATIONS NAMES AND FUNCTION COEFFICIENTS
        read(44,*,iostat=ios)(gname(n),n=1,ng)
        Write(51,*)(gname(n),n=1,ng)
        if(ios.ne.0)then
          print*,' error reading precip station name n=',n
          print*,  (sng(n),n=1,ng)
          STOP ' program aborted in ragmet.for @ 166'
        endif
        read(44,*,iostat=ios)(xsta(n),n=1,ng)
        Write(51,*)(xsta(n),n=1,ng)
        if(ios.ne.0)then
          print*,' error reading x coordinate n=',n
          print*,(xsta(n),n=1,ng)
          STOP ' program aborted in ragmet.for @ 178'
        endif
        read(44,*,iostat=ios)(ysta(n),n=1,ng)
        write(51,*)(ysta(n),n=1,ng)
        if(ios.ne.0)then
          print*,' error reading y coordinate n=',n
          print*,(ysta(n),n=1,ng)
          STOP ' program aborted in ragmet.for @ 172'
        endif
        read(44,*,iostat=ios)junk         !#
        write(51,5005)junk                   !#
        read(44,*,iostat=ios)junk         !: endHeader
        write(51,5005)junk                        !: endHeader
        if(iopt.eq.2)print*,' In ragmet @ checkpoint a2'
        write(51,*)
!       convert to grid coordinates:
        do n=1,ng
          ysta(n)=(ysta(n)-yorigin)/ydelta        
          xsta(n)=(xsta(n)-xorigin)/xdelta
        end do
      else


!       old format
!       old format
!       old format
        rewind 44

!       INPUT GRID DATA
        read(44,2001,iostat=ios)rgrd,rads,radn,radw,rade,grdn,grde
	  if(ios.ne.0)then
	    print*,' Problems reading line 1 in the tag file'
	    print*,' Possible cause: wrong formats'
          print*
          stop ' Program aborted in rdtmp @ line 55'
	  endif
	  
        read(44,2000,iostat=ios)ng,nhg,deltat,conv


c        if(iall.eq.0)then

        if(.not.allocated(xsta))then
          allocate(smcrag(ng),xsta(ng),ysta(ng),
     *       sng(ng),sngmin(ng),ewg(ng),
     *       ewgmin(ng),gname(ng),stat=iAllocateStatus)
          iall=iall+1
          if(iAllocateStatus.ne.0) STOP
     *     '***Allocation of AREA16A arrays failed in INPRAGA***'
         else

!         allocations were previously made using # of gauges in
!         fln(3) (.pdl) but there are a different # of gauges in the
!         .rag file. So re-allocation has to be done:  nk 03/12/04         
! nk - DEALLOCATION OF AREA16A ARRAYS
          deallocate(smcrag,xsta,ysta,
     *       sng,sngmin,ewg,
     *       ewgmin,gname,
     *       stat=iDeallocateStatus)
          if (iDeallocateStatus.ne.0) STOP   
     *    '**Deallocation for AREA16A arrays in tmpdsta failed**' 

          allocate(smcrag(ng),xsta(ng),ysta(ng),
     *       sng(ng),sngmin(ng),ewg(ng),
     *       ewgmin(ng),gname(ng),stat=iAllocateStatus)
          iall=iall+1
          if(iAllocateStatus.ne.0) STOP
     *     '***Allocation of AREA16A arrays failed in INPRAGA***'
        endif

! TS - ALLOCATION STATEMENT TO DYNAMICALLY ALLOCATE ARRAY FOR AREA16A
! ALLOCATIONS FOR XSTA,YSTA,EWG,SNG,GNAME OCCUR IN INPGRDA

        if(.not.allocated(rrain))then
          allocate(rrain(ng,nhg),temp(ng),smc(ng),    !w(iyr,ixr,ng),
     *    ds(ng),stat=iAllocateStatus)
          if (iAllocateStatus.ne.0) STOP 
     *    '**Allocation failed for area16 in tmpdsta l72**'
        endif

        if(ios.ne.0)then
          print*,' Problems reading line 2 in the tag file'
          print*,' Possible cause: wrong format'
          print*
          stop ' Program aborted in rdtmp @ line 63'
	  endif

        do n=1,ng
          read(44,2100,iostat=ios)sng(n),ewg(n),gname(n)
	    if(ios.ne.0)then
	      print*,' Problems reading the locations in the '
	      print*,' tag file at gauge number',n
	      print*
            stop 'Program aborted @ line 72'
	    endif

          xsta(n)=float((ewg(n)-radw+rgrd)/rgrd)
          ysta(n)=float((sng(n)-rads+rgrd)/rgrd)
  	    write(51,1001)n,ewg(n),sng(n),gname(n),xsta(n),ysta(n)
        end do
!       end of old format

      endif



! INPUT THE TEMPERATURE DATA
      do nh=deltat,nhg,deltat

!         read(44,5000,iostat=ios)(rrain(is,nh),is=1,ng)
         read(44,*,iostat=ios)(rrain(is,nh),is=1,ng)

         if(ios.ne.0)then
           print*,' Error in tag file at'
             print*,' hour= ',nh
           write(51,*)
           write(51,*)' Error in tag file at hour= ',nh
             write(51,*)
           print*
         endif
         write(51,5000)(rrain(is,nh),is=1,ng)
      end do

      junk='good'
!     Error checking
      do nh=deltat,nhg,deltat
        do is=1,ng
          if(rrain(is,nh).gt.50.0.or.rrain(is,nh).lt.-100.0)then
            print*,'At station',is,'      hour',nh
            print*,'temperature read =',rrain(is,nh)
            print*,'temp assumed wrong and set to -99.0'
            rrain(is,nh)=-99.99
            junk='bad'
          endif
        end do
      end do
c      if(junk.eq.'bad')then
c        print*
c        print*
c        pause 'Bad data found. Hit enter to continue. Ctrl c to abort'
c      endif

! CLOSE INPUT FILE
      close(unit=44)

! FORMATS:

 1001 format(' n,ewg,sng,gname,xsta,ysta/',3i5,2x,a12,2f5.0)
 2000 format(3i5,f5.2)
 2001 format(5i5,2f5.0)
 2002 format(' ',a25, 'ignored')
 2003 format(a25)
 2100 format(2i5,1x,a12)
!NOTE: THIS FORMAT IS TO COMPLY WITH THE VBASIC I/O MENUS
 5000 format(999f8.0)


 5005 format(a20,i5)


      RETURN

      END SUBROUTINE rdtmp
