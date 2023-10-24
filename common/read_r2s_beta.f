	subroutine read_r2s_beta(ir,filename1)

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
!     executable:  glb.exe

!     Purpose:
!     Read individual r2s files with CMC_glb precip & temperature 

!     readn by nk.  29/06/05	
!     rev. 9.2.02  Jun.  29/05  - NK: Added read_r2s 
!     rev. 9.8.30  Oct.  16/12  - NK: remove p(i,j)=0.0 from precip_adjust
!                               - also (likely) changed this s/r to parse the r2s file
!                               - needs a lot more work to be a full parser
 
!NOTE:
!this s/r can only read a static r2s file with one frame right now.
!Used in: bsn.f  to read the dem.r2s file
!and  precip_adjust.f to read the paf.r2s file
!so check if that works ok if changes are made.


!     input variable inarray(i,j) 

      use area_watflood
      implicit none
      save

	integer*4  :: i,j,ir,ios,iAllocate
	CHARACTER(50) :: junk
	character(30) :: filename1
	character(80) :: line,frameTime

      open(unit=ir,file=filename1,status='unknown',iostat=ios)
      if(ios.eq.0)then
        write(*,*)'opened unit =',ir,' filename =',filename1
      else
        print*,'Error opening ',filename1
        print*
        stop 'Program aborted in read_r2s @ 24'
      endif

      line='aaaaaaaaaaaaaaaaaaaaaaaa'
c	do while(line(1:10).ne.':endHeader'.or.line(1:10).ne.':EndHeader')
	do while(line(1:10).ne.':endHeader')
	  read(ir,50001)line
	  if(line(1:10).eq.':EndHeader')line=':endHeader'
50001   format(a80)
	  write(*,*)line(1:70)
	  if(line(1:11).eq.':Projection')read(line,*)junk,coordsys_temp
	  if(line(1:9) .eq.':CoordSys')read(line,*)junk,coordsys_temp
	  if(line(1:10).eq.':Ellipsoid')read(line,*)junk,datum_temp
	  if(line(1:6) .eq.':Datum')read(line,*)junk,datum_temp
	  if(line(1:5) .eq.':Zone')read(line,*)junk,zone_temp
	  if(line(1:8) .eq.':xOrigin')read(line,*)junk,xorigin_temp
	  if(line(1:8) .eq.':yOrigin')read(line,*)junk,yorigin_temp
	  if(line(1:7) .eq.':xCount')read(line,*)junk,xcount_temp
	  if(line(1:7) .eq.':yCount')read(line,*)junk,ycount_temp
	  if(line(1:7) .eq.':xDelta')read(line,*)junk,xdelta_temp
	  if(line(1:7) .eq.':yDelta')read(line,*)junk,ydelta_temp
	  if(line(1:10).eq.':FrameTime')read(line,50002)junk,frameTime
50002   format(a11,a10)	  
	end do
	
	
	
c	print*,frameTime(1:10)
c	read(frameTime,6000)year,junk,month,junk,day
c6000  format(i4,a1,i2,a1,i2)
c	print*,'reading ',year,month,day 

!     GreenKenue uses LatLong - code below uses LATLONG
      if(coordsys_temp.eq.'LatLong   ')coordsys_temp='LATLONG   '
      if(coordsys_temp.eq.'Catesian  ')coordsys_temp='CARTESIAN '

      if(coordsys_temp.eq.'UTM       ')then
	  continue
      elseif(coordsys_temp.eq.'CARTESIAN ')then
	  continue
      elseif(coordsys_temp.eq.'LATLONG   ')then
	  continue
      elseif(coordsys_temp.eq.'LatLong   ')then
	  continue
      else
        print*,' Valid format'
	  print*,'(UTM, CARTESIAN, LatLong or LATLONG) not found'
        print*,'format found is ***',coordsys_temp,'***'
        print*,'If wrong, please assign a Projection (coordsys)'
        print*,'in Green Kenue'
        print*
        stop 'Program aborted in read_r2s_beta@ 407'
      endif

      print*
      print*,'Values found:'
	print*,'projection=',coordsys_temp
	print*,'datum1=',datum_temp
	print*,'zone=',zone_temp
	print*,'xorigin=',xorigin_temp
	print*,'yorigin=',yorigin_temp
	print*,'xcount=',xcount_temp
	print*,'ycount=',ycount_temp
	print*,'xdelta=',xdelta_temp
	print*,'ydelta=',ydelta_temp
	
	
      if(.NOT.allocated(inarray))then      
        allocate(inarray(ycount_temp,xcount_temp),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *        'Error with allocation of inarray in read_r2s_beta @ 941'
        print*
      endif

c      read(ir,*)((inarray(i,j),j=1,xcount_temp),i=1,ycount_temp)
      do i=1,ycount_temp
        read(ir,*)(inarray(i,j),j=1,xcount_temp)
      end do

      print*,'last value read in filename1(1:40)'      
      print*


      close(unit=ir,status='keep')

	return


	end subroutine read_r2s_beta