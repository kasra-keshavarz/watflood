	subroutine write_r2s(ir,filename1)

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
!     writen by nk.  26/06/05	
!     rev. 9.2.01  Jun.  29/05  - NK: Added write_r2s 

!     output variable outwfo(j,i) is borrowed from the watflood.wfo code
!	in the calling code, simply define outwfo() with the values you
!     want to write to an .r2s file

      use area_watflood
      implicit none
      save

	integer*4  :: i,j,ir,ios
	CHARACTER(10) :: time
      CHARACTER(8)  :: cday
	character(30) :: filename1

      real*4, dimension(:,:),allocatable::   r2s_variable

      open(unit=ir,file=filename1,status='unknown')
      if(ios.ne.0)then    ! added Nov. 10/14  nk
        print*
        print*,'Unable to open file',filename1
        print*,'Possible cause(s):'
        print*,'file in use by another application'
        print*,'or target directory does not exist'
        stop 'Program aborted in write_r2s.f @ 27'
      endif
!	print*,'opened unit =',ir,' filename = ',filename1
      write(ir,6701)
      write(ir,6702)
      write(ir,6703)application
      write(ir,6704)author
      
	call date_and_time(cday,time)
      
	write(ir,6705)
     *     cday(1:4),cday(5:6),cday(7:8),time(1:2),time(3:4),time(5:6)
      write(ir,6702)
      write(ir,6706)coordsys1
      write(ir,6710)xorigin+xdelta/2.0
      write(ir,6711)yorigin+ydelta/2.0
      write(ir,6702)
      write(ir,6707)attribute_count
      write(ir,6708)attribute_units
      write(ir,6709)attribute_name
	write(ir,6702)
      write(ir,6712)xcount_temp
      write(ir,6713)ycount_temp
      write(ir,6714)xdelta
      write(ir,6715)ydelta
      write(ir,6716)
      write(ir,6702)
      write(ir,6717)


      do i=1,ycount_temp
	  write(ir,6000)(outarray(i,j),j=1,xcount_temp)
      end do
      close(unit=ir,status='keep')
	print*,'written to ',filename1
	print*,'closed unit =',ir,' filename = ',filename1

	return

!	don't get rid of the formats as * formats write a blank character

6000  format(<xcount>g12.3)
6701  format(':Filetype r2s ASCII for ENSIM')
6702  format('#')
6703  format(':Application        ',a60)
6704  format(':WrittenBy          ',a60)
6705  format(':Creation date      '
     *        ,a4,'/',a2,'/',a2,'  ',a2,':',a2,':',a2,2x)
6706  format(':CoordSys           ',a10)
6707  format(':AttributeCount     ',i10)
6708  format(':AttributeUnits     ',a60)
6709  format(':AttributeName      ',a60)
6710  format(':xOrigin            ',f16.7)
6711  format(':yOrigin            ',f16.7)
6712  format(':xCount             ',i10)
6713  format(':ycount             ',i10)
6714  format(':xdelta             ',f16.7)
6715  format(':ydelta             ',f16.7)
6716  format(':Angle                0.000000')
6717  format(':EndHeader')
!6718  format(':CoordSys LATLONG')

	end subroutine write_r2s