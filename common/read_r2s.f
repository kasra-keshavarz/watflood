	subroutine read_r2s(ir,filename1)

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
!     readn by nk.  29/06/05	
!     rev. 9.2.02  Jun.  29/05  - NK: Added read_r2s 

!     output variable outwfo(j,i) is borrowed from the watflood.wfo code
!	in the calling code, simply define outwfo() with the values you
!     want to read to an .r2s file

      use area_watflood
	implicit none

	integer*4  :: i,j,ir,ios,iallocate
	CHARACTER(20) :: junk
	character(30) :: filename1
	logical       :: exists

      inquire(FILE=filename1,EXIST=exists)
      if(exists)then
        open(unit=ir,file=filename1,status='unknown')
	  write(*,*)'opened unit =',ir,' filename =',filename1
        read(ir,*)junk
        write(*,*)junk
        read(ir,*)junk
        write(*,*)junk
        read(ir,*)junk,application
        write(*,*)junk,application
        read(ir,*)junk,author
        write(*,*)junk,author
	  read(ir,*)junk
        write(*,*)junk
        read(ir,*)junk
        write(*,*)junk
        read(ir,*)junk,coordsys_temp
        write(*,*)junk,coordsys_temp
        read(ir,*)junk,xorigin_temp
        write(*,*)junk,xorigin_temp
        read(ir,*)junk,yorigin_temp
        write(*,*)junk,yorigin_temp
        read(ir,*)junk
        write(*,*)junk            
    	  read(ir,*,iostat=ios)junk,attribute_count
	  write(*,*)'attribute count ',attribute_count
        read(ir,*,iostat=ios)junk,attribute_units
        write(*,*)junk,attribute_units
        read(ir,*,iostat=ios)junk,attribute_name
        write(*,*)junk,attribute_name
        read(ir,*)junk
	  write(*,*)junk
        read(ir,*)junk,xcount_temp
	  write(*,*)junk,xcount_temp
        read(ir,*)junk,ycount_temp
	  write(*,*)junk,ycount_temp
        read(ir,*)junk,xdelta_temp
	  write(*,*)junk,xdelta_temp
        read(ir,*)junk,ydelta_temp
	  write(*,*)junk,ydelta_temp
        read(ir,*)junk                      !,angle_temp
        write(*,*)junk                      !,angle_temp
        read(ir,*)junk
        write(*,*)junk
        read(ir,*)junk
        write(*,*)junk
      
      
        allocate(inarray(ycount_temp,xcount_temp),stat=iAllocate)
            if(iAllocate.ne.0) STOP
     *      'Error with allocation of inarray in read_r2s'      
      
        do i=1,ycount_temp
          read(ir,*)(inarray(i,j),j=1,xcount_temp)
        end do
      else ! file does not exist
        print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'  
        print*,'WARNING'
        print*,'trying to open file ',filename1(1:30)
        print*,'but the file was not found'
        print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'  
      endif

      close(unit=ir,status='keep')
	return


	end subroutine read_r2s