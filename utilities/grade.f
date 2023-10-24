      SUBROUTINE grade(slopemin,no_errors)

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
!*********************************************************************
! WRITTEN IN 1972!
!
! GRADEA - THIS SUBROUTINE CALCULATES THE slope OF THE CHANNEL AT THE
!          OUTLET OF EACH SQUARE ELEMENT
!
! Modified by Tricia Stadnyk - September 2000
! Converted common blocks to modules and added dynamically allocated
! run-time arrays as part of Fortran 90 conversion.
!
!
!*********************************************************************

c      USE area1
c      USE area2
c	USE area3
	USE area17

      use area_watflood

      real  ::   slopemin
	integer :: no_errors,n,i,j,ind

! FOR METER CONTOURS elvCONV WILL BE 1.0
! FOR FT CONTOURS elvCONV WILL BE 0.305


!     This section code upgrade NK 11/09/04 

      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        ind=s_2d(i,j)
        if(s_2d(i,j).eq.1)then
          next_2d(i,j)=rank_2d(i+1,j+1)
          ch_length_2d(i,j)=al*1.4142
          slope_2d(i,j)=-(elv_2d(i+1,j+1)-elv_2d(i,j))/ch_length_2d(i,j)
!          slope(i,j)=-slope(i,j)
        elseif(s_2d(i,j).eq.2)then
          next_2d(i,j)=rank_2d(i,j+1)
          ch_length_2d(i,j)=al
          slope_2d(i,j)=-(elv_2d(i  ,j+1)-elv_2d(i,j))/ch_length_2d(i,j)
!         slope(i,j)=-slope(i,j)
        elseif(s_2d(i,j).eq.3)then
          next_2d(i,j)=rank_2d(i-1,j+1)
          ch_length_2d(i,j)=al*1.4142
          slope_2d(i,j)=-(elv_2d(i-1,j+1)-elv_2d(i,j))/ch_length_2d(i,j)
!          slope(i,j)=-slope(i,j)
        elseif(s_2d(i,j).eq.4)then
          next_2d(i,j)=rank_2d(i-1,j)
          ch_length_2d(i,j)=al
          slope_2d(i,j)=-(elv_2d(i-1,j  )-elv_2d(i,j))/ch_length_2d(i,j)
!          slope(i,j)=-slope(i,j)
        elseif(s_2d(i,j).eq.5)then
          next_2d(i,j)=rank_2d(i-1,j-1)
          ch_length_2d(i,j)=al*1.4142
          slope_2d(i,j)=-(elv_2d(i-1,j-1)-elv_2d(i,j))/ch_length_2d(i,j)
!          slope(i,j)=-slope(i,j)
        elseif(s_2d(i,j).eq.6)then
          next_2d(i,j)=rank_2d(i,j-1)
          ch_length_2d(i,j)=al
          slope_2d(i,j)=-(elv_2d(i  ,j-1)-elv_2d(i,j))/ch_length_2d(i,j)
!         slope(i,j)=-slope(i,j)
        elseif(s_2d(i,j).eq.7)then
          next_2d(i,j)=rank_2d(i+1,j-1)
          ch_length_2d(i,j)=al*1.4142
          slope_2d(i,j)=-(elv_2d(i+1,j-1)-elv_2d(i,j))/ch_length_2d(i,j)
!          slope(i,j)=-slope(i,j)
        elseif(s_2d(i,j).eq.8)then
          next_2d(i,j)=rank_2d(i+1,j)
          ch_length_2d(i,j)=al
          slope_2d(i,j)=-(elv_2d(i+1,j  )-elv_2d(i,j))/ch_length_2d(i,j)
!          slope(i,j)=-slope(i,j)
        endif

      
        slope_2d(i,j)=amax1(slopemin,slope_2d(i,j),0.0000001)


        if(next_2d(i,j).le.n)then
	     no_errors=no_errors+1
!         write(*,7771)n,i,j,next(i,j)
          if(no_errors.eq.1)then
            write(39,7771)n,i,j,next_2d(i,j)
          else
            write(39,7772)n,i,j,next_2d(i,j)
          endif
          next_2d(i,j)=n+1
	    slope_2d(i,j)=-slope_2d(i,j)
       endif  
      
      
      write(1004,*)n,next_2d(i,j)
      

      end do
      open(unit=1002,file='next.xyz')
      do i=1,ycount
        do j=1,xcount
          write(1002,*)xorigin+xdelta/2.0+(j-1)*xdelta,
     *                 yorigin+ydelta/4.0+(i-1)*ydelta,
     *                 rank_2d(i,j),next_2d(i,j),order_2d(i,j)
        end do
      end do
      close(unit=1002)
      
      open(unit=1005,file='rank_n2.xyz')
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        write(1005,*)xorigin+xdelta/2.0+(j-1)*xdelta,
     *                 yorigin+ydelta/4.0+(i-1)*ydelta,
     *                 rank_2d(i,j),next_2d(i,j),order_2d(i,j)
      end do
      close(unit=1005)
      


! FORMATS:

 7771 format(' ERROR: grid does not drain to lower grid'/
     *  ' next(i,j) set to n+1 <<<<<<FIX FIX FIX'/ 
     *  ' please fix - location n,row,col,next(i,j)/',4i5)
 7772 format(' please fix - location n,row,col,next(i,j)/',4i5) 


      RETURN

      END SUBROUTINE grade
