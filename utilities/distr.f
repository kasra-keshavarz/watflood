      SUBROUTINE distr(nh,ixr,iyr,ng)

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
!      include 'DEBUG.for'

!***********************************************************************
! - THIS SUBROUTINE DISTRIBUTES RAINFALL ACCORDING TO WEIGHTS.

!***********************************************************************

c      USE area16
	use area_watflood
     
      do i=1,iyr
        do j=1,ixr
           do is=1,ng
                 p(i,j)= p(i,j)+w(i,j,is)*rrain(is,nh)
           end do
        end do
      end do

c      write(6,*)
c      write(6,6002)ng,nh,(rrain(is,nh),is=1,ng)
c6002      format('+','NG,NH,VARBL/',2i5,999f5.1)
c      pause
c      do i=iyr,1,-1
c      write(6,6000)(p(i,j),j=1,ixr)
c6000      format(' ',999f5.1)
c      end do
c      pause

      return

      END SUBROUTINE distr
