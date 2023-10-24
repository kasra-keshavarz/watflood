!***********************************************************************
      subroutine date_time(cday,time)

      
!***********************************************************************
!    Copyright (C) 2003 by Nicholas Kouwen 
        
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
      

!	dummmy s/r for unix
	implicit none
      save

      CHARACTER(10) :: time
      CHARACTER(8)  :: cday
      call date_and_time(cday,time)

      RETURN
c
      END subroutine date_time
