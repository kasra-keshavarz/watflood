      function julian_day(year,month,day)
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
      
      implicit none
      save
      
      
      INTEGER        :: ju_mon(24),year,month,day
      integer     :: julian_day
      
      DATA ju_mon/1, 32, 60, 91,121,152,182,213,244,274,305,335,
     *            1, 32, 60, 91,121,152,182,213,244,274,305,335/

         julian_day=ju_mon(month)+day-1    
         julian_day=max(julian_day,1)
         if(mod(year,4).eq.0.and.month.ge.3)julian_day=julian_day+1
         
      return
      
      end function      