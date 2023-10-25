
      MODULE area17

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
!     area17 is for the bsn.f program

      real, dimension(:),   allocatable :: distance
      real, dimension(:,:),   allocatable :: elv_2d,slope_2d,
     *                        frac_2d,frac_nca_2d,da_2d,da_nca_2d,
     *                        bnkfll_2d,ch_length_2d,
     *                        dummy,areamet,nca

      real, dimension(:,:,:), allocatable :: aclass_3d,fetch_2d

      integer, dimension(:,:), allocatable :: idummy,ijtemp

      integer, dimension(:),  allocatable :: distflg,channel,
     *         new,newxxx,newyyy,newrank,newnext

      integer, dimension(:,:),allocatable :: rank_2d,s_2d,next_2d,
     *                   last_2d,order_2d,
     *                   ielv,iak,irough_2d,ichnl_2d,
     *                   ireach_2d,ireach_2d_temp

      integer  :: split_no  
      real     :: elevconv

      END MODULE area17



