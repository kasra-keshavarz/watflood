     SUBROUTINE read_nca

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
!*****************************************************************************
!  READ_FLOW_EF - written Aug/06 by Dave Watson, CHC
!	- Derived from rdflow written by Nick Kouwen
!	This subroutine reads streamflow hydrograph (STR) file 
!	(tb0 format)
!*****************************************************************************

!     rev. 10.5.17 Aug.  10/23  = NK Added call_read_nca - can be used to decrease contributing grid area

    
      use area_watflood
      USE EF_Module
      
      Integer          :: rank,i,j,n,ios,iallocate
      real(4)          :: value,ylat,xlong      
      integer,  dimension(:),   allocatable :: nnnnn
!      real(4),  dimension(:),   allocatable :: nca




      allocate(nnnnn(na),stat=iAllocate)
      if(.not.allocated(nca_1d))then
           allocate(nca_1d(na),stat=iAllocate)
      endif
        
      nnnnn=0 
      nca_1d=0.0
      ios=0

      open(unit=99,file='basin\nca.xyz',iostat=ios,status='unknown')
    
      do while(ios.eq.0)
            read(99,*,iostat=ios)xlong,ylat,value
            if(ios.ne.0)exit
!			convert iy and jx to local coordinates
!			Jul. 21/04  nk
            i=int((ylat-yorigin)/ydelta)+1
            j=int((xlong-xorigin)/xdelta)+1

            n=s(i,j)
            if(n.ne.0)then
                nca_1d(n)=nca_1d(n)+value
                nnnnn(n)=nnnnn(n)+1
            endif
      end do
      close(unit=99)
      do n=1,naa
          if(nnnnn(n).gt.0)then
              nca_1d(n)=nca_1d(n)/float(nnnnn(n))
!              print*,'xxyyzz',n,nca_1d(n),nnnnn(n)
          endif
      end do
            
            
      end subroutine read_nca