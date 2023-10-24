      SUBROUTINE arrange()


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
!*************************************************************************
! WRITTEN IN 1972
!
! ARRANGE - ASSIGNMENT OF ELEMENT YYY, XXX AND CALCULATION OF DRAINAGE 
!            AREAS.  THE ELEMENT WITH THE HIGHEST ELEVATION WILL HAVE 
!            ORDER=1 AND THE ATIONS WILL OCCUR TO SUCCESSIVELY LOWER 
!            ELVATIONS HAVING HIGHER
!
! Modified by Tricia Stadnyk - September 2000
! Converted common blocks to modules and added dynamically allocated
! run-time arrays as part of Fortran 90 conversion.
!
!
!*************************************************************************

      USE area17

      use area_watflood

!     REAL    :: dummy(999,999)  moved to area17 Apr. 10/02 nk.
      REAL(4) :: maxelv
	integer :: i,j,n,max_order,order_start,order_end
	logical :: not_done
	character*1 :: answer


! TO CALCULATE STORAGE/DISCHARGE,MUST FIND DRAINAGE AREA ABOVE A
! DEFINE DA(I,J) AS INTEGER NO. OF SQUARE KILOMETERS
      do i=1,ycount
         do j=1,xcount
            dummy(i,j)=elv_2d(i,j)*1.0
         end do
      end do

      do i=1,ycount
         do j=1,xcount
            rank_2d(i,j)=0
            next_2d(i,j)=0
            last_2d(i,j)=0
            order_2d(i,j)=0
         end do
      end do


! ELEVATIONS ARE ENTERED OAND STRED 
! THE FOLLOWING SEQUENCE WILL ORDER ARRAY  SUBROUTINE AMAXO FIND
! T INTEGER IN AN ARRAY OF INTEGERS
      n=0
      no=0

   12 maxelv=0.0

!       Find the highest elevation in the domain
        do i=1,ycount
           do j=1,xcount
              maxelv=max(dummy(i,j),maxelv)
           end do
        end do

        if(maxelv.eq.0.0) GO TO 13    !  <<<<<<<<<<<<<<<<<<<<<<<<go to!

        do i=1,ycount
          do j=1,xcount
            if(dummy(i,j).eq.maxelv)then
    5         n=n+1
!             add to total no of grids if there is an elv.
              yyy(n)=i
              xxx(n)=j
              rank_2d(i,j)=n
              dummy(i,j)=0
!             find the number of outlets
!             if a grid has zero area but has an elv, it is part of a watershed
              if(frac_2d(i,j).le.0.0)then
                  no=no+1
                  write(*,6000)'outlet grid rank,row,col,frac',
     *               n,i,j,frac_2d(i,j) 
6000              format(a30,4i5)                  
              endif
            endif
          end do
        end do
   10   CONTINUE


      GO TO 12              ! <<<<<<<<<<<<<<<<<<<<<<<<go to!

! INITIAL VALUE IS SET EQUAL ZERO,THEN START AT ORDER=1 AND SUM
! AREA GOING TO SUCCESSIVE LOWER ELVATIONS IN THE WATERSHED
! DA(I,J) IS THE DRAINAGE AREA ABOVE THE OUTLET OF ELEMENT(I,J)
! DA IS SQ. KM AND CAN BE  < 1.0
! DEFINE DA AND S =0 FOR ELEMENTS OUTSIDE THE WATERSHED
! S IS AN INDICATOR NE=1,E=2,SE=3,S=4,SW=5,W=6,NW=7,N=8. S IS TH
! ION OF FLOW OUTOF THE SQUARE GRID ELEMENT IN 8 POSSIBLE DIRECT


   13 na=n
      naa=n-no
      
      print*,'We have found ',no,'outlets'
      print*,'If this is not correct, check that all receiving grids'
      print*,'have the same elv &  > 0.0  (e.g.0.001) and frac = 0'
      pause 'Hit enter to continue or Ctrl C to quit'
      

!     REORDER
!     find all first order grids (i.e. no upstream contribution
!     this means no drainage direction into this grid
!  Version 11.0 - Jul. 2016 - Reordering of grids for // processing 
      
      not_done=.true.
      
      print*
      print*,'Arrange grids for // computing y/n ?'
      print*
      read*,answer
      if(answer.eq.'y')then
      

      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        order_2d(i,j)=1   !Assume order > 0
!       check 8 directions if an arrow points in and if so, 
!       make it 0 again
!       what's left is order 1
!       avoid looking outside the edges
        if(i-1.gt.1.and.j-1.gt.1)then
          if(s_2d(i-1,j-1).eq.1)order_2d(i,j)=0
        endif
        if(j-1.gt.1)then
          if(s_2d(i  ,j-1).eq.2)order_2d(i,j)=0
        endif
        if(i+1.lt.ycount.and.j-1.gt.1)then
          if(s_2d(i+1,j-1).eq.3)order_2d(i,j)=0
        endif
        if(i+1.lt.ycount)then
          if(s_2d(i+1,j  ).eq.4)order_2d(i,j)=0
        endif
        if(i+1.lt.ycount.and.j+1.lt.xcount)then
          if(s_2d(i+1,j+1).eq.5)order_2d(i,j)=0
        endif
        if(j+1.lt.xcount)then
          if(s_2d(i  ,j+1).eq.6)order_2d(i,j)=0
        endif
        if(i-1.gt.1.and.j+1.lt.xcount)then
          if(s_2d(i-1,j+1).eq.7)order_2d(i,j)=0
        endif
        if(i-1.gt.1)then
          if(s_2d(i-1,j  ).eq.8)order_2d(i,j)=0
        endif
      end do
      do i=ycount,1,-1
        write(1000,1000)(order_2d(i,j),j=1,xcount)
1000  format(<xcount>i4)
      end do   
      
!     Work downstream and increase order order sequentially
!     from each 1st order grid            
      max_order=0
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
!       check 8 directions if an arrow points in and if so, 
!       take the max of this order or u/s order +1
        if(i-1.gt.1.and.j-1.gt.1)then
          if(s_2d(i-1,j-1).eq.1)
     *     order_2d(i,j)=max0(order_2d(i,j),order_2d(i-1,j-1)+1)
        endif
        if(j-1.gt.1)then
          if(s_2d(i  ,j-1).eq.2)
     *     order_2d(i,j)=max0(order_2d(i,j),order_2d(i  ,j-1)+1)
        endif
        if(i+1.lt.ycount.and.j-1.gt.1)then
          if(s_2d(i+1,j-1).eq.3)
     *     order_2d(i,j)=max0(order_2d(i,j),order_2d(i+1,j-1)+1)
        endif
        if(i+1.lt.ycount)then
          if(s_2d(i+1,j  ).eq.4)
     *     order_2d(i,j)=max0(order_2d(i,j),order_2d(i+1,j  )+1)
        endif
        if(i+1.lt.ycount.and.j+1.lt.xcount)then
          if(s_2d(i+1,j+1).eq.5)
     *     order_2d(i,j)=max0(order_2d(i,j),order_2d(i+1,j+1)+1)
        endif
        if(j+1.lt.xcount)then
          if(s_2d(i  ,j+1).eq.6)
     *     order_2d(i,j)=max0(order_2d(i,j),order_2d(i  ,j+1)+1)
        endif
        if(i-1.gt.1.and.j+1.lt.xcount)then
          if(s_2d(i-1,j+1).eq.7)
     *     order_2d(i,j)=max0(order_2d(i,j),order_2d(i-1,j+1)+1)
        endif
        if(i-1.gt.1)then
          if(s_2d(i-1,j  ).eq.8)
     *      order_2d(i,j)=max0(order_2d(i,j),order_2d(i-1,j  )+1)
        endif
      end do
      
      do i=1,ycount
        do j=1,xcount
          max_order=amax0(max_order,order_2d(i,j))
        end do
      end do
      
      do n=naa+1,na
        i=yyy(n)
        j=xxx(n)
        order_2d(i,j)=max_order+1
      end do

      do i=ycount,1,-1
        write(1000,10011)(order_2d(i,j),j=1,xcount)
      end do   

 
!     Now the hard part. Re-order the grids 
c        max_order=order_2d(yyy(naa),xxx(naa))
        print*,na,order_2d(yyy(naa),xxx(naa))
        
!     The original order from the first part        
      do i=ycount,1,-1
        write(1000,10011)(rank_2d(i,j),j=1,xcount)
10011  format(<xcount>i5)
      end do   

!     destroy the old ranks - just in case
      do i=1,ycount
        do j=1,xcount
          rank_2d(i,j)=0
        end do
      end do

!     Re-number the grids - new ranks 
      open(unit=99,file='LoopsByOrder.txt',status='unknown') 
      write(99,*)'      GridOrder     Start'
      n=0
      order_start=1
      order_end=1
      order_last=0
      do k=1,max_order+1
        do i=1,ycount
          do j=1,xcount
            if(order_2d(i,j).eq.k)then
              n=n+1
              yyy(n)=i
              xxx(n)=j
              rank_2d(i,j)=n
!             find the number of outlets (no)
!             if a grid has zero area but has an elv, it is part of a watershed
              write(1003,*)order_2d(i,j),rank_2d(i,j)  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
              
              if(k.ne.order_last)then
                write(99,*)k,n
                order_last=k
              endif

            endif
          end do
        end do
      end do
      close(unit=99,status='keep')
!  END      
!  Version 11.0 - Jul. 2016 - Reordering of grids for // processing 
      
      
      end if    !  if answer // 
      
      
!     compute the drainage area for each grid:
      do n=1,naa
         i=yyy(n)
         j=xxx(n)
c         da_2d(i,j)=step2*frac_2d(i,j)
!        frac is now the gridarea in m**2
!        da_2d is the drainage area in km**2
         da_2d(i,j)=frac_2d(i,j)*1.0e-06
c         if(da_2d(i,j).le.0.0.or.frac_2d(i,j).le.0.0)then
!           write(*,1012)i,j
!           write(39,1012)i,j
!         endif        
!         write(*,1011)n,i,j,da(i,j),i-1,i+1,j-1,j+1
        if(i.eq.1.or.j.eq.1.or.i.eq.ycount.or.j.eq.xcount)then
!         this can happen if people did NOT leave a blank edge!!!
          print*,'WARNING: missing blank edge(s) - crash possible'
        else
          if(s_2d(i-1,j-1).eq.1)da_2d(i,j)=da_2d(i,j)+da_2d(i-1,j-1)
          if(s_2d(i  ,j-1).eq.2)da_2d(i,j)=da_2d(i,j)+da_2d(i  ,j-1)
          if(s_2d(i+1,j-1).eq.3)da_2d(i,j)=da_2d(i,j)+da_2d(i+1,j-1)
          if(s_2d(i+1,j  ).eq.4)da_2d(i,j)=da_2d(i,j)+da_2d(i+1,j  )
          if(s_2d(i+1,j+1).eq.5)da_2d(i,j)=da_2d(i,j)+da_2d(i+1,j+1)
          if(s_2d(i  ,j+1).eq.6)da_2d(i,j)=da_2d(i,j)+da_2d(i  ,j+1)
          if(s_2d(i-1,j+1).eq.7)da_2d(i,j)=da_2d(i,j)+da_2d(i-1,j+1)
          if(s_2d(i-1,j  ).eq.8)da_2d(i,j)=da_2d(i,j)+da_2d(i-1,j  )
        endif
      end do

!     compute the non-contributing area above each grid
      if(non_ca_flg)then
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
!         frac is now the gridarea in m**2
!         da_nca_2d is the drainage area in km**2
          da_nca_2d(i,j)=frac_nca_2d(i,j)*1.0e-06
          if(i.eq.1.or.j.eq.1.or.i.eq.ycount.or.j.eq.xcount)then
!         this can happen if people did NOT leave a blank edge!!!
            print*,'WARNING: missing blank edge(s) - crash possible'
          else
            if(s_2d(i-1,j-1).eq.1)da_nca_2d(i,j)=
     *                             da_nca_2d(i,j)+da_nca_2d(i-1,j-1)
            if(s_2d(i  ,j-1).eq.2)da_nca_2d(i,j)=
     *                             da_nca_2d(i,j)+da_nca_2d(i  ,j-1)
            if(s_2d(i+1,j-1).eq.3)da_nca_2d(i,j)=
     *                             da_nca_2d(i,j)+da_nca_2d(i+1,j-1)
            if(s_2d(i+1,j  ).eq.4)da_nca_2d(i,j)=
     *                             da_nca_2d(i,j)+da_nca_2d(i+1,j  )
            if(s_2d(i+1,j+1).eq.5)da_nca_2d(i,j)=
     *                             da_nca_2d(i,j)+da_nca_2d(i+1,j+1)
            if(s_2d(i  ,j+1).eq.6)da_nca_2d(i,j)=
     *                             da_nca_2d(i,j)+da_nca_2d(i  ,j+1)
            if(s_2d(i-1,j+1).eq.7)da_nca_2d(i,j)=
     *                             da_nca_2d(i,j)+da_nca_2d(i-1,j+1)
            if(s_2d(i-1,j  ).eq.8)da_nca_2d(i,j)=
     *                             da_nca_2d(i,j)+da_nca_2d(i-1,j  )
          endif
        end do
      endif

! FORMATS:

 1001 format(20I4)
!1107 format(/' sampled elevations of the channel bottoms')
!1010 format(' MAXELV = ',I5/)
 1011 format(' ',3I5,F7.2,4I5)
 1012 format(' da or frac <= 0.0 in i,j=',2i5,'- please fix')


      RETURN

      END SUBROUTINE arrange
