      SUBROUTINE ftch()


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
! WRITTEN IN 2011 by Nick Kouwen
! Resvised Oct. 2013 NK
!
!!
!
!*************************************************************************

      USE area17

      use area_watflood
      implicit none

!     REAL    :: dummy(999,999)  moved to area17 Apr. 10/02 nk.
	integer :: i,j,n,reachcount,ii,jj,iAllocateStatus

!     count the number of reaches
      reachcount=0
      do i=1,ycount
        do j=1,xcount
	    if(ireach_2d(i,j).ne.0)then
	      reachcount=max(reachcount,ireach_2d(i,j))
	    endif
	  end do
	end do
	print*,'reachcount=',reachcount,' in ftch'
	  
c      allocate(fetch_2d(ycount,xcount,reachcount),stat=iAllocateStatus)
      allocate(fetch_2d(ycount,xcount,8),stat=iAllocateStatus)
      if (iAllocateStatus .ne. 0) STOP 
     *    '**Allocation failed  in fetch_2d @41**' 

      do i=1,ycount
        do j=1,xcount
          do jj=1,8
            fetch_2d(i,j,jj)=0.0
          end do
        end do
      end do

!     fetch to the north east      
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
!       are we in a reach?  
	  if(ireach_2d(i,j).ne.0)then
!         find fetch to the east
          fetch_2d(i,j,1)=0.7
          ii=i+1
          jj=j+1
	    do while(ii.lt.ycount.and.jj.lt.xcount)
	      if(ireach_2d(ii,jj).eq.ireach_2d(i,j))then
              fetch_2d(i,j,1)=sqrt(float((jj-j)**2+(ii-i)**2))+0.7
            else
!             have come to a shore so terminate the loop            
              ii=ycount
              jj=xcount
            endif
c         print*,ireach_2d(i,j),i,j,jj,int(fetch_2d(i,j,2))
            ii=ii+1
            jj=jj+1
          end do
        endif
      end do
d      write(51,*)'to the north east'
d      do i=ycount,1,-1
d        write(51,69200)(int(fetch_2d(i,j,1)),j=1,xcount)
d      end do
69200 format(999i5)      

!     fetch to the east      
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
!       are we in a reach?  
	  if(ireach_2d(i,j).ne.0)then
!         find fetch to the east
          fetch_2d(i,j,2)=0.5
          jj=j+1
	    do while(jj.lt.xcount)
	      if(ireach_2d(i,jj).eq.ireach_2d(i,j))then
              fetch_2d(i,j,2)=float((jj-j))+0.5
            else
              jj=xcount
c         print*,ireach_2d(i,j),i,j,jj,int(fetch_2d(i,j,2))
            endif
            jj=jj+1
          end do
        endif
      end do
d      write(51,*)'to the east'
d      do i=ycount,1,-1
d        write(51,69200)(int(fetch_2d(i,j,2)),j=1,xcount)
d      end do

!     fetch_2d to the south east      
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
!       are we in a reach?  
	  if(ireach_2d(i,j).ne.0)then
!         find fetch_2d to the s east
          fetch_2d(i,j,3)=0.7
          ii=i-1
          jj=j+1
	    do while(ii.gt.1.and.jj.gt.1)
	      if(ireach_2d(ii,jj).eq.ireach_2d(i,j))then
              fetch_2d(i,j,3)=sqrt(float((jj-j)**2+(ii-i)**2))+0.7
            else
!             have come to a shore so terminate the loop            
              ii=-1
              jj=xcount
            endif
c         print*,ireach_2d(i,j),i,j,ii,jj,int(fetch_2d(i,j,3))
            ii=ii-1
            jj=jj+1
          end do
        endif
      end do
d      write(51,*)'to the south east'
d      do i=ycount,1,-1
d        write(51,69200)(int(fetch_2d(i,j,3)),j=1,xcount)
d      end do

!     fetch_2d to the south      
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
!       are we in a reach?  
	  if(ireach_2d(i,j).ne.0)then
!         find fetch_2d to the south
          fetch_2d(i,j,4)=0.5
c         do ii=1,i
          ii=i
	    do while(ii.gt.1)
	      if(ireach_2d(ii,j).eq.ireach_2d(i,j))then
              fetch_2d(i,j,4)=float((i-ii))+0.5
            else
              ii=-1
c         print*,ireach_2d(i,j),i,j,ii,int(fetch_2d(i,j,4))
            endif
            ii=ii-1
          end do
        endif
      end do
d      write(51,*)'to the south'
d      do i=ycount,1,-1
d        write(51,69200)(int(fetch_2d(i,j,4)),j=1,xcount)
d      end do

!     fetch_2d to the south west      
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
!       are we in a reach?  
	  if(ireach_2d(i,j).ne.0)then
!         find fetch_2d to the sw
          fetch_2d(i,j,5)=0.7
          ii=i-1
          jj=j-1
	    do while(ii.gt.1.and.jj.gt.1)
	      if(ireach_2d(ii,jj).eq.ireach_2d(i,j))then
              fetch_2d(i,j,5)=sqrt(float((j-jj)**2+(i-ii)**2))+0.7
            else
!             have come to a shore so terminate the loop            
              ii=-1
              jj=-1
            endif
c         print*,ireach_2d(i,j),i,j,ii,jj,int(fetch_2d(i,j,5))
            ii=ii-1
            jj=jj-1
          end do
        endif
      end do
d      write(51,*)'to the south west'
d      do i=ycount,1,-1
d        write(51,69200)(int(fetch_2d(i,j,5)),j=1,xcount)
d      end do

!     fetch_2d to the west      
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
!       are we in a reach?  
	  if(ireach_2d(i,j).ne.0)then
!         find fetch_2d to the west
c	    do jj=1,xcount-j
c	    do jj=1,j
          fetch_2d(i,j,6)=0.5
          jj=j
	    do while(jj.gt.1)
	      if(ireach_2d(i,jj).eq.ireach_2d(i,j))then
              fetch_2d(i,j,6)=float((j-jj))+0.5
            else
              jj=-1
c         print*,ireach_2d(i,j),i,j,jj,int(fetch_2d(i,j,6))
            endif
            jj=jj-1
          end do
        endif
      end do
d      write(51,*)'to the west'
d      do i=ycount,1,-1
d        write(51,69200)(int(fetch_2d(i,j,6)),j=1,xcount)
d      end do

!     fetch_2d to the north west      
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
!       are we in a reach?  
	  if(ireach_2d(i,j).ne.0)then
!         find fetch_2d to the nw
          fetch_2d(i,j,7)=0.7
          ii=i+1
          jj=j-1
	    do while(ii.gt.1.and.jj.gt.1)
	      if(ireach_2d(ii,jj).eq.ireach_2d(i,j))then
              fetch_2d(i,j,7)=sqrt(float((j-jj)**2+(i-ii)**2)+0.7)
            else
!             have come to a shore so terminate the loop            
              ii=ycount
              jj=-1
            endif
c         print*,ireach_2d(i,j),i,j,ii,jj,int(fetch_2d(i,j,7))
            ii=ii+1
            jj=jj-1
          end do
        endif
      end do
d      write(51,*)'to the north west'
d      do i=ycount,1,-1
d        write(51,69200)(int(fetch_2d(i,j,7)),j=1,xcount)
d      end do

!     fetch_2d to the north      
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
!       are we in a reach?  
	  if(ireach_2d(i,j).ne.0)then
!         find fetch_2d to the north
          fetch_2d(i,j,8)=0.5
          ii=i+1
	    do while(ii.lt.ycount)
	      if(ireach_2d(ii,j).eq.ireach_2d(i,j))then
              fetch_2d(i,j,8)=float((ii-i))+.5
            else
              ii=ycount
c         print*,ireach_2d(i,j),i,j,ii,int(fetch_2d(i,j,8))
            endif
            ii=ii+1
          end do
        endif
      end do
d      write(51,*)'to the north'
d      do i=ycount,1,-1
d        write(51,69200)(int(fetch_2d(i,j,8)),j=1,xcount)
d      end do


      RETURN

      END SUBROUTINE ftch
