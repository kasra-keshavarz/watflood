      SUBROUTINE weight(nh,ng,ixr,iyr,cut)

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
!      include 'debug.for'

!***********************************************************************
! - THIS SUBROUTINE CALCULATES WEIGHTS FOR EACH GRID SQUARE FOR
!   EACH RAIN GAUGE.

!***********************************************************************

	use area_watflood
c      USE area16
      implicit none

      save

      integer(4)   :: nh,ng,is,i,j,iyr,ixr,ii,iallocate  !          iopt
      character*1  :: firstpass
	real*4       :: xa,yb,sumv,sumw,cut,old_radinfl
      integer      :: iAllocateStatus,ng_old
      real,    dimension(:,:),   allocatable :: weight_check
      real,    dimension(:,:),   allocatable :: separation
      real,    dimension(:),     allocatable :: nearest_distance
	logical, dimension(:,:,:), allocatable :: sta_relv
!     NW IS THE dstANCE WEIGHTING FACTOR
!!!!!!!!      nw=2    program argument   eg.   calmet 2

!              print*,' Calculating weights '

      DATA firstpass/'y'/
      DATA old_radinfl/-1.0/

      if(firstpass.eq.'y')then
    
c      print*,ng,'in weight'
c	pause 'in weight'
        ng_old=ng
        print*,'Allocating with ',ng,' stations'
        allocate(dst(ng),stat=iAllocate)
        if (iAllocate.ne.0) STOP 
     *    'Error with allocation of dst() in weight' 
        allocate(separation(ng,ng),stat=iAllocate)
        if (iAllocate.ne.0) STOP 
     *    'Error with allocation of separation() in weight' 
        allocate(nearest_distance(ng),stat=iAllocate)
        if (iAllocate.ne.0) STOP 
     *    'Error with allocation of nearest_distance() in weight' 
        allocate(weight_check(iyr,ixr),stat=iAllocateStatus)
        if (iAllocateStatus .ne. 0) STOP 
     *    '**Allocation failed for weight_check()' 
        allocate(sta_relv(ycount,xcount,ng),
     *        stat=iAllocateStatus)
        if (iAllocateStatus .ne. 0) STOP  
     *    '**Allocation failed for sta_relv array in weight @ 47**'  
      endif
      


      if(ng.ne.ng_old)print*,'ng,ng_old',ng,ng_old

      if(ng.gt.ng_old)then
        deallocate(w,dst,separation,nearest_distance,
     *         weight_check,sta_relv,stat=iAllocate)
        if (iAllocate.ne.0) STOP 
     *    'Warning: error with de-allocation in weight @ 58' 
        allocate(w(iyr,ixr,ng),
     *         dst(ng),separation(ng,ng),weight_check(iyr,ixr),
     *         sta_relv(ycount,xcount,ng),nearest_distance(ng),
     *                       stat=iAllocate)
        if (iAllocate.ne.0) STOP 
     *    'Error with re-allocation in weight @ 63' 
      endif


! INITIALIZE WEIGHTS

      do is=1,ng
        do i=1,iyr
          do j=1,ixr
            w(i,j,is)=0.0
	      sta_relv(i,j,is)=.true.
          end do
        end do
      end do

! IF THERE IS ONLY ONE STATION ALL WEIGHTS = 1.0

      if(ng.eq.1)then
        do i=1,iyr
          do j=1,ixr
            w(i,j,1)=1.0
	      do is=1,ng
	        sta_relv(i,j,is)=.true.
	      end do
          end do
        end do
	  firstpass='n'
        return
      endif

!	find nearest station to each
!     look at one station at a time
!     if rad_influence is given a value in the basin\ragmet.par file, 
!     it will be taken as a min value. 
!     If set as a large number, all statins will be used for all grids
!     Default is a large number if no file present
!     
	rad_influence=amax1(radinfl,0.0)

!     first, for each station find the nearest neighbour
      do is=1,ng
	  nearest_distance(is)=1.0e+32
	  do ii=1,ng     
          xa=xsta(is)-xsta(ii)
          yb=ysta(is)-ysta(ii)
          separation(is,ii)=sqrt(xa*xa+yb*yb)
	    if(is.ne.ii)nearest_distance(is)=
     *		min(nearest_distance(is),separation(is,ii))
        end do
c	   print*,is,nearest_distance(is),' cells'
	end do
!	then find the max separaton between nearest stations
!     and make this the radius of influence.
      old_radinfl=rad_influence
	do is=1,ng
	  rad_influence=amax1(rad_influence,nearest_distance(is)+2)
	end do

      if(Int(old_radinfl).ne.int(rad_influence))then
	  print*
	  print*,'NEW:'
        print*,'radius of influence ~',rad_influence*al/1000.,' km'
	  print*,'smoothing distance ~',smoothdist*al/1000.,' km'
	  print*,'---------------------------------------'
	  print*
      endif
      
! LOOP THROUGH GRID SQUARES
!     distances are in # of grids

      do i=1,iyr
         do j=1,ixr

! LOOP THROUGH  GAUGES

          ii=0
          sumw=0.0 
          do is=1,ng
            if(rrain(is,nh).ge.cut)then
              xa=xsta(is)-float(j)
              yb=ysta(is)-float(i)
!             min_distance added March 6/10 NK
!             this reduces the weight at a stations location
!             so provides smoothing
!             The larger the value, the more smoothing of the precip.
              dst(is)=sqrt(xa*xa+yb*yb)+smoothdist
!             dst is # grids
              
               if(dst(is).gt.rad_influence+2)sta_relv(i,j,is)=.false.
              
d      if(firstpass.eq.'y')then
d              write(51,6000)is,j,i,xsta(is),ysta(is),xa,yb,dst(is)
d6000          format(' IS,J,I,xsta,ysta,xa,yb,dst(IS)/',3i4,5f10.2)
d      endif

              if(sta_relv(i,j,is))then     
                if(dst(is).le.0.99)then
                  ii=ii+1
                else
                  sumw=sumw+1.0/dst(is)**nw
                endif
              endif
            else  
!             added Feb. 02/2020  NK   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  nbeeded here????????????                          
              sta_relv(i,j,is)=.false.
              
            endif
          end do

!        CUT IS THE VARIABLE CUTOFF VALUE. 
!        IF THE VALUE IS < CUT, THERE IS NO DATA

         do is=1,ng
            if(rrain(is,nh).ge.cut)then
              if(dst(is).le.0.99)then
                w(i,j,is)=1.0/float(ii)
              else
                if(sta_relv(i,j,is))then     
                  if(ii.lt.1)then
!                   THERE IS NO STATION IN THIS Grid
!                   write(6,6000)is,i,j,d(is),sumw
!6002                format(' IS,I,J,D(IS),sumw/',3i2,3f10.1)
                    w(i,j,is)=1.0/(dst(is)**nw)/sumw
                  endif
	          else
                  w(i,j,is)=0.0
	          endif
              endif
            endif
          end do

         end do
      end do

!	weight check - make sure they add up to 1.00000000000
      if(iopt.ge.1)then
        write(51,*)'weight check:'
        do i=1,iyr
	    do j=1,ixr
            weight_check(i,j)=0.0
          end do
	  end do

        do i=1,iyr
	    do j=1,ixr
            if(s(i,j).gt.0)then
!             only if it's in the basin
              do is=1,ng
                weight_check(i,j)=weight_check(i,j)+w(i,j,is)
              end do

	 	    if(weight_check(i,j).gt.1.01.or.
     *            weight_check(i,j).lt.0.99)then
d                write(*,*)'wt_check',i,j,s(i,j),weight_check(i,j)
      	    endif
	      endif
          end do
	  end do
	endif

d      if(firstpass.eq.'y'.and.iopt.ge.1)then
d        do is=1,ng
d          write(51,*)'recalculated weight for gauge no',is
d          do i=iyr,1,-1
d            write(51,6001)(w(i,j,is),j=1,ixr)
c            write(*,6001)(w(i,j,is),j=1,ixr)
d6001        format(999f5.2)
d          end do

d          write(51,*)
d        end do
d	endif

      if(iopt.ge.1)print*,'New weights hr= ',nh
      firstpass='n'

      return

      END SUBROUTINE weight




