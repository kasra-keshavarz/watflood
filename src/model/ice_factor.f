      SUBROUTINE ice_factor

!***********************************************************************
!    Copyright (C) 2016 by Nicholas Kouwen  
        
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
      
!     REV. 10.1.16 Jan.  11/16  - NK: Added subroutine ice_factor.f
!     rev. 10.5.06 Apr.  07/23  = NK Fixed ice_fctr for resume

      use area_watflood
      use area_debug
      implicit none
      save

      integer     :: i,j,l,n,jz
      logical     :: firstpass
      real(4)     :: gap(9999),gapLast(9999),ice_fctr_last(9999)

      data firstpass/.true./
      if(debug_output)print*,' checkpoint 1 in ice_factor'

!     ice_factor is called if icerivflg and/or icelakeflg = y
!     Note:
!     dd_ice  &  dd_thaw calculated in etin.f
!     dd_ice reset julian day 243 = Aug. 31   Needs to be reset before 1st frost
!     dd_thaw reset julian day 60 = March 1   Needs to be reset before 1st +ve tmp.   

!     rev. 10.1.07 Dec.  02/15  - NK: Added ice_fctr(n) to route 
!     rev. 10.1.11 Dec.  11/15  - NK: Revised ice factor initialization and calculation   
!     the equations were taken from the Lake Winnipeg model
      if(firstpass)then
!     REV. 10.1.19 Jan.  15/16  - NK: Fixed initialization of ice_factr - moved from lake_ice > runof6
c!       default ice_fctr        
!     rev. 10.3.03 Jan.  21/20  = NK added dd_ice & dd_thaw to the resume.txt file      
          if(resumflg.ne.'y')then
!             initialize the  ice factor
              do n=1,naa
                  dd_ice(n)=0.0
                  dd_thaw(n)=0.0
                  ice_fctr(n)=1.0
              end do
          endif
        firstpass=.false.
        do n=1,naa
            gapLast(n)=99999.
            ice_fctr_last(n)=ice_fctr(n)
        end do
      endif
      if(debug_output)print*,' checkpoint 2 in ice_factor'
      
!      See rating curves:
!     Lake Winnipegosis.pdf      
      
!       Calculate the ice factors for each grid

      if(debug_output)print*,' checkpoint 2a in ice_factor'
      do n=1,naa
          gap(n)=dd_thaw(n)-dd_ice(n)
          if(gap(n).lt.gapLast(n))then
!             fall freeze up & winter           
              if(dd_ice(n).gt.0.0)then
!               dd_ice (freezing degree days) must be +ve   
!               equation for Waterhen:       
!               lif=3.0981177*dd_ice(n)**-0.3324723
                ice_fctr(n)=3.7222496*dd_ice(n)**-0.4
c                if(n.eq.nnprint)
c     *   write(444,*)jul_day_now,gap(nnPrint),gapLast(nnprint),
c     *           ice_fctr(n),'freeze'
              endif
      if(n.eq.nnprint)write(678,*)jul_day_now,gap(n),'freeze'
          else  !spring thaw & summer
              if(dd_thaw(n).gt.0.0)then
!               taken from Lake Winnipeg              
c               ice_fctr(n)=1.0/(3.7222496*dd_thaw(n)**-0.2747330)
                ice_fctr(n)=1.0/(3.7222496*dd_thaw(n)**-.3)
c              else
c                ice_fctr(n)=1.0
c                if(n.eq.nnprint)
c     *   write(444,*)jul_day_now,gap(nnPrint),gapLast(nnprint),
c     *           ice_fctr(n),'thaw'
              endif
c            endif
      if(n.eq.nnprint)write(678,*)jul_day_now,gap(n),'thaw'
          endif
      end do  ! NOTE: ice_fctr written in rte.txt fln(55)
      gapLast=gap
      if(debug_output)print*,' checkpoint 3 in ice_factor'

      do n=1,naa
!              ice_fctr(n)=0.1*ice_fctr(n)+0.9*ice_fctr_last(n)
             ice_fctr(n)=0.05*ice_fctr(n)+0.95*ice_fctr_last(n)
!               in the winter is should not be less than 0.25
             ice_fctr(n)=amax1(0.50,ice_fctr(n))
!               in the summer it should not be more than 1.0
             ice_fctr(n)=amin1(1.0,ice_fctr(n))
             ice_fctr_last(n)=ice_fctr(n)
      end do
      if(debug_output)print*,' checkpoint 4 in ice_factor'

          
!     rev. 10.1.44 Dec.  02/15  - NK: Reworked icerivflg & icelakeflg
      if(.not.icefactorfile)then
        if(icelakeflg.eq.'y')then
!         ice factors for lakes if the resrl\ice_factor.tb0 
!         file does NOT exists 
          do l=1,noresv
            i=ires(l)
            j=jres(l)
            n=s(i,j)
!           lake ice factor should probably not be less than 0.50
            lake_ice_factor(l,month_now)=amax1(0.90,ice_fctr(n))
          end do
c          print*,'lake ice factors are computed'
        endif
      endif
      if(debug_output)print*,' checkpoint 5 = end in ice_factor'
      
      end subroutine ice_factor
