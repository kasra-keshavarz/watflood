      SUBROUTINE synflw(time,thr)

!***********************************************************************
!    Copyright (C) 1987 by Nicholas Kouwen 
        
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
     
!***********************************************************************

!  THIS SUBROUTINE CALCULATES THE FLOW AT THE RECORDED FLOW TIMES, 
!  THEREFORE INTERPOLATES WHEN ROUTE DOES NOT GIVE FLOWS AT THE 
!  EVEN HOURS

!     REV. 7.6  -  SEDIMENT ROUTINE INCLUDED
!     REV. 9.00    Mar.  2000 - TS: CONVERTED TO FORTRAN 90

!  id  - is the storm identification number

!***********************************************************************

      use area_watflood
	implicit none

      integer*4    :: kr,kq,l,i,j,n,k,kkt
      real*4       :: time,thr

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE
      

      kr=int(time-thr)+1
      kq=int(time)
!     save the flows at the gage locations:
      
      do l=1,no
       if(inbsnflg(l).eq.1)then
          i=iy(l)
          j=jx(l)
          n=s(i,j)
          qsyn(l,kq)=qo2(n)
        else
          qsyn(l,kq)=-1.0  !added Dec. 14/10 nk
	  endif
      end do

!      print*,inbsnflg(3),s(iy(3),jx(3)),iy(3),jx(3),qo2(s(iy(3),jx(3)))
!      pause 11111
      
      do l=1,no
       if(inbsnflg(l).eq.1)then
        i=iy(l)
        j=jx(l)
        n=s(i,j)
        if(sedflg.eq.'y')then
          do k=kr,kq
            kkt=k/kt
            if(kkt.lt.1)kkt=1
!            rev. 7.6
             sedsyn(l,kkt)=yfinal(n)
             nitsyn(l,kkt)=nfinal(n)
             phssyn(l,kkt)=pfinal(n)
!            Calculate Mass(kg): kg/m3 * m3/s * 3600 (/1000 in tons for sed)
             sedmss(l,kkt)=yfinal(n)*qsyn(l,kkt)*3.6
             nitmss(l,kkt)=nfinal(n)*qsyn(l,kkt)*3600
             phsmss(l,kkt)=pfinal(n)*qsyn(l,kkt)*3600
          end do
        endif
       endif   !inbsnflg
      end do

!     SAVE THE FLOWS AT THE DAMAGE LOCATIONS:
      do l=1,ndam
!         i=iys(l)
!         j=jxs(l)
!         n=s(i,j)
        if(iys(l).le.0.or.jxs(l).le.0)then
!         this can happen if the stream gauge is outside the waterhsed
!         as when subwaterhseds are modelled as separate watersheds
!         added Mar 14/04 nk.          
          n=0
        else
          n=s(iys(l),jxs(l))
        endif
         if(n.eq.1)then
!          only if it is in the watershed           
           do k=kr,kq
             qloc(l,k)=qo2(n)-(qo2(n)-qo1(n))*(time-float(k))/thr
           end do
         endif
         end do
         
!     rev. 10.4.32 Jan.  10/21  = NK Added diversion output
!     Save the flows at the diversion outflow lacations
      if(wrtdiverflg)then
          do l=1,nDivOut
               qdiv(l,kq) =qo2(divGrid(l))
          end do    
      endif
         
         

! FORMATS

 5678 format(' ','in synflw: element ',2i5,' is outside basin')


      RETURN

      END SUBROUTINE synflw

