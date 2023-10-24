      subroutine rt_channel(div,thr,dtmin,jz,iz,time,date,n)
      
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
     
!     rev. 10.2.08 Nov.  04/17  - NK: New rt_channel & rt-_wetland subroutines      
!     rev. 10.2.44 Jan.  21/19  - NK: Added qOld() to allow previous weighted outflow on the resume file
      
!             CHANNEL ROUTING:
!             ADD UPSTREAM CONTRIBUTIONS AND LOCAL INFLOW:
!              if(aclass(n,classcount-2).gt.0.0)

      use area_watflood
      implicit none

      SAVE          ! SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT

      INTEGER :: istate(400,400),rbin,lsta,nnn1,jz,n,lll,l,iz,jjz,&
                  i,j,ii,ic,jj,ios,itracker,unt,ln,n_dt_min,&
                  month_last,ic_max
      integer :: iallocate
      REAL    :: oldwet,convert,wth,old_at,aaa
      REAL(4) :: time,newstore,try1,try2,try3,div,thr,at,dtmin,&
                  wt,atemp,tfA,cfA,at1,dt_min_n,&
                  hold,whold,strlossvol,qtest
      integer :: frame_no1
!        REAL(4), DIMENSION(:) :: att(10000)
      character*1 :: flowsta,flwinitflg
      character*80 :: junk
!        REAL    :: qdlast(6)
      character(14) :: date
      logical  :: printmsg,firstpass,nantest
      
!     needs to fixed for dynamic mem
      DATA month_last/0/
      DATA istate /160000*0/
      DATA firstpass/.true./
      data printmsg/.false./
      data nantest/.false./
      
      REAL(4), DIMENSION(:), ALLOCATABLE :: tdum
      
      
!     rev. 10.1.98 Oct   04/17  - NK: Deal with -ve flows in route
      if(iopt99)then
        if(strlossvol*2.0.gt.store1(n)+(qi1(n)+qi2(n))*div)then
           write(98,*)'Warning:',time,'strloss set = 0.0 for grid ',n
           strloss(n)=0.0
           warning(2)=.true.
        endif
      endif
!####################################################################################      

!     rev. 9.9.59  Feb.  18/15  - NK: In route: strloss option fracflg y/n
!     rev. 10.1.62 Jan.  08/17  - NK: Checkup on strloss effect on low flows

!     rev. 9.2.11  Sep.  15/05  - NK: added Manning's n  r1n & r2n
!     rev. 9.5.06  Feb.  05/08  - NK: added pool and pool_o in rdpar & route
!             pool-riffle storage added
!             pool = % of bankfull area attributed to dead storage
!     tfA = total water cross-sectional area
!     cap = full channel volume  m**3
!     cfA = cap/rl = channel cross-section area  m**2
!     chaxa = bankfull area               
                
!     if(firstpass)then
!         open(unit=405,file='405.csv') 
!         open(unit=507,file='507.csv') 
!         open(unit=566,file='566.csv') 
!         open(unit=641,file='641.csv') 
!         open(unit=869,file='869.csv') 
!         open(unit=1835,file='1835.csv') 
!         write(611,*)'totaltime,tfA,cfa ,hfp,hcha2,hfp+hcha2,strloss,qo2,qi2,qr,qstream'
!     endif
      
      
     qi2(n)=qi2(n)+qr(n)+qstream(n)-strloss(n)
!     qi2(n)=qi2(n)+qr(n)-strloss(n)

!     qi2(n)=amax1(0.0,qi2(n))
      
      lll=next(n)
      qold(n)=qo1(n)
      hold=1.0e+25
      do ic=1,20
!       UP TO 20 ITERATIONS ARE ALLOWED
        if(abs(hold-store2(n)).gt.0.0003*hold)then
            hfp(n)=0.0
!         THIS ITERATES TO 3% OR ALLOWS UP TO 15 ITERATIONS
          if(store2(n).le.0.0)then
!           NO OUTFLOW - CHANNEL IS EMPTY
            tfA=0.0                  ! = xa
            qo2(n)=0.001
          else     ! (store2(n).le.0.0)
            over(n)=(store2(n)-cap(n))/rl(n)
            tfA=store2(n)/rl(n)      ! = xa
            if(over(n).le.0.0)then !+++++++++++++++++++++++++++++++++++++++++++++
!             CHANNEL FLOW ONLY
              cfA=cap(n)/rl(n)   ! added Jun 21/06 nk  but not ysed - only reported below
              if(tfA/chaxa(n).gt.pool(n))then
!               calculate flow for area above pool area:                  
                qo2(n)=(tfA-pool(n)*chaxa(n))**1.67*slope(n)/chawid(n)**0.667/r2n(n)*ice_fctr(n)
                NaNtest=ISNAN(qo2(n))
                if(NaNtest)then
                      write(98,*)'Error: channel routing NaN in grid ',n,'at t= ',totaltime
                      write(98,*)'Error: Check upstream for abrubt flow changes'
                      write(98,*)'Error: Eg. ponds over ~10% of grid area'
                      write(*,*)'Error: channel routing NaN in grid ',n,'at t= ',totaltime
                      write(*,*)'Error: Check upstream for abrubt flow changes'
                      write(*,*)'Error: Eg. ponds over ~10% of grid area'
                endif

!        if(n.eq.3315)write(886,887)-jz,tfA,pool(n),chaxa(n),slope(n),chawid(n),r2n(n),ice_fctr(n),qo2(n)
!        if(n.eq.3315)write(884,887)-totaltime,(tfA-pool(n)*chaxa(n))**1.67*slope(n)/chawid(n)**0.667/r2n(n)*ice_fctr(n)
  
                if(.not.qo2(n).eq.qo2(n))Print*,'a',jx,n,tfA,pool(n)
              else
                qo2(n)=0.0
              endif
              hcha2(n)=tfA/chawid(n)
!      if(n.eq.3315)write(887,887)totaltime,tfA,cfa,hfp(n),over(n),qo2(n)
887  format(99(f12.3,',  '))
      
            else     ! (over(n).le.0.0)!+++++++++++++++++++++++++++++++++++++++++++++

!             CHANNEL + FLOOD PLAIN FLOW
!     rev. 9.2.43  Jun.  21/06  - NK: fixed spikes in route
!              ax=cap(n)/rl(n)   ! added Jun 21/06 nk             WRONG!
!              cfA=cap(n)/rl(n)   ! added Jun 21/06 nk   Channel Flow Area
!     rev. 9.2.11  Sep.  15/05  - NK: added Manning's n  r1n & r2n
!     rev. 9.8.60  May   14/13  - NK: fixed ice factor for whole x-section
!             0.17 factor is based on 100:1 fp w/d ratio
!             flood plain width/depth assumes as 100
!             use quadratic equation to solve for fp. depth
!             400 from 4*100=4*FPwidth              
              hfp(n)=(-1.0+sqrt(1.+400.0*over(n)))/200.0
!             hcha2(n) is the bankfull depth here
              hcha2(n)=chaxa(n)/chawid(n)
!             cfA is the total main channel flow cross section area
              cfA=(hfp(n)+hcha2(n)*pool(n))*chawid(n)
!             Note: sq. root of slope taken earlier              
              if(over(n)-hfp(n)*chawid(n).gt.0.0)then  ! overbank
                qo2(n)=&                               ! cfA = Channel Flow Area
                   ((cfA-pool(n)*chaxa(n))**1.67*slope(n)/chawid(n)**0.667/r2n(n)&
                    +fpFactor(n)/2.6667**0.6667*hfp(n)**2.6667*slope(n)/r1n(n))*ice_fctr(n)
                    
!                   +(over(n)-hfp(n)*chawid(n))**1.33&
!                            *slope(n)*fpFactor(n)/r1n(n))*ice_fctr(n)

                if(.not.qo2(n).eq.qo2(n))Print*,'b',jx,n,tfA,pool(n)
              else
                qo2(n)=cfA**1.67*slope(n)/chawid(n)**0.667/r2n(n)*ice_fctr(n)
                if(.not.qo2(n).eq.qo2(n))Print*,'c',jx,n,tfA,pool(n)
              endif     ! (over(n)-hfp(n)*chawid(n).gt.0.0)
              NaNtest=ISNAN(qo2(n))
              if(NaNtest)then
                      write(98,*)'Error: channel routing NaN in grid ',n,'at t= ',totaltime
                      write(98,*)'Error: Check upstream for abrubt flow changes'
                      write(98,*)'Error: Eg. ponds over ~10% of grid area'
                      write(*,*)'Error: channel routing NaN in grid ',n,'at t= ',totaltime
                      write(*,*)'Error: Check upstream for abrubt flow changes'
                      write(*,*)'Error: Eg. ponds over ~10% of grid area'
                      error_msg=2
              endif
!      if(n.eq.386)write(887,887)totaltime,tfA,cfa,hfp(n),hch(n),strloss(n)
!      if(n.eq.3315)write(884,887)totaltime,&
!                  cfA**1.67*slope(n)/chawid(n)**0.667/r2n(n),&
!                  fpFactor(n)/2.6667**0.6667*hfp(n)**2.6667*slope(n)/r1n(n),&
!                 fpFactor(n)
            endif      ! (over(n).le.0.0)!+++++++++++++++++++++++++++++++++++++++++++++
!     rev. 9.8.12  Dec.  07/11  - NK: removed 30 char limit on find filetype 
            wt=amax1(0.5,float(ic)/21.)
!            wt=amax1(0.5,0.5+float(ic)/41.)
            wt=amin1(1.0,wt)
            qo2(n)=(1.0-wt)*qo2(n)+wt*qold(n)
            qold(n)=qo2(n)
            endif        ! (store2(n).le.0.0)
            hold=store2(n)
!     rev. 9.8.49  Feb.  20/13  - NK: Added n=municipal & irrigation withdrawals
            
            store2(n)=store1(n)+(qi1(n)+qi2(n)-qo1(n)-qo2(n)-2.0*qwdr(n,month_now))*div
!                     qwdr multiplied by 2 to offset /2 in div     

!     rev. 9.9.52  Jan.  14/15  - NK: Fixed bug for channel store < 0 for withdrawals
          if(store2(n).le.0.0)then
!           cut off withdrawals                  
            store2(n)=store1(n)+(qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
            if(store2(n).le.0.0)then
!             if store2 still = 0 then reduce calculated outflow
!             but this should not happen as it never did before 
!             withdrawals were added to the code
              store2(n)=0.0
!              qo2(n)=store1(n)/div+qi1(n)+qi2(n)-qo1(n)         
              qo2(n)=0.0        
            endif
          endif
        else
            if(iopt99)then   
               if(.not.qo2(n).eq.qo2(n))pause 88888
            endif
!          convergence to 3%
           GO TO 16
         endif
         
      end do    ! ic=1,20
       
     
      
16    continue
      
      
!      if(n.eq.405)write(405,887)totaltime,tfA,cfa,hfp(n),hcha2(n),hfp(n)+hcha2(n),strloss(n),qo2(n),qi2(n),qr(n),qstream(n),area(405)
!      if(n.eq.507)write(507,887)totaltime,tfA,cfa,hfp(n),hcha2(n),hfp(n)+hcha2(n),strloss(n),qo2(n),qi2(n),qr(n),qstream(n),area(507)
!      if(n.eq.566)write(566,887)totaltime,tfA,cfa,hfp(n),hcha2(n),hfp(n)+hcha2(n),strloss(n),qo2(n),qi2(n),qr(n),qstream(n),area(566)
!      if(n.eq.641)write(641,887)totaltime,tfA,cfa,hfp(n),hcha2(n),hfp(n)+hcha2(n),strloss(n),qo2(n),qi2(n),qr(n),qstream(n),area(641)
!      if(n.eq.869)write(869,887)totaltime,tfA,cfa,hfp(n),hcha2(n),hfp(n)+hcha2(n),strloss(n),qo2(n),qi2(n),qr(n),qstream(n),area(869)
!      if(n.eq.1835)write(1835,887)totaltime,tfA,cfa,hfp(n),hcha2(n),hfp(n)+hcha2(n),strloss(n),qo2(n),qi2(n),qr(n),qstream(n),area(869)

!     rev. 9.5.46  Dec.  23/08  - NK: trying to fix problem with -ve storage. Changed conditional to .lt.
!     rev. 9.6.05  Apr.  06/10  - NK: added store_error_flag for -ve storage grids
!     rev. 10.1.98 Oct   04/17  - NK: Deal with -ve flows in route
      if(store2(n).lt.0.0.or.qo2(n).lt.0.0)then
          store_error_flag(n)='true'
          qo2(n)=store1(n)/div+qi1(n)+qi2(n)-qo1(n)
          store2(n)=0.0
          dtmin=a6
!     rev. 9.5.47  Dec.  26/08  - NK: add flwinitflg to warn about initial flows
!         if we are early in the first event,possibly
!         the initial flows are too low.
!         This can happen if run is started in mid-winter
!         and the lowest downstream flow station is not at
!          outlet.
!         Might be fixed by increasing initial flow in the str file
  
          flwinitflg='y'
      endif     ! (store2(n).lt.0.0)
      
      if(qo2(n).gt.0.000001)then
!         CALCULATE THE VELOCITY THROUGH EACH SQUARE
!         CALCULATE TRAVEL TIME FOR MtfAIMUM VELOCITY.
          at=store2(n)/qo2(n)
!         SELECT MIN TRAVEL TIME FOR THE TIME STEP CALC
          dtmin=amin1(at,dtmin)
          dtmin=amax1(dtmin,a6)   ! dtmin > a6 no matter what
!     rev. 9.9.49  Jan.  06/14  - NK: Added courantflg
          if(at.lt.a6)then
              courantflg=.false.
          endif
     
!         DTMIN IS THE TIME REQUIRED TO COMPLETELY DRAIN
!         THE FASTEST EMPTYING ELEMENT

!         CALCULATE THE CHANNEL STATE FOR GRAPHICAL OUTPUT:
          i=yyy(n)
          j=xxx(n)
          atemp=qo2(n)/(0.4*bnkfll(n))+1.0

!     rev. 10.2.16 Feb.  14/18  - NK: Added bankfull flow calculation
          bankfull(i,j)=qo2(n)/bnkfll(n)*100.0
          bankfull(i,j)=amax1(0.1,bankfull(i,j))
          bankfull(i,j)=amin1(10000.,bankfull(i,j))

!         TO PREVENT INTEGER UNDERFLOW OR OVEFLOW:  
          atemp=amax1(atemp,1.0)
          atemp=amin1(atemp,99.0)

          istate(i,j)=int(atemp)
!         istate(i,j)=min(istate(i,j),99)
!         istate(i,j)=max(istate(i,j),1)
      endif
              
      firstpass=.false.
      
      if(iopt.ge.4)write(55,1002)i,j,n,istate(i,j),bnkfll(n),qo2(n)
 1002 format(' i,j,n,istate,bnkfull,qo2/',4i5,2f10.3)

      return
      
      end subroutine rt_channel
              