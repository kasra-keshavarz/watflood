      subroutine par_assign()

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
      use area_watflood
	implicit none

!     rev. 9.6.02  Mar.  15/10  - NK: add sublimation to optimization
!     rev. 9.6.06  Apr.  18/10  - NK: added glacier adjust for optimization
!     rev. 9.7.29  Jul.  07/11  - NK: Add sublim_rate to set sublimation rate/day to pat file
!     rev. 9.8.01  Jul.  21/11  - NK: added ragmet optimization to dds setup
!     rev. 9.8.08  Nov.  18/11  - NK: added fratio for interception hight optimization

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      CHARACTER(128):: qstr
      CHARACTER(72) :: junk
!      CHARACTER(14) :: date
      CHARACTER(1)  :: smok 
	INTEGER    :: iallcnt,icnt(5),ndir(5),nchr,ix,ios,icase,
     *              iallocate,igrdshft,iyshiftmin,iyshiftmax,
     *              jxshiftmin,jxshiftmax,ishift,jshift,inum,jnum,
     *              l,iw1,iw2,iv,iflg,i,n,ii,j,nhr,nhf
      REAL(4) ::   smc5(16),errold(5),err(5),chng(5),best(5)
      REAL(4)    :: optlow,e1,scale,ddtemp,cc1,cc2,crit,conv,best1
	real(4)    :: optlast
      integer*2 result1,ntest

      numa=0
      if(classcount.gt.0)then
        do iv=1,classcount
          if(akdlt(iv).gt.0.0)then
            numa=numa+1
            ak(iv)=a(numa)
          endif
        end do
        do iv=1,classcount
          if(akfsdlt(iv).gt.0.0)then
            numa=numa+1
            akfs(iv)=a(numa)
          endif
        end do
        do iv=1,classcount
          if(recdlt(iv).gt.0.0)then
            numa=numa+1
            rec(iv)=a(numa)
          endif
        end do
        do iv=1,classcount
          if(r3dlt(iv).gt.0.0)then
            numa=numa+1
            r3(iv)=a(numa)
          endif
        end do
        do iv=1,classcount
          if(fpetdlt(iv).gt.0.0)then
            numa=numa+1
            fpet(iv)=a(numa)
          endif
        end do
        do iv=1,classcount
          if(ftalldlt(iv).gt.0.0)then
            numa=numa+1
            ftall(iv)=a(numa)
          endif
        end do
!     rev. 9.8.08  Nov.  18/11  - NK: added fratio for interception hight optimization        do iv=1,classcount
        do iv=1,classcount
          if(fratiodlt(iv).gt.0.0)then
            numa=numa+1
            fratio(iv)=a(numa)
          endif
        end do
      endif

      if(classcount.gt.0)then
        do iv=1,classcount
          if(fmdlt(iv).gt.0.0)then
            numa=numa+1
            fm(iv)=a(numa)
            endif
        end do
        do iv=1,classcount
          if(basdlt(iv).gt.0.0)then
            numa=numa+1
            base(iv)=a(numa)-273.
          endif
        end do
!     rev. 9.6.02  Mar.  15/10  - NK: add sublimation to optimization
        if(sublimflg.eq.1)then
          do iv=1,classcount
            if(subdlt(iv).gt.0.0)then
              numa=numa+1
              sublim_factor(iv)=a(numa)
            endif
          end do
        else  
!     rev. 9.7.29  Jul.  07/11  - NK: Add sublim_rate to set sublimation rate/day to par file
          do iv=1,classcount
            if(subdlt(iv).gt.0.0)then
              numa=numa+1
              sublim_rate(iv)=a(numa)
            endif
          end do
        endif  
      endif

      if(classcount.gt.0)then
        do iv=1,classcount
          if(retndlt(iv).gt.0.0)then
            numa=numa+1
            retn(iv)=a(numa)
          endif
        end do
        do iv=1,classcount
          if(ak2dlt(iv).gt.0.0)then
            numa=numa+1
            ak2(iv)=a(numa)
          endif
        end do
        do iv=1,classcount
          if(ak2fsdlt(iv).gt.0.0)then
            numa=numa+1
            ak2fs(iv)=a(numa)
          endif
        end do
      endif

      if(nbsn.gt.0)then
        do iv=1,nbsn
          if(flzdlt(iv).gt.0.0)then
            numa=numa+1          
            flz_o(iv)=a(numa)
          endif
        end do
        do iv=1,nbsn
          if(pwrdlt(iv).gt.0.0)then
            numa=numa+1
            pwr_o(iv)=a(numa)
          endif
        end do
!     rev. 9.2.11  Sep.  15/05  - NK: added Manning's n  r1n & r2n
          do iv=1,nbsn
            if(r2ndlt(iv).gt.0.0)then
              numa=numa+1
              r2n_o(iv)=a(numa)
            endif
          end do

        do iv=1,nbsn
          if(thetadlt(iv).gt.0.0)then
            numa=numa+1
            theta_o(iv)=a(numa)
          endif
        end do
        do n=1,naa
          wcap(n)=rl(n)*wetwid(n)*chadep(n)*theta_o(ibn(n))
          wetxa(n)=wcap(n)/rl(n)/theta_o(ibn(n))
        end do
        do iv=1,nbsn      
          if(kconddlt(iv).gt.0.0)then
            numa=numa+1
            kcond_o(iv)=a(numa)
          endif
        end do
        do iv=1,nbsn      
          if(rlakedlt(iv).gt.0.0)then
            numa=numa+1
            rlake_o(iv)=a(numa)
          endif
        end do
      endif

!     rev. 9.4.07  May.  15/07  - NK: converted opt to gridded routing parameters
	do n=1,naa
	  flz(n)=flz_o(ibn(n))
	  pwr(n)=pwr_o(ibn(n))
	  r2n(n)=r2n_o(ibn(n))
	  theta(n)=theta_o(ibn(n))
	  kcond(n)=kcond_o(ibn(n))
	  rlake(n)=rlake_o(ibn(n))
	end do

      if(a5dlt.gt.0.0)then
        numa=numa+1
        a5    =a(numa)
      endif

!     rev. 9.6.06  Apr.  18/10  - NK: added glacier adjust for optimization
      if(fmadjustdlt.gt.0.0)then
        numa=numa+1
        fmadjust=a(numa)
      endif

      if(fmalowdlt.gt.0.0)then
        numa=numa+1
        fmalow=a(numa)
      endif

      if(fmahighdlt.gt.0.0)then
        numa=numa+1
        fmahigh=a(numa)
      endif

      if(gladjustdlt.gt.0.0)then
        numa=numa+1
        gladjust=a(numa)
      endif


!     rev. 9.8.01  Jul.  21/11  - NK: added ragmet optimization to dds setup

      if(rlapsedlt.gt.0.0)then
        numa=numa+1
        rlapse=a(numa)
      endif

      if(tlapsedlt.gt.0.0)then
        numa=numa+1
        tlapse=a(numa)
      endif

      if(rainsnowtempdlt.gt.0.0)then
        numa=numa+1
        rainsnowtemp=a(numa)
      endif

      if(radinfldlt.gt.0.0)then
        numa=numa+1
        radinfl=a(numa)
      endif

      if(smoothdistdlt.gt.0.0)then
        numa=numa+1
        smoothdist=a(numa)
      endif


! FORMATS

 1001 format
     *(' you cannot run smc optimization for more than 1 storm')
 1002 format('        you have',i3,' events')
 6000 format(' optimized smcs:'/5f10.3/)
 6500 format(6e10.3)
 6502 format(' a(',i2,')too close to upper constraint')
c 6503 format(' ddtemp    cc2       checkl    checkh    a(i)')
 6510 format(' the constraints on a(',i3,') are too close together')
 8001 format(' err=',i5,5e10.3/)
 8003 format('    smc=',5f10.5)
 8004 format('   chng=',5f10.5)
 8123 format(' ','a(1),scale,best1,optim,optlow/',5f10.5)
 9345 format(' classcount,nbsn,numa/',3i5/)
 9901 format(10i5)
 9902 format(/' grid shift index =',i2,'  0=no shift  1=shifting')

      RETURN

      END SUBROUTINE par_assign
