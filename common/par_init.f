      subroutine par_init()

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
!     rev. 9.4.07  May.  15/07  - NK: converted opt to gridded routing parameters
!     rev. 9.8.08  Nov.  18/11  - NK: added fratio for interception hight optimization

      use area_watflood
	implicit none

!     rev. 9.6.02  Mar.  15/10  - NK: add sublimation to optimization
!     rev. 9.6.06  Apr.  18/10  - NK: added glacier adjust for optimization

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      CHARACTER(128):: qstr
      CHARACTER(72) :: junk
!      CHARACTER(14) :: date
      CHARACTER(1)  :: smok 
	character(5) :: var_name(1000)
	integer       :: var_no(1000),dds_no(1000)
	INTEGER    :: iallcnt,icnt(5),ndir(5),nchr,ix,ios,icase,
     *              iallocate,igrdshft,iyshiftmin,iyshiftmax,
     *              jxshiftmin,jxshiftmax,ishift,jshift,inum,jnum,
     *              l,iw1,iw2,iv,iflg,i,n,ii,j,nhr,nhf
      REAL(4) ::   smc5(16),errold(5),err(5),chng(5),best(5)
      REAL(4)    :: optlow,e1,scale,ddtemp,cc1,cc2,crit,conv,best1
	real(4)    :: optlast
      integer*2 result1,ntest

      numa=0

!     surface model parameters (not gridded)

c      if(classcount.gt.0)then
        do iv=1,classcount
          if(akdlt(iv).gt.0.0)then
            numa=numa+1
            a(numa)     =ak(iv)
            ddelta(numa)=akdlt(iv)
            checkl(numa)=aklow(iv)
            checkh(numa)=akhgh(iv)
	      var_name(numa)='ak'
	      var_no(numa)=iv
	      dds_no(numa)=numa
          endif
        end do
        do iv=1,classcount
          if(akfsdlt(iv).gt.0.0)then
            numa=numa+1
            a(numa)     =akfs(iv)
            ddelta(numa)=akfsdlt(iv)
            checkl(numa)=akfslow(iv)
            checkh(numa)=akfshgh(iv)
	      var_name(numa)='akfs'
	      var_no(numa)=iv
	      dds_no(numa)=numa
          endif
        end do
        do iv=1,classcount
          if(recdlt(iv).gt.0.0)then
            numa=numa+1
            a     (numa)=rec(iv)
            ddelta(numa)=recdlt(iv)
            checkl(numa)=reclow(iv)
            checkh(numa)=rechgh(iv)
	      var_name(numa)='rec'
	      var_no(numa)=iv
	      dds_no(numa)=numa
          endif
        end do
        do iv=1,classcount
          if(r3dlt(iv).gt.0.0)then
            numa=numa+1
            a     (numa)=r3(iv)
            ddelta(numa)=r3dlt(iv)
            checkl(numa)=r3low(iv)
            checkh(numa)=r3hgh(iv)
	      var_name(numa)='r3'
	      var_no(numa)=iv
	      dds_no(numa)=numa
          endif
        end do
        do iv=1,classcount
          if(fpetdlt(iv).gt.0.0)then
            numa=numa+1
            a     (numa)=fpet(iv)
            ddelta(numa)=fpetdlt(iv)
            checkl(numa)=fpetlow(iv)
            checkh(numa)=fpethgh(iv)
	      var_name(numa)='fpet'
	      var_no(numa)=iv
	      dds_no(numa)=numa
          endif
        end do
        do iv=1,classcount
          if(ftalldlt(iv).gt.0.0)then
            numa=numa+1
            a     (numa)=ftall(iv)
            ddelta(numa)=ftalldlt(iv)
            checkl(numa)=ftalllow(iv)
            checkh(numa)=ftallhgh(iv)
	      var_name(numa)='ftall'
	      var_no(numa)=iv
	      dds_no(numa)=numa
          endif
        end do
!     rev. 9.8.08  Nov.  18/11  - NK: added fratio for interception hight optimization        
        do iv=1,classcount
          if(fratiodlt(iv).gt.0.0)then
            numa=numa+1
            a     (numa)=fratio(iv)
            ddelta(numa)=fratiodlt(iv)
            checkl(numa)=fratiolow(iv)
            checkh(numa)=fratiohgh(iv)
	      var_name(numa)='fratio'
	      var_no(numa)=iv
	      dds_no(numa)=numa
          endif
        end do
        do iv=1,classcount
          if(fmdlt(iv).gt.0.0)then
            numa=numa+1
            a     (numa)=fm(iv)
            ddelta(numa)=fmdlt(iv)
            checkl(numa)=fmlow(iv)
            checkh(numa)=fmhgh(iv)
	      var_name(numa)='mf'
	      var_no(numa)=iv
	      dds_no(numa)=numa
          endif
        end do
        do iv=1,classcount
          if(basdlt(iv).gt.0.0)then
            numa=numa+1
            a     (numa)=base(iv)+273.
            ddelta(numa)=basdlt(iv)
            checkl(numa)=baslow(iv)+273.
            checkh(numa)=bashgh(iv)+273.
	      var_name(numa)='base'
	      var_no(numa)=iv
	      dds_no(numa)=numa
          endif
        end do


!     rev. 9.6.02  Mar.  15/10  - NK: add sublimation to optimization
        do iv=1,classcount
          if(subdlt(iv).gt.0.0)then
            numa=numa+1
            if(sublimflg.eq.1)then          ! for sublim_factor only
              a(numa)=sublim_factor(iv)
            else
              a(numa)=sublim_rate(iv)
            endif  
            ddelta(numa)=subdlt(iv)
            checkl(numa)=sublow(iv)
            checkh(numa)=subhgh(iv)
	      var_name(numa)='sublim'
	      var_no(numa)=iv
	      dds_no(numa)=numa
          endif
        end do
        do iv=1,classcount
          if(retndlt(iv).gt.0.0)then
            numa=numa+1
            a     (numa)=retn(iv)
            ddelta(numa)=retndlt(iv)
            checkl(numa)=retnlow(iv)
            checkh(numa)=retnhgh(iv)
	      var_name(numa)='retn'
	      var_no(numa)=iv
	      dds_no(numa)=numa
          endif
        end do
        do iv=1,classcount
          if(ak2dlt(iv).gt.0.0)then
            numa=numa+1
            a     (numa)=ak2(iv)
            ddelta(numa)=ak2dlt(iv)
            checkl(numa)=ak2low(iv)
            checkh(numa)=ak2hgh(iv)
	      var_name(numa)='ak2'
	      var_no(numa)=iv
	      dds_no(numa)=numa
          endif
        end do
        do iv=1,classcount
          if(ak2fsdlt(iv).gt.0.0)then
            numa=numa+1
            a     (numa)=ak2fs(iv)
            ddelta(numa)=ak2fsdlt(iv)
            checkl(numa)=ak2fslow(iv)
            checkh(numa)=ak2fshgh(iv)
	      var_name(numa)='ak2fs'
	      var_no(numa)=iv
	      dds_no(numa)=numa
          endif
        end do

c        if(classcount.gt.classcount)then
c          write(6,*)'too many classes being optimized'
c          write(6,*)'classcount max. = classcount'
c          STOP 'reduce classcount in the spl runtime options menu'
c        endif

c        if(nbsn.gt.5)then
c          write(6,*)'too many basin roughness classes being optimized'
c          write(6,*)'nbsn max. = 5'
c          STOP 'reduce nbsn in the spl runtime options menu'
c        endif

c      endif

!     routing parameters (gridded)
!     rev. 9.4.07  May.  15/07  - NK: converted opt to gridded routing parameters

c      if(nbsn.gt.0)then
        do iv=1,nbsn
          if(flzdlt(iv).gt.0.0)then  
            numa=numa+1
            a     (numa)=flz_o(iv)
            ddelta(numa)=flzdlt(iv)
            checkl(numa)=flzlow(iv)
            checkh(numa)=flzhgh(iv)
	      var_name(numa)='flz'
	      var_no(numa)=iv
	      dds_no(numa)=numa
          endif
        end do  
        do iv=1,nbsn
          if(pwrdlt(iv).gt.0.0)then
            numa=numa+1
            a     (numa)=pwr_o(iv)
            ddelta(numa)=pwrdlt(iv)
            checkl(numa)=pwrlow(iv)
            checkh(numa)=pwrhgh(iv)
	      var_name(numa)='pwr'
	      var_no(numa)=iv
	      dds_no(numa)=numa
          endif
        end do
!     rev. 9.2.11  Sep.  15/05  - NK: added Manning's n  r1n & r2n
          do iv=1,nbsn
            if(r2ndlt(iv).gt.0.0)then
              numa=numa+1
              a     (numa)=r2n_o(iv)
              ddelta(numa)=r2ndlt(iv)
              checkl(numa)=r2nlow(iv)
              checkh(numa)=r2nhgh(iv)
	      var_name(numa)='r2n'
	      var_no(numa)=iv
	      dds_no(numa)=numa
            endif
          end do
        do iv=1,nbsn
          if(thetadlt(iv).gt.0.0)then
!           theta (wetland)
            numa=numa+1
            a     (numa)=theta_o(iv)
            ddelta(numa)=thetadlt(iv)
            checkl(numa)=thetalow(iv)
            checkh(numa)=thetahgh(iv)
	      var_name(numa)='theta'
	      var_no(numa)=iv
	      dds_no(numa)=numa
          endif
        end do
!       For what it's worth, widep should not an optimizable parameter..
!       It should be measured in the field!
        do iv=1,nbsn
          if(kconddlt(iv).gt.0.0)then
!           kcond (wetland)
            numa=numa+1
            a     (numa)=kcond_o(iv)
            ddelta(numa)=kconddlt(iv)
            checkl(numa)=kcondlow(iv)
            checkh(numa)=kcondhgh(iv)
	      var_name(numa)='kcond'
	      var_no(numa)=iv
	      dds_no(numa)=numa
          endif
        end do
        do iv=1,nbsn
          if(rlakedlt(iv).gt.0.0)then
!           rlake - increase n with water area (lakes)
            numa=numa+1
            a     (numa)=rlake_o(iv)
            ddelta(numa)=rlakedlt(iv)
            checkl(numa)=rlakelow(iv)
            checkh(numa)=rlakehgh(iv)
	      var_name(numa)='rlake'
	      dds_no(numa)=numa
	      var_no(numa)=iv
          endif
        end do
c      endif

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
        a     (numa)=a5
        ddelta(numa)=a5dlt
        checkl(numa)=a5low
        checkh(numa)=a5hgh
	  var_name(numa)='a5'
	  dds_no(numa)=numa
	  var_no(numa)=-1
      endif

!     rev. 9.6.06  Apr.  18/10  - NK: added glacier adjust for optimization
      if(fmadjustdlt.gt.0.0)then
        numa=numa+1
        a     (numa)=fmadjust
        ddelta(numa)=fmadjustdlt
        checkl(numa)=fmadjustlow
        checkh(numa)=fmadjusthgh
	  var_name(numa)='fmadjust'
	  dds_no(numa)=numa
	  var_no(numa)=-1
      endif

      if(fmalowdlt.gt.0.0)then
        numa=numa+1
        a     (numa)=fmalow
        ddelta(numa)=fmalowdlt
        checkl(numa)=fmalowlow
        checkh(numa)=fmalowhgh
	  var_name(numa)='fmalow'
	  dds_no(numa)=numa
	  var_no(numa)=-1
      endif

      if(fmahighdlt.gt.0.0)then
        numa=numa+1
        a     (numa)=fmahigh
        ddelta(numa)=fmahighdlt
        checkl(numa)=fmahighlow
        checkh(numa)=fmahighhgh
	  var_name(numa)='fmahigh'
	  dds_no(numa)=numa
	  var_no(numa)=-1
      endif

      if(gladjustdlt.gt.0.0)then
        numa=numa+1
        a     (numa)=gladjust
        ddelta(numa)=gladjustdlt
        checkl(numa)=gladjustlow
        checkh(numa)=gladjusthgh
	  var_name(numa)='gladjust'
	  dds_no(numa)=numa
	  var_no(numa)=-1
      endif


!     rev. 9.8.01  Jul.  21/11  - NK: added ragmet optimization to dds setup

      if(rlapsedlt.gt.0.0)then
        numa=numa+1
        a     (numa)=rlapse
        ddelta(numa)=rlapsedlt
        checkl(numa)=rlapselow
        checkh(numa)=rlapsehgh
	  var_name(numa)='rlapse'
	  dds_no(numa)=numa
	  var_no(numa)=-1
      endif

      if(tlapsedlt.gt.0.0)then
        numa=numa+1
        a     (numa)=tlapse
        ddelta(numa)=tlapsedlt
        checkl(numa)=tlapselow
        checkh(numa)=tlapsehgh
	  var_name(numa)='tlapse'
	  dds_no(numa)=numa
	  var_no(numa)=-1
      endif
      
!     rev. 10.4.31 Dec.  22/20  = NK fix climate end of year
      if(rainsnowtempdlt.gt.0.0)then
        numa=numa+1
        a     (numa)=rainsnowtemp
        ddelta(numa)=rainsnowtempdlt
        checkl(numa)=rainsnowtemplow
        checkh(numa)=rainsnowtemphgh
	  var_name(numa)='rainsnowtemp'
	  dds_no(numa)=numa
	  var_no(numa)=-1
      endif
      
      if(radinfldlt.gt.0.0)then
        numa=numa+1
        a     (numa)=radinfl
        ddelta(numa)=radinfldlt
        checkl(numa)=radinfllow
        checkh(numa)=radinflhgh
	  var_name(numa)='radinfl'
	  dds_no(numa)=numa
	  var_no(numa)=-1
      endif

      if(smoothdistdlt.gt.0.0)then
        numa=numa+1
        a     (numa)=smoothdist
        ddelta(numa)=smoothdistdlt
        checkl(numa)=smoothdistlow
        checkh(numa)=smoothdisthgh
	  var_name(numa)='smoothdist'
	  dds_no(numa)=numa
	  var_no(numa)=-1
      endif

!     rev. 10.1.92 May   25/17  - NK: Changed to max 200 dds variables
      if(numa.gt.200)then
	  print*,'in par_init: too many optimization params'
	  print*,'You seem to have ',numa
	  print*,'max allowed = 200'
        STOP 
	endif

!     REV. 7.52 - CHECK CONSTRAINTS
      iflg=0

      iflg=0
!      write(52,6503)
      if(dds_flag.ne.1)then
        do i=1,numa
!         for the pattern search:
          if(nper.le.0)then
            ddtemp=abs(ddelta(i))
          else
            ddtemp=abs(ddelta(i)*a(i))
          endif

          if(checkh(i)-checkl(i).le.3.0*ddtemp)then
!           THE CONSTRAINTS ARE TO CLOSE TOGETHER
            write(6,6510)i
            write(51,6510)i
            iflg=1
          endif

          cc1=a(i)-1.01*ddtemp
          if(cc1.le.checkl(i))then
            print*,' ddtemp, cc1, checkl, checkh, a'
            write(6,6501)i
            write(6,6500)ddtemp,cc1,checkl(i),checkh(i),a(i)
            write(51,6501)i
            write(51,6500)ddtemp,cc1,checkl(i),checkh(i),a(i)
            iflg=-1
          endif

          cc2=a(i)+1.01*ddtemp
          if(cc2.ge.checkh(i))then
            print*,' ddtemp, cc1, checkl, checkh, a'
            write(6,6502)i
            write(6,6500)ddtemp,cc2,checkl(i),checkh(i),a(i)
            write(51,6502)i
            write(51,6500)ddtemp,cc2,checkl(i),checkh(i),a(i)
            iflg=-1
          endif
        end do
	endif
c	else            ! added Feb 20/11 nk
!       for dds:

	  print*,'!!!!!!!!!!dds_init_value=',dds_init_value
	  open(unit=99,file='dds_init_value.txt',status='unknown',iostat=ios)
        if(ios.ne.0)then    ! added Nov. 10/14  nk
          print*
          print*,'Unable to open file dds_init_value.txt'
          print*,'Possible cause(s):'
          print*,'file in use by another application'
          stop 'Program aborted in par_init @ 472'
        endif
	  if(dds_init_value.lt.-990.0)then    ! check only with first trial
	    print*,'No of variables flaggged for DDS =',numa
          print*,'check lower constraints:'
	    print*,'    dds# variable   class #   low_limit     init_value'
	    write(99,*)'No of variables flaggged for DDS =',numa
          write(99,*)'check lower constraints:'
	    write(99,*)'    dds# variable   class #   low_limit     init_value'
          do i=1,numa
            if(a(i).le.checkl(i))then
              print*,dds_no(i),var_name(i),var_no(i),checkl(i),a(i)
              write(99,*)dds_no(i),var_name(i),var_no(i),checkl(i),a(i)
              iflg=-1
            endif
          end do

          print*,'check upper constraints:'
	    print*,'    dds# variable   class #  high_limit     init_value'
          write(99,*)'check upper constraints:'
	    write(99,*)'    dds# variable   class #  high_limit     init_value'
          do i=1,numa   
            if(a(i).ge.checkh(i))then
c             if(var_name(i).eq.'base')then
c                 print*,dds_no(i),var_name(i),var_no(i),
c    *             checkh(i)-273.,a(i)-273.
c                 write(99,*)dds_no(i),var_name(i),var_no(i),
c    *             checkh(i)-273.,a(i)-273.
c             elseif(var_name(i).eq.'subli')then
c                 print*,dds_no(i),var_name(i),var_no(i),
c    *             checkh(i)*24.,a(i)*24.
c                 write(99,*)dds_no(i),var_name(i),var_no(i),
c    *             checkh(i)*24.,a(i)*24.
c             else
                  print*,dds_no(i),var_name(i),var_no(i),
     *             checkh(i),a(i)
                  write(99,*)dds_no(i),var_name(i),var_no(i),
     *             checkh(i),a(i)
c              endif
              iflg=-1
            endif
          end do
	  endif
	  print*
	  write(99,*)

!       write the parameter table:
	  print*,
     *' dds# variable class# low_limit init_value high_limit'
	  write(99,*)
     *' dds# variable class# low_limit init_value high_limit'
	  do i=1,numa
c          if(var_name(i).eq.'base')then
c              write(*,9903)dds_no(i),var_name(i),var_no(i),
c     *               	checkl(i)-273.,a(i)-273.,checkh(i)-273.
c              write(99,9903)dds_no(i),var_name(i),var_no(i),
c     *               	checkl(i)-273.,a(i)-173.,checkh(i)-273.
c          elseif(var_name(i).eq.'subli')then
c              write(*,9903)dds_no(i),var_name(i),var_no(i),
c     *               	checkl(i)*24.,a(i)*24.,checkh(i)*24.
c              write(99,9903)dds_no(i),var_name(i),var_no(i),
c     *               	checkl(i)*24.,a(i)*24.,checkh(i)*24.
c          else
              write(*,9903)dds_no(i),var_name(i),var_no(i),
     *               	checkl(i),a(i),checkh(i)
              write(99,9903)dds_no(i),var_name(i),var_no(i),
     *               	checkl(i),a(i),checkh(i)
c          endif
9903      format(i5,a10,i5,3g12.3)
        end do
c	endif  ! dds_flag
        close(unit=99,status='keep')

c      write(52,6503)
      if(iflg.ne.0.and.dds_flag.eq.1)then   ! added Feb 20/11 nk
	  if(dds_init_value.lt.-990.0)then    ! check only with first trial
          close(unit=52,status='keep')
          print*
	    print*,'Please fix the variables & constraints'
          print*,'listed above'
	    print*,'DDS will not execute properly with these entries' 
	    print*,'Hit return to continue & then abort the batch file'
	    print*,'with Ctrl C'
	    print*
          print*
	    print*,'Please fix the variables & constraints'
          print*,'listed above'
	    print*,'DDS will not execute properly with these entries' 
	    print*,'Hit return to continue & then abort the batch file'
	    print*,'with Ctrl C'
	    print*
          pause  'Program paused in par_init @ 410'
!         return
          stop
	  endif
      endif

c      write(52,6503)
      if(iflg.ne.0.and.dds_flag.eq.0)then
        close(unit=52,status='keep')
        print*
	  print*,'Format for the opt part of the par file has changed to'
        print*,'free format'
	  print*,'Please ensure there are blanks between each entry in the'
	  print*,'columns in part two'
	  print*
        STOP 'please change init values or constraints or format'
      endif

      write(6,9345)classcount,nbsn,numa

!	write header for the summary.txt file
!     added Dec. 1/11  nk
      open(unit=9, file='summary_header.txt',status='unknown',
     *                                 iostat=ios)
      if(ios.ne.0)then    ! added Nov. 10/14  nk
        print*
        print*,'Unable to open file  summary_header.txt'
        print*,'Possible cause(s):'
        print*,'file in use by another application'
        stop 'Program aborted in par_init.f @ 566'
      endif
      write(9,90013)'      variable',(var_name(i),i=1,numa)
      write(9,90014)'        obj_fn',(var_no(i),i=1,numa)
90013 format(a14,<numa>a14)
90014 format(a14,<numa>i14)
      close(unit=9,status='keep')

 1001 format
     *(' you cannot run smc optimization for more than 1 storm')
 1002 format('        you have',i3,' events')
 6000 format(' optimized smcs:'/5f10.3/)
 6500 format(6e10.3)
 6501 format(' a(',i2,')too close to lower constraint')
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

      return

      END SUBROUTINE par_init

