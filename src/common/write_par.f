      SUBROUTINE write_par(ir,iverflg)

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
!***********************************************************************
!  rev. 9.2.09  Sep.  11/05  - NK: removed write_par.for from rdpar.for
!  rev. 9.5.43  Oct.  27/08  - NK: changed bottom part of par file to be free format
!  rev. 9.5.45  Dec.  16/08  - NK: added various error calculations - user's choice with errflg
!  rev. 9.6.01  Mar.  01/10  - NK: DDS capability added
!  rev. 9.6.06  Apr.  18/10  - NK: added glacier adjust for optimization
!  rev. 9.7.14  Nov.  22/10  - NK: Allow 30 land cover classes
!
!  s/r created Sept. 11/05 to separate writing from reading 
!  in the rdpar.for subroutine. 
!
!  version numbers are added for version 7.0 and later. this allows
!  parameter files to backward compatible when more parameters are 
!  added to the file in the future. in version 7.0 some headers are 
!  added as well for readability.
!
!
!  iprtflg - if eq 2 write a new parameter file after optimization
!  r1      - the roughness of the flood plain
!  r2      - the roughness of the largest part of the channel
!  r1      - factor for raising r2 ie if r1=2 then f.p. roughness
!            is 2 times channel roughness
!  zz      - is an exponent in calculating the mannings  n
!  h       - crop height
!  ix      - exponent in defining porosity
!  e1       - void ratio
!  ak      - permeability
!  sm      -soil moisture(average for month)
!  nbsn    - no of basins with different r2 to be optimized
!          - must be smaller than nrvr
!  classcount  - no of classes to be optimized max=4
!  errflg  - picks the type of error calculation 
!  a1 .. a12 parameters variously used mostly in runoff
!  a7      - weighting factor for old vs. new sca value default=.90
!  a8      - time offset to check for min temp on rdtemp
!  a9      - heat deficit to swe ratio  default=0.333
!  a10     - power on interflow discharge function default = 1.0
!  a11     - equivalent vegetation height for bare ground evaporation
!  a12     - min precip rate for smearing used in rain.for default=0.0
!
!***********************************************************************

      use area_watflood
	implicit none

      DIMENSION :: ajunk(17)

      CHARACTER(10) :: time
      CHARACTER(8)  :: cday
      CHARACTER(128):: qstr
      CHARACTER(60) :: junk
      CHARACTER(30) :: filename1
      CHARACTER(1)  :: errorflg,answer,linetype,tempflg
      INTEGER(kind=2) :: result1
	INTEGER       :: nrvr1,nchr,iprtflg,linecount,i,
     *                 ios,iverflg,ix,iallocate,ijunko,ii,j,ir,numb
      real*4        :: hmax,e1,ajunk,kcondflg

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE


!     ADDED THIS FROM TODD'S ETPAR.FOR - FRANK S: NOV/97 
      alamb=2.478
      den=999.
      alpha=1.35

! READ SNOW COVER PARAMETERS:   
!  - CALLED ON FIRST TIME STEP FOR SIMULATION
!  - CALLED TO WRITE TO newsnow.par DURING OPTIMIZATION 
!    IF OPTIMIZATION PARAMETERS ARE UPDATED
!               -- first the snow output flags
!               -- melt factors
!               -- now mbase
!               -- negative melt factors
!               -- wind function
!               -- ati decay
!               -- conversion snow density depth to we
!               -- liquid water holding capacity
!               -- daily ground melt
!---------------------------------------------------------

      itype=int(type1)

!     rev. 9.5.06  Feb.  05/08  - NK: added pool and pool_o in rdpar & route

      if(ir.ne.51)then

!       if iverflg=1 then it is called from spl

!       if iverflg=0 it is called from connector_dds_watflood
!          and the file is already opened in the calling program

        ver=9.50 
	  
        if(program_name(1:3).eq.'bsn')then
          open(ir,file='new.par',status='unknown',iostat=ios)
          if(ios.ne.0)then
            print*,' Error opening basin\new.par file'
	      print*,'on unit=',ir
            print*
            stop ' Program aborted in write_par @ ~104'
          endif
	  else
          if(iverflg.eq.1)then
            open(ir,file='basin\new.par',status='unknown',iostat=ios)
            if(ios.ne.0)then
              print*,' Error opening basin\new.par file'
	        print*,'on unit=',ir
              print*
              stop ' Program aborted in write_par @ ~111'
            endif
	    endif
	  endif

	else
!       here for either spl or bsn 
	  write(ir,*)
	  write(ir,*)'echo the parameter file fln=',fln(2)
	  write(ir,*)
      endif

      numb=0  ! so new.par can be copied to bsnm.par for batch job

!       copied from readpar.f
        title(71)='ver'
        title(72)='iopt'
        title(73)='itype'
        title(74)='numa'
        title(75)='nper'
        title(76)='kc'
        title(77)='maxn'
        title(78)='ddsfl'
        title(79)='itrce'
        title(80)='errfl'
        title(81)='typeo'
        title(82)='nbsn'
        title(83)='mndr'
        title(84)='aa2'
        title(85)='aa3'
        title(86)='aa4'
        title(87)='theta'
        title(88)='widep'
        title(89)='kcond'
	  title(90)='pool'
	  title(70)='Rlake'
        title(91)='a1'
        title(92)='a2'
        title(93)='a3'
        title(94)='a4'
        title(95)='a5'
        title(96)='a6'
        title(97)='a7'
        title(98)='a8'
        title(99)='a9'
        title(67)='a10'
        title(68)='a11'
        title(69)='a12'
        notes(71)='parameter file version number'
        notes(72)='debug level'
        notes(73)='Not used'
        notes(74)='optimization 0=no 1=yes'
        notes(75)='opt delta 0=absolute 1=ratio'
        notes(76)='no of times delta halved'
        notes(77)='max no of trials'
        notes(78)='0=single run  1=DDS'
        notes(79)='Tracer choice'
        notes(80)='1=wMSE 2=SSE 3=wSSE 4=Dv 5=wDv 6=wSAE 7=Nash'
        notes(81)='no of land classes optimized(part 2)'
        notes(82)='no of river classes optimized (part 2)'
        notes(91)='ice cover weighting factor'
        notes(92)='Manning`s n correction for instream lakes'
        notes(93)='error penalty coefficient'
        notes(94)='error penalty threshold'
        notes(95)='API coefficient'
        notes(96)='Minimum routing time step in seconds'
        notes(97)='weighting factor - old vs. new sca value'
        notes(98)='min temperature time offset'
        notes(99)='max heat deficit /swe ratio'
        notes(67)='exponent on uz discharce function'
        notes(68)='bare ground equivalent veg height for evap`n'
        notes(69)='min precip rate for smearing'


      call date_and_time(cday,time)
      write(ir,6011) time(1:2),time(3:4),time(5:6)
      write(ir,6012) cday(1:4),cday(5:6),cday(7:8)
      write(ir,9807)title(71),ver,notes(71)
      write(ir,9805)title(72),iopt,notes(72)
      write(ir,9805)title(73),itype,notes(73)
      write(ir,9805)title(74),numb,notes(74)
      write(ir,9805)title(75),nper,notes(75)
      write(ir,9805)title(76),kc,notes(76)
      write(ir,9805)title(77),maxn,notes(77)
      write(ir,9805)title(78),dds_flag,notes(78)
      write(ir,9805)title(79),itrace,notes(79)
      write(ir,9805)title(80),errflg,notes(80)
      write(ir,9805)title(81),classcount,notes(81)
      write(ir,9805)title(82),nbsn,notes(82)
      write(ir,9807)title(91),a1,notes(91) 
      write(ir,9807)title(92),a2,notes(92)  
      write(ir,9807)title(93),a3,notes(93)
      write(ir,9807)title(94),a4,notes(94)
      write(ir,9807)title(95),a5,notes(95)
      write(ir,9807)title(96),a6,notes(96)
      write(ir,9807)title(97),a7,notes(97)
      write(ir,9807)title(98),a8,notes(98)
      write(ir,9807)title(99),a9,notes(99)
      write(ir,9807)title(67),a10,notes(67)
      write(ir,9807)title(68),a11,notes(68)
      write(ir,9807)title(69),a12,notes(69)
 
!     rev. 9.4.07  May.  15/07  - NK: converted opt to gridded routing parameters
      write(ir,5009)(rivtype(i),i=1,nrvr)
      write(ir,1001)title(14),(flz_o(i),i=1,nrvr)   
      write(ir,1001)title(15),(pwr_o(i),i=1,nrvr)   
      write(ir,1001)title(16),(r1n_o(i),i=1,nrvr)    
      write(ir,1001)title(17),(r2n_o(i),i=1,nrvr) 
      write(ir,1001)title(83),(mndr_o(i),i=1,nrvr)
      write(ir,1001)title(84),(aa2_o(i),i=1,nrvr)
      write(ir,1001)title(85),(aa3_o(i),i=1,nrvr)
      write(ir,1001)title(86),(aa4_o(i),i=1,nrvr)
      write(ir,1001)title(87),(theta_o(i),i=1,nrvr)
      write(ir,1001)title(88),(widep_o(i),i=1,nrvr)
      write(ir,1001)title(89),(kcond_o(i),i=1,nrvr)
!     rev. 9.5.06  Feb.  05/08  - NK: added pool and pool_o in rdpar & route
      write(ir,1001)title(90),(pool_o(i),i=1,nrvr)
      write(ir,1001)title(70),(rlake_o(i),i=1,nrvr)

      write(ir,5009)(nclass(i),i=1,classcount)
      write(ir,1001)title(6),(ds(i),i=1,classcount)    
      write(ir,1001)title(7),(dsfs(i),i=1,classcount)    
      write(ir,1001)title(8),(rec(i),i=1,classcount)   
      write(ir,1001)title(9),(ak(i),i=1,classcount)    
      write(ir,1001)title(10),(akfs(i),i=1,classcount)    
      write(ir,1001)title(11),(retn(i),i=1,classcount)   
      write(ir,1001)title(12),(ak2(i),i=1,classcount)   
      write(ir,1001)title(13),(ak2fs(i),i=1,classcount)   
      write(ir,1001)title(18),(r3(i),i=1,classcount)    
      write(ir,1001)title(19),(r3fs(i),i=1,classcount)    
      write(ir,1001)title(20),(r4(i),i=1,classcount)  
        
      write(ir,1001)title(21),(chnl(i),i=1,5)  

      write(ir,1001)title(22),(fm(ii),ii=1,classcount)
      write(ir,1001)title(23),(base(ii),ii=1,classcount)
      write(ir,1001)title(24),(fmn(ii),ii=1,classcount)
      write(ir,1001)title(25),(uadj(ii),ii=1,classcount)
      write(ir,1001)title(26),(tipm(ii),ii=1,classcount)
      write(ir,1001)title(27),(rho(ii),ii=1,classcount)
      write(ir,1001)title(28),(whcl(ii),ii=1,classcount)
      write(ir,5004)'fmadj',fmadjust
      write(ir,5004)'fmlow',fmalow
      write(ir,5004)'fmhgh',fmahigh
      write(ir,5004)'gladj',gladjust
      write(ir,5004)'rlaps',rlapse
      write(ir,5004)'elvrf',elvref
	junk=' 1 = pan; 2 = Hargreaves; 3 = Priestley-Taylor'
      write(ir,9804)title(31),flgevp2,junk
      write(ir,9802)title(32),albe
      write(ir,9802)title(33),(alb(i),i=1,classcount)
      write(ir,9802)title(34),(fpet(i),i=1,classcount)
      write(ir,9802)title(35),(ftall(i),i=1,classcount)
      write(ir,9801)title(36),(flint(i),i=1,classcount)
      write(ir,9802)title(37),(fcap(i),i=1,classcount)
      write(ir,9802)title(38),(ffcap(i),i=1,classcount)
      write(ir,9802)title(39),(spore(i),i=1,classcount)
      write(ir,9802)title(40),(sublim_factor(i),i=1,classcount)
      write(ir,9801)title(41),tempa2
      write(ir,9801)title(42),tempa3
      write(ir,9801)title(43),tton
      write(ir,9801)title(44),lat
      write(ir,9803)title(45),(diff(i),i=1,12)
      write(ir,9803)title(46),(hu(i),i=1,12)
      write(ir,9803)title(47),(pres(i),i=1,12)

      write(ir,5005) heading(2)

      if(iopt.eq.2)print*, 'param 14'

      do ii=1,classcount
         write(ir,1003)title(50+ii),(h(i,ii),i=1,12)
      end do

      write(ir,5005)heading(3)

      do i=1,classcount
        write(ir,1002) 
     *  title(9),akdlt(i),aklow(i),akhgh(i),ak(i)
	end do
      do i=1,classcount
        write(ir,1002) 
     *  title(10),akfsdlt(i),akfslow(i),akfshgh(i),akfs(i)
	end do
      do i=1,classcount
        write(ir,1002) 
     *  title(8),recdlt(i),reclow(i),rechgh(i),rec(i)
	end do
      do i=1,classcount
        write(ir,1002) 
     *  title(18),r3dlt(i),r3low(i),r3hgh(i),r3(i)
	end do
      do i=1,classcount
        write(ir,1002)
     *  title(34),fpetdlt(i),fpetlow(i),fpethgh(i),fpet(i)
	end do
      do i=1,classcount
        write(ir,1002)
     *  title(35),ftalldlt(i),ftalllow(i),ftallhgh(i),ftall(i)
	end do
      do i=1,classcount
        write(ir,1002)
     *  title(22),fmdlt(i),fmlow(i),fmhgh(i),fm(i)   
	end do
      do i=1,classcount
        write(ir,1002)
     *  title(23),basdlt(i),baslow(i),bashgh(i),base(i)
	end do
      do i=1,classcount
        write(ir,1002)
     *  title(40),subdlt(i),sublow(i),subhgh(i),sublim_factor(i)
	end do
      do i=1,classcount
        write(ir,5003)
     *  title(11),retndlt(i),retnlow(i),retnhgh(i),retn(i)
	end do
      do i=1,classcount
        write(ir,5003)
     *  title(12),ak2dlt(i),ak2low(i),ak2hgh(i),ak2(i)
	end do
      do i=1,classcount
        write(ir,5003)
     *  title(13),ak2fsdlt(i),ak2fslow(i),ak2fshgh(i),ak2fs(i)
!     rev. 9.4.07  May.  15/07  - NK: converted opt to gridded routing parameters
	end do
      do i=1,nbsn
        write(ir,5003)
     *  title(14),flzdlt(i),flzlow(i),flzhgh(i),flz_o(i)
	end do
      do i=1,nbsn
        write(ir,5003)
     *  title(15),pwrdlt(i),pwrlow(i),pwrhgh(i),pwr_o(i)
	end do
      do i=1,nbsn
        write(ir,1002)'r2n  ', 
     *  r2ndlt(i),r2nlow(i),r2nhgh(i),r2n_o(i)
	end do
      do i=1,nbsn
        write(ir,5003)'theta',
     *  thetadlt(i),thetalow(i),thetahgh(i),theta_o(i)
	end do
      do i=1,nbsn
        write(ir,5003)'kcond',
     *  kconddlt(i),kcondlow(i),kcondhgh(i),kcond_o(i)
	end do
      do i=1,nbsn
        write(ir,5003)'rlake',
     *  rlakedlt(i),rlakelow(i),rlakehgh(i),rlake_o(i)
	end do
      write(ir,5003)'a5   ',a5dlt,a5low,a5hgh,a5
!     rev. 9.6.06  Apr.  18/10  - NK: added glacier adjust for optimization
	write(ir,5003)'fmadj',fmadjustdlt,fmadjustlow,fmadjusthgh,fmadjust
	write(ir,5003)'fmlow',fmalowdlt,fmalowlow,fmalowhgh,fmalow
	write(ir,5003)'fmhgh',fmahighdlt,fmahighlow,fmahighhgh,fmahigh
	write(ir,5003)'gladj',gladjustdlt,gladjustlow,gladjusthgh,gladjust


      if(ir.ne.51)print*,' new.par file written'

      if(ir.eq.51)then     ! write a few blank lines for readability
	  write(ir,*)  
	  write(ir,*)
        return
	endif
      
      close(ir,status='keep')

d	print*,'Closed unit=',ir
d	print*,'file name =','basin\new.par'

      if(iopt.eq.2)print*, 'param before return @ 9989'

 9989 RETURN

! FORMATS

 1000 format(a5,g10.3,2i5,a65)
 1001 format(a5,99e10.3)
 1002 format(a5,4(e12.3))
 1003 format(a5,12f5.2)
 1004 format(a5,99e10.3)
 1005 format(a5,99e10.3)
 1006 format(a5,99e10.3)
 1007 format(a5,99e10.3)
 1008 format(a5,99e10.3)
 1009 format('     e       ix')
 1010 format(' maximun interception storage in mm')
 1011 format(' depression storage ds(i)')
 1012 format(' permeabilities ak(i)')
 1013 format(' roughness r3(i)  -  pervious area')
 1014 format(' velocity factor for channels/square   chnl(i)')
 1015 format(' roughness r4(i)  -  impervious area')
 1019 format(' a1 ... a12')
 1022 format(' interflow recession constant rec(i)')
 1023 format(' ','ti            delta        low         high   param') 
 1025 format(/' ','flgevp2 set to',f5.1,' in rdpar.for'/)

 1031 format(' lower zone discharge function flz(i)')
 1032 format(' lower zone power pwr(i)')
 1033 format(' flood plain roughness multiplier r1(i)')
 1034 format(' river roughness r2(i)')
 1035 format(' meander length multiplier mndr(i)')
 1036 format(' bankfull constant aa2(i)')
 1037 format(' bankfull constant aa3(i)')
 1038 format(' bankfull constant aa4(i)')
 1039 format(' wetland porosity theta(i)')
 1040 format(' wetland width/depth ratio widep(i)')
 1041 format(' wetland conductivity kcond(i)')

 3000 format(a5,6i5,25x,f10.0) 
 3010 format(a1,a79)
 3101 format(12f5.0)
 3203 format(' Warning: basin/evap.dat table incomplete'/
     *         '          zero values are inserted for evap.dat'/)
 5000 format(a5,2i5,a75)
 5001 format(a5,17f10.3)
 5002 format(17f10.3)
 5003 format(a5,4(e12.3))
 5004 format(a5,6f10.3)
 5005 format(a80)
 5006 format(' warning: default values for retn, ak2, flz & pwr'/
     *'          are used. par file is out of date'/)
 5007 format(' warning: errflg in d2 =',i5,' reset to 1 in rdpar')
 5009 format(5x,99a10)
 5010 format(a1)
 5100 format(' in rdpar - problem opening basin/evap.dat file'/
     *         '          zero values are inserted for evap.dat'/)
5110  format(//,' Value for a9 heat deficit to swe outside range',/
     *     f10.3,' assumed')
5111  format(//,' Value for a10 uz exponent outside range',/
     *     f10.3,' assumed')
 5999 format(' mflag=',i5,'   errflg=',i5)
 6003 format(a5,17e10.3)

 6011 format('# runtime    ',2(a2,':'),a2)
 6012 format('# rundate  ',a4,'-',a2,'-',a2)


 8429 format(//' in runof5,flz/',5f10.6/
     * ' wrong parameter value - change flz(i) to +ve value'/)
 9801 format(a5,99f10.0)
 9802 format(a5,99f10.2)
 9803 format(a5,12f5.1)
 9804 format(a5,f10.2,a60)
 9805 format(a5,i10,5x,a40)
 9806 format(a5,f10.0,5x,a40)
 9807 format(a5,f10.3,5x,a40)
 9808 format(//' Parameters used for this run:'/)
 9809 format('A new par file called basin\new.par has been created')

      RETURN

      END SUBROUTINE write_par

