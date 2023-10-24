      SUBROUTINE write_par_10(unitNum,flnNum)

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
!***********************************************************************
!  rev. 9.2.09  Sep.  11/05  - NK: removed write_par.for from rdpar.for
!  rev. 9.5.43  Oct.  27/08  - NK: changed bottom part of par file to be free format
!  rev. 9.5.45  Dec.  16/08  - NK: added various error calculations - user's choice with errflg
!  rev. 9.6.01  Mar.  01/10  - NK: DDS capability added
!  rev. 9.6.06  Apr.  18/10  - NK: added glacier adjust for optimization
!  rev. 9.8.01  Jul.  21/11  - NK: added ragmet optimization to dds setup
!  rev. 9.8.04  Sep.  02/11  - NK: Fix bug in write)par_10 when reading old par file
!  rev. 9.8.08  Nov.  18/11  - NK: added fratio for interception hight optimization!
!  rev. 9.8.08  Nov.  18/11  - NK: ver=10.1
!  rev. 9.8.17  Apr.  24/11  - NK: Moved dds flags to top of par file
!  rev. 9.8.40  Nov.  26/12  - NK: convert interception cap: h(,)*fratio()
!     REV. 10.1.13 Dec.  28/15  - NK: Rearranged the par file blocks & contents
!     rev. 10.4.60 Oct.  10/22  = NK Reordered lines in write_par_10 to accomodate OstrichMPI files 
!     rev. 10.5.15 June  11/23  = NK Added sdcd, r1n, aa2, aa3, widep, pool to Ostin.txt list
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
	INTEGER       :: nrvr1,nchr,iprtflg,linecount,i,unitNum,flnNum,
     *                 ios,iverflg,ix,iallocate,ijunko,ii,j,numb,nnts
      real*4        :: hmax,e1,ajunk,kcondflg
      logical       :: exists

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

!      if(unitNum.ne.51)then

!       unitNum51 is spl.txt and is already open

!       if iverflg=1 then it is called from spl

!       if iverflg=0 it is called from connector_dds_watflood
!          and the file is already opened in the calling program

        ver=10.1
	  
c        if(program_name(1:3).eq.'bsn')then

!       if filename(0) the data is written to spl.txt
        if(flnNum.ne.0)then
          open(unitNum,file=filename(flnNum),
     *        	status='unknown',iostat=ios)
          if(ios.ne.0)then
            print*,' Error opening',filename(flnNum)
	      print*,'on unit=',unitNum
            print*
            print*,'When this happens, output files ended up read only'
            print*,'Delete all files in the debug directory & '
            print*,'delete all most recent output files in the '
            print*,'results directory and try again'
            print*
            stop ' Program aborted in write_par_10 @ ~104'
	    endif
        endif

!	else
!       here for either spl or bsn 
!	  write(unitNum,*)
!	  write(unitNum,*)'echo the parameter file fln=',fln(2)
!	  write(unitNum,*)
!      endif

      numb=0  ! so new.par can be copied to bsnm.par for batch job

      write(unitNum,97001,iostat=ios)':FileType, ','WatfloodParameter'
     *                 ,ver,',# parameter file version number'
97001 format(a11,a19,f10.2,a31)     
      if(ios.ne.0)then
          print*,'ios = ',ios
          print*,'error writing ',filename(flnnum)
          print*,'Check if files in results are write protect'
          stop 'program aborted in write_par.f @ 160'
      endif

      call date_and_time(cday,time)
      write(unitNum,98002)cday(1:4),cday(5:6),cday(7:8),
     *                    time(1:2),time(3:4),time(5:6)
98002 format(':CreationDate ,',a4,'-',a2,'-',a2,'  ',2(a2,':'),a2)

!     REV. 10.1.29 May   04/16  - NK: Added parfile comments
      if(ParFileCommentNumber.gt.0)then
        do i=1,parFileCommentNumber
          write(unitNum,98007)ParFileComment(i)
98007     format(a72)          
        end do
      endif

      write(unitNum,98003)
98003 format(':GlobalParameters')


c      write(unitNum,98004)'# Global parameters'
c98004 format(':ver,            ',f7.2,',# parameter file version')

!     rev. 10.4.27 Sep.  28/20  = NK Write parfile.scv during dds run
c      write(unitNum,98005)iopt
      write(unitNum,98005)1
98005 format(':iopt,           ',i7,',# debug level')

      write(unitNum,98006)itype
98006 format(':itype,          ',i7,',# channel type - floodplain/no')

      write(unitNum,98013)itrace
98013 format(':itrace,         ',i7,',# Tracer choice')

c      write(unitNum,98014)classcount
c98014 format(':classcount,',i7,',# no of land classes optimized(part 2)')

c      write(unitNum,98016)nbsn
c98016 format(':nbsn,  ',i7,',# no of river classes optimized (part 2)')

      write(unitNum,98017)a1
98017 format(':a1,          ',f10.3,',# ice cover weighting factor')

      write(unitNum,98018)a2
98018 format(':a2,          ',f10.3,',# swe correction threshold')

      write(unitNum,98019)a3
98019 format(':a3,          ',f10.3,',# error penalty coefficient')

      write(unitNum,98020)a4
98020 format(':a4,          ',f10.3,',# error penalty threshold')

      write(unitNum,98021)a5
98021 format(':a5,          ',f10.3,',# API coefficien')

      write(unitNum,98022)a6
98022 format(':a6,          ',f10.3,
     *                      ',# Minimum routing time step in seconds')

      write(unitNum,98023)a7
98023 format(':a7,          ',f10.3,
     *                      ',# weighting - old vs. new sca value')
  
      write(unitNum,98024)a8
98024 format(':a8,          ',f10.3,',# min temperature time offset')

      write(unitNum,98025)a9
98025 format(':a9,          ',f10.3,',# max heat deficit /swe ratio')

      write(unitNum,98026)a10
98026 format(':a10,         ',f10.3,
     *                     ',# exponent on uz discharce function')

      write(unitNum,98027)a11
98027 format(':a11,         ',f10.3,
     *                     ',# bare ground equiv. veg height for ev')

      write(unitNum,98028)a12
98028 format(':a12,         ',f10.3,',# min precip rate for smearing')
 
!     rev. 10.4.30 Dec.  18/20  = NK Added a13 parameter for rain/snow tmp > base tmp
      write(unitNum,98050)a13
98050 format(':a13,         ',f10.3,',# ')
 
      write(unitNum,98029)fmadjust
98029 format(':fmadjust,    ',f10.3,',# snowmelt ripening rate')
 
      write(unitNum,98030)fmalow
98030 format(':fmalow,      ',f10.3,',# min melt factor multiplier')
 
      write(unitNum,98031)fmahigh
98031 format(':fmahigh,     ',f10.3,',# max melt factor multiplier')
 
      write(unitNum,98032)gladjust
98032 format(':gladjust,    ',f10.3,',# glacier melt factor multiplier')
 
      write(unitNum,98033)rlapse
98033 format(':rlapse,      ',f10.6,',# precip lapse rate mm/m')

      write(unitNum,98034)tlapse
98034 format(':tlapse,      ',f10.6,',# temperature lapse rate dC/m')
 
c      if(rainsnow
      write(unitNum,98036)rainSnowTemp
98036 format(':rainsnowtemp,',f10.3,',# rain/snow temperature')
 
      write(unitNum,98037)radinfl
98037 format(':radiusinflce,',f10.3,',# radius of influence km')
 
      write(unitNum,98038)smoothDist
98038 format(':smoothdist,  ',f10.3,',# smoothing distance km')
 
      write(unitNum,98035)elvref
98035 format(':elvref,      ',f10.3,',# reference elevation')
 
      write(unitNum,98039)flgevp2
98039 format(':flgevp2  ,   ',f10.3,
     *                ',# 1=pan;4=Hargreaves;3=Priestley-Taylor')
 
      write(unitNum,98040)albe
98040 format(':albe  ,      ',f10.3,',# albedo????')

      write(unitNum,98041)tempa2
98041 format(':tempa2,      ',f10.3,',# ')

      write(unitNum,98042)tempa3
98042 format(':tempa3,      ',f10.3,',# ')

      write(unitNum,98043)tton
98043 format(':tton  ,      ',f10.3,',# ')

      write(unitNum,98044)lat
98044 format(':lat   ,      ',f10.3,',# latitude')

      write(unitNum,98045)chnl(1)
98045 format(':chnl(1),     ',f10.3,',# manning`s n multiplier')
 
      write(unitNum,98046)chnl(2)
98046 format(':chnl(2),     ',f10.3,',# manning`s n multiplier')
 
      write(unitNum,98047)chnl(3)
98047 format(':chnl(3),     ',f10.3,',# manning`s n multiplier')
 
      write(unitNum,98048)chnl(4)
98048 format(':chnl(4),     ',f10.3,',# manning`s n multiplier')
 
      write(unitNum,98049)chnl(5)
98049 format(':chnl(5),     ',f10.3,',# manning`s n multiplier')
 
      write(unitNum,98099)
98099 format(':EndGlobalParameters')
      write(unitNum,98098)
98098 format('#')

!     rev. 9.8.17  Apr.  24/11  - NK: Moved dds flags to top of par file
      write(unitNum,98500)
98500 format(':OptimizationSwitches')
	write(unitNum,98507)numa
98507 format(':numa,  ',i7,',# PS optimization 1=yes 0=no')
      write(unitNum,98508)nper
98508 format(':nper,  ',i7,',# opt 1=delta 0=absolute')
      write(unitNum,98509)kc
98509 format(':kc,    ',i7,',# no of times delta halved')
      write(unitNum,98510)maxn
98510 format(':maxn,  ',i7,',# max no of trials')
!     rev. 10.4.27 Sep.  28/20  = NK Write parfile.csv during dds run
      if(flnNum.eq.25)write(unitNum,98511)1    ! for dds
      if(flnNum.eq.26)write(unitNum,98511)0    ! for results\parfile.csv & charm singlerun
98511 format(':ddsflg,',i7,',# 0=single run  1=DDS ')
      write(unitNum,98512)errflg
98512 format(':errflg,',i7,',# 1=wMSE 2=SSE 3=wSSE 4=VOL ')
      write(unitNum,98513)
98513 format(':EndOptimizationSwitches')
      write(unitNum,98098)


!     *************************************************************** 
      write(unitNum,98100)
98100 format(':RoutingParameters')
      write(unitNum,98198)nrvr
98198 format(':RiverClasses,',i12)
      write(unitNum,98101)(rivtype(i),i=1,nrvr)
98101 format(':RiverClassName,  ',<nrvr>(a10,'  ,'))
c     *       '# in channel lake retardation coefficient')	  
      write(unitNum,98102)(flz_o(i),i=1,nrvr) 
98102 format(':flz,             ',<nrvr>(g12.3,','),
     *       '# lower zone oefficient')	  
      write(unitNum,98103)(pwr_o(i),i=1,nrvr)   
98103 format(':pwr,             ',<nrvr>(g12.3,',')
     *       '# lower zone exponent')	  
      write(unitNum,98105)(r2n_o(i),i=1,nrvr) 
98105 format(':r2n,             ',<nrvr>(g12.3,',')
     *       '# channel Manning`s n')	 
      write(unitNum,98110)(theta_o(i),i=1,nrvr)
98110 format(':theta,           ',<nrvr>(g12.3,',')
     *       '# wetland or bank porosity')	  
      write(unitNum,98112)(kcond_o(i),i=1,nrvr)
98112 format(':kcond,           ',<nrvr>(g12.3,',')
     *       '# wetland/bank lateral conductivity')	  
      write(unitNum,98114)(rlake_o(i),i=1,nrvr)
98114 format(':rlake,           ',<nrvr>(g12.3,',')
     *       '# in channel lake retardation coefficient')	  
      write(unitNum,98104)(r1n_o(i),i=1,nrvr)    
98104 format(':r1n,             ',<nrvr>(g12.3,',')
     *       '# overbank Manning`s n')	  
      write(unitNum,98107)(aa2_o(i),i=1,nrvr)
98107 format(':aa2,             ',<nrvr>(g12.3,',')
     *       '# channel area intercept = min channel xsect area')	  
      write(unitNum,98108)(aa3_o(i),i=1,nrvr)
98108 format(':aa3,             ',<nrvr>(g12.3,',')
     *       '# channel area coefficient')	  
      write(unitNum,98111)(widep_o(i),i=1,nrvr)
98111 format(':widep,           ',<nrvr>(g12.3,',')
     *       '# channel width to depth ratio')	  
      write(unitNum,98113)(pool_o(i),i=1,nrvr)
98113 format(':pool,            ',<nrvr>(g12.3,',')
     *       '# average area of zero flow pools')	  
!
! The pars below not optimized      
      write(unitNum,98106)(mndr_o(i),i=1,nrvr)
98106 format(':mndr,            ',<nrvr>(g12.3,',')
     *       '# meander channel length multiplier')	  
      write(unitNum,98109)(aa4_o(i),i=1,nrvr)
98109 format(':aa4,             ',<nrvr>(g12.3,',')
     *       '# channel area exponent')	  
      write(unitNum,98115)(fpfactor_o(i),i=1,nrvr)
98115 format(':fpfactor,        ',<nrvr>(g12.3,',')
     *       '# average area of zero flow pools')	  
      write(unitNum,98199)
98199 format(':EndRoutingParameters')
c98199 format(':EndRiverClasses')
      write(unitNum,98098)


!     *************************************************************** 
      write(unitNum,98200)
98200 format(':HydrologicalParameters')
      write(unitNum,98201)classcount
98201 format(':LandCoverClasses,',i12)
      write(unitNum,98202)(nclass(ii),ii=1,classcount)
98202 format(':ClassName       ,',<classcount>(a10,'  ,')	  
     *       '# class name')	
      WRITE(UnitNum,*)'#Vegetationparameters'
      write(unitNum,98222)(fpet(ii),ii=1,classcount)
98222 format(':fpet,            ',<classcount>(g12.3,',')	  
     *       '# interception evaporation factor * pet')	  
      write(unitNum,98223)(ftall(ii),ii=1,classcount)
98223 format(':ftall,           ',<classcount>(g12.3,',')	  
     *       '# reduction in PET for tall vegetation')	  
!     rev. 9.8.40  Jan.  14/13  - NK: convert interception cap: h(,)*fratio()
      if(flnNum.eq.26)then
        write(unitNum,98296)(1.0,ii=1,classcount)
      else
        write(unitNum,98296)(fratio(ii),ii=1,classcount)
      endif
!     rev. 9.8.40  end
98296 format(':fratio,          ',<classcount>(g12.3,',')	  
     *       '# int. capacity multiplier')	 
      write(unitNum,*)'#SoilParameters'      
      write(unitNum,98205)(rec(ii),ii=1,classcount)   
98205 format(':rec,             ',<classcount>(g12.3,',')	  
     *       '# interflow coefficient')	  
      write(unitNum,98206)(ak(ii),ii=1,classcount)    
98206 format(':ak,              ',<classcount>(g12.3,',')	  
     *       '# infiltration coefficient bare ground')	  
      write(unitNum,98207)(akfs(ii),ii=1,classcount)    
98207 format(':akfs,            ',<classcount>(g12.3,',')	  
     *       '# infiltration coefficient snow covered ground')	  
      write(unitNum,98208)(retn(ii),ii=1,classcount)   
98208 format(':retn,            ',<classcount>(g12.3,',')	  
     *       '# upper zone retention mm')	  
      write(unitNum,98209)(ak2(ii),ii=1,classcount)   
98209 format(':ak2,             ',<classcount>(g12.3,',')	  
     *       '# recharge coefficient bare ground')	  
      write(unitNum,98210)(ak2fs(ii),ii=1,classcount)   
98210 format(':ak2fs,           ',<classcount>(g12.3,',')	  
     *       '# recharge coefficient snow covered ground')	  
      write(unitNum,98211)(r3(ii),ii=1,classcount)    
98211 format(':r3,              ',<classcount>(g12.3,',')	  
     *       '# overland flow roughness coefficient bare ground')	  
      
! The pars below not optimized      
      write(unitNum,98203)(ds(ii),ii=1,classcount)    
98203 format(':ds,              ',<classcount>(g12.3,',')	  
     *       '# depression storage bare ground mm')	  
      write(unitNum,98204)(dsfs(ii),ii=1,classcount)    
98204 format(':dsfs,            ',<classcount>(g12.3,',')	  
     *       '# depression storage snow covered area mm')	  
      write(unitNum,98212)(r3fs(ii),ii=1,classcount)    
98212 format(':r3fs,            ',<classcount>(g12.3,',')	  
     *       '# overland flow roughness coefficient snow covered grnd')	  
      write(unitNum,98213)(r4(ii),ii=1,classcount)  
98213 format(':r4,              ',<classcount>(g12.3,',')	  
     *       '# overland flow roughness coefficient impervious area')	  
      write(unitNum,98224)(flint(ii),ii=1,classcount)
98224 format(':flint,           ',<classcount>(g12.3,',')	  
     *       '# interception flag  1=on  <1=off')	  
      write(unitNum,98225)(fcap(ii),ii=1,classcount)
98225 format(':fcap,            ',<classcount>(g12.3,',')	  
     *       '# not used - replaced by retn (retention)')	  
      write(unitNum,98226)(ffcap(ii),ii=1,classcount)
98226 format(':ffcap,           ',<classcount>(g12.3,',')	  
     *       '# wilting point - mm of water in uzs')	  
      write(unitNum,98227)(spore(ii),ii=1,classcount)
98227 format(':spore,           ',<classcount>(g12.3,',')	  
     *       '# soil porosity')	 
      write(unitNum,98298) 
98298 format(':EndHydrologicalParameters')
      write(unitNum,98098)

!     *************************************************************** 
      write(unitNum,98299)
98299 format(':SnowParameters')
      write(unitNum,98214)(fm(ii),ii=1,classcount)
98214 format(':fm,              ',<classcount>(f12.3,',')	  
     *       '# melt factor mm/dC/hour')	  
      write(unitNum,98215)(base(ii),ii=1,classcount)
98215 format(':base,            ',<classcount>(f12.3,',')	  
     *       '# base temperature dC')	  
!     rev. 9.7.29  Jul.  07/11  - NK: Add sublim_rate to set sublimation rate/day to par file
      if(sublimflg.eq.1)then
        write(unitNum,98228)(sublim_factor(ii),ii=1,classcount)
98228   format(':sublim_factor,   ',<classcount>(f12.3,',')
     *       '# sublimation factor ratio')	
      else  
!       sublim_rate changed back to 24 hour value for the par file
        write(unitNum,98229)(sublim_rate(i)*24.0,i=1,classcount)
98229   format(':sublim_rate,     ',<classcount>(f12.3,',')
     *       '# sublimation rate mm/day')	
      endif
      write(unitNum,98234)(sdcd(2,ii),  ii=1,classcount)
98234 format(':sdcd,            ',<classcount>(f12.3,',')
     *       '# swe for 100% snow covered area')	  
      write(unitNum,98216)(fmn(ii),ii=1,classcount)
98216 format(':fmn,             ',<classcount>(f12.3,',')	  
     *       '# -ve melt factor')	  
      write(unitNum,98217)(uadj(ii),ii=1,classcount)
98217 format(':uadj,            ',<classcount>(f12.3,',')	  
     *       '# not used')	  
      write(unitNum,98218)(tipm(ii),ii=1,classcount)
98218 format(':tipm,            ',<classcount>(f12.3,',')	  
     *       '# coefficient for ati')	  
      write(unitNum,98219)(rho(ii),ii=1,classcount)
98219 format(':rho,             ',<classcount>(f12.3,',')	  
     *       '# snow density')	  
      write(unitNum,98220)(whcl(ii),ii=1,classcount)
98220 format(':whcl,            ',<classcount>(f12.3,',')	  
     *       '# fraction of swe as water in ripe snow')	 
      write(unitNum,98221)(alb(ii),ii=1,classcount)
98221 format(':alb,             ',<classcount>(f12.3,',')	  
     *       '# albedo')	 
      write(unitNum,98231)(idump(ii),   ii=1,classcount)
98231 format(':idump,           ',<classcount>(i12,  ',')
     *       '# receiving class for snow redistribution')	  
      write(unitNum,98232)(snocap(ii),  ii=1,classcount)
98232 format(':snocap,          ',<classcount>(f12.3,',')
     *       '# max swe before redistribution')	  
      write(unitNum,98230)(nsdc(ii),    ii=1,classcount)
98230 format(':nsdc,            ',<classcount>(i12,  ',')
     *       '# no of points on scd curve - only 1 allowed')	  
      write(unitNum,98233)(sdcsca(2,ii),ii=1,classcount)
98233 format(':sdcsca,          ',<classcount>(f12.3,',')
     *       '# snow covered area - ratio=1.0')	  
      write(unitNum,98297)
98297 format(':EndSnowParameters')
      write(unitNum,98098)
      

c      write(unitNum,98229)title(45),(diff(i),i=1,12)
c98229 format(':diff,            ',12(f12.3,',')	  
c      write(unitNum,98230)title(46),(hu(i),i=1,12)
c98230 format(':hu,              ',12(f12.3,',')	  
c      write(unitNum,98231)title(47),(pres(i),i=1,12)
c98231 format(':pres,            ',12(f12.3,',')	  

      write(unitNum,98300)
98300 format(':InterceptionCapacityTable ')
!     rev. 9.8.40  Jan.  14/13  - NK: convert interception cap: h(,)*fratio()
      if(flnNum.ne.26)then
!       write the interception capacity as in the bsnm_par.csv file      
        write(unitNum,98301)(h(1,ii),ii=1,classcount)
98301   format(':IntCap_Jan,      ',<classcount>(f12.3,',')	  
     *       '# interception capacity jan mm')	  
        write(unitNum,98302)(h(2,ii),ii=1,classcount)
98302   format(':IntCap_Feb,      ',<classcount>(f12.3,',')	  
     *       '# interception capacity feb mm')	  
        write(unitNum,98303)(h(3,ii),ii=1,classcount)
98303   format(':IntCap_Mar,      ',<classcount>(f12.3,',')	  
     *       '# interception capacity mar mm')	  
        write(unitNum,98304)(h(4,ii),ii=1,classcount)
98304   format(':IntCap_Apr,      ',<classcount>(f12.3,',')	  
     *       '# interception capacity apr mm')	  
        write(unitNum,98305)(h(5,ii),ii=1,classcount)
98305   format(':IntCap_May,      ',<classcount>(f12.3,',')	  
     *       '# interception capacity may mm')	  
        write(unitNum,98306)(h(6,ii),ii=1,classcount)
98306   format(':IntCap_Jun,      ',<classcount>(f12.3,',')	  
     *       '# interception capacity jun mm')	  
        write(unitNum,98307)(h(7,ii),ii=1,classcount)
98307   format(':IntCap_Jul,      ',<classcount>(f12.3,',')	  
     *       '# interception capacity jul mm')	  
        write(unitNum,98308)(h(8,ii),ii=1,classcount)
98308   format(':IntCap_Aug,      ',<classcount>(f12.3,',')	  
     *       '# interception capacity aug mm')	  
        write(unitNum,98309)(h(9,ii),ii=1,classcount)
98309   format(':IntCap_Sep,      ',<classcount>(f12.3,',')	  
     *       '# interception capacity sep mm')	  
        write(unitNum,98310)(h(10,ii),ii=1,classcount)
98310   format(':IntCap_Oct,      ',<classcount>(f12.3,',')	  
     *       '# interception capacity oct mm')	  
        write(unitNum,98311)(h(11,ii),ii=1,classcount)
98311   format(':IntCap_Nov,      ',<classcount>(f12.3,',')	  
     *       '# interception capacity nov mm')	  
        write(unitNum,98312)(h(12,ii),ii=1,classcount)
98312   format(':IntCap_Dec,      ',<classcount>(f12.3,',')	  
     *       '# interception capacity dec mm')	
      else
!       write results\parfile.csv with the converted values:
        write(unitNum,98301)(h(1,ii)*fratio(ii),ii=1,classcount)
        write(unitNum,98302)(h(2,ii)*fratio(ii),ii=1,classcount)
        write(unitNum,98303)(h(3,ii)*fratio(ii),ii=1,classcount)
        write(unitNum,98304)(h(4,ii)*fratio(ii),ii=1,classcount)
        write(unitNum,98305)(h(5,ii)*fratio(ii),ii=1,classcount)
        write(unitNum,98306)(h(6,ii)*fratio(ii),ii=1,classcount)
        write(unitNum,98307)(h(7,ii)*fratio(ii),ii=1,classcount)
        write(unitNum,98308)(h(8,ii)*fratio(ii),ii=1,classcount)
        write(unitNum,98309)(h(9,ii)*fratio(ii),ii=1,classcount)
        write(unitNum,98310)(h(10,ii)*fratio(ii),ii=1,classcount)
        write(unitNum,98311)(h(11,ii)*fratio(ii),ii=1,classcount)
        write(unitNum,98312)(h(12,ii)*fratio(ii),ii=1,classcount)
      endif
      write(unitNum,98399)
98399 format(':EndInterceptionCapacityTable')        
      write(unitNum,98098)
      


!     Default monthly ET values - to be replaced by input
      do j=1,12
	  do ii=1,classcount
	    evap(ii,j)=0.0
	  end do
      end do

      
!      write(unitNum,98400)
98400 format(':MonthlyEvapotranspirationTable ')

!      write(unitNum,98401)(evap(ii,1),ii=1,classcount)
98401 format(':Montly_ET_Jan,   ',<classcount>(f12.1,',')	  
     *       '# monthly evapotranspiration jan mm')	  
!      write(unitNum,98402)(evap(ii,2),ii=1,classcount)
98402 format(':Montly_ET_Feb,   ',<classcount>(f12.1,',')	  
     *       '# monthly evapotranspiration feb mm')	  
!      write(unitNum,98403)(evap(ii,3),ii=1,classcount)
98403 format(':Montly_ET_Mar,   ',<classcount>(f12.1,',')	  
     *       '# monthly evapotranspiration mar mm')	  
!      write(unitNum,98404)(evap(ii,4),ii=1,classcount)
98404 format(':Montly_ET_Apr,   ',<classcount>(f12.1,',')	  
     *       '# monthly evapotranspiration apr mm')	  
!      write(unitNum,98405)(evap(ii,5),ii=1,classcount)
98405 format(':Montly_ET_May,   ',<classcount>(f12.1,',')	  
     *       '# monthly evapotranspiration may mm')	  
!      write(unitNum,98406)(evap(ii,6),ii=1,classcount)
98406 format(':Montly_ET_Jun,   ',<classcount>(f12.1,',')	  
     *       '# monthly evapotranspiration jun mm')	  
!      write(unitNum,98407)(evap(ii,7),ii=1,classcount)
98407 format(':Montly_ET_Jul,   ',<classcount>(f12.1,',')	  
     *       '# monthly evapotranspiration jul mm')	  
!      write(unitNum,98408)(evap(ii,8),ii=1,classcount)
98408 format(':Montly_ET_Aug,   ',<classcount>(f12.1,',')	  
     *       '# monthly evapotranspiration aug mm')	  
!      write(unitNum,98409)(evap(ii,9),ii=1,classcount)
98409 format(':Montly_ET_Sep,   ',<classcount>(f12.1,',')	  
     *       '# monthly evapotranspiration sep mm')	  
!      write(unitNum,98410)(evap(ii,10),ii=1,classcount)
98410 format(':Montly_ET_Oct,   ',<classcount>(f12.1,',')	  
     *       '# monthly evapotranspiration oct mm')	  
!      write(unitNum,98411)(evap(ii,11),ii=1,classcount)
98411 format(':Montly_ET_Nov,   ',<classcount>(f12.1,',')	  
     *       '# monthly evapotranspiration nov mm')	  
!      write(unitNum,98412)(evap(ii,12),ii=1,classcount)
98412 format(':Montly_ET_Dec,   ',<classcount>(f12.1,',')	  
     *       '# monthly evapotranspiration dec mm')	  
!      write(unitNum,98499)
98499 format(':EndMonthlyEvapotranspirationTable')        
!      write(unitNum,98098)

!     END PARAMETERS
!     END PARAMETERS


!     rev. 9.8.01  Jul.  21/11  - NK: added ragmet optimization to dds setup
!     snow parameters
      write(unitNum,98910)
98910 format(':GlobalSnowParLimits')
      write(unitNum,98914)
98914 format('# snowmelt ripening rate')
      write(unitNum,98911)fmadjustdlt
98911 format(':fmadjustdlt,       ',g12.3)	  
      write(unitNum,98912)fmadjustlow
98912 format(':fmadjustlow,       ',g12.3)	  
      write(unitNum,98913)fmadjusthgh
98913 format(':fmadjusthgh,       ',g12.3)	  

      write(unitNum,98924)
98924 format('# min melt factor multiplier')
      write(unitNum,98921)fmalowdlt
98921 format(':fmalowdlt,         ',g12.3)	  
      write(unitNum,98922)fmalowlow
98922 format(':fmalowlow,         ',g12.3)	  
      write(unitNum,98923)fmalowhgh
98923 format(':fmalowhgh,         ',g12.3)	  

      write(unitNum,98929)
98929 format('# max melt factor multiplier')
      write(unitNum,98926)fmahighdlt
98926 format(':fmahighdlt,        ',g12.3)	  
      write(unitNum,98927)fmahighlow
98927 format(':fmahighlow,        ',g12.3)	  
      write(unitNum,98928)fmahighhgh
98928 format(':fmahighhgh,        ',g12.3)	  

      write(unitNum,98935)
98935 format('# glacier melt factor multiplier')
      write(unitNum,98931)gladjustdlt
98931 format(':gladjustdlt,       ',g12.3)	  
      write(unitNum,98932)gladjustlow
98932 format(':gladjustlow,       ',g12.3)	  
      write(unitNum,98933)gladjusthgh
98933 format(':gladjusthgh,       ',g12.3)	  
      write(unitNum,98934)
98934 format(':EndGlobalSnowParLimits')
      write(unitNum,98098)

      write(unitNum,99940)
99940 format(':GlobalParLimits')
      write(unitNum,99944)
99944 format('# precip lapse rate')
      write(unitNum,99941)rlapsedlt
99941 format(':rlapsedlt,       ',g12.3)	  
      write(unitNum,99942)rlapselow
99942 format(':rlapselow,       ',g12.3)	  
      write(unitNum,99943)rlapsehgh
99943 format(':rlapsehgh,       ',g12.3)	  

      write(unitNum,99954)
99954 format('# temperature lapse rate')
      write(unitNum,99951)tlapsedlt
99951 format(':tlapsedlt,       ',g12.3)	  
      write(unitNum,99952)tlapselow
99952 format(':tlapselow,       ',g12.3)	  
      write(unitNum,99953)tlapsehgh
99953 format(':tlapsehgh,       ',g12.3)	  

      write(unitNum,99984)
99984 format('# rain/snow temperature')
      write(unitNum,99986)rainsnowtempdlt
99986 format(':rainsnowtempdlt, ',g12.3)	  
      write(unitNum,99987)rainsnowtemplow
99987 format(':rainsnowtemplow, ',g12.3)	  
      write(unitNum,99988)rainsnowtemphgh
99988 format(':rainsnowtemphgh, ',g12.3)	  

      write(unitNum,99969)
99969 format('# radius of influence')
      write(unitNum,99966)radinfldlt
99966 format(':radinfldlt,      ',g12.3)	  
      write(unitNum,99967)radinfllow
99967 format(':radinfllow,      ',g12.3)	  
      write(unitNum,99968)radinflhgh
99968 format(':radinflhgh,      ',g12.3)	  

      write(unitNum,99975)
99975 format('# smoothing distance')
      write(unitNum,99971)smoothdistdlt
99971 format(':smoothdisdlt,    ',g12.3)	  
      write(unitNum,99972)smoothdistlow
99972 format(':smoothdislow,    ',g12.3)	  
      write(unitNum,99973)smoothdisthgh
99973 format(':smoothdishgh,    ',g12.3)	  
      write(unitNum,99974)
99974 format(':EndGlobalParLimits')
      write(unitNum,98098)

!     API FLAGS & LIMITS
      write(unitNum,98900)
98900 format(':APILimits')
      write(unitNum,98901)a5dlt
98901 format(':a5dlt,             ',g12.3)	  
      write(unitNum,98902)a5low
98902 format(':a5low,             ',g12.3)	  
      write(unitNum,98903)a5hgh
98903 format(':a5hgh,             ',g12.3)	  
      write(unitNum,98904)
98904 format(':EndAPILimits')
      write(unitNum,98098)


!     Routing parameters FLAGS & LIMITS
      write(unitNum,98800)
98800 Format(':RoutingParLimits')
      write(unitNum,98101)(rivtype(i),i=1,nrvr)
      write(unitNum,98801)
98801 format('# lower zone oefficient')
      write(unitNum,98802)(flzdlt(i),i=1,nrvr)
98802 format(':flzdlt,          ',<nrvr>(g12.3,','))	  
      write(unitNum,98803)(flzlow(i),i=1,nrvr)
98803 format(':flzlow,          ',<nrvr>(g12.3,','))	  
      write(unitNum,98804)(flzhgh(i),i=1,nrvr)
98804 format(':flzhgh,          ',<nrvr>(g12.3,','))	  

      write(unitNum,98805)
98805 format('# lower zone exponent')
      write(unitNum,98806)(pwrdlt(i),i=1,nrvr)
98806 format(':pwrdlt,          ',<nrvr>(g12.3,','))	  
      write(unitNum,98807)(pwrlow(i),i=1,nrvr)
98807 format(':pwrlow,          ',<nrvr>(g12.3,','))	  
      write(unitNum,98808)(pwrhgh(i),i=1,nrvr)
98808 format(':pwrhgh,          ',<nrvr>(g12.3,','))	  

      write(unitNum,98810)
98810 format('# channel Manning`s n')
      write(unitNum,98811)(r2ndlt(i),i=1,nrvr)
98811 format(':r2ndlt,          ',<nrvr>(g12.3,','))	  
      write(unitNum,98812)(r2nlow(i),i=1,nrvr)
98812 format(':r2nlow,          ',<nrvr>(g12.3,','))	  
      write(unitNum,98813)(r2nhgh(i),i=1,nrvr)
98813 format(':r2nhgh,          ',<nrvr>(g12.3,','))	  

      write(unitNum,98820)
98820 format('# wetland or bank porosity')
      write(unitNum,98821)(thetadlt(i),i=1,nrvr)
98821 format(':thetadlt,        ',<nrvr>(g12.3,','))	  
      write(unitNum,98822)(thetalow(i),i=1,nrvr)
98822 format(':thetalow,        ',<nrvr>(g12.3,','))	  
      write(unitNum,98823)(thetahgh(i),i=1,nrvr)
98823 format(':thetahgh,        ',<nrvr>(g12.3,','))	  

      write(unitNum,98825)
98825 format('# wetland/bank lateral conductivity')
      write(unitNum,98826)(kconddlt(i),i=1,nrvr)
98826 format(':kconddlt,        ',<nrvr>(g12.3,','))	  
      write(unitNum,98827)(kcondlow(i),i=1,nrvr)
98827 format(':kcondlow,        ',<nrvr>(g12.3,','))	  
      write(unitNum,98828)(kcondhgh(i),i=1,nrvr)
98828 format(':kcondhgh,        ',<nrvr>(g12.3,','))	  

      write(unitNum,98830)
98830 format('# in channel lake retardation coefficient')
      write(unitNum,98831)(rlakedlt(i),i=1,nrvr)
98831 format(':rlakedlt,        ',<nrvr>(g12.3,','))	  
      write(unitNum,98832)(rlakelow(i),i=1,nrvr)
98832 format(':rlakelow,        ',<nrvr>(g12.3,','))	  
      write(unitNum,98833)(rlakehgh(i),i=1,nrvr)
98833 format(':rlakehgh,        ',<nrvr>(g12.3,','))	  

      write(unitNum,97230)
97230 format('# overbank retardation coefficient')
      write(unitNum,97231)(r1ndlt(i),i=1,nrvr)
97231 format(':r1ndlt,        ',<nrvr>(g12.3,','))	  
      write(unitNum,97232)(r1nlow(i),i=1,nrvr)
97232 format(':r1nlow,        ',<nrvr>(g12.3,','))	  
      write(unitNum,97233)(r1nhgh(i),i=1,nrvr)
97233 format(':r1nhgh,        ',<nrvr>(g12.3,','))	  

      write(unitNum,97330)
97330 format('# channel area constant')
      write(unitNum,97331)(aa2dlt(i),i=1,nrvr)
97331 format(':aa2dlt,        ',<nrvr>(g12.3,','))	  
      write(unitNum,97332)(aa2low(i),i=1,nrvr)
97332 format(':aa2low,        ',<nrvr>(g12.3,','))	  
      write(unitNum,97333)(aa2hgh(i),i=1,nrvr)
97333 format(':aa2hgh,        ',<nrvr>(g12.3,','))	  

      write(unitNum,97430)
97430 format('# channel area coefficient')
      write(unitNum,97431)(aa3dlt(i),i=1,nrvr)
97431 format(':aa3dlt,        ',<nrvr>(g12.3,','))	  
      write(unitNum,97432)(aa3low(i),i=1,nrvr)
97432 format(':aa3low,        ',<nrvr>(g12.3,','))	  
      write(unitNum,97433)(aa3hgh(i),i=1,nrvr)
97433 format(':aa3hgh,        ',<nrvr>(g12.3,','))	  

      write(unitNum,97530)
97530 format('# width / depth ratio')
      write(unitNum,97531)(widepdlt(i),i=1,nrvr)
97531 format(':widepdlt,        ',<nrvr>(g12.3,','))	  
      write(unitNum,97532)(wideplow(i),i=1,nrvr)
97532 format(':wideplow,        ',<nrvr>(g12.3,','))	  
      write(unitNum,97533)(widephgh(i),i=1,nrvr)
97533 format(':widephgh,        ',<nrvr>(g12.3,','))	  

      write(unitNum,97630)
97630 format('# pool (fill & spill) channel fraction')
      write(unitNum,97631)(pooldlt(i),i=1,nrvr)
97631 format(':pooldlt,        ',<nrvr>(g12.3,','))	  
      write(unitNum,97632)(poollow(i),i=1,nrvr)
97632 format(':poollow,        ',<nrvr>(g12.3,','))	  
      write(unitNum,97633)(poolhgh(i),i=1,nrvr)
97633 format(':poolhgh,        ',<nrvr>(g12.3,','))	  

      write(unitNum,98834)
98834 Format(':EndRoutingParLimits')
      write(unitNum,98098)


      
!HYDROLOGICAL PARAMETER FLAGS & LIMITS
      write(unitNum,98600)
98600 format(':HydrologicalParLimits')
      write(unitNum,98620)
      write(unitNum,98202)(nclass(ii),ii=1,classcount)
      write(unitNum,*)'#VegetationParameters'
      write(unitNum,98640)
98640 format('# interception evaporation factor * pet')
      write(unitNum,98641)(fpetdlt(ii),ii=1,classcount)
98641 format(':fpetdlt,         ',<classcount>(g12.3,','))	  
      write(unitNum,98642)(fpetlow(ii),ii=1,classcount)
98642 format(':fpetlow,         ',<classcount>(g12.3,','))	  
      write(unitNum,98643)(fpethgh(ii),ii=1,classcount)
98643 format(':fpethgh,         ',<classcount>(g12.3,','))	  

      write(unitNum,98650)
98650 format('# reduction in PET for tall vegetation')
      write(unitNum,98651)(ftalldlt(ii),ii=1,classcount)
98651 format(':ftalldlt,        ',<classcount>(g12.3,','))	  
      write(unitNum,98652)(ftalllow(ii),ii=1,classcount)
98652 format(':ftalllow,        ',<classcount>(g12.3,','))	  
      write(unitNum,98653)(ftallhgh(ii),ii=1,classcount)
98653 format(':ftallhgh,        ',<classcount>(g12.3,','))	  

      write(unitNum,98750)
98750 format('# multiplier for interception capacity')
      write(unitNum,98751)(fratiodlt(ii),ii=1,classcount)
98751 format(':fratiodlt,        ',<classcount>(g12.3,','))	  
      write(unitNum,98752)(fratiolow(ii),ii=1,classcount)
98752 format(':fratiolow,        ',<classcount>(g12.3,','))	  
      write(unitNum,98753)(fratiohgh(ii),ii=1,classcount)
98753 format(':fratiohgh,        ',<classcount>(g12.3,','))	  

      write(unitNum,*)'#SoilParameters'
98620 format('# interflow coefficient')
      write(unitNum,98621)(recdlt(ii),ii=1,classcount)
98621 format(':recdlt,          ',<classcount>(g12.3,','))	  
      write(unitNum,98622)(reclow(ii),ii=1,classcount)
98622 format(':reclow,          ',<classcount>(g12.3,','))	  
      write(unitNum,98623)(rechgh(ii),ii=1,classcount)
98623 format(':rechgh,          ',<classcount>(g12.3,','))	  

      write(unitNum,98601)
98601 format('# infiltration coefficient bare ground')
      write(unitNum,98602)(akdlt(ii),ii=1,classcount)
98602 format(':akdlt,           ',<classcount>(f12.3,','))	  
      write(unitNum,98603)(aklow(ii),ii=1,classcount)
98603 format(':aklow,           ',<classcount>(f12.3,','))	  
      write(unitNum,98604)(akhgh(ii),ii=1,classcount)
98604 format(':akhgh,           ',<classcount>(f12.3,','))	  

      write(unitNum,98610)
98610 format('# infiltration coefficient snow covered ground')
      write(unitNum,98611)(akfsdlt(ii),ii=1,classcount)
98611 format(':akfsdlt,         ',<classcount>(f12.3,','))	  
      write(unitNum,98612)(akfslow(ii),ii=1,classcount)
98612 format(':akfslow,         ',<classcount>(f12.3,','))	  
      write(unitNum,98613)(akfshgh(ii),ii=1,classcount)
98613 format(':akfshgh,         ',<classcount>(f12.3,','))	  

      write(unitNum,98870)
98870 format('# upper zone retention mm')
      write(unitNum,98871)(retndlt(ii),ii=1,classcount)
98871 format(':retndlt,         ',<classcount>(g12.3,','))	  
      write(unitNum,98872)(retnlow(ii),ii=1,classcount)
98872 format(':retnlow,         ',<classcount>(g12.3,','))	  
      write(unitNum,98873)(retnhgh(ii),ii=1,classcount)
98873 format(':retnhgh,         ',<classcount>(g12.3,','))	  

      write(unitNum,98970)
98970 format('# recharge coefficient bare ground')
      write(unitNum,98971)(ak2dlt(ii),ii=1,classcount)
98971 format(':ak2dlt,          ',<classcount>(g12.3,','))	  
      write(unitNum,98972)(ak2low(ii),ii=1,classcount)
98972 format(':ak2low,          ',<classcount>(g12.3,','))	  
      write(unitNum,98973)(ak2hgh(ii),ii=1,classcount)
98973 format(':ak2hgh,          ',<classcount>(g12.3,','))	  
c      write(unitNum,98974)(ak2(ii),ii=1,classcount)
c98974 format(':ak2,             ',<classcount>(g12.3,','))	  

      write(unitNum,98975)
98975 format('# recharge coefficient snow covered ground')
      write(unitNum,98976)(ak2fsdlt(ii),ii=1,classcount)
98976 format(':ak2fsdlt,          ',<classcount>(g12.3,','))	  
      write(unitNum,98977)(ak2fslow(ii),ii=1,classcount)
98977 format(':ak2fslow,          ',<classcount>(g12.3,','))	  
      write(unitNum,98978)(ak2fshgh(ii),ii=1,classcount)
98978 format(':ak2fshgh,          ',<classcount>(g12.3,','))	  
      
      write(unitNum,98630)
98630 format('# overland flow roughness coeff bare ground')
      write(unitNum,98631)(r3dlt(ii),ii=1,classcount)
98631 format(':r3dlt,           ',<classcount>(g12.3,','))	  
      write(unitNum,98632)(r3low(ii),ii=1,classcount)
98632 format(':r3low,           ',<classcount>(g12.3,','))	  
      write(unitNum,98633)(r3hgh(ii),ii=1,classcount)
98633 format(':r3hgh,           ',<classcount>(g12.3,','))	  

      write(unitNum,98979)
98979 format(':EndHydrologicalParLimits')
      write(unitNum,98098)


!     snow parameters
      write(unitNum,98660)
98660 Format(':SnowParLimits')
      write(unitNum,98202)(nclass(ii),ii=1,classcount)
      write(unitNum,98661)
98661 Format('# melt factor mm/dC/hour')
      write(unitNum,98662)(fmdlt(ii),ii=1,classcount)
98662 format(':fmdlt,           ',<classcount>(g12.3,','))	  
      write(unitNum,98663)(fmlow(ii),ii=1,classcount)
98663 format(':fmlow,           ',<classcount>(g12.3,','))	  
      write(unitNum,98664)(fmhgh(ii),ii=1,classcount)
98664 format(':fmhgh,           ',<classcount>(g12.3,','))	  

      write(unitNum,98670)
98670 format('# base temperature dC')
      write(unitNum,98671)(basdlt(ii),ii=1,classcount)
98671 format(':basedlt,         ',<classcount>(g12.3,','))	  
      write(unitNum,98672)(baslow(ii),ii=1,classcount)
98672 format(':baselow,         ',<classcount>(g12.3,','))	  
      write(unitNum,98673)(bashgh(ii),ii=1,classcount)
98673 format(':basehgh,         ',<classcount>(g12.3,','))	  

      write(unitNum,98770)
98770 format('# sublimation factor OR ratio')
      write(unitNum,98771)(subdlt(ii),ii=1,classcount)
98771 format(':subdlt,          ',<classcount>(g12.3,','))	  
      if(sublimflg.eq.1)then
        write(unitNum,98772)(sublow(ii),ii=1,classcount)
98772   format(':sublow,          ',<classcount>(g12.3,','))	  
        write(unitNum,98773)(subhgh(ii),ii=1,classcount)
98773   format(':subhgh,          ',<classcount>(g12.3,','))	  
      else
        write(unitNum,98772)(sublow(ii)*24,ii=1,classcount)
        write(unitNum,98773)(subhgh(ii)*24,ii=1,classcount)
      endif

      write(unitNum,97170)
97170 format('# sdcd')
      write(unitNum,97171)(sdcddlt(ii),ii=1,classcount)
97171 format(':sdcddlt,         ',<classcount>(g12.3,','))	  
      write(unitNum,97172)(sdcdlow(ii),ii=1,classcount)
97172 format(':sdcdlow,         ',<classcount>(g12.3,','))	  
      write(unitNum,97173)(sdcdhgh(ii),ii=1,classcount)
97173 format(':sdcdhgh,         ',<classcount>(g12.3,','))	  

      
      write(unitNum,98674)
98674 Format(':EndSnowParLimits')
      write(unitNum,98098)
      

c      print*,'flnNum=',flnNum
      
      close(unitNum,status='keep')
c	print*,'Closed unit=',unitNum
	
      Write(*,*)filename(flnNum)(1:50),'written'

      
      if(iopt.eq.2)print*, 'param before return @ 9989'

 9989 RETURN

! FORMATS

 1000 format(a5,g10.3,2i5,a65)
 1001 format(a5,17e10.3)
 1002 format(a5,4(e12.3))
 1003 format(a5,12f5.2)
 1004 format(a5,17e10.3)
 1005 format(a5,17e10.3)
 1006 format(a5,17e10.3)
 1007 format(a5,17e10.3)
 1008 format(a5,17e10.3)
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
 5009 format(5x,17a10)
 5010 format(a1)
 5100 format(' in rdpar - problem opening basin/evap.dat file'/
     *         '          zero values are inserted for evap.dat'/)
5110  format(//,' Value for a9 heat deficit to swe outside range',/
     *     f10.3,' assumed')
5111  format(//,' Value for a10 uz exponent outside range',/
     *     f10.3,' assumed')
 5999 format(' mflag=',i5,'   errflg=',i5)
 6003 format(a5,17e10.3)



 8429 format(//' in runof5,flz/',5f10.6/
     * ' wrong parameter value - change flz(i) to +ve value'/)
 9801 format(a5,17f10.0)
 9802 format(a5,17f10.2)
 9803 format(a5,12f5.1)
 9804 format(a5,f10.2,a60)
 9805 format(a5,i10,5x,a40)
 9806 format(a5,f10.0,5x,a40)
 9807 format(a5,f10.3,5x,a40)
 9808 format(//' Parameters used for this run:'/)
 9809 format('A new par file called basin\new.par has been created')

      RETURN

      END SUBROUTINE write_par_10

