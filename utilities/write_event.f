      SUBROUTINE write_event(date,conv,scale,smc5,mo,nhr,nhf,
     *        yy1,mm1,dd1,hh1,evtname,evtlength,evtType)


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
c      USE area2
c      USE area12


      use area_watflood

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

!      include 'debug.for'
	
      CHARACTER(12) :: date
      DIMENSION     :: smc5(5)
      integer  :: ios,nhr,nhf,evtlength
      real*4   :: conv,scale,smc5
      character*2   :: mm1,dd1,hh1
      character*4   :: yy1
	character*8   :: evtname,evtType

! WRITE A NEW EVENT.evt FILE:
! SEE INPEVTA FOR VARIABLE DEFS

      open(unit=99,file=fln(99),status='unknown',iostat=ios)
      if(ios.ne.0)then
        print*,'problem opening',fln(99),'in write_event.for @ line 15'
        print*,'ios=',ios,' possible problem: lile is locked'
        stop 'program aborted in write_event'
      endif

      write(99,99000)'#'
      write(99,99001)':filetype                     .evt '
      write(99,99001)':fileversionno                9.9  '
      write(99,99002)':year                         ',yy1
      write(99,99003)':month                        ',mm1
      write(99,99003)':day                          ',dd1
      write(99,99003)':hour                         ',hh1
      write(99,99000)'#'
      write(99,99006)':snwflg                       ',snwflg
      write(99,99001)':sedflg                       n    '
      write(99,99006)':vapflg                       ',vapflg
      write(99,99006)':smrflg                       ',smrflg
      write(99,99006)':resinflg                     ',resinflg
      write(99,99001)':tbcflg                       n    '
      write(99,99001)':resumflg                     n    '
      write(99,99001)':contflg                      n    '
      write(99,99001)':routeflg                     n    '
      write(99,99001)':crseflg                      n    '
      write(99,99001)':kenueflg                     n    '
      write(99,99001)':picflg                       n    '
      write(99,99006)':wetflg                       ',wetflg
      write(99,99001)':modelflg                     n    '
      write(99,99001)':shdflg                       n    '
      write(99,99001)':trcflg                       n    '
      write(99,99001)':frcflg                       n    '
      write(99,99001)':initflg                      n    '
      write(99,99001)':hdrflg                       n    '
      write(99,99001)':grdflg                       n    '
      write(99,99001)':ntrlflg                      n    '
      write(99,99001)':nudgeflg                     n    '
      write(99,99001)':resetflg                     n    '
      write(99,99006)':divertflg                    ',divertflg
      write(99,99001)':pafflg                       n    '
      write(99,99001)':fliflg                       n    '
      write(99,99006)':lakeflg                      ',lakeflg
      write(99,99006)':iceflg                       ',iceflg
      write(99,99001)':FEWSflg                      n    '
      write(99,99001)':netCDFoutflg                 n    '

      write(99,99000)'#'
      write(99,99004)
     *    ':intsoilmoisture              0.25 0.25 0.25 0.25 0.25'
      write(99,99001)':rainconvfactor                1.00'
      write(99,99001)':eventprecipscalefactor        1.00'
      write(99,99001)':precipscalefactor             0.00'
      write(99,99001)':eventsnowscalefactor          0.00'
      write(99,99001)':snowscalefactor               0.00'
      write(99,99001)':eventtempscalefactor          0.00'
      write(99,99001)':tempscalefactor               0.00'
	write(99,99011)':disaggregate                  ',smearfactor

      write(99,99000)'#'
      write(99,99010)':hoursraindata                 ',evtlength
      write(99,99010)':hoursflowdata                 ',evtlength
      write(99,99010)':deltat_report                 ',deltat_report
      write(99,99010)':spinupevents                  ',0

      write(99,99000)'#'
      write(99,99005)':basinfilename                 ',fln(1)
      write(99,99005)':parfilename                   ',fln(2) 
      write(99,99005)':channelparfile                ',fln(3)
      write(99,99005)':pointdatalocations            ',fln(4)
c      write(99,99005)':snowcoverdepletioncurve       ',fln(5)
c      write(99,99005)':waterqualitydatafile          ',fln(6)
c     write(99,99005)':climateNormalsDiff            ',fln(7)

      write(99,99000)'# station climate data'
      write(99,997)
     * ':pointsoilmoisture            moist\',evtname,'_psm.pt2'
      write(99,997)
     * ':pointprecip                  raing\',evtname,'_rag.tb0'
      write(99,997)
     * ':pointHourlyPrecip            hrlpp\',evtname,'_rag.tb0'
      write(99,997)
     * ':pointDailyPrecip             dlypp\',evtname,'_rag.tb0'
      write(99,997)
     * ':pointtemps                   tempg\',evtname,'_tag.tb0'
c      write(99,997)':pointnetradiation                  '
c      write(99,997)
c     * ':pointhumidity                humid\',evtname,'_hum.tb0'
c      write(99,997)
c     * ':pointwindspd                 winds\',evtname,'_spd.tb0'
c      write(99,997)
c     * ':pointwinddir                 winds\',evtname,'_dir.tb0'
c      write(99,996)':pointlongwave                      '
c      write(99,996)':pointshortwave                     '
c      write(99,996)':pointatmpressure                   '
!     next 3 added for ver = 9.8  nk  oct 31/07
      write(99,997)
     * ':pointsnow                    snowg\',evtname,'_snw.tb0'
      write(99,997)
     * ':pointdrain                   drain\',evtname,'_drn.tb0'
      write(99,997)
     * ':pointdsnow                   dsnow\',evtname,'_dsn.tb0'

      write(99,99000)'# hydrometric station data'
      write(99,997)
     * ':streamflowdatafile           strfw\',evtname,'_str.tb0'
      write(99,997)
     * ':reservoirreleasefile         resrl\',evtname,'_rel.tb0'
      write(99,997)
     * ':reservoirinflowfile          resrl\',evtname,'_rin.tb0'
      write(99,997)
     * ':diversionflowfile            diver\',evtname,'_div.tb0'
      write(99,997)
     * ':snowcoursefile               snow1\',evtname,'_crs.pt2'
      write(99,997)
     * ':initlakelevel                level\',evtname,'_ill.pt2'
      write(99,997)
     * ':observedlakelevel            level\',evtname,'_lvl.tb0'
      write(99,997)
c     * ':flowinit                     strfw\',evtname,'_fli.tb0'
      write(99,99000)'# gridded data'

      if(evtType(1:2).eq.'na')then
        write(99,90001)evtname
90001   format(':griddedmodelfile             model\',a8,'_capa.r2c')
      else
        strLength = LEN_TRIM(evtType)
        write(99,90002)evtname,evtType
90002   format(':griddedmodelfile             model\',
     *                                   a8,'_',a<strlength>,'_pcp.r2c')
      endif  
        
        
      write(99,997)
     * ':griddedinitsnowweq           snow1\',evtname,'_swe.r2c'
      write(99,997)
     * ':griddedinitsoilmoisture      moist\',evtname,'_gsm.r2c'
      write(99,996)':griddedinitlzs                     '
      
     
c      write(99,997)
c    * ':griddedrainfile              radcl\',evtname,'_met.r2c'
      if(evtType(1:2).eq.'na')then
        write(99,90003)evtname
90003   format(':griddedrainfile              radcl\',a8,'_met.r2c')
      else
        strLength = LEN_TRIM(evtType)
        write(99,90004)evtname,evtType
90004   format(':griddedrainfile              radcl\',
     *                                   a8,'_',a<strlength>,'_pcp.r2c')
      endif  
       
       write(99,997)
     * ':griddedHourlyPrecip          hrlgp\',evtname,'_met.r2c'
       write(99,997)
     * ':griddedDailyPrecip           dlygp\',evtname,'_met.r2c'
      write(99,996)':griddedsnowfile                    '
      
      
c      write(99,997)
c     * ':griddedtemperaturefile       tempr\',evtname,'_tem.r2c'
      if(evtType(1:2).eq.'na')then
        write(99,90005)evtname
90005   format(':griddedtemperaturefile       tempr\',a8,'_tem.r2c')
      else
        strLength = LEN_TRIM(evtType)
        write(99,90006)evtname,evtType
90006   format(':griddedtemperaturefile       tempr\',
     *                                   a8,'_',a<strlength>,'_tmp.r2c')
      endif  
      
c      write(99,997)
c     * ':griddeddailydifference       tempr\',evtname,'_dif.r2c'
      if(evtType(1:2).eq.'na')then
        write(99,90007)evtname
90007   format(':griddeddailydifference       tempr\',a8,'_dif.r2c')
      else
        strLength = LEN_TRIM(evtType)
        write(99,90008)evtname,evtType
90008   format(':griddeddailydifference       tempr\',
     *                                   a8,'_',a<strlength>,'_dif.r2c')
      endif  
      
      
      
c      write(99,996)':griddednetradiation                '
      write(99,997)
     * ':griddedhumidity              humid\',evtname,'_hum.r2c'
      write(99,997)
     * ':griddedwindspd               winds\',evtname,'_spd.r2c'
      write(99,997)
     * ':griddedwinddir               winds\',evtname,'_dir.r2c'
c      write(99,996)':griddedlongwave                    '
c      write(99,996)':griddedshortwave                   '
c      write(99,996)':griddedatmpressure                 '

!     next 3 added for ver = 9.8  nk  oct 31/07
      write(99,997)
     * ':griddedsnow                  snowg\',evtname,'_snw.r2c'
      write(99,997)
     * ':griddeddrain                 drain\',evtname,'_drn.r2c'
      write(99,997)
     * ':griddeddsnow                 dsnow\',evtname,'_dsn.r2c'

      write(99,997)
     * ':griddedrunoff                runof\',evtname,'_rff.r2c'
      write(99,997)
     * ':griddedrecharge              rchrg\',evtname,'_rch.r2c'
      write(99,997)
     * ':griddedleakage               lkage\',evtname,'_lkg.r2c'
      write(99,99000)'#'
      write(99,996)':noeventstofollow                 00'
      write(99,99000)'#'
      write(99,99007)'eof'

      close (unit=99)
99000 format(a1)
99001 format(a35)
99002 format(a30,a4)
99003 format(a30,a2)
99004 format(a54)
99005 format(a30,a30)
99006 format(a30,a1)
99007 format(a3)
99008 format(a30,a40)
99009 format(a25,i5)
99010 format(a30,i5)
99011 format(a30,f5.2)
996   format(a36,a8,a4)
997   format(a36,a8,a8)
998   format(a36,a8,a9)
999   format(a36,a8,a10)
      RETURN

      END SUBROUTINE write_event
