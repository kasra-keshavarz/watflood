      module area_watflood

      
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
      
      
      
ccccc      MODULE area1

c!DEC$ ATTRIBUTES DLLEXPORT :: area_watflood

!     entrees rearranged in alphabetical order.  nk  Mar. 7/07
c     rev. 9.5.20  Mar.  06/08  - NK: added resvstore for iso mosed
!     rev. 9.5.48  Dec.  26/08  - NK: added event_fln() to allow unlimited events
!     rev. 9.5.58  Apr.  16/09  - NK: added nudgeflg for forcing gauge flows
!     rev. 9.5.62  Sep.  04/09  - NK: new tb0 file for DW routing
!     rev. 9.5.68  Oct.  07/09  - NK: debugged read_resvin_ef.f
!     rev. 9.6.02  Mar.  15/10  - NK: add sublimation to optimization
!     rev. 9.6.06  Apr.  18/10  - NK: added glacier adjust for optimization
!     rev. 9.8.21  Jun.  18/12  - NK: Assed swe observed date & report
!     rev. 9.8.22  Jul.  17/12  - NK: Added resetflg to reset cumm. precip Sept.1
!     rev. 9.8.24  Aug.  07/12  - NK: Added reading yyyymmdd_lvl.tb0 for lake levels
!     rev. 9.8.51  Mar.  11/13  - NK: Link skiphours in s/r stats to value1 in the str file
!     rev. 9.8.53  Mar.  20/13  - NK: Add Lake St. Joseph diversion algorithm to REROUT.f
!     rev. 9.8.57  Apr.  12/13  - NK: Added lakeEflg to stop lake evaporation whan levels very low
!     rev. 9.8.59  May   14/13  - NK: REmoved psmear & punused from the program
!     rev. 9.8.90  Oct.  30/13  - NK: Added fetch to the shd file 
!     rev. 9.9.00  Dec.  08/13  - NK: Added Lake Evaporation model
!     rev. 9.9.20  Jul.  24/14  - NK: Added dead storage for lakes "stroe_dead"
!     rev. 9.9.38  Nov.  12/14  - NK: Added LKdepth to ill file
!     rev. 9.9.49  Jan.  06/14  - NK: Added courantflg
!     rev. 9.9.62  Mar.  21/15  - NK: Change zone from character to integer
!                                    also zone1. zone2, zone3 & zone_temp
!     rev. 9.9.65  Apr.  03/15  - NK: Added rule s/r; resrl\rules.txt & ruleflg
!     REV. 10.1.41 Oct   11/16  - NK: Added tb0flg to write lake_*.tb0 files
!     rev. 10.1.64 Jan.  26/17  - NK: Added XML output file 
!     rev. 10.2.13 Jan.  31/18  - NK: Re-wrote rules.f to mimic stop log operations    -> rules_sl.f
!     rev. 10.2.14 Jan.  31/18  - NK: Renamed rules.f to rules_tl.f - for use with target levels
!     rev. 10.2.20 Apr.  06/18  - NK: Added res_next for NBS
!     rev. 10.2.27 Jul.  08/18  - NK: replaced lake_inflow_sum with temp_flow_sum


      character*10         :: program_name,start_date
      character*12         :: yyyymmdd12(366) ! used for writing the xml file
      integer*4      ::       maxr,nl_max,nrel_max,no_data_lines_read
      integer*4      ::       skiplines,error_msg,lowgrid,lowid
      real(4)        ::       min_flow_cutoff
      logical        ::       skipflg,courantflg 
!     rev. 10.1.98 Oct   04/17  - NK: Deal with -ve flows in route
      logical        ::       warning(100)
      real*4         ::       por,sqlz,slzinflw,leakage,dacheck,
     *                        totaltime,qlow,lowtime,
     *                        par_file_version,par_file_version_latest
      character*10   ::       program_version,program_date

      character*1, dimension(:),allocatable :: glacier_flag

      integer*4,dimension(:), allocatable   :: next,ibn,ichnl,irough,
     *                        ireach,nreach,res,res_next,xxx,yyy
     *                                       
      logical,dimension(:),   allocatable      ::lakeEflg
      logical,dimension(:),   allocatable      ::flag1
!     rev. 9.9.65  Apr.  03/15  - NK: Added rule s/r; resrl\rules.txt & ruleflg
      logical,dimension(:),   allocatable      ::resRuleFlg
!     rev. 10.2.13 Jan.  31/18  - NK: Re-wrote rules.f to mimic stop log operations    -> rules_sl.f
!     rev. 10.2.14 Jan.  31/18  - NK: Renamed rules.f to rules_tl.f - for use with target levels
!     rev. 10.2.35 Oct.  08/18  - NK: Moved logical def. to area_watflood
      logical, dimension(:), allocatable  :: printwarning
      logical           ::    ruleflg
      character*12      ::    ruletype
      

!     rev. 9.9.38  Nov.  12/14  - NK: Added LKdepth to ill file
      real*4, dimension(:),   allocatable :: 
     *                        att,bnkfll,bin_precip,ch_length,da,d2,
     *                        cap,fake,fakefs,frac,
     *                        evap_convert,grid_area,channel_area,
     *                        nca_1d,lake_area,lake_tmp, 
     *                        lzs,over,pot,potfs,      !psmear,punused,
     *                        qi1,qi2,qo1,qo2,qda,qr,qlz,qbase,qmax,
!     rev. 10.2.44 Jan.  21/19  - NK: Added qOld() to allow previous weighted outflow on the resume file
     *                        qqlow,qOld,     
     *                        qiob1,qiob2,qoob1,qoob2,      !rev 9.5
     *                        qUS1,qUS2,     !     REV. 10.1.24 Jan.  30/16  - NK: Added qUS1 & qUS2 for watbal
     *                        qstream,qstrm,qdrng,qdrngfs,
     *                        rechrg,rl,
     *                        sl1,slope,sl2,sq1,sq1fs,sqint,sqintfs,
     *                        store1,store2,storinit,sr,strloss,
     *                        sdrng,sdrngfs,sexcess,
     *                        sump,sumrff,sumqint,sumqintfs,
     *                        sumq1,sumq1fs,sumrechrg,
!     rev. 10.4.21 Apr.  21/20  = NK Add UZS deficit to wfo file = UZS(class=classcount)
     *                        totd1,totuzs,totuzsdef,totsnw,totchnl,
     *                        totgrid,x4,x5,resvstore1,resvstore2,
     *                        sum_grid_inflow,totwetl,
     *                        netinflow,netoutflow,
     *                        LKdepth,LKinvert,
     *                        ice_fctr,ice_fctr_min,ice_fctr_max
!     rev. 10.1.11 Dec.  11/15  - NK: Revised ice factor initialization and calculation   
     
!     rev. 9.9.65  Apr.  03/15  - NK: Added rule s/r; resrl\rules.txt & ruleflg
!     rev. 9.9.73  Aug.  31/15  - NK: Finshed rules s/r - ready for beta testing
      integer*4                           :: nrules
      integer, dimension(:),  allocatable :: ruleNo,resvNo
!     rev. 10.2.37 Oct.  18/18  - NK: Changed target levels to real*8 
      real*8,  dimension(:,:),allocatable :: upper_range,lower_range
      real*4,  dimension(:),  allocatable :: safe_max,max_range,qmin
      real*4,  dimension(:),  allocatable :: DecayT
      REAL*4,  DIMENSION(:),  ALLOCATABLE :: QQsum,QRsum
                   
      real*4,  dimension(:),  allocatable :: sta_area,areasum,area_error
      real*4,  dimension(:,:),allocatable :: areaclass                 !,nca
      real*4                              :: effective_nca

!     note:  frac is to be replaced by grid_area


!     rev. 9.1.30  Nov.  08/02  - added q1, qint & drng to the wfo file
!     REV. 10.1.38 Jul   28/16  - NK: Added noDataValue to WFO & tb0 files
      real*4, dimension(:,:),  allocatable::
     *                         aclass,d1,d1fs,effpor,drng,drngfs,
     *                         df,dffs,q1,q1fs,qint,qintfs,
     *                         snow,sumf,sumffs,r,v,uzs,uzsfs,
     *                         qdrng2,qdrngfs2,qdf,qdffs,dprecip
!    *                                       intcap

!     rev. 9.4.02  Apr.  18/07  - NK: moved rf, rffs from areawq to area1
      real*4, dimension(:,:),  allocatable:: rf,rffs
      real*4                   :: noDataValue,wfo_spec_version_number
!     rev. 10.4.33 May.  26/21  = NK Added totUZS, API & heatdef to the wfo file
      integer  npick
	logical, dimension(:),  allocatable:: store_error_flag
	logical, dimension(:),  allocatable:: wetland_flag
!     rev. 9.8.29  Oct.  15/12  - NK: added wetland_flag to speed up route.f

cccc      END MODULE area1


! NOTES:
!
! NOV. 2000, TS: Added QLZ variable during WATFLOOD swamp routing changes 
!                b/c needed to access it from routea (not just runof5a).
! AUG 27/03  TS: Added sumq1,sumq1fs,sumqint,and sumqintfs,qdrngfs arrays and made
!                qdrng and qdrngfs arrays.  FOR TRACER S/R
   



ccccc      MODULE area2


!     general program variables

!       flag for calling program
        character*10    calling_program_flg

        integer*1,      dimension(:), allocatable :: wfo_pick
        character*50,   dimension(:), allocatable :: wfo_attributes
        character*10,   dimension(:), allocatable :: column_units
        character*10,   dimension(:), allocatable :: column_type
        integer      :: iopt,itype,iymin,iymax,jxmin,jxmax,
     *                  ib,it,no,ndam,ni,id,nl,mhtot,mhrd,kt,
     *                  irads,iradn,jradw,jrade,iyoffset,jxoffset,
     *                  iyshift,jxshift,istep,nblock,
     *                  xcount_precip,ycount_precip,
     *                  ireport,ireport_start,ireport_end,
     *                  ioflg,ichsm,nnprint,iiprint,iopt_start,
     *                  ktri,glacier_class_number
!     REV. 10.1.29 May   04/16  - NK: Added parfile comments
      integer*4         ParFileCommentNumber
      character(72)     ParFileComment(100)
      logical        :: iopt1,iopt2,iopt3,iopt4,iopt99


!     rev. 9.8.24  Aug.  07/12  - NK: Added reading yyyymmdd_lvl.tb0 for lake levels
        logical      :: lvlflg
        integer      :: ktlvl,nolvl,nlvl
        integer,        dimension(:),   allocatable :: ewg_lvl,sng_lvl
        integer,        dimension(:),   allocatable :: lvl_reach,noptlvl
        real*4,         dimension(:),   allocatable :: xlvl,ylvl
        real*4,         dimension(:,:), allocatable :: lvl_obs,lvl_calc
	  character(12),  dimension(:),   allocatable :: gname_lvl

!     rev. 9.8.49  Feb.  20/13  - NK: Added n=municipal & irrigation withdrawals
!       withdrawal amounts will be done with each time step on a monthly basis        
        logical      :: withdrawflg
        real*4,         dimension(:,:), allocatable :: qwdr
     
        integer*4    :: rgrd,rads,radw,rade,radn

        real*4       :: rgrdn,rgrde     !changed to real Jul. 27/04 nk
        real*4       :: evt_version

        integer*4    :: mm,yy,iall,nclt,nch,nnch,irdt,ii_water,ii_wetl
        integer(4)   :: year,month,day,hrs,mins,secs,hsecs,hh,dd
!        integer*4, parameter :: flen=9000
!        integer*4, parameter :: flen='16384'
        integer*4, parameter :: flen='none'
c      real(4)  :: scaleallsnw
        
        real*4       :: ver,grdn,grde,al,astep,step2,readscale,
     *                     scalesnw,readscalesnw,scaleallsnw,scaleall,
     *                     scaletem,scalealltem,sstep,rdt,flowunitconv,
     *                     smearfactor
	  logical(1)   :: keyflg,precflg,found_data_end,rd_evp_flg,optflg
	  logical(1)   :: newparflg,sensitivityflg,debugflg
        CHARACTER(5) :: title(200)
        character(40) :: notes(100)
        CHARACTER(80) :: heading(10)
!       flags flags flags flags flags flags flags & more flags         
!     REV. 10.1.41 Oct   11/16  - NK: Added tb0flg to write lake_*.tb0 files
!     rev. 10.1.64 Jan.  26/17  - NK: Added XML output file 
!     rev. 10.2.28 Jul.  08/18  - NK: Revised for TSW = NBS without lake precip/evap'n
!     rev. 10.4.42 Sep.  29/21  = NK Added flowresetflg to the event file
        character(1) :: snwflg,sedflg,vapflg,smrflg,resinflg,
     *                  resumflg,tbcflg,contflg,routeflg,crseflg,
     *                  ensimflg,leapflg,llflg,picflg,wetflg,
     *                  modelflg,shdflg,wfo_open_flg,trcflg,frcflg,
     *                  newevtflg,manningflg,translateflg,flowfillflg,
     *                  outfileflg,initflg,frc_file_flg,cum_precip_flg,
     *                  grdflg,ntrlflg,nudgeflg,nudgeflg1,divertflg,
     *                  pafflg,fliflg,trcOffFlg,xmlflg,nbsflg,
     *                  fewsflg,netCDFoutflg,flowresetflg,
     *                  lakeflg,iceflg,icerivflg,icelakeflg,tb0flg
        logical         resetflg,hdrflg0,dlyflg,icefactorfile,netCDFflg
        logical     ::  writeflg(999)


!     rev. 10.1.88 May   23/17  - NK: Fixed Julian_day problems for iso R/W
        logical         leapyear
        character(1) :: ssmc_firstpass
!        character(1), dimension(:), allocatable :: glacier_flag
!       these things taken out of some argument lists  nk 05/10/04
        integer      :: year1,mo1,day1,hour1,jul_day1
	  integer      :: year_now,month_now,day_now,hour_now,jul_day_now
        integer      :: min_now,sec_now,mins_precip_start,mins_tmp_start
	  integer      :: year_fc,month_fc,day_fc,hour_fc,jul_day_fc
	  integer      :: year_now1,month_now1,day_now1,hour_now1,jul_day_now1
	  integer      :: year_now2,month_now2,day_now2,hour_now2,jul_day_now2
        integer      :: hour_now_fews
!     rev. 10.1.01 Oct.  05/15  - NK: Isotope update: added 2H
!     rev. 10.1.04 Oct.  10/15  - NK: Added year_last variable for use in reading isotope data
	  integer      :: year_last  ! needed for reading annual recorded isotope data
!       *start used to save the very first date for later use in the program        
	  integer(4)   :: year_start,mo_start,day_start,hour_start
        integer      :: epoch_min_start
!       *swe is the date & time for a swe update in FEWS
!       update time for snow1\swe.nc  time series file
!       update (mult or add) from snow1\swe_update.xml
        integer      :: yy_swe,mm_swe,dd_swe,hh_swe
        integer      :: yy_uzs,mm_uzs,dd_uzs,hh_uzs
        real(4)      :: swe_add,swe_mult,uzs_add,uzs_mult
        logical      :: swe_update,uzs_update,use_nc_swe_update
        logical      :: swe_update_done
	  character(2) :: yy2,mm2,dd2,hh2
c        character(4) :: yyyy4

        character(5) :: source,rdr,data_source
        character(80):: querystring
!     rev. 10.1.90 Jul.  27/17  - NK: Added date_now for i/o files
        character(10)   :: date_now
        character(8)    :: time_now
!     rev. 9.1.55  Jun.  12/04  - NK: write new str files to strfw\newfmt folder.
        character(10)   :: coordsys1,datum1
	  character(40)   :: attribute_name,attribute_units,attribute_type,
     *                      application,author,char_block,name
!       2 attribute names here:
!         one for writing r2c files with just 1 attribute    attribute_name
!         another for files where there are any number of attributes  
	  character(40),  dimension(:),   allocatable :: attributeName
 	  character(40),  dimension(:),   allocatable :: attributeType
        character(30)   :: source_file_name
	  integer         :: attribute_count,zone1
        real*4          :: init_heat_deficit,unit_conversion
!       this set of variables needed when reading files other than the 
!       shed file so coincidence of data can be checked
!       fix fix   do this in all rd**** files someday

!       for the gridded precip file in rdrain
        character(10) :: coordsys2,datum2
        real          :: xorigin2,yorigin2,xdelta2,ydelta2,
     *                     dtrain,convrain
        integer       :: xcount2,ycount2,deltat2,nhrain,zone2

!       for the gridded temperature file in rdtemp
        character(10) :: coordsys3,datum3
        real          :: xorigin3,yorigin3,xdelta3,ydelta3,
     *                     dttemp,convtemp
        integer       :: xcount3,ycount3,deltat3,nhtemp,zone3

!       for the generic read & write modules (e.g write_r2c
        character(10) :: coordsys_temp,datum_temp
        real          :: xorigin_temp,yorigin_temp,
     *                     xdelta_temp,ydelta_temp,UnitConv_temp

        integer       :: xcount_temp,ycount_temp,deltat_temp,
     *	               deltat_report,no_hdrcomments,zone_temp

	  character(40), dimension(:), allocatable:: hdrcomment
	  real*4,        dimension(:), allocatable:: hdrcomment_value
!     rev. 9.9.76  Sep.  11/15  - NK: Added recorded isotope concentrations
	  real*4,        dimension(:), allocatable:: x_temp,y_temp       !     rev. 9.9.76


!       please note that coordsys and datum are used elsewhere.

!      integer :: jan,ii,n,i,j,i3,ii1,ii2
!      real ::    aintvl


ccccc      END MODULE area2

!      llflg  - when 'y', coordinates for .str .snw etc in lat-long
!      ioflg  - if .ge. 1 read the outfiles (note: integer)
!       leapflg  - to indicate a leap year - not an input

!	snwflg - whether there is snow to melt
!	sedflg - whether the sediment routine is used y or n
!	vapflg - turn on evap routine (Todd)
!	smrflg - turns on smearing (smear precip data over data dt)
!	resinflg - will use resin record for comparison
!	tbcflg - read resume.txt file for run init values
!	resumflg - resume.txt file written at end of run (mem dump)
!     contflg  - for continuing statistics upon resume = input
!     routeflg - output qr grids for routing program 
!     crseflg  - read snow course data to replace resume file data
!     ensimflg - write the wfo file for ENSIM
!     ensimflg1 = 'y' for first time needed, else = 'n' (inpevt)
!     wfoflg - set to 'y' initially. Changed to 'n' once wfo header=written
!     picflg   - write the simout/pic.txt file for mapper
!     wetid1flg- run the wetland routing module - read in first event only
!     modelflg  - pick model
!     shedflg  - replace the watershed file basin\bsnm.shd

!	source - what is the data source - radar, mc2, erf, rag, etc.




ccccc      MODULE area3

        integer*4 :: na,naa,ij,ji,ls,ks,js,ijk,ih,ipr,jpr,iw,ntype,
     *             nbsn,nrvr,nsnow,nlz,net,nastart,naend

        integer*4 :: classcount,classcount_temp

ccccc      END MODULE area3


ccccc      MODULE area4
!     rev. 9.5.06  Feb.  05/08  - NK: added pool and pool_o in rdpar & route
!     REV. 10.1.17 Jan.  11/16  - NK: Added fpetLakeOverride factor
!     rev. 10.5.09 May   11/23  = NK Added fpFactor - flood plain width/depth .17 = 100/1

      real*4,    dimension(:), allocatable ::  r1,r4,ds,dsfs,chnl,
     *                         r2,r3,r3fs,rec,ak,akfs,r2low,r3low,
     *                         r3fslow,reclow,aklow,akfslow,ak2fslow,
     *                         r2hgh,r3hgh,r3fshgh,rechgh,akhgh,akfshgh,
     *                         ak2fshgh,r2dlt,r3dlt,r3fsdlt,recdlt,
     *                         akdlt,akfsdlt,ak2fsdlt,retn,ak2,flz,flz2,
     *                         pwr,pwr2,retnlow,ak2low,flzlow,pwrlow,
     *                         retnhgh,ak2hgh,flzhgh,pwrhgh,retndlt,
     *                         ak2dlt,flzdlt,pwrdlt,retfs, ak2fs,
     *            fpet,fpetdlt,fpetlow,fpethgh,
     *            ftall,ftalldlt,ftalllow,ftallhgh,
     *            fratio,fratiodlt,fratiolow,fratiohgh,
     *            mndr,aa2,aa3,aa4,pool,
     *            thetadlt,thetalow,thetahgh,widepdlt,wideplow,widephgh,
     *            kconddlt,kcondlow,kcondhgh,
     *            flz_o,pwr_o,r1n_o,r2n_o,rlake_o,mndr_o,
     *            aa2_o,aa3_o,aa4_o,fpFactor_o,
     *            theta_o,widep_o,kcond_o,pool_o,
     *            fpetLakeOverride
      logical fpetLakeOverrideFlag
      logical fratioflg

      real*4, dimension(:,:), allocatable :: nrvr_array,class_array

!     rev. 9.2.11  Sep.  11/05  - NK: added Manning's n  r1n & r2n
!     rev. 10.5.15 June  11/23  = NK Added sdcd, r1n, aa2, aa3, widep, pool to Ostin.txt list
      real*4,    dimension(:), allocatable :: r1n,r2n,rlake,fpFactor, 
     *                           r2nlow,r2nhgh,r2ndlt,
     *                           rlakelow,rlakehgh,rlakedlt,
     *                           poollow,poolhgh,pooldlt,
     *                           r1nlow,r1nhgh,r1ndlt,
     *                           aa2low,aa2hgh,aa2dlt,
     *                           aa3low,aa3hgh,aa3dlt,
     *                           fpflow,fpfhgh,fpfdlt
!     rev. 10.2.07 Nov.  03/17  - NK: New rt_pond subroutine for channel pond routing      
      logical,   dimension(:), allocatable :: pondflg
     
!     REV. 10.1.43 Oct   21/16  - NK: lake_ice_facter changed from : to :,:  
      real*4,    dimension(:,:), allocatable :: lake_ice_factor

      real*4,    dimension(:,:), allocatable :: h,fpetmo 

      integer*2,   dimension(:), allocatable :: iiclass

      character*10, dimension(:), allocatable :: nclass,rivtype

!     rev. 10.4.30 Dec.  18/20  = NK Added a13 parameter for rain/snow tmp > base tmp
      real*4    :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,type1

ccccc      END MODULE area4

ccccc      MODULE area5

      real*4,  dimension(:,:), allocatable :: qrel,qinfl,qdwpr,
     *                     lake_stor,lake_outflow,lake_inflow,
     *                     obs_elv,lake_elv,lake_release,
     *                     del_stor,net_lake_inflow,net_lake_outflow
     *                     ,qstream_sum,strloss_sum,
     *                     net_basin_supply


      real*4,    dimension(:), allocatable :: temp_flow_sum,q24,
     *                                        qdwprmd

!     rev. 9.9.20  Jul.  24/14  - NK: Added dead storage for lakes "store_dead"
!     rev. 10.1.46 Nov.  08/16  - NK: Changed B1 - 5 to real*8
      real*8,    dimension(:), allocatable :: b1,b2,b3,b4,b5
      real*4,    dimension(:), allocatable ::   ppsum,
     *                                          yres,xres,
     *                                          yresin,xresin,
     *                                          b6,b7,fpet_lake,
     *                                          store_dead

      integer*4, dimension(:), allocatable :: ires,jres,nnsum,nopti
      integer*4, dimension(:), allocatable :: iresin,jresin,
     *                                          resin_reach

      integer*4  noresv,local,mhto,index,nrel,noresvi,mhtoi,nreli,ktr

	character*12, dimension(:), allocatable :: resname,resnamei
      character(10) :: starttime1,startdate1
      character(8)  :: starttimeXML
      character(10) :: startdateXML

      character*1 :: poliflg


!     rev. 9.5.52  Jan.  20/09  - NK: added reading yyyymmdd_div.pt2 for diversions
      integer*4 :: nodivert,ndiv,nDivOut
      integer*4, dimension(:),   allocatable  :: gridrecieve,gridsource
      integer*4, dimension(:),   allocatable  :: itake,jtake,igive,jgive
      real*4,    dimension(:),   allocatable  :: ysource,xsource
      real*4,    dimension(:),   allocatable  :: yrecieve,xrecieve
	real*4,    dimension(:),   allocatable  :: upstrda
	real*4,    dimension(:,:), allocatable  :: qdivert,qdiv  !qdiv added for diversion output
	character*12, dimension(:), allocatable :: divertname
	logical     , dimension(:), allocatable :: diverflg,divert_inbasin_flg
	logical                                 :: diversion,wrtdiverflg
	
!     rev. 9.9.16  Jun.  06/14  - NK: Added location file for Root R. diversion
	integer  :: divX(4),divY(4),divGridNo(4)
	real*4   :: divlon(4),divlat(4),qtweak
	character*20 :: cjunk(20)

!     for diversion output      
      character(64) :: outputDirectory
      integer,  dimension(:), allocatable :: divRow,divCol,divGrid
      character(16),  dimension(:), allocatable ::  divOutName,divInName
      
	

ccccc      END MODULE area5

ccccc      MODULE area6
!     rev. 10.1.12 Dec.  12/15  - NK: Added Nash Efficiency nasheff.r2c file unit-66
      real*4,    dimension(:,:), allocatable :: sn1,ssmc,api,
     *                                          precadj,basinerr,
     *                                          nasheff

      real*4,    dimension(:,:), allocatable :: sum_precip,last_precip

      integer*4, dimension(:,:), allocatable :: nhyd,nbasin

      real*4,      dimension(:), allocatable :: statnerr,suberr,
     *                                        subarea,statnash 

      integer*4,   dimension(:), allocatable :: nxtbasin,jl,jr,
     *                                        iflowgrid,inbsnflg

      integer(kind=2), dimension(:,:), allocatable :: s

ccccc      END MODULE area6

! TS - MARCH 2000
! Removed p( , ) array decclaration because  it's decclared in area16a
!     rev. 9.1.18  Jun.  03/02  - Added sub-watershed modelling capability inbsnflg




ccccc      MODULE area7
c      integer*4 :: lc,itt,izy,ncoun,icoun,ifirs,ldelt,nsave,ll
c      real*4    :: ys1,yy1,yx1
ccccc      END MODULE area7

ccccc      MODULE area8
      real*4, dimension(:), allocatable :: a,ddelta,checkl,checkh,ssave,
     *                                   les,ba,b,iclosl,iclosh,nsign,
     *                                   odelta
ccccc      END MODULE area8




ccccc      MODULE area9

        real(4) :: optim

        integer*4      :: numa,nstart,nper,kc,maxn,nnn,dds_flag,pre_flag
	  real*4         :: dds_init_value,pre_emption_value
	  logical        :: dds
        REAL*4, DIMENSION(:), ALLOCATABLE :: sta_weight,mse
        real*4, dimension(:), allocatable :: penalty
        real*4, dimension(:), allocatable :: pre_sse,pre_sum,pre_value
!       sse - sum squared error for dds

ccccc      END MODULE area9




ccccc      MODULE area10

	
      real,   dimension(:,:), allocatable :: qhyd,qsyn,qloc,
     *                              delta,frc,frcs,
     *                              qhyd_dly,qsyn_dly,qhyd_mly,qsyn_mly,
     *                              qhyd_mean_evt,qsyn_mean_evt

	logical, dimension(:), allocatable :: flowflag
	integer, dimension(:), allocatable :: us_sta,ds_sta

      real,   dimension(:), allocatable :: tmp,wl,qbar,
     *                     datum,area,qpeakh,qpeaks,dpeakh,
     *                     dpeaks,volsyn,volhyd,
     *                     hydfctr,synfctr,ystr,xstr,
     *                     qhyd_dly_sum,qsyn_dly_sum,
     *                     qhyd_mly_sum,qsyn_mly_sum,
     *                     qhyd_sum,qsyn_sum,
     *                     qhyd_mean,qsyn_mean,
     *                     sum_num,sum_den,
     *                     mean_flow_error
!     rev. 10.2.15 Feb.  05/18  - NK: Added 'results\monthly_peaks'
      real,   dimension(:,:), allocatable :: qp_month_hyd,qp_month_syn
      
	integer, dimension(:), allocatable :: num_obs

c                                 vel,   already in areatr
      
      integer, dimension(:), allocatable :: iy,jx,iys,jxs,nlow,
     *                                              nopt,nopt1



!         these will replace iys & jxs respectively
      real,    dimension(:),allocatable :: xdamg_locn,ydamg_locn
      
!     rev. 9.9.66  Apr.  `3/15  - NK: Added options to write_tb0 str files
	logical                           :: wrtstrflg
      real,    dimension(:),allocatable :: coeff1,coeff2,coeff3,coeff4
      real,    dimension(:),allocatable :: coeff5,value1


!        integer,       dimension(:), allocatable :: iylatdeg,iylatmin,
!     *                                              jxlondeg,jxlonmin


        character(12), dimension(:), allocatable :: gage,damage
        character(80), dimension(:), allocatable :: note

        integer, parameter :: nnote=100
        integer :: ns,nr,mo,itogo,month1,month2,nnnote
        real    :: arain,tm

ccccc      END MODULE area10


!	qhyd( , ) are the observedlows
!	qsyn( , ) are the computed flows at the gage locations
!	qdam( , ) are the computed flows at the damage sites
!	iy, jx are the streamgauge locations
!	these are read from the str file
!	iys, jxs are the damage site locations
!	these are read from the basin\xxxx.str file is s/r strfw
!	frc( , ) are the polinimial coefficients used for the rating 
!	curves: first at the gages and then at the damage sites.
!	gage are the gauge names - up to 12 characters
!     qpeakh( , ) peak instantaneous recorded
!     qpeaks( , ) peak instantaneous computed
!     dpeakh( , ) peak daily mean recorded
!     dpeaks( , ) peak daily mean computed

ccccc      MODULE area11

        real, dimension(:,:), allocatable :: qrgrid

ccccc      END MODULE area11

! qrgrid - flow exchange grid



ccccc      MODULE area12

c        character(31), dimension(:), allocatable :: fln
c        character(30), dimension(:), allocatable :: filename
c        character(30), dimension(:), allocatable :: outfln
        character(256), dimension(:), allocatable :: fln
        character(256), dimension(:), allocatable :: filename
        character(256), dimension(:), allocatable :: outfln
        character(256), dimension(:), allocatable :: event_fln
        character(256)                            :: par_fln
        character(256)                            :: shd_fln
        character(3)                             :: filetype
        character(3)                             :: FLtype(999)
        character(1)                             :: chars(999) 
        integer(4)                               :: filename_length

ccccc      END MODULE area12





ccccc      MODULE AREA14
c        real, dimension(:,:,:), allocatable :: sn
ccccc      END MODULE area14


ccccc      MODULE area15
c        real, dimension(:), allocatable :: bsmc
ccccc      END MODULE area15

! bsmc is to use soil init moisture basin by basin
! and assures compatibility with old files


ccccc      MODULE area16

      real, dimension(:,:,:), allocatable :: distance  !added Feb. 20/09 nk
      real, dimension(:),     allocatable :: dst  
      real, dimension(:),     allocatable :: radius_influence
!             dst replaces d (below) in ragmet & tmp

      real, dimension(:),     allocatable :: smc,smcrag,xsta,ysta,
c     *                                         temp,d,rsum,gsum,ratio,
     *                                         temp,rsum,gsum,ratio,
     *                                         axx,ayy,xsta1,ysta1,
     *                                         sta_elv,                 ! added Sep. 25/09 nk
     *                                         flow_xsta,flow_ysta,      ! for netCDF
     *                                         divOut_xsta,divOut_ysta,  ! for diversion output
     *                                         divIn_xsta,divIn_ysta    ! for diversion outpu
      real, dimension(:,:),   allocatable :: p,psum,pmin,rrain,radp,
     *                                         sumr,cltr,f,f1,sto,g,sum,
     *	                                     smc_class
      real, dimension(:,:,:), allocatable :: w
      integer, dimension(:),  allocatable :: ntogo,ewg,sng,sngmin,
     *                                         ewgmin
      integer                             :: nw,frameDayLast
      real                                :: lapse_precip,lapse_temp
	real             :: rad_influence_par,min_distance_par
	real             :: rad_influence,min_distance
      character(12), dimension(:),allocatable :: gname,gname1
      character(8),  dimension(:),allocatable :: flow_sta_name
	character(12), dimension(:),allocatable :: gname_temp
	real             ::  monthly_temperature_delta(12)
	real             ::  monthly_precipitation_delta(12)
	logical          ::  climate_delta_flg,non_ca_flg
	logical          ::  new_precip_flg,precip_flg
	logical          ::  new_temp_flg,temp_flg
	logical          ::  new_precip_grid_flg,new_temp_grid_flg
	character*80     ::  new_precip_frame_header,precip_frame_header
	character*80     ::  new_temp_frame_header,temp_frame_header

ccccc      END MODULE area16




ccccc      MODULE areacg

! COMMON BLOCK FOR CRAIG & GORDON ISOTOPE FRACTIONATION ROUTINE     

ccccc      END MODULE areacg



c      REAL*4 :: pint(17500,99),deficit(99)



ccccc      MODULE areaet
!     rev. 9.1.80  Mar.  31/05  - NK: added sublimation   (sublim)

        real*4, dimension(:,:), allocatable :: vo,evap,intev,ev,
     *                        pet,x2,x3,intevt,evt,ssumr,v1,rad,
     *                        sublim,sum_sublim,sum_et,sum_pet,pint
c     *                       , strloss_sum,qstream_sum

        real*4,   dimension(:), allocatable :: 
     *                        fpet2,flint,tto,ttomin,ttomax,
     *                        dd_ice,dd_ice_min,dd_ice_sum,dd_thaw,  !added june 2011 nk
     *                        eloss,diff,hu,ffcap,fcap,totint,
     *                        spore,alb,uzsinit,radv,pres,flgtemp,rh,
     *                        sublim_factor,sublim_rate,
     *                        deficit
        real*4, dimension(:),   allocatable :: sinlat,coslat,tanlat
        real*4, dimension(:),   allocatable :: dly_diff
        real*4, dimension(:,:), allocatable :: mean_dly_diff

        real*4    :: lat,flgevp2,tempa1,tempa2,tempa3,uzsid,
     *                    albe,alamb,den,alpha,tton,akt,flgtmp1
c        real*4    ::  so  ! removed for etharh rewrite
     
        integer*4 :: ngauges,ntemps,ittoflg,sublimflg

        logical  diffflg

ccccc      END MODULE areaet

ccccc      MODULE areaindx
c        integer :: n,ii,i,j
ccccc      END MODULE areaindx


ccccc      MODULE areamelt

        real*4, dimension(:,:),allocatable :: excess1,snowf,snowt,raint,
     *                                 ta,deld,dsno,top,bot,water,wlmax

!     rev. 9.8.21  Jun.  18/12  - NK: Assed swe observed date & report
        real*4,  dimension(:,:),allocatable :: snw,snowc,dsn,ttemp,tmx,
     *                         tmn,sca,oldsca,fexcess,snowcmin,wcl,sdcd,
     *                         sdcsca,ati,def,el,course_obs,course_calc

        real*4,    dimension(:), allocatable :: dsnow,tempv,tempvmin,
     *                         tmax,tmin,
     *                         tmin1,tmin2,base,fm,fmn,whcl,snocap,
     *                         tipm,uadj,elev,rho,qtot,robg,rosn,qnet,
     *                         smelt,excess,extra,qrain,qsnow,qrn,
     *                         qsn,glmelt,qe,qh,qn,qp,refrz,fmadj,
     *                         xswe,yswe
      integer*4                                :: nswe
      integer*4,    dimension(:), allocatable  ::  ewg_swe,sng_swe   
      character*12, dimension(:), allocatable  ::  gname_swe        


        real*4    :: conv31,conv32,tlst31,tlst32,conv33,tlst33,daygm,
     *           elvref,
     *	         fmadjust,fmadjustdlt,fmadjustlow,fmadjusthgh,
     *	         fmalow,fmalowdlt,fmalowlow,fmalowhgh,
     *	         fmahigh,fmahighdlt,fmahighlow,fmahighhgh,
     *           gladjust,gladjustdlt,gladjustlow,gladjusthgh
!     rev. 9.6.06  Apr.  18/10  - NK: added glacier adjust for optimization


!     rev. 9.8.01  Jul.  21/11  - NK: added ragmet optimization to dds setup
!     REV. 10.1.32 May   18/16  - NK: Separate radinfl for precip & temperature
        real*4    ::rlapse,rlapsedlt,rlapselow,rlapsehgh,
     *           tlapse,tlapsedlt,tlapselow,tlapsehgh,
     *           radinfl,radinfldlt,radinfllow,radinflhgh,
     *           rainsnowtemp,rainsnowtempdlt,
     *           rainsnowtemplow,rainsnowtemphgh,
     *           pradinfl,pradinfldlt,pradinfllow,pradinflhgh,
     *           tradinfl,tradinfldlt,tradinfllow,tradinflhgh,
     *           smoothdist,smoothdistdlt,smoothdistlow,smoothdisthgh

        integer*4, dimension(:), allocatable :: nsdc,idump

        integer*4 :: idt31,idt32,mltpivot,mflag,errflg

ccccc      END MODULE areamelt

ccccc      MODULE areanash

        real*4,    dimension(:), allocatable :: aa,bb,cc,ashnum,ashden,
     *                        rsquare,qbarobs,qbarcalc,mean_observed

        integer*4, dimension(:), allocatable :: nq,nqc

ccccc      END MODULE areanash



    
ccccc      MODULE areaopts
        real*4, dimension(:), allocatable :: fmdlt,fmlow,fmhgh,fmndlt,
     *                         fmnlow,fmnhgh,uajdlt,uajlow,uajhgh,
     *                         basdlt,baslow,bashgh,
     *                         subdlt,sublow,subhgh,
     *                         sdcddlt,sdcdlow,sdcdhgh
        integer*4,dimension(:), allocatable :: mbsdlt,mbslow,mbshgh
	  integer*4 :: idskip
        real*4 :: a5dlt,a5low,a5hgh,a8dlt,a8low,a8hgh,a9dlt,a9low,a9hgh,
     *          a10dlt,a10low,a10hgh,a11dlt,a11low,a11hgh,
     *          a12dlt,a12low,a12hgh
c        integer*4 :: nrec       >> used in wqnut only - see below
ccccc      END MODULE areaopts

!  snow optimization parameters



ccccc      MODULE areatrc

! COMMON BLOCK FOR TRACER/ISOTOPE ROUTINE     

      INTEGER*4 :: nofer
      INTEGER   :: itrace,    ! 0=Sub-basin   1=Glacier  2=Landcover
!                             3=Rain-on-stream 4=Flow-type   5=Snowmelt
!                             100=Original GW Tracer (NK)  
!                             101=Wetland Tracer
     *             icount
      REAL*4    :: ncrn,nscn,ncpw,nrec,nlec,pscn,pcpw,prec,plec,ndec,
     *             snwwt,isowt,oldiso,oldiso1,oldiso2,oldiso3,oldiso4,
     *             oldiso5,oldiso6,oldiso7,oldiso8,wetinit
      REAL*4    :: pdec,sdep,navr,navs,ndmv,nrmv,pdmv,prmv,nrnc
	REAL*8    :: sqerr
      REAL*4, DIMENSION(:), ALLOCATABLE :: 
     * massin,massout,masstore,ISOdelta,wmassin,wmassout,wmasstore,
     * wISOdelta,tt,vel,disp,coeff
	INTEGER, DIMENSION(:), ALLOCATABLE :: nn
      
!      DATA itrace/5/
	DATA wetinit/0.75/
	DATA icount/0/


!     SUB-BASIN TRACERS = TRACER 0
      REAL*4, DIMENSION(:), ALLOCATABLE :: 
     *isoin1IBN,isoin2IBN,isoout1IBN,isoout2IBN,
     *isoconcIBN,isostore1IBN,isostore2IBN,
     *isosumQ,isosumQfs,
!     rev. 10.5.08 May   03/23  = NK Base flow Index BFI  added
     *isoconcGWsum,qo2Sum

!     GLACIER TRACERS = TRACER 1


!     LANDCOVER TRACERS = TRACER 2
      REAL*4, DIMENSION(:,:,:), ALLOCATABLE :: 
     *isoin1LC,isoin2LC,isoout1LC,isoout2LC,
     *isoconcLC,isostore1LC,isostore2LC,isoBLZS1,isoBLZS2,
     *isoin1BLZS,isoin2BLZS,isoout1BLZS,isoout2BLZS,isoconcBLZS
      REAL*4, DIMENSION(:,:), ALLOCATABLE :: isosum

!     RAIN TRACERS = TRACER 3
      REAL*4, DIMENSION(:,:), ALLOCATABLE :: 
     *isoin1P,isoin2P,isoout1P,isoout2P,
     *isoconcP,isostore1P,isostore2P

!     FLOW-TYPE TRACERS = TRACER 4
      REAL*4, DIMENSION(:,:), ALLOCATABLE :: 
     *isoin1GW,isoin2GW,isoout1GW,isoout2GW,
     *isoconcGW,isostore1GW,isostore2GW,
     *isoin1SW,isoin2SW,isoout1SW,isoout2SW,
     *isoconcSW,isostore1SW,isostore2SW,
     *isoin1IF,isoin2IF,isoout1IF,isoout2IF,
     *isoconcIF,isostore1IF,isostore2IF

ccccc      MODULE area2


!     SNOWMELT TRACERS = TRACER 5
      REAL*4, DIMENSION(:,:), ALLOCATABLE :: 
     *isoin1SWfs,isoin2SWfs,isoout1SWfs,isoout2SWfs,
     *isoconcSWfs,isostore1SWfs,isostore2SWfs,
     *isoin1IFfs,isoin2IFfs,isoout1IFfs,isoout2IFfs,
     *isoconcIFfs,isostore1IFfs,isostore2IFfs,
     *isoin1GWfs,isoin2GWfs,isoout1GWfs,isoout2GWfs,
     *isoconcGWfs,isostore1GWfs,isostore2GWfs,
     *isoLZS1fs,isoLZS2fs,isoin1LZSfs,isoin2LZSfs,isoout1LZSfs,
     *isoout2LZSfs,isoconcLZSfs,isoLZS1,isoLZS2,
     *isoin1LZS,isoin2LZS,isoout1LZS,isoout2LZS,isoconcLZS

!     WETLAND TRACERS = TRACER 101
      REAL*4, DIMENSION(:,:), ALLOCATABLE :: 
     *isoin1QWET,isoin2QWET,isoout1QWET,isoout2QWET,isoconcQWET,
     *isostore1QWET,isostore2QWET

!     LAKE TRACER INFLOWS
      REAL*4, DIMENSION(:), ALLOCATABLE :: isolakeGW,isolakeIF,
     *isolakeSW,isolakeGWst,isolakest

ccccc      END MODULE areatrc


ccccc      MODULE areawet


      real*4, dimension(:), allocatable :: wetwid,chawid,chadep,qswrain,
     *             chflxa,obflxa,totwid,chflfrac,obflfrac,
     *             wstore1,wstore2,wcap,flowxa,
     *             chaxa,satxa,wetxa,qin,hcha1,
     *             hcha2,hwet1,hwet2,qiwet1,hfp,
     *             qiwet2,qowet1,qowet2,qswevp,
     *             kcond,theta,widep,chaarea,wetarea,wsat,wetfrac,
     *             qiwetsum,qowetsum,wstoreinit	   ! added Oct. 15/07 nk	
      
      REAL :: isowetwt,isowetold,oldisow,isooutwet
      REAL :: isowetold1,isowetold2,isowetold3,isowetold4,isowetold5,
     *        isowetold6,oldisow1,oldisow2,oldisow3,oldisow4,oldisow5,
     *        oldisow6,isooutwet1,isooutwet2,isooutwet3,isooutwet4,
     *        isooutwet5,isooutwet6

	real*4, dimension(:,:), allocatable :: isoin1wet,isoin2wet,
     *             isoout1wet,isoout2wet,isoconcwet,isowstore1,
     *             isowstore2
	real*4, dimension(:,:), allocatable :: isoin1SWwet,isoin2SWwet,
     *             isoout1SWwet,isoout2SWwet,isoconcSWwet,isowstore1SW,
     *             isowstore2SW,isoin1IFwet,isoin2IFwet,
     *             isoout1IFwet,isoout2IFwet,isoconcIFwet,isowstore1IF,
     *             isowstore2IF
	real*4, dimension(:,:), allocatable :: isoin1fswet,isoin2fswet,
     *             isoout1fswet,isoout2fswet,isoconcfswet,isowstore1fs,
     *             isowstore2fs,isoin1SWfswet,isoin2SWfswet,
     *      isoout1SWfswet,isoout2SWfswet,isoconcSWfswet,isowstore1SWfs,
     *      isowstore2SWfs,isoin1IFfswet,isoin2IFfswet,
     *      isoout1IFfswet,isoout2IFfswet,isoconcIFfswet,isowstore1IFfs,
     *      isowstore2IFfs



ccccc      END MODULE areawet
  

! TS - Nov 20/03: Added isotope tracer parameters to wetland module
! TS - Jan 30/04: Added more isotope tracer parameters to wetland module
! TS - Apr 25/05: Added more isotope tracer parameters (for snowmelt)

!     rev. 9.8.90  Oct.  30/13  - NK: Added fetch to the shd file 
      real*4,  dimension(:,:), allocatable :: fetch
      real*4,  dimension(:),   allocatable :: WindSpd
      integer, dimension(:),   allocatable :: WindDir
      logical  windflg,windsnwflg

! PARAMETER LIST:
!
! wetwid(n)  - width of wetland coverage to one side of stream channel
! chanwid(n) - width of the stream channel
! chandep(n) - depth of the full stream channel
! wstore1(n) - wetland storage at beginning of time step
! wstore2(n) - wetland storage at end of time step
! wcap(n)    - wetland capacity (maximum storage)
! chflxa(n)  - cross-sectional flow area in the channel
! obflxa(n)  - cross-sectional flow area in the overbank
! flowxa(n)  - cross-sectional flow area of the depth of flow in the channel
! chaxa(n)   - cross-sectional bankfull area of the channel
! satxa(n)   - cross-sectional area of the saturation depth in the wetland
! wetxa(n)   - cross-sectional area of the wetland
! hcha(n)    - height of water in the stream channel
! hwet(n)    - height of water in the wetland
! wetarea(n) - plan area of wetlands in m^2
! channel_area(n) - plan area of channel in m^2
! widep      - width-to-depth ratio (a11 in PAR file)
! chawid     - channel width
! wetwid     - overbank or wetland width
! chflfrac   - channel fraction of total flow width
! obflfrac   - overbank fraction of total flow width (=0 for channel flow only)
! theta      - wetland soil porosity (a9 in PAR file)
! kcond      - soil conductivity (a10 in PAR file)




ccccc      MODULE areawfo


!     note that outwfo is (column,row) & outarray is (row,column)
!     in the wind program, there are 2 input arrays and 2 output arrays
!     the rest have just one of each

      real*4, dimension(:,:),    allocatable:: outwfo
      real*4, dimension(:,:,:),  allocatable:: inarray3
      real*4, dimension(:,:),    allocatable:: outarray1,inarray1
      real*4, dimension(:,:) ,   allocatable:: outarray2,inarray2
      integer,dimension(:,:),    allocatable:: outarrayi
      real*4, dimension(:,:),    allocatable:: outarray,inarray

      real*4, dimension(:),   allocatable::wfo_sum_p,wfo_cum_p 
!     REV. 10.1.34 Jul   05/16  - NK: Added Obs. & Model mean flows to wfo file
      real*4, dimension(:),   allocatable::wfo_sum_qhyd,wfo_sum_qsyn
      real*4, dimension(:),   allocatable::wfo_qsyn 
      real*4, dimension(:,:), allocatable::wfo_qhyd 
      real*4, dimension(:,:), allocatable::bankfull 



      real         :: xorigin,yorigin,xdelta,ydelta,angle
      integer      :: xcount,ycount,deltat
      character(10) ::  starttime,startdate

      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: attname
      CHARACTER(32), DIMENSION(:), ALLOCATABLE :: attunits

      real*4       :: east_boundary,north_boundary
      logical,       dimension(:), allocatable :: indomainflg
      logical                                  :: courseflg

ccccc      END MODULE areawfo
! common for sediment and nutrient subroutines:      

ccccc      module areawq

c      real*4 :: gamma,ro,viskin,grav,a,b
      real*4 :: gamma,ro,viskin,grav,a_wq,b_wq
!       changed a & b to a_wq and b_wq

!     rev. 9.4.02  Apr.  18/07  - NK: moved rf, rffs from areawq to area1
c      real*4, dimension(:,:),  allocatable::
c     *                         rf,rffs


      Real*4 :: d50(17500),spg(17500),
     *erod(17500),ernfl(17500),y(17500,16),kf(16),gc(16),cf(16),
     *yrot(17500),qs(17500,16),hsed(17500,16),ql(17500,16),
     *diam(400,400),spew(400,400),erodi(400,400),
     *nfer(17500),nfa(17500),pfer(17500),pfa(17500),
     *cronrot(17500),croprot(17500),
     *cron(17500,16),crop(17500,16),
     *mnfer(400,400),mnfa(400,400),mpfer(400,400),mpfa(400,400),
     *WI1(17500),WI2(17500),WO1(17500),WO2(17500),yfinal(17500),
     *WI1n(17500),WI2n(17500),WO1n(17500),WO2n(17500),nfinal(17500),
     *WI1p(17500),WI2p(17500),WO1p(17500),WO2p(17500),pfinal(17500),
     *ss1(17500),ss2(17500),SEDHYD(60,8784),SEDSYN(60,8784),
     *ss1n(17500),ss2n(17500),NITHYD(60,8784),NITSYN(60,8784),
     *ss1p(17500),ss2p(17500),PHSHYD(60,8784),PHSSYN(60,8784),
     *sedmss(60,8784),nitmss(60,8784),phsmss(60,8784),
     *psat(17500),sat(400,400)
      

c          these are already in areatr
c      real*4 :: ncrn,nscn,ncpw,nrec,nlec,pscn,pcpw,prec,plec,ndec
c      real*4 :: pdec,sdep,navr,navs,ndmv,nrmv,pdmv,prmv,nrnc
c      integer*4 :: nofer

      real*4 ::    flow_shear,rey,crit_shear,c_val,y_crit,
     *           coef_cover,grf,ero,gro,y_pot

c           phi, already in areatr


!     rev. 9.9.26  Sep.  16/14  - NK: Added precip adjust for forecast & fcstflg
      logical :: fcst_exists,fcst_mode ! declare  in header of 
      integer :: fcst_yr,fcst_hr0,fcst_days2avg
      real*4  :: fcst_snow_adj,fcst_rain_adj




ccccc      end module areawq

!     NOTE:  THIS MODULE CANNOT BE USED IN CONJUNCTION WITH AREA1A
!            DUE TO DUPLICATE ARRAY NAMES WITH DIFFERENT SIZES

!     variables added since this file was setup initially










      end module area_watflood
