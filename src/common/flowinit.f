      SUBROUTINE flowinit()
      
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
!
!     REV. 7.72 - feb. 04/96 -  took flowinit.for from sub.for
!     REV. 8.25 - May. 22/97 -  fixed allocating the basin # in flowinit
!     REV. 8.73 - Mar.  1/98 -  changed mhrd to mhtot in flowinitFstore2(n)=

!     REV. 8.77 - June  1/98 -  added sub-basin error calculation
!
!     REV. 9.00    Mar.  2000 - TS: CONVERTED TO FORTRAN 90 store2(n)=
      
!     REV. 9.03    Nov.  2000 - TS: ADDED WATFLOOD SWAMP ROUTING
!     rev  9.1.03  July  24/01  - added polinomial to reservoir routing
!     rev. 9.1.18  Jun.  03/02  - Added sub-watershed modelling capability
!     rev. 9.1.36  Jan.  28/03  - Fixed wetland init condition in flowinit
!     rev. 9.1.59  Jul.  15/04  - NK: split rerout into two parts: rdresv & rerout
!     rev. 9.2.12  Sep.  15/05  - NK: added EXCEL eqn to flowinit
!     rev. 9.2.15  Sep.  30/05  - NK: Fixed bug for opt in flowinit
!     rev. 9.2.17  Oct.  11/05  - NK: Fixed bug for .str bounds in route
!     rev. 9.2.18  Oct.  27/05  - NK: Fixed bug in flowinit (init spike)
!     rev. 9.2.32  Feb.  10/06  - NK: Added area_check.csv to output
!     rev. 9.3.05  Nov.  13/06  - NK: adder write_flowinit.for to flowinit.for
!     rev. 9.3.10  Jan.  29/07  - NK: routing pars changed to gridded values
!     rev. 9.4.09  Jun.  19/07  - NK: added lake_area as a variable for iso
!     rev. 9.4.10  Jun.  19/07  - NK: adjusted frac for channel water area
!     rev. 9.4.13  Jul.  09/07  - NK: modified lzs to account for lake area (flowinit) 
!     rev. 9.5     Sep.  07/07  - NK: changed wetland/channel routing 
!     rev. 9.5.25  Mar.  20/08  - NK: fixed lake initiation - moved code route -> flowinit
!     rev. 9.5.28  Apr.  15/08  - NK: fixed allocation for inbsnflg in flowinit
!     rev. 9.5.32  Jun.  04/08  - NK: compute reservoir levels
!     rev. 9.5.34  Sep.  17/08  - NK: fixed lake area in flowinit
!     rev. 9.5.35  Sep.  22/08  - NK: moved flow_sta_location to flowinit
!     rev. 9.5.42  Oct.  22/08  - NK: added b7() as the initial lake surface elevation
!     rev. 9.5.55  Feb.  11/09  - NK: Correct R2n for instream lakes
!     rev. 9.5.65  Sep.  23/09  - NK: change class frac to whole basin values
!     rev. 9.5.66  Oct.  06/09  - NK: fixed bug in flowinit for init flows < 1.0
!     rev. 9.5.71  Oct.  12/09  - NK: fixed bug in lst for setting value for nhyd(,)
!     rev. 9.5.76  Oct.  26/09  - NK: fixed basin exclusion for opt if resin present
!     rev. 9.5.81  Jan.  16/09  - NK: allow reservoirs outside watershed in resv file
!     rev. 9.6.01  Mar.  01/10  - NK: rlake parameter added for Manning n correction
!     rev. 9.7.16  Jan.  05/11  - NK: Fixed init flows outside sub-basin
!     rev. 9.8.06  Oct.  18/11  - NK: Added check for `water` class name
!     rev. 9.8.07  Oct.  10/11  - NK: area_check - removed unused stations
!     rev. 9.8.39  Nov.  26/12  - NK: added check for flow stations in lakes
!     rev. 9.8.50  Feb.  27/13  - NK: Initialize store1&2() for zero lake outflow
!     rev. 9.9.22  Jul.  29/14  - NK: Fixed basin no assignment in flowinit.f
!     rev. 9.9.31  Oct.  13/14  - NK: Changed flow initialization RE: zero init flows
!     rev. 9.9.50  Jan.  07/14  - NK: Added zero - initial flow warning
!     rev. 10.1.06 Nov.  19/15  - NK: Added area_check with can_discharge_sites.xyz 
!     rev. 10.2.16 Feb.  14/18  - NK: Added bankfull flow calculation
!!
!     changes made to include c&g model stuff  nk  April. 26/07
!
!***********************************************************************

      USE area_watflood
      use area_debug
      implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      real*4,     dimension(:),    allocatable :: qinit,datemp
      real*4,     dimension(:,:),  allocatable :: qdagrd
      integer*4,  dimension(:),    allocatable :: iset
      character(12), dimension(:), allocatable :: sta_id
      character(80), dimension(:), allocatable :: sta_name,temp_sta_name
      real*4,        dimension(:), allocatable :: sta_lat
      real*4,        dimension(:), allocatable :: sta_long
c      real*4,        dimension(:), allocatable :: sta_area
      real*4,        dimension(:), allocatable :: temp_area
      logical    :: exists
      real*4     :: da1,qda1,qch,obdepth,tdum,trialq,flow_max,convert,
     *              class_sum,awr,try1,xxyyzz
      integer    :: n,i,j,istep2,l,ktt,k,ktemp,jj,nnn1,lll,
     *              iasdf,ios,noread,nnx,inx,jnx,ii,nresv,nnext,
     *              iallocate      
	integer    :: ll,n1,n2,noo,flowtype,ll1,ll2
!     rev. 10.1.06 Nov.  19/15  - NK: Added area_check with can_discharge_sites.xyz 
	logical    :: found_sta,fatal,area_check1,area_check2
      CHARACTER(30) :: msg
      CHARACTER(10) :: time
      CHARACTER(8)  :: cday
      character(2)  :: colum1(1000),colum2(1000),colum3(1000)
      character(256) :: line

      CHARACTER(1)  :: msg1,errflag,firstmsg    !,smok
      CHARACTER(1)  :: col0(13),col1(13),col2(26)
      character(2)  :: colum(1000),col20(26),col21(26)
      LOGICAL       :: resflag,firstpass,messages,msgflg

      data firstmsg/'y'/firstpass/.true./
      data msgflg/.true./
      DATA col1/'b','d','f','h','j','l','n','p','r','t','v','x','z'/
      DATA col0/'c','e','g','i','k','m','o','q','s','u','w','y','a'/
      DATA col2/' ','a','b','c','d','e','f','g','h','i','j','k','l'
     *         ,'m','n','o','p','q','r','s','t','u','v','w','x','y'/

      DATA col20/'a','b','c','d','e','f','g','h','i','j','k','l','m'
     *          ,'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      DATA col21/'a','b','c','d','e','f','g','h','i','j','k','l','m'
     *          ,'n','o','p','q','r','s','t','u','v','w','x','y','z'/

     
     
     
      open(unit=99,file='junk99.txt')
     
     
     
     
     
c      if(nnn.le.0)then
!       pass by when optimizing after first iteration
!       NK - ALLOCATIONS    Jun. 03/02
      if(.not.allocated(sta_id))then
        allocate(sta_id(no),sta_name(no),sta_lat(no),sta_long(no),
     *                 temp_sta_name(no),   
     *         sta_area(no),temp_area(no),area_error(no),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *   'Error with allocation  in flowinit @ 73'

        allocate(datemp(na),qdagrd(ycount,xcount),
     *           iset(na),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *   'Error with allocation  in flowinit @ 78'

        allocate(qinit(noresv),nreach(noresv),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *   'Error with allocation of nreach in flowinit @ 97'
      endif
     
!     REV. 10.1.43 Oct   21/16  - NK: lake_ice_facter changed from : to :,:  
!       to accomodate 12 monthly values if read from a file
      if(.not.allocated(lake_ice_factor))then
        allocate(lake_ice_factor(noresv,12),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *     'Error with allocation of lake_ice_factor in flowinit @ 100'
        do l=1,noresv
          do i=1,12
            lake_ice_factor(l,i)=1.0     ! initialze default
          end do
        end do
      endif
      
!     rev. 9.9.39  Nov.  14/14  - NK: Modifications for watroute
      if(modelflg.ne.'n')then
        if(.not.allocated(p))then
          allocate(p(ycount,xcount),stat=iAllocate)
          if(iAllocate.ne.0) STOP
     *     'Error with allocation of p in flowinit @ 117'
        endif
        if(.not.allocated(ttemp))then
          allocate(ttemp(ycount,xcount),stat=iAllocate)
          if(iAllocate.ne.0) STOP
     *     'Error with allocation of ttemp in flowinit @ 125'
        endif
      endif  

! * * * * * * * * * * *
d      if(iopt.eq.2)print*,'In flowinit @ 87'

!     wetland geometry initialization (moved from spl May 25/05)

!     this section moved from rdshed   28/12/04  nk
!     S(I,J) IS AN ARRAY OF ORDER NUMBER FOR THE WATERSHED IN  
!     PROGRAM SIMPLE

      if(.not.allocated(nclass))then
        allocate(nclass(classcount),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *     'Error with allocation of nclass array in flowinit'
	endif

!     rev. 9.4.09  Jun.  19/07  - NK: added lake_area as a variable for iso
      if(.not.allocated(lake_area))then
        allocate(lake_area(noresv),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *     'Error with allocation of lake_area vector in flowinit'
        do l=1,noresv
          lake_area(l)=0.0
	  end do
	endif

!     rev. 9.5.28  Apr.  15/08  - NK: fixed allocation for in flowinit
!     these are normally done in read_resv_ef but when there is no
!     rel file, these allocations are not made so need to be done here.
!     rev. 9.5.68  Oct.  07/09  - NK: debugged read_resvin_ef.f
!     rev. 10.2.15 Feb.  05/18  - NK: Added 'results\monthly_peaks'
      if(.not.allocated(inbsnflg))then
        allocate(inbsnflg(no+noresvi),delta(ni,no+noresvi),
     *	  volhyd(no+noresvi),volsyn(no+noresvi),
     *	  qpeakh(no+noresvi),dpeakh(no+noresvi),
     *    qpeaks(no+noresvi),dpeaks(no+noresvi),
     *   stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *     'Error with allocation of inbsnflg vector in flowinit'
  	endif

!         TS - ALLOCATION OF AREANASHA ARRAYS
      if(.not.allocated(aa))then
	  allocate(aa(no+noresvi),bb(no+noresvi),cc(no+noresvi),
     *         ashnum(no+noresvi),ashden(no+noresvi),
     *		 rsquare(no+noresvi),nq(no+noresvi),nqc(no+noresvi),
     *         qbarobs(no+noresvi),stat=iAllocate)
	  if(iAllocate.ne.0) STOP
     *		'Error with allocation of areanasha arrays in flowinit'
      endif

	messages=.false.

!     initialize sums for stats calculation
      do l=1,no+noresvi  
        nq(l)=0
        nqc(l)=0
        aa(l)=0.
        bb(l)=0.
        cc(l)=0.
        ashnum(l)=0.0
        ashden(l)=0.0
        qpeakh(l)=0.1e-10
        dpeakh(l)=0.1e-10
        qpeaks(l)=0.0e-10
        dpeaks(l)=0.0e-10
      end do
        
!     what class is the water class?

d     write(63,*)(nclass(i),i=1,classcount)   	  
d	write(63,*)'------',nclass(classcount-1),'------'

!     rev. 9.8.83  Sep.  10/13  - NK: Set classcount=0 for fli.exe program only
      if(program_name(1:3).ne.'fli')then
!     rev. 9.8.06  Oct.  18/11  - NK: Added check for `water` class name
        ii_water=0
        do ii=2,classcount
          if(nclass(ii)(1:5).eq.'water')ii_water=ii
	  end do
	  if(ii_water.eq.0)then
	    print*
	    print*,'The keyword `water` is not found in the par file'
	    print*
	    stop 'Progam aborted in flowinit @ 185'
	  endif
	  if(ii_water.ne.classcount-1)then
	    print*
	    print*,'The water class is in the wrong order in the par file'
	    print*,'ìi_water=',ii_water
	    print*,'classcount-1=',classcount-1
	    print*,'Check the map file also'
	    print*,'It must be the 2nd last class before impervious'
	    print*
	    stop 'program aborted in flowinit @ 193'
	  endif
	else  ! for program fli.exe
!       this added Sept. 10/13 nk	
	  ii_water=classcount-1  ! no checl made as name is not read in
	endif
	if(wetflg.eq.'y')then
	  ii_wetl=classcount-2
	else
	  ii_wetl=0
	endif

      if(iopt99)then
        write(55,*)'In flowinit.for'
        write(55,*)'~~~~~~~~~~~~~~~'
      endif

      write(68,*)'    n,grid_area(n),aclass(n,classcount-1),
     * chaarea(n)    ,chawid(n),    chadep(n),      cap(n),  pool(n)'

      
      convert=al*al*1000.0  ! converts vol in m^3 to mm on the unit grid
      flow_max=0.0

!     Determine how many reaches have been specified
!     usually = no of lakes but can be more.
      maxr=0
      
	do n=1,naa
        if(ireach(n).gt.0)then
          maxr=max0(maxr,ireach(n))
        endif
      end do
d     write(63,*)
d	write(63,*)'highest reach no found (by flowinit)    =',maxr
d     write(63,*)'no of reservoirs found in the rel files =',noresv
d     write(63,*)
	if(maxr.lt.noresv.and.dds_flag.eq.0)then
	  print*
	  print*,'WARNING:'
	  print*,'Reach numbers for one or more reservoirs are'
	  print*,'not specified in the bsnm_shd.r2c file'
	  print*,'please assign reach numbers to grids in reservoirs'
	  print*,'and lakes in the map file & rerun bsn.exe'
	  print*,'Possibly ok if you are running a sub-basin'
d       pause 'Hit enter to continue (flowinit @ 141)'
	endif
	if(maxr.gt.noresv.and.dds_flag.eq.0)then
	  print*,'No of resrvoirs or lakes in the yyyymmdd_rel.tbo file'
	  print*,'are less than specified in the bsnm_shd.r2c file'
	  print*,'please fix the yyyymmdd_rel.tb0 file'
	  print*,'For more info, please run with iopt = 1'
	  print*
	  stop 'Program aborted in flowinit @ 165'
	endif

	do n=1,naa
!     rev. 9.2.15  Sep.  30/05  - NK: Fixed bug for opt in flowinit
!        removed slope(n)=..... and elev(n)=.....  already done in rdshed
!        rl(n)=al*mndr(ii)    !  moved to spl after reading the par & shd files
!     rev. 9.1.60  Jul.  27/04  - NK: reversed definitions for sl1 & sl2 Int. Slope
!        sl2(n)=sqrt(sl1(n))  !moved to read_shed_ef  Jun. 5/07 nk

! * * * TS * * * 
!       CAP IS THE VOLUME OF WATER IN A REACH FOR THE MEAN ANNUAL FLO
!        widep=a11

!     calculate the channel cross sectional area
      if(aa4(n).gt.0.0)then
          chaxa(n)=(aa2(n)+aa3(n)*da(n)**aa4(n))
        else
!         rev. 9.2.12  Sep.  15/05  - NK: added EXCEL eqn to flowinit
!         EXCEL compatible equation. aa4 must be -ve in the par file
          chaxa(n)=10.0**(aa2(n)*alog10(da(n))+aa3(n))
!         had to put a lower bound on channel xsec area to avoid NaN in resume file
!         NK  Oct. 5/05
          chaxa(n)=amax1(1.0,chaxa(n))
        endif

!     rev. 10.1.59 Dec.  18/16  - NK: Fixed missing # channel correction chnl(1-5)
!     rev. 10.2.19 Mar.  13/18  - NK: Fixed array fault read_divert.f
c        print*,'n= ',n
c        print*,'ichnl(n)=',ichnl(n)
        if(ichnl(n).le.0)then
d            print*,'# channels in gid misssing in the shd file'
d            print*,'for grid # ',n,' Set to 1'
            ichnl(n)=1
d            pause 'In flowinit @ 329'
        endif

        r2n(n)=r2n(n)/chnl(ichnl(n))

        ice_fctr(n)=1.0
        

        
        
        
!     rev. 10.5.17 July  06/23  = NK fill & spill
       if(pool(n).le.0.000001)then        
!     rev. 10.1.98 Oct   04/17  - NK: Deal with -ve flows in route
            cap(n)=chaxa(n)*rl(n)    ! vol of channel @ bankfull m**3
            chadep(n)=SQRT(chaxa(n)/widep(n))
            chawid(n)=chaxa(n)/chadep(n)
            chaarea(n)=chawid(n)*rl(n)
        else
!           For fill & spill      
!           for pool > 0.0        
            chaarea(n)=grid_area(n)*aclass(n,classcount-1)
            chawid(n)=chaarea(n)/rl(n)        ! channel width
            chadep(n)=1.0                     ! channel depth -
            cap(n)=chawid(n)*chadep(n)*rl(n)  ! total volume of water at bankfull
            chaxa(n)=chadep(n)*chawid(n)      !  reset this here -bankfull channel cross section area
            widep(n)=chawid(n)/chadep(n)      !  width-depth ratio - now based on total pond area
        endif
        
      write(68,887)n,grid_area(n),aclass(n,classcount-1),
     *                   chaarea(n),chawid(n),chadep(n),cap(n),pool(n)  
887   format(i10,8f15.3)
        
!     rev. 10.1.59 Dec.  18/16  - NK: Fixed missing # channel correction chnl(1-5)

!       bankfull flow:      
!     rev. 10.2.16 Feb.  14/18  - NK: Added bankfull flow calculation
!            note:  slope = sqrt of slope in read_shed
        bnkfll(n)=chaxa(n)**1.67*slope(n)*slope(n)/
     *                            chawid(n)**0.667/r2n(n)
        
          
c          write(717,71700)n,ichnl(n),chaxa(n),slope(n),chawid(n),
c     *          chnl(ichnl(n)),          
c     *          r2n(n),da(n),aa2(n),aa3(n),aa4(n),bnkfll(n)          
c71700     format(2i5,99f15.3)
        
!     rev. 9.4.10  Jun.  19/07  - NK: adjusted frac for channel water area
!       adjust class areas for channel area
!       Fixed Nov. 20/07 nk

        if(aclass(n,classcount-1).lt.chaarea(n)/grid_area(n))then

!         water area in the shd file is less than the area calculated
!         from the channel parameters and the reach length
!         If added to the water class it
!         will make the sum of the class fractions > 0.0
!         so take the area away from the other classes
	    aclass(n,classcount-1)=chaarea(n)/grid_area(n)
!         Adjust the classed for the added channel area:
          class_sum=0.0
	    do ii=1,classcount
            class_sum=class_sum+aclass(n,ii)
	    end do
	    awr=(class_sum-aclass(n,classcount-1))/class_sum
	    do ii=1,classcount
            aclass(n,ii)=aclass(n,ii)*awr
          end do
        endif
        
!       NOTE:
!       if rlake > 1.0 Manning's n will be corrected for ponding
!       if rlake < 1.0e-04  Mannings n will be left but pond routing used
!     rev. 10.2.70 Nov.  11/19  - NK Lowered acceptable rvalue for Manning n correction 1.0 > 0.1
c        if(rlake(n).lt.1.0.and.rlake(n).gt.1.0e-04)then
        if(rlake(n).lt.0.1.and.rlake(n).gt.1.0e-04)then
            print*,'rlake(',ibn(n),') has the value of ',rlake(n)
            print*,'For pond routing,' 
            print*,'rlake should be approx. 1.0e-12 - 1.0e-11'
            print*,'For the old way of adjusting Manning`s n, the '
            print*,'value of rlake should be >1.0'
            print*,'altough a value > 0.1 is accepted'
            stop 'Program aborted in flowinit @ 361'
        endif

!     rev. 10.2.09 Nov.  04/17  - NK: Reinstated old Manning's n correction for legacy files      
!       From before pond routing - keep new program compatible with old method
c        if(rlake(n).ge.1.0)then
        if(rlake(n).ge.0.1)then
!         Special case
!         Water area in the shed file > than default from equations
!         reservoirs have ireach > 0
!         grids with reservoirs are not touched
!         the water area is set by the shd file input  OR
!         could be set as a river class with high roughness  OR
!         We can adjust the manning's n proportional to the water area:
          if(a2.gt.0.0)then
!           when the water area is larger than the default water area
!           width * reach length
!     rev. 9.6.01  Mar.  01/10  - NK: rlake parameter added for Manning n correction
            r2n(n)=amax1(r2n(n),                              ! Manning n correction
     *        r2n(n)*(grid_area(n)*aclass(n,classcount-1))/   ! Manning n correction
     *                                  chaarea(n)*rlake(n))  ! Manning n correction                   
            if(debug_output)
     *         write(51,4950)n,r2n(n),aclass(n,classcount-1)
4950        format('Grid #, Adjusted R2n, water frac:',i5,2f10.3)
          endif
!         modify water surface areas from the shd file:
!         but it can't be less than calculated above
!     rev. 9.9.51  Jan.  13/15  - NK: Added min channel area in flowinit
          chaarea(n)=
     *           amax1(chaarea(n),grid_area(n)*aclass(n,classcount-1)) 
!         channel width from area/length
          chawid(n)=chaarea(n)/rl(n)
!         same channel depth is used as for channels with no area from shed
!         calculated above
          cap(n)=chaarea(n)*chadep(n)
          chaxa(n)=cap(n)/rl(n)        !redefined here
!         chadep(n) is left the same regardless of chaarea
        else
!         pond routing
!                  qo2(n)=rlake(n)*store2(n)**1.75
          store1(n)=(qo1(n)/rlake(n))**0.57143
          store2(n)=store1(n)
          qo2(n)=qo1(n)
        endif

!     rev. 10.2.07 Nov.  03/17  - NK: New rt_pond subroutine for channel pond routing    
!       use pond routing when the water area is  greater
!       than the channel area calculation and we're not in a 
!       lake or wetland grid            
!       Not to be used for grids with wetlands or coded lakes            
        pondflg(n)=.false.
        if(.not.wetland_flag(n).and.
     *            ireach(n).eq.0.and.
     *            pool(n).eq.0.0)then     ! added pool  Jul/23
            if(aclass(n,classcount-1).gt.chaarea(n)/grid_area(n).and.
     *                                        rlake(n).lt.1.0e-4.and.       
     *                                        rlake(n).gt.0.0)then
              pondflg(n)=.true.
            endif
        endif
        

        class_sum=0.0
	  do ii=1,classcount
          class_sum=class_sum+aclass(n,ii)
	  end do

        if(iopt99)then
          if(n.eq.1)write(55,55001)
     *'  n     frac(n)  chawid(n)  rl(n)   chaarea(n)
     * aclass(n,ii),ii=1,classcount)     class_sum   ' 
55001     format(a90)   
          write(55,55500)n,frac(n),chawid(n),rl(n),chaarea(n),
     *            (aclass(n,ii),ii=1,classcount),class_sum
55500     format(i5,4f10.0,99f10.5)
        endif


!     rev. 9.4.10  
! TS: moved from in 'wetflg.eq.y' statement since required w/o wetlands too, Jul 26/06.
        if(wetflg.eq.'y')then
*         theta=a9

c          wetwid(n)=(al*rl(n))*frac(n)*aclass(n,classcount-2)/rl(n)
!         for wetlands, the floodplain width is fixed
          wetwid(n)=grid_area(n)*aclass(n,classcount-2)/rl(n)
          totwid(n)=chawid(n)+wetwid(n)
          chflfrac(n)=chawid(n)/totwid(n)               ! not used anywhere
          obflfrac(n)=wetwid(n)/totwid(n)               ! not used anywhere
          wcap(n)=rl(n)*wetwid(n)*chadep(n)*abs(theta(n))
          wetxa(n)=wcap(n)/rl(n)/abs(theta(n))
!         next 2 added by NK Dec. 15/02
          wetarea(n)=wetwid(n)*rl(n)   ! OR grid_area(n)*aclass(n,classcount-2)
!         added by TAS Mar 28/07 - for isotopes
c          wetfrac(n)=wetarea(n)/(al*al*frac(n))
          wetfrac(n)=wetarea(n)/(grid_area(n)*frac(n))  ! not used anywhere
	  else
          wetwid(n)=0.0
!         for wetlands, the floodplain width is fixed
          totwid(n)=chawid(n)+wetwid(n)
          chflfrac(n)=chawid(n)/totwid(n)
          obflfrac(n)=0.0
          wcap(n)=0.0
          wetxa(n)=0.0
!         next 2 added by NK Dec. 15/02
          wetarea(n)=0.0   ! OR grid_area(n)*aclass(n,classcount-2)
!         added by TAS Mar 28/07 - for isotopes
c          wetfrac(n)=wetarea(n)/(al*al*frac(n))
          wetfrac(n)=0.0  ! not used anywhere

        endif
!       added the wetflg check Jul. 13/05  nk
        if(iopt99.and.wetflg.eq.'y')then
          if(n.eq.1)write(51,51003)
          if(n.eq.1)write(51,51002)
          write(51,51001)n,da(n),cap(n),chaxa(n),chadep(n),chawid(n),
     *        wetwid(n),wcap(n),wetxa(n),aclass(n,classcount-2),sl2(n)
        endif
      end do

      if(iopt99)then         !added Nov. 10/14  nk
	  open(unit=99,file='debug\chaarea.xyz')
        write(99,*)'Water/channel surface area in km**2'
  	  do n=1,naa
          i=yyy(n)
          j=xxx(n)
          write(99,*)(j-0.5)*xdelta+xorigin,(i-0.5)*ydelta+yorigin,
     *	      chaarea(n)/1000000.
        end do
	  close(unit=99)
	endif

d      if(iopt.eq.2)print*,'In flowinit @ 169'

!     FLOW INITIALIZATION SECTION
!     THIS SECTION COMPUTES THE INITIAL FLOW AND STORAGE IN EACH RIVER 
!     SEGMENT USING THE FIRST DOWNSTREAM GAUGE AND PRORATES THE FLOW 
!     ACCORDING TO DRAINAGE AREA.F

!     STEP 1: RESERVOIR RELEASES ARE SUBTRACTED

!     QDAGRD IS QDA IN A GRID ARRAY USED FOR PRINTING IN THIS S/R 
!     AND USED TO FILL IN 

c      do n=1,naa
      do n=1,na
         i=yyy(n)
         j=xxx(n)
         iset(n)=0
         qda(n)=0.0
         qbase(n)=0.0
         datemp(n)=0.0
      end do

	do i=1,ycount
	  do j=1,xcount
	    qdagrd(i,j)=0.0
          nhyd(i,j)=0  
          bankfull(i,j)=-1.0
	  end do
	end do

      istep2=int(al/1000.)**2

      if(numa.eq.0) write(51,6000)

d      if(iopt.eq.2)print*,'In flowinit @ 221'

!     ERROR CHECKING:
!     check to see flow station is in a grid
      if(iopt99)then
         write(55,*)
         write(55,*)'Checking to see if flow statin is in a grid:'
         write(55,*)'If it shows, it`s in a grid'
      endif
      
      errflag='n'
      do l=1,no
        i=iy(l)
        j=jx(l)
        inbsnflg(l)=1    !  assume flow station is in a grid
!        if(iy(l).le.0.or.iy(l).gt.ycount.or.jx(l).le.0.or.jx(l).gt.xmax)then


!     rev. 9.2.17  Oct.  11/05  - NK: Fixed bug for .str bounds in route
        if(iy(l).le.0.or.iy(l).gt.ycount
     *             .or.jx(l).le.0.or.jx(l).gt.xcount)then

c	    print*,iy(l),ycount,jx(l),xcount
!         this can happen if the stream gauge is outside the waterhsed
!         as when subwaterhseds are modelled as separate watersheds
!         added Mar 14/04 nk.          
          n=0
!         nhyd(i,j) is already 0 from above         
        else
          n=s(i,j)
!         Step 1: assign the basin no at the flow statio.          
          nhyd(i,j)=l
          if(iopt99)write(55,1777)l,i,j,s(i,j),nhyd(i,j)
        endif

        if(n.eq.0.and.nnn.le.0.and.dds_flag.ne.1)then   ! print first trial only
!         station is not in a grid  -  don't mess with this
!         write just once
          messages=.true.
          if(iopt99)then
              write(98,1778)l
              write(55,1778)l
          endif
        endif

        if(n.eq.0)then   
!         station is not in a grid  -  don't mess with this
          inbsnflg(l)=0       ! flow station is not in a grid
	    nopt(l)=0           ! can't use it for optimization either 
          errflag='y'
        endif
      end do
      
!     output added Jul. 29/14  NK
c      if(iopt.ge.2)then
c        open(unit=99,file='nhyd_step_1.xyz',status='unknown')
c        do i=1,ycount
c          do j=1,xcount
c            write(99,*)float(j)*xdelta+xorigin-0.5*xdelta,
c     *                 float(i)*ydelta+yorigin-0.5*ydelta,
c     *                 nhyd(i,j)             
c          end do
c        end do
c        close(unit=99,status='keep')
c      endif

      if(iopt99)then
	  if(messages)then
          write(*,1779)
          write(55,1779)
          write(98,1779)
          write(*,*)
          write(55,*)
          write(98,*)
	  endif

        write(55,*)'    grid no ',
     *     '  reach #      ireach       nreach     inbsnflg    maxr'
	endif
      
      do n=1,naa
!       CONVERSIONS &
!       STORE THE GRIDS CONTRIBUTING TO A DWOPER REACH
        if(ireach(n).ne.0)then
	    l=ireach(n)
          nreach(l)=n
          if(iopt99)write(55,*)n,l,
     *                ireach(n),nreach(l),inbsnflg(l),maxr
        endif
	end do

!     calculate the lake areas:
      if(iopt99)then
        write(55,*)
     *'    lake no      grid no   frac water     grid area    lake_area'
        write(53,*)'      l       inbsnflg(l)'
      endif
      do n=1,naa
!       CONVERSIONS &
!       STORE THE GRIDS CONTRIBUTING TO A DWOPER REACH
          if(ireach(n).ne.0)then
	      l=ireach(n)
c	        if(inbsnflg(l).eq.1)then  ! added Dec. 25/10 nk
	        if(l.gt.0)then  ! added Dec. 25/10 nk  Fixed May 15/12
!     rev. 9.5.34  Sep.  17/08  - NK: fixed lake area in flowinit
!     rev. 9.8.20  MAy.  15/12  - NK: fixed lake area in flowinit 9.5.34
	        lake_area(l)=lake_area(l)+aclass(n,ii_water)*grid_area(n)
c             write(55,*)n,frac(ii_water),grid_area(n),lake_area(l)
              if(iopt99)write(55,*)l,n,
     *                aclass(n,ii_water),grid_area(n),lake_area(l)
            endif
          endif
 	end do

      if(iopt99)then
        write(53,*)
        write(53,*)'    lake no  lake_area in m**2'
        do l=1,noresv
          write(53,*)l,lake_area(l)
        end do
        write(53,*)
        
        write(55,*)'   n,   chawid,      rl,    chaarea,
     *	 aclass(1-classcount),     last col=sum_class'
	endif

!     rev. 9.5.81  Jan.  16/09  - NK: allow reservoirs outside watershed in resv file
!     check to see reservoir is in a grid
      errflag='n'
      do l=1,noresv
        i=ires(l)
        j=jres(l)
        inbsnflg(no+l)=1    !  assume reservoir is in a grid

        if(ires(l).le.0.or.ires(l).gt.ycount
     *             .or.jres(l).le.0.or.jres(l).gt.xcount)then
!         this can happen if the reservoir is outside the waterhsed
!         as when subwaterhseds are modelled as separate watersheds
!         added Jan. 18/09 nk.          
          n=0
        else
          n=s(i,j)
!         n will be 0 if not in the watershed.
        endif
 
        if(n.eq.0.and.nnn.le.0.and.dds_flag.ne.1)then   ! print first trial only
!         reservoir is not in a grid  -  don't mess with this
	    print*,'Resv no',l,' is outside basin & ignored'
        endif

        if(n.eq.0)then   
!         reservoir is not in a grid  -  don't mess with this
          inbsnflg(no+l)=0       ! flow station is not in a grid
          errflag='y'
        endif
       
        if(noresvi.gt.0)then
          if(n.eq.0)then   
!           reservoir is not in a grid  -  don't mess with this
            nopti(l)=0
          endif
        endif

c	print*,no,l,i,j,inbsnflg(no+l)
      end do


d      if(iopt.eq.2)print*,'In flowinit @ 258'

!     SET UP THE INITIAL CHANNEL CONDITION FOR THE FIRST HYDROGRAPH:
c      if(iopt99)then
c        write(55,*)'  id    l    i    j    n   ',
c     *         'nlow(l) qda(n) area(l)'
c      endif

!     *****************************************************************
!     STEP 1: FIND THE FIRST +VE FLOW FOR THE FIRST HYDROGRAPH:
!     STEP 1: FIND THE FIRST +VE FLOW FOR THE FIRST HYDROGRAPH:
!     STEP 1: FIND THE FIRST +VE FLOW FOR THE FIRST HYDROGRAPH:
      if(iopt99)then
         write(55,*)
         write(55,*)'STEP 1: '
         write(55,*)'FIND THE FIRST +VE FLOW FOR THE FIRST HYDROGRAPH'
      endif
      do l=1,no
        area(l)=-1.0  ! default  nk Nov. 29/10
        if(inbsnflg(l).eq.1)then	  
          i=iy(l)
          j=jx(l)
          n=s(i,j)


!         for watroute nr needs to be min(nh,nhtot)
!         for watroute use the flag to set this up

          if(iopt99)write(55,*)nr,nl,mhtot

           ktt=nl    ! no of hours of data in the _str.r2c file

!         WHEN THERE IS NO RAINFALL AT ALL, NR = 0
!         THIS MIGHT HAPPEN WHEN THERE IS SNOWCOURSE DATA AND MELT
!         OCCURS WITHOUT RAIN
          ktt=max(kt,ktt)
          flowflag(l)=.false.


          do k=kt,ktt,kt
!            if(qhyd(l,k).ge.0.0)then
            if(qhyd(l,k).gt.0.0)then
!              THE FIRST +VE FLOW IS FOUND
               flowflag(l)=.true.
               if(n.ge.1)then
                 qda(n)=qhyd(l,k)
                     
!     rev. 9.7.24  Apr.  20/11  - NK: Added diverflg to indicate if a diversion is in grid
c                 if(diverflg(n))qda(n)=0.00001
c                 if(diverflg(n))qda(n)=0.00001
                     
                 iset(n)=1   !```````````````````````````````````````
                 datemp(n)=da(n)
               endif
               nlow(l)=k
               ktemp=k
c               if(iopt99)write(55,*)'*',n,qda(n),da(n),datemp(n)
               if(iopt99)write(55,*)'*',n,qda(n),da(n),datemp(n)
               GO TO 15
            endif
          end do
          if(iopt99)write(55,*)flowflag(l),kt,ktt,l,n,i,j

 15       if(.not.flowflag(l))then
!           ALL FLOWS AT THIS STATION ARE -1 (i.e. NO DATA)
!           SO SET AN INIT FLOW = 0.1 AT KTT=KT

!           DONE JAN. 14/00 IN ZURICH: LET INIT FLOW = DA/1000
!           FOR LACK OF ANYTHING BETTER <<<<<<<<<<<<<<<<<<<<<<<<<<<

            qda(n)=0.001*da(n)
            nlow(l)=kt
            ktt=kt

!           revised Feb. 18/02
            ktemp=kt
!     if there is no +ve flow, abort the program and ask user to estimate
!     an initial flow.
!     an initial flow.
            if(resumflg.ne.'y'.and.iopt99)then
!             Only print these warnings when resume file is not used
              if(firstmsg.eq.'y')then
                print*
                print*,' WARNING:'
                print*,' No flows found @ station',l,' for first event' 
                print*,' Need init. flows to initialize routing'
                print*,' and lower zone storage in each grid.'
	          print*,' Default = DA/1000 is assumed'
                print*,' Please edit .str file - insert est. init flow'
                print*,' at each station.        18/02/02 nk'
	          firstmsg='n'
!               stop 'program aborted in flowinit @ 164'
              else
                if(firstpass)print*,' also station',l
	        endif
             endif
	     endif
	     
           if(iopt99.and.l.eq.1)then
             write(55,*)'  id    l    i    j    n ',
     *         'nlow(l)   qda(n)     area(l)'
           endif

		   area(l)=da(n)

           if(iopt99)write(55,5009)id,l,i,j,n,nlow(l),qda(n),area(l)

c!          ASSIGN QDA AT EACH STREAMFLOW STATION
c!          FIND LOWEST PRE-RISE FLOW:
c!          TIS SECTION IF ACTIVATED WILL INITIATE FLOWS AT EACH GAUGE  
c!          THE LOWEST PRE-RISE FLOWS
c           do k=ktemp+kt,ktt,kt
c             if(qhyd(l,k).lt.0.0) GO TO 10
c!            THIS MEANS WE'VE REACHED THE END OF THE FLOW DATA
c             if(qhyd(l,k).le.qda(n))then
c!              THIS MEANS WE'RE ON A RECESSION CURVE - NO GOOD
c!              ESTIMATE OF INITIAL FLOWS IN BASIN IS TOO HIGH
c!              WE'RE LOOKING FOR THE LOWEST FLOW PRECEDING THE RISE
c!              BUT CUT OFF THE SEARCH AT THE END OF THE RAINFALL
c               qda(n)=qhyd(l,k)
c               nlow(l)=k
c             elseif(qhyd(l,k).gt.qda(n).and.qda(n).gt.0.0)then
c!              WE'RE OUT OF THE RECESSION CURVE
c               GO TO 10	
c             endif
c           end do




   10      if(iopt99)then
             write(55,5009)id,l,i,j,n,nlow(l),qda(n),area(l)
          endif
        endif   ! inbsnflg.eq.?
      end do

!     rev. 9.7.16  Jan.  05/11  - NK: Fixed init flows outside sub-basin
!     for cases where sub-basins are created, flows need to be initialized
!     outside the watershed to prevent infinities etc.
c      do n=1,naa
c	  if(da(next(n)).le.0.0.and.next(n).gt.0)then    !?????????????????
c          do i=1,naa
c	      da(next(n))=da(next(n))+da(i)
c	      qda(next(n))=qda(next(n))+qda(i)
c	    end do
c       endif
c	end do


!     rev. 9.5.35  Sep.  22/08  - NK: moved flow_sta_location to flowinit
!     Changed for FEWS to always write the flow station location as FEWS does weird things
      if(FLtype(6).eq.'nc')then
        open(unit=99,file='debug\flow_station_location.xyz',
     *            status='unknown',iostat=ios)
!       rev. 9.5.33  Sep.  12/08  - NK: added column labels for grapher in flow_station_location.xyz
        i=0
        
        do l=1,no
          i=l
          j=1+i/13
          k=j
          if(l.gt.13)i=l-(i/13)*13
          if(i.eq.0)i=13
          if(mod(i,13).eq.0)j=j-1
          write(99,6012)xstr(l),ystr(l),l,flow_sta_name(l),
     *                  col2(j),col1(i),col2(k),col0(i)
        end do
        close(unit=99,status='keep')
      elseif(iopt99)then
        open(unit=99,file='flow_station_location.xyz',
     *            status='unknown',iostat=ios)
!       rev. 9.5.33  Sep.  12/08  - NK: added column labels for grapher in flow_station_location.xyz
        i=0
        do l=1,no
          i=l
          j=1+i/13
          k=j
          if(l.gt.13)i=l-(i/13)*13
          if(i.eq.0)i=13
          if(mod(i,13).eq.0)j=j-1
          write(99,6012)xstr(l),ystr(l),l,gage(l),
     *                  col2(j),col1(i),col2(k),col0(i),area(l)
        end do
        close(unit=99,status='keep')
      endif

      

!     rev. 9.5.35  Sep.  22/08  - NK: moved flow_sta_location to flowinit
!     rev. 9.5.33  Sep.  12/08  - NK: added column labels for grapher in flow_station_location.xyz
      i=0
      open(unit=99,file='temp_junk.txt',status='unknown')
!     only one character for the first 26      
      do l=1,26
        write(99,62000)' ',col20(l)
62000   format(2a1,3i10) 
      end do
!     then 2 characters      
      do l=27,no*3
        if(mod(l,26).ne.0)then
          i=l/26
          write(99,62000)col20(i),col20(mod(l,26)),l,i,mod(l,26)
        else 
          i=l/26-1
          write(99,62000)col20(i),col20(mod(l,26)+26),l,i,mod(l,26)
        endif  
      end do
      close(unit=99,status='keep')
      open(unit=99,file='temp_junk.txt',status='unknown')
      do l=1,no*3
        read(99,62001)colum1(l)
62001   format(a2)
      end do
      close(unit=99,status='delete')
      
      if(iopt99)then
      open(unit=99,file='tracer_column_entrees.txt',status='unknown')
      do i=1,no*3,3   ! the first column is bypassed
        l=(i-1)/3+1
        write(99,62002),xstr(l),ystr(l),l,gage(l),
     *                    colum1(i+1),colum1(i+2),colum1(i+3)
62002   format(2f15.3,i10,5x,a10,3a5)     
      end do
      close(unit=99,status='keep')
      endif  ! iopt99

!     Compare computed to actual drainage areas
      area_check1=.false.
      area_check2=.false.

!     rev. 10.1.06 Nov.  19/15  - NK: Added area_check with can_discharge_sites.xyz 
!     this is a new format based on the WSC fav sites file as generated by HYDAT
!     It needs a bit of reformatting to make it FORTRAN readable - i.e. 
!     station name is one chracter string
!     This file can also be used to automatically load station titles into GRAPHER
!     Get rid of spaces in names and also brackets or other characters
      if(iopt99)then
          
!     rev. 10.3.01 Jan.  05/20  = Use nudge_flaf.xyz file for station area check      
      INQUIRE(FILE='strfw\nudge_flags.xyz',EXIST=exists)
      IF(exists)THEN
        open(unit=99,file='strfw\nudge_flags.xyz',
     *                              status='unknown',iostat=ios)
        if(ios.ne.0)then  ! added Nov. 10.14  nk
          print*
          print*,'Error opening strfw\nudge_flags.xyz'
          print*,'ios=',ios
          stop 'Program aborted in flowinit @ 791'
        endif
c          read(99,*,iostat=ios)  ! header line - ignore
        do l=1,no
          read(99,*,iostat=ios)sta_long(l),sta_lat(l),i,sta_id(l),
     *                        temp_sta_name(l),temp_area(l)
        end do

        do i=1,no    ! nudge file
              do l=1,no! flow file
                  ll1 = LEN_TRIM(sta_id(i))
                  ll2 = LEN_TRIM(gage(l))
                  ll1=min(ll1,ll2)
                  if(sta_id(i)(1:ll1).eq.gage(l)(1:ll1))then
                          sta_area(l)=temp_area(i)
                          sta_name(l)=temp_sta_name(i)
                  endif
              end do
         end do
       
        close(unit=99,status='keep')
        area_check1=.true.
      endif  
      INQUIRE(FILE='strfw\can_discharge_sites.csv',EXIST=exists)

!     Added back in Jan. 14, 2023  NK
      if(exists)then
        IF(exists)THEN
        open(unit=99,file='strfw\can_discharge_sites.csv',
     *                              status='unknown',iostat=ios)
        if(ios.ne.0)then  ! added Nov. 10.14  nk
          print*
          print*,'Error opening strfw\can_discharge_sites.csv'
          print*,'ios=',ios
          stop 'Program aborted in flowinit @ 1007'
        endif
c          read(99,*,iostat=ios)  ! header line - ignore
        do l=1,no
          read(99,*,iostat=ios)sta_long(l),sta_lat(l),i,sta_id(l),
     *                        temp_area(l),temp_sta_name(l)
d          write(*,*,iostat=ios)sta_long(l),sta_lat(l),i,sta_id(l),
d    *                        temp_area(l),temp_sta_name(l)
        end do
        do i=1,no    ! nudge file
              do l=1,no! flow file
                  ll1 = LEN_TRIM(sta_id(i))
                  ll2 = LEN_TRIM(gage(l))
                  ll1=min(ll1,ll2)
                  if(sta_id(i)(1:ll1).eq.gage(l)(1:ll1))then
                          sta_area(l)=temp_area(i)
                          sta_name(l)=temp_sta_name(i)
                  endif
              end do
         end do
       
          close(unit=99,status='keep')
          area_check2=.true.
        endif
      else
        INQUIRE(FILE='strfw\WSC_data\discharge_sites.xyz',EXIST=exists)
        IF(exists)THEN
          print*,'Please note:'
          print*,'strfw\WSC_data\discharge_sites.xyz  no longer used'
          print*,'Replaced by nudge_flags.xyz  that has same info.'
          print*
          print*,'Please copy strfw\wsc_data\nudge_flags.new to'
          print*,'the working directory as nudge_flags.xyz'
          print*
          write(68,68001)
          write(68,68001)
          write(68,68001)
          write(68,68001)
68001     format('Warning: strfw\WSC_data\discharge_sites.xyz  ')
68002     format('Warning: Replaced by nudge_flags.xyz  w/same info.')
68003     format('Warning: Copy strfw\wsc_data\nudge_flags.new to')
68004     format('Warning: the working directory as nudge_flags.xyz')
!          pause 'Hit enter to continue without a check'
        endif
      endif
      endif  !iopt99

      if(iopt99)then
      if(area_check1.or.area_check2)then
        open(unit=99,file='area_check.xyz',status='unknown')
        write(99,99003)'   longitude    latitude    #  station',
     *                 '     station name                          ',     
     *          '                 actual     model    %difference'
99003   format(a38,a38,a53)
        do l=1,no
          If(sta_area(l).gt.0)then
            area_error(l)=(area(l)-sta_area(l))/sta_area(l)*100.
          else
            area_error(l)=99999.0
          endif

          if(area(l).gt.0)then
            if(area_check1)then
              write(99,99001)xstr(l),ystr(l),l,gage(l),
     *          sta_name(l)(1:50),sta_area(l),area(l),area_error(l)
            else
              write(99,99002)xstr(l),ystr(l),l,sta_id(l),
     *          sta_area(l),area(l),area_error(l)
            endif
          endif
        end do
        close(unit=99,status='keep')
      else
!     rev. 9.5.65  Sep.  23/09  - NK: change class frac to whole basin values
!       initialize sta_area
        do l=1,no
	    sta_area(l)=-1.0
	  end do
        open(unit=99,file='area_check.xyz',status='unknown')
	  Write(99,*)'file  strfw\nudge_flags.xyz  not found'
	  write(99,*)'flow station areas can not be compared'
        write(99,*)'  longitude    latitude    #  station   model area'
        do l=1,no
          write(99,99002)xstr(l),ystr(l),l,sta_id(l),area(l)
        end do
        close(unit=99,status='keep')
      end if
      endif   !iopt99

      if(iopt99)then
        iasdf=1
        msg=' init flow at gauges'
        write(55,5551)iasdf,msg
c        do n=1,naa
        do n=1,na
          i=yyy(n)
          j=xxx(n)
          qdagrd(i,j)=qda(n)
	    flow_max=amax1(flow_max,qda(n))
        end do
	  if(flow_max.gt.99.9)then
          do i=ycount,1,-1
            write(55,5554)(qdagrd(i,j),j=1,xcount)
          end do
	  else
          do i=ycount,1,-1
            write(55,5555)(qdagrd(i,j),j=1,xcount)
          end do
	  endif
      endif

d 	if(iopt.eq.2)print*, 'in flowinit at 300'

!     *****************************************************************
!     STEP 2: SUBTRACT RESERVOIR RELEASES FROM D/S GAUGE FLOWS
!     STEP 2: SUBTRACT RESERVOIR RELEASES FROM D/S GAUGE FLOWS
!     STEP 2: SUBTRACT RESERVOIR RELEASES FROM D/S GAUGE FLOWS
      if(iopt99)then
        write(55,*)
        write(55,*)'STEP 2: '
        write(55,*)'SUBTRACT RESERVOIR RELEASES FROM D/S GAUGE FLOWS'
	  write(55,*)'Basin no assignment at each flow station:'
	  write(55,*)
      endif

!     rev. 9.1.59  Jul.  15/04  - NK: split rerout into two parts: rdresv & rerout
!     A whole section that read the .rel file is gone from here.

!     reservoir data and releases nor entered in rdresv
!     INITIALIZE THE INIT RES RELEASE 
      noread=noresv

      if(noread.gt.0)then   !~~~~~~~~~~~~~~~~~~~~~~~~~~start
        do i=1,noread
          qinit(i)=0.0
        end do

!     CHECK FOR 1ST +'VE FLOW AND SAVE IT
d 	  if(iopt.eq.2)print*, 'in flowinit at 360'
        do j=ktr,nrel,ktr
          do i=1,noread
            if(qinit(i).le.0.0.and.qrel(i,j).ge.0.0)then 
              qinit(i)=qrel(i,j)
            endif
          end do
        end do
d 	if(iopt.eq.2)print*, 'in flowinit at 400'

!     SUBTRACT RESERVOIR RELEASES FROM MEASURED FLOWS:
!     LOOP THRU ALL ELEMENTS AND SUBTRACT UPSTREAM RESERVOIR RELEASE
!     FLOW FROM DOWNSTREAM GAUGE FLOWS.

!     subtract the reservoir releases from the recorded flows at gauges
        do k=1,noread
          resflag=.false.
          if(inbsnflg(no+k).eq.1)then
            i=ires(k)
            j=jres(k)		    
            n=s(i,j)
!           keep going until we leave the (sub)watershed
            do while(.not.resflag.and.n.le.naa
     *		        .and.n.gt.0.and.next(n).gt.0)
               if(qda(n).gt.0.0)then	

!               WE'RE AT A GAUGE AND WE'LL SUBTRACT OUT THE RELEASE
!               RELEASE CAN'T BE GREATER THAN THE GAUGE FLOW
!               nothing is taken out if flow = natural

!     rev. 9.5.66  Oct.  06/09  - NK: fixed bug in flowinit for init flows < 1.0
                qda(n)=amax1(qda(n)-qinit(k),0.0)
!               CHECK TO SEE IF WE'VE RUN INTO ANOTHER RESERVOIR
!               WE HAVE TO CHECK THEM ALL
                do mm=1,noread
                  if(yyy(n).eq.ires(mm).and.xxx(n).eq.jres(mm))then
!                    WE'VE FOUND ANOTHER RESERVOIR
                     resflag=.true.
                  endif
                end do
              endif
              n=next(n)
            end do
          endif     !inbsnflg
        end do

      endif   !~~~~~~~~~~~~~~~~~~~~~~~~~~end

      if(iopt99)then
         write(55,6105)
         write(55,6100)  ! [in flowinit]
c         write(55,6102)
c     *        (n,yyy(n),xxx(n),iset(n),da(n),qda(n),qbase(n),n=1,naa)
         write(55,6102)
         do n=1,na
           IF(qda(n).ne.0.0.or.qbase(n).ne.0.0)then
!            find the gauge number in this grid           
             Do l=1,no
c               if(xxx(n).eq.int((xstr(l)-xorigin)/xdelta+1).and.
c     *            yyy(n).eq.int((ystr(l)-yorigin)/ydelta+1))i=l    
               if(xxx(n).eq.jx(l).and.yyy(n).eq.iy(l))i=l    
             end do
	       write(55,6101)
     *        n,yyy(n),xxx(n),iset(n),da(n),qda(n),qbase(n),i
c	       write(55,6101)
	     endif
	   end do
	   write(55,*)'In station order:'
	   write(55,*)'Note stations witn no initial flows!!'
!     rev. 9.9.50  Jan.  07/14  - NK: Added zero - initial flow warning
         do l=1,no
           if(iy(l).ge.1.and.jx(l).ge.1.and.
     *                iy(l).le.ycount.and.jx(l).le.xcount)then    ! added Jan 23/15  NK   
             n=s(iy(l),jx(l))
             if(n.gt.0.and.n.le.naa)then  ! added Jan 23/15  NK
	         write(55,6101)
     *          n,jx(l),iy(l),iset(n),da(n),qda(n),qbase(n),l
               if(qda(n).eq.0.0.and.iopt99)then
                print*,'Warning: zero init. flow at gauge #',l
               endif
             endif
           endif
	   end do
         do i=ycount,1,-1
            write(55,6104)(nhyd(i,j),j=1,xcount)
         end do
         
        iasdf=3
        msg=' gauge flows - reservoir flows'
        write(55,5551)iasdf,msg

c        do n=1,naa
        do n=1,na
          i=yyy(n)
          j=xxx(n)
          qdagrd(i,j)=qda(n)
	    flow_max=amax1(flow_max,qda(n))
        end do
	  if(flow_max.gt.99.9)then
          do i=ycount,1,-1
            write(55,5554)(qdagrd(i,j),j=1,xcount)
          end do
	  else
          do i=ycount,1,-1
            write(55,5555)(qdagrd(i,j),j=1,xcount)
          end do
	  endif
      endif

d 	if(iopt.eq.2)print*, 'in flowinit at 500'

!     *****************************************************************
!     REV. 8.25 - May.  22/97 -  FIXED ALLOCATING THE BASIN # IN FLOWINIT
!     rev. 9.9.22  Jul.  29/14  - NK: Fixed basin no assignment in flowinit.f
!     ASSIGN THE GRID TO A DRAINAGE BASIN FOR CALC OF SUM PRECIP 
!     Step 3
!     WORK UPSTREAM:
!     Assign basin numbers to grids u/s of gauging stations
      if(iopt99)then
      Write(55,*)
      write(55,*)'STEP 3:'
      write(55,*)'Assign u/s grids with basin numbers - d/s gauge no.'
      write(55,*)'Assign init flow to u/s stations - prorated with d.a.'
      endif
      do n=naa,1,-1
!        CURRENT LOCATION:
         i=yyy(n)
         j=xxx(n)
!        LOCATION OF DOWNSTREAM ELEMENT:
         nnx=next(n)
         if(nnx.le.0)then
           if(dds_flag.eq.0)then  ! added Feb. 8/11 nk
	       print*,'WARNING:'
             print*,' grid number',n,' does not have a recieving'
             print*,' within the grid limits'
             print*,' Possible problem: no blank grid around the'
             print*,' watershed'
             print*
c            STOP ' Program aborted in flowinit at line ~339'
           endif
         endif
         if(nnx.gt.0)then
           inx=yyy(nnx)
           jnx=xxx(nnx)
!           WE ARE NOT IN ONE OF THE OUTLET ELEMENTS
            if(nhyd(i,j).le.0)then
!              WE ARE NOT ENTERING ANOTHER SUB-AREA
!              and we are not at a flow station already assigend
               if(da(nnx).gt.0.0)then
!                 WE'RE ABOVE A GAUGE AND WE CAN ALLCOATE THE TRIBUTARY
!                 AREA TO THE DOWNSTREAM GAUGE:
                  nhyd(i,j)=nhyd(inx,jnx)
               endif
            endif
         endif
      end do

!     output added Jul. 29/14  NK
c      if(iopt99)then
c        open(unit=99,file='nhyd_step_2.xyz',status='unknown')
c        do i=1,ycount
c          do j=1,xcount
c            write(99,*)float(j)*xdelta+xorigin-0.5*xdelta,
c     *                 float(i)*ydelta+yorigin-0.5*ydelta,
c     *                 nhyd(i,j)             
c          end do
c        end do
c        close(unit=99,status='keep')
c      endif

      do n=1,no
         nxtbasin(n)=0      !added (n) nk
      end do

d 	if(iopt.eq.2)print*, 'in flowinit at 600'

      if(iopt99)write(55,6044)
      do n=1,naa
!        CURRENT LOCATION:
         i=yyy(n)
         j=xxx(n)
!        LOCATION OF DOWNSTREAM ELEMENT:
         nnx=next(n)
         if(nnx.gt.0)then
            inx=yyy(nnx)
            jnx=xxx(nnx)
!           WE ARE NOT IN ONE OF THE OUTLET ELEMENTS
            if(nhyd(i,j).ne.nhyd(inx,jnx))then
!              WE ARE ENTERING ANOTHER SUB-AREA
               if(nhyd(inx,jnx).gt.0)then 
!                 THIS SUB AREA IS A TRIBUTARY AREA TO THE NEXT
!                 SUB-BASIN :
                  if(nhyd(i,j).le.0)then
                     write(*,6046)n,inx,jnx,nhyd(i,j)
                     write(98,6046)n,inx,jnx

!     rev. 9.04    Jan    16/01 - fixed grid diagnosis in flowinit
!!1111111111111                     STOP 'program stopped in flowinit @65'
                  endif
                  if(nhyd(i,j).gt.0)then
                    nxtbasin(nhyd(i,j))=nhyd(inx,jnx)
                  endif
                  if(iopt99)then
                    if(nhyd(i,j).gt.0)then
                      write(55,6043)n,i,j,nhyd(i,j),nxtbasin(nhyd(i,j))
                    else
                      write(55,6043)n,i,j,nhyd(i,j),-999
                    endif
                 endif
               endif
            endif
         endif
      end do
      if(iopt99)write(55,6045)(l,nxtbasin(l),l=1,no)

d 	if(iopt.eq.2)print*, 'in flowinit at 700'

!     output added Jul. 29/14  NK
      if(iopt99)then
	   write(55,*)' basin number allocations after working upstream'
         do i=ycount,1,-1
            write(55,6104)(nhyd(i,j),j=1,xcount)
         end do
      endif
      
!     CALCULATE BASEFLOW FOR EACH ELEMENT:

!     WORK UPSTREAM 
!     THE FOLLOWING INITIALIZES THE BASE FLOWS FOR EACH GRID POINT
!     ACCORDING TO THE NEAREST DOWNSTREAM GAUGE          
!     RECORDED FLOW AT THE BEGINNING OF THE SIMULATION
!     RESERVOIR RELEASES ARE TAKEN OFF ABOVE

      do n=naa,1,-1
         if(qda(n).le.0.0)then
!           THIS MEANS WE ARE NOT AT A GAUGE (WITH FLOWS)
!           WE'LL ONLY ASSIGN A FLOW IF IT HAS NOT DONE BEFORE
!           THIS KEEPS FROM GOING PAST A GAUGE
            nnx=next(n)
            i=yyy(n)
            j=xxx(n)
            if(nnx.gt.0)then
            inx=yyy(nnx)
            jnx=xxx(nnx)
!              WHEN NNX = 0 WE ARE IN ONE OF THE OUTLET ELEMENTS
               if(da(nnx).gt.0.0)then
!                 WE'RE ABOVE A GAUGE AND WE CAN CALC BASE FLOW
!                 AND ALLOCATE THE TRIBUTARY AREA
!                 FIXED MAR. 23/97
!                 WE HAVE TO ALSO PASS THROUGH THE LAKES WITH 
!                 NATURAL CONTROLS 
!!                  if(ireach(n).eq.0.or.ireach(n).gt.noread)then
!                           taken out Dec. 21/05
                     qda(n)=qda(nnx)*da(n)/da(nnx)
                     datemp(n)=da(n)
!                    IF ISET=1 THEN DON'T TOUCH THIS BASE FLOW AGAIN!
                     iset(n)=1
!!                  endif    ! taken out Dec. 21/05
!                 QBASE IS USED IN RUNOF5 TO INITIALIZE LZS FOR BASEFLOW
                  qbase(n)=qda(n)
!                 JUST TO MAKE SURE:
                  if(qda(n).le.0.0)iset(n)=0
               endif
            endif
         endif
      end do

d      if(iopt.ge.2.and.iopt.le.10)then
d         write(55,6108)
d         write(55,6100)
d         write(55,6102)
d    *       (n,yyy(n),xxx(n),iset(n),da(n),qda(n),qbase(n),n=1,naa)
d      endif

      if(iopt99)then
        iasdf=4
        msg=' distributed base flows'
        write(55,5551)iasdf,msg
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          qdagrd(i,j)=qda(n)
	    flow_max=amax1(flow_max,qda(n))
        end do
	  if(flow_max.gt.99.9)then
          do i=ycount,1,-1
            write(55,5554)(qdagrd(i,j),j=1,xcount)
          end do
	  else
          do i=ycount,1,-1
            write(55,5555)(qdagrd(i,j),j=1,xcount)
          end do
	  endif
      endif

d 	if(iopt.eq.2)print*, 'in flowinit at 800'


!     *****************************************************************
!     STEP 4
!     FILL IN REMAINING GRIDS
!     WORK DOWNSTREAM
      Write(55,*)
      write(55,*)'STEP 4:'
      write(55,*)'Assign d/s grids with basin numbers '
      write(55,*)'Assign init flow to d/s stations - prorated with d.a.'

!     rev. 9.9.22  Jul.  29/14  - NK: Fixed basin no assignment in flowinit.f
      do n=1,naa
!        this added Nov. 21/07  nk
!        to ensure that all grids belong to a basin even if below 
!        a stream flow gauge.
!        CURRENT LOCATION:
!        LOCATION OF DOWNSTREAM ELEMENT:
         nnx=next(n)
!     REV. 10.1.22 Jan.  25/16  - NK: Fixed flowinit for partioa basins
         if(nnx.gt.0)then
!          this does not work for partial watersheds so this is bypassed for this case         
           j=xxx(nnx)
           i=yyy(nnx)
           if(nhyd(i,j).eq.0)then
	       nhyd(i,j)=nhyd(yyy(n),xxx(n))
	     endif
	   endif
	end do
!     output added Oct. 14/14  NK
      if(iopt99)then
	   write(55,*)' basin number allocations after working downstream'
         do i=ycount,1,-1
            write(55,6104)(nhyd(i,j),j=1,xcount)
         end do
      endif
!     output added Jul. 29/14  NK
c      if(iopt99)then
c        open(unit=99,file='nhyd_step_3.xyz',status='unknown')
c        do i=1,ycount
c          do j=1,xcount
c            write(99,*)float(j)*xdelta+xorigin-0.5*xdelta,
c     *                 float(i)*ydelta+yorigin-0.5*ydelta,
c     *                 nhyd(i,j)             
c          end do
c        end do
c        close(unit=99,status='keep')
c      endif
      
!     rev. 9.9.22  Jul.  29/14  - NK: Fixed basin no assignment in flowinit.f
!     WORK UPSTREAM again to fill in below last gauges:
      do n=naa,1,-1
!        CURRENT LOCATION:
         i=yyy(n)
         j=xxx(n)
!        LOCATION OF DOWNSTREAM ELEMENT:
         nnx=next(n)
         if(nnx.le.0)then
           if(dds_flag.eq.0)then  ! added Feb. 8/11 nk
	       print*,'WARNING:'
             print*,' grid number',n,' does not have a recieving'
             print*,' within the grid limits'
             print*,' Possible problem: no blank grid around the'
             print*,' watershed'
             print*
c            STOP ' Program aborted in flowinit at line ~339'
           endif
         endif
         if(nnx.gt.0)then
           inx=yyy(nnx)  ! y location of d/s grid
           jnx=xxx(nnx)  ! x location of d/s grid
!           WE ARE NOT IN ONE OF THE OUTLET grids
            if(nhyd(i,j).eq.0)then
!              WE ARE NOT ENTERING ANOTHER SUB-AREA already assigend
!              and we are not at a flow station already assigend
               if(da(nnx).gt.0.0)then
!                 WE'RE ABOVE A GAUGE AND WE CAN ALLCOATE THE TRIBUTARY
!                 AREA TO an upstream GAUGE:
                  nhyd(i,j)=nhyd(inx,jnx)
               endif
            endif
         endif
      end do
!     output added Oct. 14/14  NK
      if(iopt99)then
	   write(55,*)' basin number allocations after working upstream'
	   write(55,*)' to fill in remaining holes'
         do i=ycount,1,-1
            write(55,6104)(nhyd(i,j),j=1,xcount)
         end do
      endif
!     output added Jul. 29/14  NK
c      if(iopt99)then
c        open(unit=99,file='nhyd_step_4.xyz',status='unknown')
c        do i=1,ycount
c          do j=1,xcount
c            write(99,*)float(j)*xdelta+xorigin-0.5*xdelta,
c     *                 float(i)*ydelta+yorigin-0.5*ydelta,
c     *                 nhyd(i,j)             
c          end do
c        end do
c        close(unit=99,status='keep')
c      endif

!     MAKESURE GAUGES ARE MARKED AND NOT PASSED
      do l=1,no
	  if(inbsnflg(l).eq.1)then
         i=iy(l)
         j=jx(l)
         n=s(i,j)
         iset(n)=1
        endif
      end do

!     ROUTE ALL GAUGE initial FLOWS DOWNSTREAM WHEN THERE'S
!     NO DOWNSTREAM GAUGE


	write(55,*)
	write(55,*)'start of ds flow assignments'

      if(iopt99)write(55,55099)
     *          '      l      i      j     n  next(n)   iset(n) 
     *  datemp        qda      datemp        qda'
55099 format(a90)
      do l=1,no
	  if(inbsnflg(l).eq.1)then
         i=iy(l)
         j=jx(l)
         n=s(i,j)
!        taken out to fix problem oct 16/07 for canag  nk
         datemp(n)=0.0
         da1=da(n)
         qda1=qda(n)
!        START AT EACH GAUGE AND WORK DOWNSTREAM
!        IF THE NEXT ELEMENT HAS A FLOW, 
!        FLOW MAY HAVE BEEN ROUTED DOWN A TRIBUTARY OR
!        THERE IS A GAUGE DOWNSTREAM
         nnx=next(n)
         if(n.eq.0)then
           write(55,*)'Gauge #',l,' is not in the watershed'
           write(55,*)'row & column:',i,j
           write(*,*)'Gauge #',l,' is not in the watershed'
           write(*,*)'row & column:',i,j
         elseif(nnx.eq.0)then
           print*
           print*,'ERROR:'
           print*,'In grid #',n
           print*,'row & column:',i,j
           print*,'next grid',nnx,' is not in the watershed'
           print*,'Please fix the map & shd files'
           print*
           stop 'Program aborted in flowinit @ 1153'           
         else
           msg1='a'
d          if(iopt.ge.2.and.iopt.le.10)then
d            write(55,*)
d    *          '      l    i    j    n next(n) iset(n) datemp',
d    *          '         qda      datemp        qda'
d            write(55,5556)msg1,l,i,j,n,next(n),iset(n),datemp(n),
d    *          qda(n),datemp(nnx),qda(nnx)
d          endif
         endif
!        this loop revised Dec. 22/05  nk
c         do while(n.lt.naa.and.iset(nnx).ne.1)

!        check added Apr. 26/12 nk
c         if(nnx.le.0)then
         if(n.le.0)then
           print*,'Rank of the next downstrean cell =',nnx
           print*,'in cell with rank =',n
           print*,'row',yyy(n),' col',xxx(n)
           print*
           print*,'Probably cause:'
           print*,'There can not be a flow station in a receiving grid'
           print*,'I.e. a flow station in the grid downstream from'
           print*,'a grid you chose as an outlet grid.'
           print*
           stop 'Program aborted in flowinit @ 1127'
         endif
         do while(n.lt.na.and.iset(nnx).ne.1
     *		        .and.n.gt.0.and.next(n).gt.0)    !??????????????????????
!          DOWNSTREAM FLOW HAS NOT BEEN SET BY A GAUGE
!          IT COULD HAVE BEEN SET BY A GAUGE UP ANOTHER TRIBUTARY 
!          IF ISET=1 THEN FLOW HAS BEEN SET BY A DOWNSTREAM GAUGE
!          AND WE DON'T TOUCH IT.
           nnx=next(n)    ! leave in - we're in an n loop!

c	write(55,*)n,qda(n),nnx,qda(nnx)

           if(qda(nnx).le.0.0)then
!            FIRST TIME IN THIS ELEMENT
!            KEEP TRACK OF DRAINAGE AREA CONTRIBUTING TO THE FLOW
             datemp(nnx)=da1
             qda(nnx)=qda(n)     !  *da(nnx)/da(n)
             iset(nnx)=2
!            IF ISET=2 flow can be modified by tributary
             msg1='b'
           else
!            FLOW HAS TO BE COMBINED WITH FLOW FROM ANOTHER TRIBUTARY
!            SUM DRAINAGE AREAS AND FLOWS
             datemp(nnx)=datemp(nnx)+da1
             qda(nnx)=qda(nnx)+qda1
             msg1='c'
             iset(nnx)=2
           endif
           if(iopt99)then
            write(55,5556)msg1,nhyd(i,j),xxx(n),yyy(n),n,next(n),
     *             iset(n),datemp(n),qda(n),datemp(nnx),qda(nnx)
           endif
           n=next(n)
         end do
	 endif   ! inbsn
      end do
      

!     output added Jul. 29/14  NK
      if(iopt99)then

c        MOVED TO LST WITH ADDITIONAL DATA          
!     rev. 10.4.59 aUG.  30/22  = NK Added class ET to basin_no.xyz for iopt>=1 
c        open(unit=99,file='debug\basin_no.xyz',status='unknown')
c        do i=1,ycount
c          do j=1,xcount
c            write(99,*)float(j)*xdelta+xorigin-0.5*xdelta,
c     *               float(i)*ydelta+yorigin-0.5*ydelta,
c     *               nhyd(i,j)             
c          end do
c        end do
c        close(unit=99,status='keep')
      
!     output added Jul. 29/14  NK
        do i=1,ycount
          do j=1,xcount
            outarray(i,j)=nhyd(i,j)             
          end do
        end do
!     rev. 9.9.75  Sep.  11/15  - NK: Added basin_no.r2c output to flowinit.f
        author='flowinit'     
        name='basin number'
        coordsys_temp=coordsys1
!       GreenKenue uses LatLong - code below uses LATLONG
        if(coordsys_temp.eq.'LatLong   ')coordsys_temp='LATLONG   '
        if(coordsys_temp.eq.'Cartesian ')coordsys_temp='CARTESIAN '
        zone_temp=zone1
        datum_temp=datum1
        xorigin_temp=xorigin
        yorigin_temp=yorigin
        xcount_temp=xcount
        ycount_temp=ycount
        xdelta_temp=xdelta
        ydelta_temp=ydelta
        attribute_name='basin_number'
        attribute_units='1' 
        attribute_type='integer'
	    startdate='unknown   '
	    starttime='unknown   '
        source_file_name='unknown'     
!       write the header for gridded withdrawal flows
        fln(99)='results\basin_no.r2c'
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call write_r2c(99,99,0,1,0,0,1)   
        call write_r2c(99,99,0,1,0,1,1)   
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      endif
     
	write(55,*)'end of ds assignments'
	write(55,*)
	write(55,*)

d 	if(iopt.eq.2)print*, 'in flowinit at 801'

!     AFTER WORKING UPSTREAM, STORE THE LOCATION FOR ERROR CORRECTION
!     NBASIN() WILL BE USED IN LST.FOR
!        moved from above Nov. 21/07
      do i=1,ycount
         do j=1,xcount
            nbasin(i,j)=nhyd(i,j)
         end do
      end do
      
!     rev. 9.9.25  Sep.  02/14  - NK: Finally fixed the error when nbasin=0
!     check to see that all cells have been assigned a basin number
!     i.e. a relevant flow gauge.
      IF(iopt99)THEN
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          if(nbasin(i,j).eq.0)then
!     rev. 9.9.25  Sep.  02/14  - NK: Finally fixed the error when nbasin=0
            if(msgflg)then
              print*,'WARNING'
              print*,'     Cell #',n,' at row',i,' col',j
              print*,'has not been assigned a sub-basin #'
              print*,'GW contribution can not be written to'
              print*,'the watflood.wfo file'
              print*,'This can happen if a watershed'
              print*,'does not have a flow station assigned to it'
              msgflg=.false.
            else
              print*,'also Cell #',n,' at row',i,' col',j
            endif
          endif
        end do
      endif

!     prorate flow wrt. drainage areas  nk  Dec. 22/05
      if(iopt99)write(55,*)'Prorating the flows downstream:'
      do n=1,naa
        if(iset(n).eq.2)then
!         IF ISET=2 flow can be modified by tributary
          if(iopt99)write(55,*)n,qda(n),da(n),datemp(n)
C          qda(n)=qda(n)*da(n)/datemp(n)
!     REV. 10.1.22 Jan.  25/16  - NK: Fixed flowinit for partial basins
          if(next(n).gt.0)then
            if(iopt99)then        ! added Feb. 2/18 NK
              fatal=.false.
              if(da(next(n)).le.0.0)then
                write(98,98001)'Warning: Downstream drainage area ='
     *           ,da(next(n)),
     *           ' in grid #, x,y:',n,xxx(n),yyy(n),' ok if outlet grid'
98001            format(a35,i10,a16,3i8,a18)  
                fatal=.true.
              endif
            endif
!     rev. 9.9.31  Oct.  13/14  - NK: Changed flow initialization RE: zero init flows
            if(da(next(n)).gt.0.0)then
c             qda(n)=qda(n)*da(n)/da(next(n))  !<<<<<<<<<<<<<<<<check check check check 
!             fixed Jan. 4/15  nk
              qda(next(n))=qda(n)*da(n)/da(next(n))  !<<<<<<<<<<<<<<<<check check check check 
            else
!             this should only happen at the outlets where da is set = 0.0            
             qda(next(n))=qda(n)
            endif
	      if(iopt99)write(55,*)n,qda(n),da(n),da(next(n))
	      if(iopt99)write(55,*)
            iset(n)=1
          endif
        endif    !     REV. 10.1.22
      end do
      if(fatal.and.iopt99)then
        print*,'OK if these locations are receiving grids'
        print*,'OK if these locations are receiving grids'
        print*
c        print*,'Paused in flowinit @ 1479'
c        pause 'Hit enter to continue'
      endif
      
      if(iopt99)then
         iasdf=5
         msg=' base flows routed downstream'
         write(55,5551)iasdf,msg
      endif
!     FIXED A BUG HERE ON APRIL 16/96
!     THIS DO EXECTUED ONLY FOR iopt99  -  NO GOOD I THINK

      if(iopt99)write(55,5557)
      do n=1,naa
         i=yyy(n)
         j=xxx(n)
         qdagrd(i,j)=qda(n)
         qbase(n)=qda(n)
         if(iopt99)then
            write(55,5556)msg1,nhyd(i,j),xxx(n),yyy(n),n,next(n),
     *                    iset(n),da(n),qda(n),datemp(n)
         endif
      end do
d 	if(iopt.eq.2)print*, 'in flowinit at 802'
      if(iopt99)then
        iasdf=5
        msg=' base flows have been routed downstream'
        write(55,5551)iasdf,msg
	  if(flow_max.gt.99.9)then
          do i=ycount,1,-1
            write(55,5554)(qdagrd(i,j),j=1,xcount)
          end do
	  else
          do i=ycount,1,-1
            write(55,5555)(qdagrd(i,j),j=1,xcount)
          end do
	  endif
      endif
d 	if(iopt.eq.2)print*, 'in flowinit at 802'

!     *****************************************************************
!     NOW WORK BACK UPSTREAM:
!     HAVE TO HIT ALL grids BUT HAVE TO MAKE SURE THAT WE 
!     DON'T CHANGE ANY FLOWS ALREADY ASSIGNED
      write(55,*)
      write(55,*)'STEP 5'
      write(55,*)'Assign flows working u/s'

!     this loop revised Dec. 22/05  nk
      do n=naa,1,-1
        if(iset(n).eq.0)then
!         THIS MEANS PROPER BASE FLOWS HAVEN'T BEEN ASSIGNED YET
          nnx=next(n)
!         GET THE BASE FLOW FROM THE DOWNSTREAM ELEMENT
!         WHICH SHOULD HAVE BEEN DEFINED BY NOW IF THERE IS AT 
!         LEAST ONE GAUGE IN THE WATERSHED
!         FIRST TIME IN THIS ELEMENT
          if(qda(nnx).gt.0.0)then
!     rev. 9.2.18  Oct.  27/05  - NK: Fixed bug in flowinit (init spike)
!           qda(n)=qda(nnx)*da(n)/datemp(nnx)
            qda(n)=qda(nnx)*da(n)/da(nnx)
          else
            qda(n)=0.1     ! ok Jul. 11/02
          endif
          qbase(n)=qda(n)
        endif
      end do

      if(iopt99)then
         iasdf=6
         msg=' distributed base flows'
         write(55,5551)iasdf,msg
         msg1='f'
        write(55,5557) 
        do n=1,naa
            i=yyy(n)
            j=xxx(n)
            write(55,5556)msg1,nhyd(i,j),xxx(n),yyy(n),n,next(n),
     *                    iset(n),da(n),qda(n),datemp(n)
         end do
      endif
d 	if(iopt.eq.2)print*, 'in flowinit at 803'

!     FIX ON APRIL 17/97    
      do n=1,naa
         i=yyy(n)
         j=xxx(n)
         qdagrd(i,j)=qda(n)
         qbase(n)=qda(n)
      end do

      if(iopt99)then 
        iasdf=6
        msg=' distributed base flows'
        write(55,5551)iasdf,msg
	  if(flow_max.gt.99.9)then
          do i=ycount,1,-1
            write(55,5554)(qdagrd(i,j),j=1,xcount)
          end do
	  else
          do i=ycount,1,-1
            write(55,5555)(qdagrd(i,j),j=1,xcount)
          end do
	  endif
      endif
d 	if(iopt.eq.2)print*, 'in flowinit at 804'

!     NICK

!     *****************************************************************
!     STEP 6
!     ADD RELEASES BACK IN:

!     NOW WE HAVE TO ADD THE RELEASES BACK ON FOR THE ROUTING INIT
!     LOOP THRU ALL RESERVOIRS
!     GO DOWNSTREAM AND ADD QINIT TO THE REACHES BELOW THE DAMS
!     TO SET PROPER INITIAL RIVER FLOWS
!     WE'VE GOT TO TRACK THE RIVER TO THE OUTLET
!     OR TO THE NEXT RESERVOIR
      write(55,*)
      write(55,*)'STEP 6'
      write(55,*)'Add the releases back in'

      do k=1,noread
	  if(inbsnflg(no+k).eq.1)then
         i=ires(k)
         j=jres(k)
         n=s(i,j)
         resflag=.false.
         do while(.not.resflag.and.n.le.naa
     *		        .and.n.gt.0.and.next(n).gt.0)     ! ???????????????
            qda(n)=qda(n)+qinit(k)		
            n=next(n)
!           CHECK TO SEE IF WE'VE RUN INTO ANOTHER RESERVOIR
!           WE HAVE TO CHECK THEM ALL
            do mm=1,noresv
               if(yyy(n).eq.ires(mm).and.xxx(n).eq.jres(mm))then
!                 WE'VE FOUND ANOTHER RESERVOIR
                  resflag=.true.
               endif
            end do
         end do
        endif   ! inbsnflg	 
      end do

      if(iopt99) write(55,6109)

d      if(iopt.ge.2.and.iopt.le.10)then
d         write(55,6100)
d         write(55,6102)
d    *        (n,yyy(n),xxx(n),iset(n),da(n),qda(n),qbase(n),n=1,naa)
d      endif
d 	if(iopt.eq.2)print*, 'in flowinit at 805'


!     check it out!!!!!!!!!!!!!!!!!1
!     iopt in the wrong place??????  changed jan 29/09 nk
c      if(iopt99)then



      if(iopt99)then
         iasdf=7
         msg=' reservoir flows added back in'
         write(55,5551)iasdf,msg
	endif

      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        qdagrd(i,j)=qda(n)
      end do

      if(iopt99)then
	  if(flow_max.gt.99.9)then
          do i=ycount,1,-1
            write(55,5554)(qdagrd(i,j),j=1,xcount)
          end do
	  else
          do i=ycount,1,-1
            write(55,5555)(qdagrd(i,j),j=1,xcount)
          end do
	  endif
      endif
d 	if(iopt.eq.2)print*, 'in flowinit at 806'

!     REPLACE THE PROPER QDA VALUES FOR LZS INITIATION IN S/R RUNOF5
      do n=1,naa
         qbase(n)=abs(qbase(n))
      end do

      if(iopt99)then
         write(55,6100)
         write(55,6102)
     *        (n,yyy(n),xxx(n),iset(n),da(n),qda(n),qbase(n),n=1,naa)
      endif

!     FIX ON June 19/06  nk    
      do n=1,naa
         i=yyy(n)
         j=xxx(n)
         qdagrd(i,j)=qda(n)
      end do

      if(iopt99)then 
         iasdf=8
         msg=' final init flows'
         write(55,5551)iasdf,msg
	  if(flow_max.gt.99.9)then
          do i=ycount,1,-1
            write(55,5554)(qdagrd(i,j),j=1,xcount)
          end do
	  else
          do i=ycount,1,-1
            write(55,5555)(qdagrd(i,j),j=1,xcount)
          end do
	  endif
      endif

      if(iopt99)then 
        iasdf=9
        msg=' write nhyd(i,j) '
	  write(55,6110)
        write(55,5551)iasdf,msg
        do i=ycount,1,-1
          write(55,5553)(nhyd(i,j),j=1,xcount)
        end do
      endif
d 	if(iopt.eq.2)print*, 'in flowinit at 808'

!     THIS SECTION CAME FROM RESET.FOR

!     rev  9.1.02  July  12/01  - put in dacheck in flowinit for flag
!         dacheck=al*al/1000000.*1000.
      dacheck=1.0e+10
      if(numa.eq.0)write(51,9801)dacheck
!     assign all flow state variables with initial conditions
      do n=1,naa
         i=yyy(n)
         j=xxx(n)
         ii=ibn(n)
         qi1(n)=qda(n)
         qi2(n)=qda(n)
         qo1(n)=qda(n)
         qo2(n)=qda(n)
!     rev. 9.5     Sep.  07/07  - NK: changed wetland/channel routing 
	   if(wetflg.eq.'y')then
	     qiob1(n)=0.0
           qiob2(n)=0.0
           qoob1(n)=0.0
           qoob2(n)=0.0
         endif

!     rev. 9.1.36  Jan.  28/03  - Fixed wetland init condition in flowinit

!        check for r2 value made in param.for
! * * * * * * * * TS - ADDED WATFLOOD ROUTING INITIALIZATION * * * * * *

         if(manningflg.eq.'y')then
           qch=(cap(n)/rl(n))**1.67*slope(n)/
     *                           chawid(n)**0.667/r2n(n)
         else
           qch=(cap(n)/rl(n))**1.33*slope(n)/r2(n)
         endif
      flowtype=0   
!        WHEN THE SLOPE IS .LE. 0.0 THE ELEMENT IS NOT IN THE BASIN
!        BUT it IS A RECEIVING ELEMENT
         if(slope(n).gt.0.0)then
            if(qda(n).le.qch)then
!             no overbank flow initially
               if(aclass(n,classcount-2).gt.0.0.and.wetflg.eq.'y')then
!                WETLAND FLOW INITIALIZATION
                 over(n)=0.0 
                 if(manningflg.eq.'y')then
                   store1(n)=rl(n)*
     *              (qo2(n)*chawid(n)**0.667*r2n(n)/slope(n))**.60
c     *              (qo2(n)*chawid(n)**0.667*r2n(ii)/slope(n))**.60
c                   if(store1(n).lt.0.0)then
c                     print*,'In flowinit - store2',n,'<0.0'
c                   endif
	           else
c                   store1(n)=rl(n)*(qo2(n)*r2(ii)/slope(n))**.75
                   store1(n)=rl(n)*(qo2(n)*r2(n)/slope(n))**.75
                 endif

                 store2(n)=store1(n)
                 flowxa(n)=store1(n)/rl(n)
	           hcha1(n)=amax1(0.0,flowxa(n)/chawid(n))!  mar 07/09 mk
	           hcha2(n)=hcha1(n)
	           hwet1(n)=hcha1(n)
	           hwet2(n)=hwet1(n)
!                assumes inbank flow only
c	           wstore1(n)=wetwid(n)*rl(n)*hwet1(n)*abs(theta(ii))
	           wstore1(n)=wetwid(n)*rl(n)*hwet1(n)*abs(theta(n))
	           wstore2(n)=wstore1(n)
c                 satxa(n)=wstore1(n)/rl(n)/abs(theta(ii))
                 satxa(n)=wstore1(n)/rl(n)/abs(theta(n))
	         else
!                CHANNEL FLOW INITIALIZATION
!                 if(aclass(n,classcount-2).gt.0.0) 
!     *              PAUSE 'initializing wetland flows as channel flows'
                 over(n)=0.0
                 if(manningflg.eq.'y')then
                   store1(n)=rl(n)*
     *              (qo2(n)*chawid(n)**0.667*r2n(n)/slope(n))**.60
c     *              (qo2(n)*chawid(n)**0.667*r2n(ii)/slope(n))**.60
	           else
c                   store1(n)=rl(n)*(qo2(n)*r2(ii)/slope(n))**.75
                   store1(n)=rl(n)*(qo2(n)*r2(n)/slope(n))**.75
                 endif
                 store2(n)=store1(n)

	           wstore1(n)=0.0
	           wstore2(n)=0.0
	           
	           flowxa(n)=store1(n)/rl(n)
	           hcha1(n)=amax1(0.0,flowxa(n)/chawid(n))!  mar 07/09 mk
	           hcha2(n)=hcha1(n)
	           hwet1(n)=0.0
	           hwet2(n)=0.0
	         endif
            else           ! overbank flow, wetland saturated
	         if(aclass(n,classcount-2).gt.0.0.and.wetflg.eq.'y')then
!                WETLAND+OVERBANK FLOW INITIALIZATION
                 if(manningflg.eq.'y')then
c                   over(n)=((qo2(n)-qch)*r1n(ii)*6.0/slope(n))**.75
                   over(n)=((qo2(n)-qch)*r1n(n)*6.0/slope(n))**.75
!                  the factor of 6.0 is incorportated in r1
                 else
c                   over(n)=((qo2(n)-qch)*r1(ii)/slope(n))**.75
                   over(n)=((qo2(n)-qch)*r1(n)/slope(n))**.75
                 endif
                 store1(n)=rl(n)*(cap(n)/rl(n)+over(n))
                 store2(n)=store1(n)
                 flowxa(n)=store1(n)/rl(n)
                 obdepth=over(n)/(wetwid(n)+chawid(n))
!	           hcha1(n)=flowxa(n)/chawid(n)
!                for initial overbank flow, 
!                all depths are assumed to be bankfull depth nk 28/01/03
	           hcha1(n)=chadep(n)+obdepth
	           hcha2(n)=hcha1(n)
!                make initial wetland depth = to 80% of channel depth
!                to ensure +ve qowet initially
                 hwet1(n)=0.80*chadep(n)
	           hwet2(n)=hwet1(n)
! WHAT IS THE POINT OF THIS ? ANSWER: needed for watbal
!	           wstore1(n)=rl(n)*(wcap(n)/rl(n)+over(n))   nk 28/01/03
	           wstore1(n)=rl(n)*wcap(n)/rl(n)
	           wstore2(n)=wstore1(n)
                 satxa(n)=wstore1(n)/rl(n)/abs(theta(n))
	         else
!                CHANNEL+OVERBANK FLOW INITIALIZATION
!                 if(aclass(n,classcount-2).gt.0.0) 
!     *              PAUSE 'initializing wetland flows as channel flows'
                 if(manningflg.eq.'y')then
c                   over(n)=((qo2(n)-qch)*r1n(ii)*6.0/slope(n))**.75
                   over(n)=((qo2(n)-qch)*r1n(n)*6.0/slope(n))**.75
!                  the factor of 6.0 is incorportated in r1
                 else
c                   over(n)=((qo2(n)-qch)*r1(ii)/slope(n))**.75
                   over(n)=((qo2(n)-qch)*r1(n)/slope(n))**.75
                 endif
                 store1(n)=rl(n)*(cap(n)/rl(n)+over(n))
                 store2(n)=store1(n)
	           wstore1(n)=0.0
	           wstore2(n)=0.0
	           flowxa(n)=store1(n)/rl(n)
	           hcha1(n)=amax1(0.0,flowxa(n)/chawid(n))!  mar 07/09 mk
	           hcha2(n)=hcha1(n)
	           hwet1(n)=0.0
	           hwet2(n)=0.0
	         endif
            endif
         endif

!     
!        IF THE ELEMENT IS A 1ST ORDER ELEMENT, THEN SET THE INFLOW 
!        FROM UPSTREAM = 0.0
         qmax(n)=0.0
c         p(i,j)=0.0            
         do ii=1,classcount
            r(n,ii)=0.0
         end do
         qr(n)=0.

      end do

d 	if(iopt.eq.2)print*, 'in flowinit at 1123'

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     STEP 7
!     assign initial LZ conditions:  qlz, lzs &
!     initial coupled wetland conditions
      write(55,*)
      write(55,*)'STEP 7'
      write(55,*)'Initialize LZ & wetland state variables'


!           QBASE(n) NOW BECOMES THE BASEFLOW CONTRIBUTION FROM EACH 
!           SQUARE.  ITS FIRST USE IS TO DETERMINE INITIAL FLOW IN THE 
!           RIVERS.  NOTE: IT IS MULTIPLIED BY FRACT LATER SO THIS IS 
!           FOR FULL ELEMENTS

!     tdum converts mm/hour to m^3/sec for a grid with frac=1.0
!     step2 is the nominal grid area - i.e. the centre grid area in km**2
      tdum=1000.*step2/3600.

!     This section moved from soilinit   June 13/03
!     Why was it there anyway???????
      do n=1,naa
        if(da(n).gt.0.0)then
          qbase(n)=qbase(n)/da(n)*step2*frac(n)
! * * * * * * *  TS  - added initializations * * * * * * * * 
          if(wetflg.eq.'y')then
	      qiwet1(n)=qbase(n)
	      qiwet2(n)=qbase(n)
	      qowet1(n)=qbase(n)
	      qowet2(n)=qbase(n)
	    endif
        else
          write(98,*)
     *     'Warning: Zero or negative drainage area found in grid # ',n
          write(98,*)
     *     'Warning:To accept this, hit return - but check outcome'
          write(98,*)'Warning:Program may crash'

          qbase(n)=0.0
! * * * * * * *  TS  - added initializations * * * * * * * * 
          if(wetflg.eq.'y')then
	      qiwet1(n)=qbase(n)
	      qiwet2(n)=qbase(n)
	      qowet1(n)=qbase(n)
	      qowet2(n)=qbase(n)
          endif
c          pause 'Paused in flowinit @ 984' 
        endif
          flz2(n)=1.0-(1.0-flz(n))
          pwr2(n)=1.0/pwr(n)       ! used in soilinit 
        if(qbase(n).gt.0.0)then
!         WITHOUT THIS, N IS SET TO 0 !!!!!  FAR OUT & VERY WIERD
c          jj=ibn(n)
c         Copied over from runof6.for (thr=1):  AKB July 11, 2002
          lzs(n)=(qbase(n)/flz2(n)/tdum/frac(n))**pwr2(n)

!     rev. 9.4.13  Jul.  09/07  - NK: modified lzs to account for lake area (flowinit) 
!         There is no lzs under laks & reservoirs so adjust:
          lzs(n)=lzs(n)*(1.0-aclass(n,classcount-1))

!         print*,n,qbase(n),flz2(jj),frac(n),pwr2(jj),lzs(n)
        else
!         print*,n,qbase(n),flz2(jj),frac(n),pwr2(jj)
          lzs(n)=1.0
!         DEFAULT VALUE
        endif
d       if(iopt.eq.2.and.n.eq.1)write(55,*)
d    *  '           n        flz         flz2       pwr         pwr2'  
d	  if(iopt.eq.2)write(55,*)n,flz(n),flz2(n),pwr(n),pwr2(n)
      end do
d 	if(iopt.eq.2)print*, 'in flowinit at 1173'

      if(iopt99)then
        do i=1,ycount
          do j=1,xcount
            qdagrd(i,j)=-1.0
          end do
        end do
        
        iasdf=10
        msg=' initial lzs in mm'
        write(55,5551)iasdf,msg
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          qdagrd(i,j)=lzs(n)
        end do
	  do i=ycount,1,-1
          write(55,5559)(qdagrd(i,j),j=1,xcount)
        end do
      endif

      qi1(na)=0.0
      qi2(na)=0.0
!      m=1
!      jan=1
!      tot1=0.
      qi2(n)=0.0

d      if(iopt.ge.2.and.iopt.le.10)then
d         write(55,6010)
d         write(55,6011)(n,i,j,qda(n),qi1(n),qi2(n),qo1(n),qo2(n),
d    *              store1(n),store2(n),cap(n),over(n),n=1,na)
d      endif

!     initialize reservoir storage
      if(iopt99)then
        write(53,*)'In flowinit'
        write(53,*)'Note: set iopt = 0 to not write this.'
        write(53,*)'Initial reservoir storage & outflow:'
        write(53,*)'     grid#        res#        try1       storage',
     *             '       outflow'
	endif


      if(numa.eq.0)write(51,*)'Flow initializing completed'

d 	if(iopt.eq.2)print*, 'in flowinit at 1223'


!     rev. 9.5.25  Mar.  20/08  - NK: fixed lake initiation - moved code route -> flowinit
!       this bunch of code was moved from route where it was repeated if timestep < 1 hour
!       This one kept me up until 1 am!

        do n=1,naa
          qowetsum(n)=0.0
          qiwetsum(n)=0.0
          wstoreinit(n)=wstore2(n)  
        end do

c        nnn11=na
        convert=al*al*1000.0  ! converts vol in m^3 to mm on the unit grid

        do n=1,naa
          res(n)=0
        end do
        
        
!       Make a table for the next downstream lake        
!           LOOK TO SEE IF THERE IS A RESERVOIR IN THIS GRID:
c            if(yyy(n).eq.ires(l).and.xxx(n).eq.jres(l))then
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
        do n=1,naa
          do l=1,noresv
            lll=next(n)
!           LOOK TO SEE IF THERE IS A RESERVOIR IN THIS GRID:
            if(yyy(n).eq.ires(l).and.xxx(n).eq.jres(l))then
              res(n)=l
              
!             IF THERE IS, SET INITIAL RESERVOIR STORAGE:
!              if(b1(l).gt.0.0.and.b2(l).gt.0.0)then
c              if(b1(l).ne.0.0.or.b2(l).ne.0.0)then
              if(b1(l).ne.0.0.and.b2(l).ne.0.0)then
!               WE HAVE A NATURAL CONTROL & WE NEED INITIAL
!               RESERVOIR STORAGE
!               INITIAL FLOWS ARE CALC IN SUB SO LEAVE THEM ALONE 
c                if(id.eq.1)then
!     rev. 9.1.07  Jan.   3/02  - check that outlet is in a lake
!                 Jan. 03/02  NK
                if(ireach(n).eq.0)then
                  print*,' in flowinit in grid number',n
                  print*,' reach number =',ireach(n)
                  print*,' i.e. lake or reservoir #',l,' outlet' 
                  print*,' is not in a lake'
                  print*,' please make grid part of a lake or  '
                  print*,' move outlet to a grid in a lake'
!                 note that dam locations with releases do not have to be in a lake!! 
                  print*
                  stop 'Program aborted in flowinit @ 170'
                endif
!               rev. 8.99mm Dec. 13/2001-     added check for <= 0 init res flow
!               check for +ve flow for res. initializations
!               probable don't need this if we are using a resume file
                if(qo1(n).le.0.0)then
                  print*,' Initial flow at reservoir',l
                  print*,' is .le. 0.0 and can not be initialized'
                  print*,' please ensure there is a downstream '
                  print*,' flow station with a flow > 0.0 '
                  print*,'We`re in row ',yyy(n),' column ',xxx(n)
                  print*,'grid no ',n,' with init flow =',qo1(n) 
                  print*
 !     rev. 9.8.19  May.  10/11  - NK: Added check on mising init flow for lakes
                  if(routeflg.ne.'q')then
                    pause 'in flowinit @ 1722 - hit enter to abort'
                    stop 'Program aborted in flowinit @ 1748'
                  endif
                endif

               if(iopt99)then
                 write(53,*)
                 write(53,*)'Reservoir or lake no:',l
                 write(53,53999)'b`s',
     *                    b1(l),b2(l),b3(l),b4(l),b5(l),b6(l),b7(l)
53999            format(a5,5g12.3,2f12.3)    
                endif

                if(b1(l).gt.0.0)then
!                 we have a rating curve and calculate initial storage                
                  if(b3(l).eq.0.0)then
!                   use 2 coefficients  
c                    if(resumflg.ne.'y')then
c                      store1(n)=(qo1(lll)/b1(l))**(1.0/b2(l))
c                      if(b3(l).eq.0.0)then                      
c                        if(b6(l).le.b7(l))then
                         if(iopt99)write(53,*)l,'power function'
                          store1(n)=(qo2(n)/b1(l))**(1/b2(l))
                          store2(n)=store1(n)
                          try1=qo1(n)
                          if(iopt99)
     *                      write(53,*)n,l,try1,qo1(n),store1(n),'done'
                          
c                        else
c!     rev. 9.5.42  Oct.  22/08  - NK: added b7() as the initial lake surface elevation
c                          write(53,*)l,'initial lake level given'
c                          store1(n)=(b6(l)-b7(l))*lake_area(l)
c                          store2(n)=store1(n)
c                          write(53,*)n,l,b6(l),b7(l),
c     *				         	store1(n)/lake_area(l)+b7(l)
c                        endif
c                      endif
c                    endif
                  else
                    if(iopt99)write(53,*)l,'bisection method'
!                   use bisection to get init flow
!                   actually, I made the int. storage a little larger
!                   needs to be fixed
c                    if(resumflg.ne.'y')then
                      store1(n)=1000.
                      try1=0.0
!                     use 2-5 coefficients
                      if(iopt99)write(53,*)'         n           l',    
     *                     '          try1         qo1         store1'
                      do while(try1.lt.qo1(lll))
!                       keep doubling the res. storage until the corresponding
!                       flow is larger than the initialized streamflow.
                        store1(n)=2.0*store1(n)                      
!     rev. 9.8.89  Oct.  27/13  - NK: Fixed undefined (NAN) problem in flowint
c                        if(b5(l).gt.0.00000)then
                        if(b5(l).ne.0.00000)then
                          try1=b1(l)*store1(n)+
     *                         b2(l)*store1(n)**2+
     *                         b3(l)*store1(n)**3+
     *                         b4(l)*store1(n)**4+
     *                         b5(l)*store1(n)**5
c                        elseif(b5(l).le.0.0000.and.b4(l).gt.0.00000)then
                        elseif(b5(l).eq.0.0000.and.b4(l).ne.0.00000)then
                          try1=b1(l)*store1(n)+
     *                         b2(l)*store1(n)**2+
     *                         b3(l)*store1(n)**3+
     *                         b4(l)*store1(n)**4
                        else
                          try1=b1(l)*store1(n)+
     *                         b2(l)*store1(n)**2+
     *                         b3(l)*store1(n)**3
                        endif
d      write(53,*)b1(l)*store1(n),b2(l)*store1(n)**2,b3(l)*store1(n)**3
d      write(53,*),b4(l)*store1(n)**4,b5(l)*store1(n)**5
                        if(iopt99)
     *                     write(53,*)n,l,try1,qo1(n),store1(n),'poli'
                        if(try1.lt.0.0)then
                          print*,' trial value for flow =',try1
                          print*,' trial reservoir outflow < 0.0'
                          print*,' polinimial functins not'
                          print*,' monotonically increasing'
                          print*,' for reservoir',l
                          print*,' Please fix the function'
                          print*,' Program may excecute but results' 
                          print*,' are approximate'
                          print*,'                   @ 167/flowinit '
                          print*
                          if(iopt99)
     *                        pause 'Hit enter to continue. '
                        endif
                      end do
c                    endif
                  endif

                  if(iopt99)write(53,*)n,l,try1,qo1(n),store1(n),'done'
                  if(iopt99)write(53,*)
                else  ! b1(l)=0.0
!     rev. 9.8.50  Feb.  27/13  - NK: Initialize store1&2() for zero lake outflow
!                 use a default storage - no outflow anyways             
                  store1(n)=100.0
                  store2(n)=100.0

                endif   !if(b1(l).ne.0.0.and.b2(l).ne.0.0)then
d 	if(iopt.eq.2)print*, 'in flowinit at 1224'

                if(resumflg.ne.'y')then
                  store1(n)=max(100.0,store1(n))
                  store2(n)=store1(n)
                  if(iopt99)then
                    write(51,8492)
                    write(51,8493)n,l,b1(l),b2(l),qda(n),store1(n)
                    write(53,8492)
                    write(53,8493)n,l,b1(l),b2(l),qda(n),store1(n)
                  endif
                endif
              endif     !if(b1(l).ne.0.0.and.b2(l).ne.0.0)then

!     rev. 9.9.34  Oct.  17/14  - NK: Added re-compute of lake storage re: new lake levels
!             this section duplicated in sub with this revision #
	        if(b6(l).ne.0.0.or.b7(l).ne.0.0)then
!               initial lake levels can be lower than the datum	        
c	          if(b7(l).lt.b6(l))then
                if(iopt99)then
                  write(53,*)
                  write(53,*)'Override flow based initial storage:'
                  write(53,*)'Override flow based initial storage:'
                  write(53,*)'Override flow based initial storage:'
                endif
!               Note: these are over written again in rerout for special cases
!               Note: these are over written again in rerout for special cases
!               Note: these are over written again in rerout for special cases
                
!     rev. 10.4.48 Jan.  09/22  = NK fixec init lake level for continuing flow or based on flow 
                if(b6(l).gt.0.0)then
!                   start off with lake  init level from the *_ill.pt2
                    store1(n)=(b6(l)-b7(l))*lake_area(l)
                    store2(n)=store1(n)
                else
!                   calculate the init lake level from the initial flow data
!                   For this, set the init eld in the pt2 file  -999.000
                    b6(l)=store1(n)/lake_area(l)+b7(l)
                endif
                
               if(iopt99)write(53,*)
     *            '*****',b6(l),b7(l),lake_area(l),store1(n)
                
                if(store2(n).ge.0.0)then
                  if(b3(l).eq.0.0)then
!     REV. 10.1.21 Jan.  23/16  - NK: Fixed lake init flow bug in flowinit
c                    qo2(n)=b1(l)*store2(n)**b2(l)
                    qo1(n)=b1(l)*store2(n)**b2(l)
                  else
                    qo1(n)=store2(n)*(b1(l)+store2(n)*(b2(l)+store2(n)*
     *                (b3(l)+store2(n)*(b4(l)+b5(l)*store2(n)))))
                  endif   
                  qo2(n)=qo1(n)
                else
                  qo2(n)=0.0
                endif
c              else            !if(b1(l).ne.0.0.or.b2(l).ne.0.0)
c                write(53,*)'In a reservoir but nothing initialized'
c                write(53,*)'Not good really!' 
              endif

!     rev. 9.9.46  Dec.  10/14  - NK: Added check on initial lake outflow
              if(qo2(n)/da(n).gt.0.10)then
                write(98,98002)
     *           'Warning: Initial calculated outflow for lake',l,
     *           ' is ',qo2(n),' CMS which seems high ',
     *           'Check  initial lake level or rating curve'
98002         format(a44,i5,a4,f10.0,a22,a41)          
              endif
d 	if(iopt.eq.2)print*, 'in flowinit at 1225'
                 
!     rev. 9.9.20  Jul.  24/14  - NK: Added dead storage for lakes "store_dead"
!     rev. 9.9.38  Nov.  12/14  - NK: Added LKdepth to ill file
c              if(b1(l).eq.0.0)then
c!               Add dead storage:
!             -----------------                
!               here dead storage is added to the live storage     
!               Use only for reservoir with releases
!               b6-b7 = current live storage
!               lake depth is the total lake depth used in the lake evap routine
!               store1 & store2 are TOTAL storage of the lake
                store_dead(l)=((b7(l)-LKinvert(l)))*lake_area(l)
                store_dead(l)=amax1(0.0,store_dead(l))
                LKdepth(l)=b6(l)-LKinvert(l)
                store1(n)=store1(n)+store_dead(l)
                store2(n)=store1(n)
c              endif
              if(iopt99)then              
              write(53,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
              write(53,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
              write(53,*)'WARNING: these will be overwritten with'
              write(53,*)'resumeflg = `y`'              
              write(53,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
              write(53,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
              write(53,*)'Lake (reach) #       ',l
              write(53,*)'Lake area m^2        ',lake_area(l)
              write(53,*)'Lake invert m        ',LKinvert(l)
              write(53,*)'Datum         m      ',b7(l)
              write(53,*)'Initial Elv.  m      ',b6(l)
              write(53,*)'Initial depth m      ',store2(n)/lake_area(l)
              write(53,*)'Initial storage m^3  ',store2(n)
!     rev. 9.9.74  Sep.  11/15  - NK: Added output to unit 53 in flowinit
              write(53,*)'Dead storage m^3     ',store_dead(l)
              write(53,*)'Live storage m^3     ',store2(n)-store_dead(l)
              write(53,*)'Lake depth m         ',LKdepth(l)
              write(53,*)'Initial outflow m^3/s',qo2(n)
              write(53,*)'WARNING: these will be overwritten with'
              write(53,*)'resumeflg = `y`'              
              write(53,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
              endif
            endif
         end do
      end do  
        
      do n=1,naa  
        if(res(n).ne.0)then
        if(rlake(n).lt.1.0.and.rlake(n).gt.0.1)then
            write(98,98003)
     *         'Warning: rlake(',ibn(n),') has the value of ',rlake(n),
     *                'Warning: Manning`s n will be decreased'
98003       format(a20,i5,a15,f10.3,a30)            
        endif
        endif
      end do
                        
!     rev. 10.2.20 Apr.  06/18  - NK: Added res_next for NBS
!     This section needs to be here because res(n) needs to be defines as above
      if(iopt99)then              
         write(53,*)
         write(53,*)'      Derived in FlowInit @ 2515'
         write(53,*)
         write(53,*)'      Lake #           Next Lake #'
      endif
!        initialize values         
         do l=1,noresv
             res_next(l)=0
         end do
         do l=1,noresv
              i=ires(l)
              j=jres(l)
!             Check that we are in the basin - needed when dealing with
!             sub-watershed models              
              if(i.le.ycount.and.j.le.xcount.and.
     *                           i.gt.0.and.j.gt.0)then
                  n=s(i,j)
              else
                  n=0
              endif
              if(n.gt.0)then
                  nnext=next(n)
              else
                  nnext=0
              endif
              if(nnext.gt.0)then
              nresv=res(nnext)
!              print*,'a',l,n,nnext,nresv
              do while(nresv.eq.0)
                  if(nresv.le.0)then
                      nnext=next(nnext)
                      if(next(nnext).gt.0)then
!                         We've found the next lake
                          nresv=res(nnext)
                          res_next(l)=res(nnext)
!                          print*,'c',nnext,next(nnext),nresv
                          if(nresv.gt.0.and.iopt99)write(53,53000)
     *                         resname(l)(1:12),l, 
     *                         res_next(l),resname(res_next(l))(1:12)
53000                     format(a12,2i5,5x,a12)                          
                      else
!                         Haven't found next lake so keep looking
                          nresv=-1
                  endif
                  endif
                  end do
              endif
              end do    

         if(iopt99)then
           write(63,*)
           write(63,*)'NEW  Nov. 2014<<<<<'
           write(63,*)'dead storage has been added to the live storage'
           write(63,*)'previously, only live storage was used'
           write(63,*)
           write(63,*)'Lake/reservoir storage initialization completed'
           write(63,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         endif    
 8492 format(' ','initialize resvr flow & storage - in flowinit'/
     *            'n,l,b1,b2,qda,store1')
 8493 format(' in route/flowinit:',2i5,4e12.3)
     
d 	if(iopt.eq.2)print*, 'in flowinit at 1226'
      
!     rev. 10.2.07 Nov.  03/17  - NK: New rt_pond subroutine for channel pond routing   
      if(iopt99)write(51,*)'Debug grid # ',nnprint
      if(iopt99)write(51,4951)
4951  format(
     *    '          #          n   waterArea(frac)   gridArea   ',
     *    ' channelArea     cA/gA     width     depth ',
     *    ' iak ireach wetflg pondflg  ')      
      do n=1,naa
          if(pondflg(n))then
             store1(n)=(qo2(n)/rlake(n))**(1.0/1.75)
c             store1(n)=(qo2(n)/1.0e-12)**(1.0/1.75)
             store2(n)=store1(n)
          endif
      end do
      if(iopt99)then
          do n=1,naa
              if(nnprint==n)write(51,*)
              write(51,4952)n,r2n(n),aclass(n,classcount-1),
     *                   grid_area(n),chaarea(n),
     *                   chaarea(n)/grid_area(n),
     *                   chawid(n),chadep(n),
     *                   ibn(n),ireach(n),
     *                   wetland_flag(n),pondflg(n)
              if(nnprint==n)write(51,*)
4952          format('Grid # :',i5,2f12.3,2f15.3,3f10.3,2i5,5x,l1,5x,l1)
          end do
      endif
!     rev. 9.7.01  Jun.  09/10  - NK: fixed error.xyz & error.r2s

!     create flow station heirarchy
d	write(63,*)'**********************'
d	write(63,*)'**********************'
d     write(63,*)'In flowinit:'     
d	write(63,*)'     flow station    downstream flow station'

c      if(resumflg.ne.'y')then    !added mar17/11 nk
      do l=1,no
        if(inbsnflg(l).eq.1)then    ! added nk Nov. 29/10
	    ds_sta(l)=-1
c	    print*
c	    print*
!        what grid no?
c          print*,'l=',l
          i=iy(l)
          j=jx(l)
	    n=s(i,j)
c	    print*,'i,j,n=',i,j,n,next(s(i,j))
!         find ds station - start in the next grid
          n2=n
          found_sta=.false.
          do while(.not.found_sta.and.n2.lt.naa
     *		        .and.n2.gt.0.and.next(n).gt.0)         ! ???????????????
	      n2=next(n2)    ! cell # we are in
!           see if there's a staton in this grid
!           loop thru stations
	      ll=0
c	      print*,'no,n2=',no,n2
	      do while(ll.lt.no.and.n2.gt.0)  ! ??????????????????
              ll=ll+1
	        if(l.ne.ll)then
                if(inbsnflg(ll).eq.1)then    ! added nk Jan. 05/11
!                 skip current station
!                 find sta in this grid & get its number
	            ii=iy(ll)
                  jj=jx(ll)
	            nnn1=s(ii,jj)  
c	            print*,n2,ll,ii,jj,nnn1
	            if(nnn1.eq.n2)then  
!                   found a station
	              ds_sta(l)=ll
	              us_sta(ds_sta(l))=l
d	              write(63,*)'sta',l,'ds_sta(l)=',ll,flowflag(l)
	              found_sta=.true.
	            endif
	          endif
	        endif
	      end do
	    end do
	  endif  ! inbsnflg = 1
	end do
c	endif    ! resumeflg
c	do l=1,no
c	  print*,l,ds_sta(l),us_sta(l)
c	end do
!     end
!     rev. 9.7.01  Jun.  09/10  - NK: fixed error.xyz & error.r2s

c      if(iopt99)then
c        write(55,*)
c        write(55,*)'nhyd - basin number'
c        write(55,*)
c	  do i=ycount,1,-1
c	    write(55,55999)(nhyd(i,j),j=1,xcount)
c	  end do
55999   format(999i3)
c      endif
d 	if(iopt.eq.2)print*, 'in flowinit at 1227'

      if(numa.ge.1.and.nnn1.eq.0)then
	    write(55,*)'   station #    nopt'
        do l=1,no
	    write(55,*)l,nopt(l)
	  end do

        write(55,*)
        write(55,*)'nhyd - grids used for opt'
        write(55,*)
	  do i=ycount,1,-1
	    write(55,55999)(nhyd(i,j),j=1,xcount)
	  end do
	endif
      
!     rev. 9.8.39  Nov.  26/12  - NK: added check for flow stations in lakes
!     Check that no flow stations are in lakse or reservoirs
!     and only in outlets if they are.
      fatal=.false.
      do l=1,no
!       rev. 9.9.05  Jan.  02/14  - NK: Add check if in-basin in flowinit
        if(inbsnflg(l).ge.1)then
          n=s(iy(l),jx(l))
          if(ireach(n).gt.0.and.res(n).eq.0)then
            if(nopt(l).eq.2)then
              fatal=.true.
              print*
              print*,'Fatal ERROR'
            endif
            write(98,98000)'Warning:Flow station',l, 
     *    'is located in a lake or reservoir that is NOT an outlet'
98000 format(a20,i5,a56)           
          endif
        endif
      end do
      print*
      if(routeflg.ne.'q')then
        if(fatal)then
            write(98,*)'Program aborted in flowinit @ 1987'
            stop 'Program aborted in flowinit @ 1987'
        endif
      ENDIF

      firstpass=.false.

      if(iopt99)then
       write(98,*)'Info: Flow initialization completed'      
       write(63,*)'Flow initialization completed'      
        write(55,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'      
        write(55,*)'Flow initialization completed'      
        write(55,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'  
      endif
            
      
      
      close(unit=99)
      
      
      
      RETURN

! FORMATS

  500 format(8f10.3)
  501 format(3i5,4x,a1)
 1011 format(' ',3x,
     *  '  i  ires(i) jres(i)    b1(i)     b2(i)    b3(i)     b4(i)')
 1013 format(' ',3x,i3,2i8,5f12.5,a12/)
 1761 format(' Error: prob. cause - gauge not in watershed'/
     *          ' gauge no. ',i3,' colum no.',i3,' row no.',i3/) 
 1777 format(' ','in flowinit: l,i,j,n,nhyd/',5i5)
 1778 format(' stream gauge ',i5,' is outside the watershed')
 1779 format('Check the gauge locations in basin/xxxx.str and'/ 
     *'in strfw/xxxx.str  Station is ignored!!!'/
     *'This is OK if you are using the sub-watershed option'/
     *'run with iopt = 1 and check the spl.err file for info'/
     *'              ')
 5003 format(2i5,4g10.3,5x,a12)
 5004 format(2i5,5g10.3,5x,a12)
 5008 format(' sub'/
     *'   id    l    i    j    n nlow    qda       area - 1st&low'/)
 5009 format(6i5,2f12.3)
 5551 format('  write #',i2,a30)
 5553 format(999i5)
 5554 format(999f8.0)
 5555 format(999f5.1)
 5556 format(' ',a1,6i7,4f12.3)
 5557 format(' msg   bsn    col   row     n  next  n    iset     da',
     *    '       qda         datemp')
 5558 format('    n  row   col       iset     da',
     *    '     qda       datemp')
 5559 format(999f6.0)
 6000 format(' ','Flowinit called to initialize flows and lzs')
 6001 format(5x,' no of storms =',i5)
 6002 format(35x,' * * * time =',f7.2,' hours * * *')
 6004 format('+',i5,' hrs rainfall available, ',i5,' hrs used')
 6005 format(' Warning: no streamflow stations selected for '/
     *' error calculation. Enter data in 1st line of .str file.')
 6007 format(' no,nl,mhtot,kt/',4i5)
 6008 format(' id,nr,tj1,mo,conv/',2i5,f5.2,i5,f5.2)
 6010 format(' reset:'/
     *  '    n    i    j         qda         qi1          qi2',                
     *  '        qo1        qo2      store1      store2         cap',
     *  '         over')
 6011 format(3i5,9f12.3)
 6043 format(5I12)
 6044 format('          n          i           j       nhyd   nxtbasin')     
 6045 format(' l,nxtbasin/',2i5)
 6046 format('n=',i5,' nhyd(',2i5,')=',i10,'  check if in watershed')
 6100 format(' [in flowinit]'/
     *'     n    i    j iset         da         qda       qbase   gage')
 6101 format(' ',4i5,3f12.3,i5)
 6102 format(' ',4i5,3f12.3)
 6103 format(' element #',i5,' base flow = 1 cms assumed')
 6104 format(999i4)
 6105 format(' initial flows modified by reservoir releases') 
 6106 format(' i,qinit(i)/',i10,f12.3)
 6107 format(' n,yyy(n),xxx(n),i,ires(i),jres(i)/',6i3)
 6108 format(' initial base flows prorated upstream')
 6109 format(' reservoir releases added back in') 
 6110 format(' basin number allocations')
 6012 format(' ',2f15.3,i5,1x,a12,3x,2a1,3x,2a1,f12.0)
 9005 format(' iymin,iymax,jxmin,jxmax/',4i10)
 9801 format(//'WARNING: DACHECK set to ',g12.3,' in flowint '//)
 5500 format(' n,i,ireach(n),nreach(i),maxr/',5i5)
51001 format(' ',i5,2g12.3,4f12.3,f16.3,3f12.3)
51002 format('     n          da         cap       chaxa      chadep',
     *'      chawid      wetwid        wcap       wetxa  %wet class',
     *'         sl2')     
51003 format(' Channel properties:' )   
99001 format(2f12.3,i5,2x,a12,a50,3f10.0,' %')
99002 format(2f12.3,i5,2x,a12,3f10.0,' %')

      END SUBROUTINE flowinit





