      SUBROUTINE flowreset()
      
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
!     rev. 10.4.41 Sep.  28/21  = NK Added flowreset at end of spinup

      USE area_watflood
      USE area_debug
      implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      real*4,        dimension(:),    allocatable :: qinit,datemp
      real*4,        dimension(:,:),  allocatable :: qdagrd
      integer*4,     dimension(:),    allocatable :: iset
      
      character(12), dimension(:), allocatable :: sta_id
      character(80), dimension(:), allocatable :: sta_name,temp_sta_name
      real*4,        dimension(:), allocatable :: sta_lat
      real*4,        dimension(:), allocatable :: sta_long
      real*4,        dimension(:), allocatable :: temp_area
      logical    :: exists
      real*4     :: da1,qda1,qch,obdepth,tdum,trialq,flow_max,convert,&
                   class_sum,awr,try1,xxyyzz
      integer    :: n,i,j,istep2,l,ktt,k,ktemp,jj,nnn1,lll,&
                   iasdf,ios,noread,nnx,inx,jnx,ii,nresv,nnext,&
                   iallocate      
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
      DATA col2/' ','a','b','c','d','e','f','g','h','i','j','k','l'&
              ,'m','n','o','p','q','r','s','t','u','v','w','x','y'/

      DATA col20/'a','b','c','d','e','f','g','h','i','j','k','l','m'&
               ,'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      DATA col21/'a','b','c','d','e','f','g','h','i','j','k','l','m'&
               ,'n','o','p','q','r','s','t','u','v','w','x','y','z'/


! * * * * * * * * * * *
      if(debug_output)print*,'In flowreset @ 87'

!     wetland geometry initialization (moved from spl May 25/05)

	messages=.false.

      if(iopt.ge.1)then
        write(55,*)'In flowreset.for'
        write(55,*)'~~~~~~~~~~~~~~~'
      endif

!      if(.not.allocated(sta_id))then
      allocate(sta_id(no),sta_name(no),sta_lat(no),sta_long(no),&
              temp_sta_name(no),&   
              temp_area(no),stat=iAllocate)
      if(iAllocate.ne.0) STOP  'Error with allocation  in flowreset @ 73'

      allocate(datemp(na),qdagrd(ycount,xcount),iset(na),stat=iAllocate)
      if(iAllocate.ne.0) STOP 'Error with allocation  in flowreset @ 78'
      
      allocate(qinit(noresv),stat=iAllocate)
      if(iAllocate.ne.0) STOP'Error with allocation of nreach in flowreset @ 97'

      
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
      if(debug_output)write(63,*)
      if(debug_output)write(63,*)'highest reach no found (by flowreset)    =',maxr
      if(debug_output)write(63,*)'no of reservoirs found in the rel files =',noresv
      if(debug_output)write(63,*)
	
	do n=1,naa
!     rev. 10.2.09 Nov.  04/17  - NK: Reinstated old Manning's n correction for legacy files      
!       From before pond routing - keep new program compatible with old method
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
            r2n(n)=amax1(r2n(n),&                              ! Manning n correction
             r2n(n)*(grid_area(n)*aclass(n,classcount-1))/&   ! Manning n correction
                                       chaarea(n)*rlake(n))  ! Manning n correction                   
            write(51,4950)n,r2n(n),aclass(n,classcount-1)
4950        format('Grid #, Adjusted R2n, water frac:',i5,2f10.3)
          endif
!         modify water surface areas from the shd file:
!         but it can't be less than calculated above
!     rev. 9.9.51  Jan.  13/15  - NK: Added min channel area in flowreset
          chaarea(n)=amax1(chaarea(n),grid_area(n)*aclass(n,classcount-1)) 
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
      end do

      if(debug_output)print*,'In flowreset @ 169'

!     FLOW INITIALIZATION SECTION
!     THIS SECTION COMPUTES THE INITIAL FLOW AND STORAGE IN EACH RIVER 
!     SEGMENT USING THE FIRST DOWNSTREAM GAUGE AND PRORATES THE FLOW 
!     ACCORDING TO DRAINAGE AREA.F

!     STEP 1: RESERVOIR RELEASES ARE SUBTRACTED

!     QDAGRD IS QDA IN A GRID ARRAY USED FOR PRINTING IN THIS S/R 
!     AND USED TO FILL IN 

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

      do n=1,naa
!       CONVERSIONS &
!       STORE THE GRIDS CONTRIBUTING TO A DWOPER REACH
        if(ireach(n).ne.0)then
	    l=ireach(n)
          nreach(l)=n
          if(iopt.ge.1)write(55,*)n,l,ireach(n),nreach(l),inbsnflg(l),maxr
        endif
	  end do




      if(debug_output)print*,'In flowreset @ 258'


!     *****************************************************************
!     STEP 1: FIND THE LAST +VE FLOW FOR THE FIRST HYDROGRAPH:
!     STEP 1: FIND THE LAST +VE FLOW FOR THE FIRST HYDROGRAPH:
!     STEP 1: FIND THE LAST +VE FLOW FOR THE FIRST HYDROGRAPH:
!     STEP 1: FIND THE LAST +VE FLOW FOR THE FIRST HYDROGRAPH:
      write(53,*)
      write(53,*)'Resetting flow '
      write(53,*)'Resetting flow '
      write(53,*)'Resetting flow '
      write(53,*)'Resetting flow '
      write(53,*)'Resetting flow '
      write(53,*)'Resetting flow '
      write(53,*)'STEP 1: '
      write(53,*)'FIND THE last +VE FLOW '
      do l=1,no
        area(l)=-1.0  ! default  nk Nov. 29/10
        if(inbsnflg(l).eq.1)then	  
          i=iy(l)
          j=jx(l)
          n=s(i,j)
          if(n.ge.1)then
!           for watroute nr needs to be min(nh,nhtot)
!           for watroute use the flag to set this up
            if(iopt.ge.1)write(55,*)nr,nl,mhtot
!           WHEN THERE IS NO RAINFALL AT ALL, NR = 0
!           THIS MIGHT HAPPEN WHEN THERE IS SNOWCOURSE DATA AND MELT
!           OCCURS WITHOUT RAIN
            flowflag(l)=.false.
!            do k=nl,nl-100,-kt             !go back 12 hours to make sure we hit WSC data
            write(53,*)
            do k=nl,nl-24,-1             !go back 24 hours to make sure we hit WSC data
!              if(qhyd(l,k).gt.0.0.and.ice_fctr(n).gt.0.99)then
              if(qhyd(l,k).gt.0.0.and.ice_fctr(n).gt.0.50)then
                  
                  
!                THE  +VE FLOW IS FOUND
                 flowflag(l)=.true.
                 qda(n)=qhyd(l,k)
                 iset(n)=1   !```````````````````````````````````````
                 datemp(n)=da(n)
                 nlow(l)=k
                 ktemp=k
                 write(53,*)'+',l,n,da(n),qda(n),qhyd(l,k)
                 exit
              else   
!     rev. 10.4.44 Oct.  10/21  = NK Fixed nudging for locations with no observed flows
!                NO  +VE FLOW IS FOUND or THERE'S ICE
!                In this case we continue on with the computed flow
!                so there should be no changes in the state variables
!                eventhough we are recomputing them
                 flowflag(l)=.true.
                 qda(n)=qo2(n)  ! continue with the computed outflow
                 iset(n)=1   !```````````````````````````````````````
                 datemp(n)=da(n)
                 nlow(l)=k
                 ktemp=k
                 write(53,*)'-',l,n,da(n),qda(n),qhyd(l,k)
                 write(*,*)'-',l,n,da(n),qda(n),qhyd(l,k)
                 
              endif
            end do
          endif
!          if(iopt.ge.1)write(55,*)flowflag(l),kt,ktt,l,n,i,j

          if(iopt.ge.1.and.l.eq.1)then
             write(55,*)'  id    l    i    j    n ','nlow(l)   qda(n)     area(l)'
          endif
          if(debug_output)write(55,5009)id,l,i,j,n,nlow(l),qda(n),area(l)
          if(debug_output)then
             write(55,5009)id,l,i,j,n,nlow(l),qda(n),area(l)
          endif
        endif   ! inbsnflg.eq.?
      end do
      
      do n=1,naa
        if(res(n).eq.0)then
!           pond routing
!                  qo2(n)=rlake(n)*store2(n)**1.75
            store1(n)=(qo1(n)/rlake(n))**0.57143
            store2(n)=store1(n)
            qo2(n)=qo1(n)
        endif
      end do
      
      
      if(iopt.ge.1)then
        iasdf=1
        msg=' reset flow at gauges'
        write(55,5551)iasdf,msg
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


      if(debug_output)print*, 'in flowreset at 300'

!     *****************************************************************
!     STEP 2: SUBTRACT RESERVOIR RELEASES FROM D/S GAUGE FLOWS
!     STEP 2: SUBTRACT RESERVOIR RELEASES FROM D/S GAUGE FLOWS
!     STEP 2: SUBTRACT RESERVOIR RELEASES FROM D/S GAUGE FLOWS
      write(55,*)
      write(55,*)'STEP 2: '
      write(55,*)'SUBTRACT RESERVOIR RELEASES FROM D/S GAUGE FLOWS'
	write(55,*)'Basin no assignment at each flow station:'
	write(55,*)

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
      if(debug_output)print*, 'in flowreset at 360'
        do j=ktr,nrel,ktr
          do i=1,noread
            if(qinit(i).le.0.0.and.qrel(i,j).ge.0.0)then 
              qinit(i)=qrel(i,j)
            endif
          end do
        end do
      if(debug_output)print*, 'in flowreset at 400'

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
            do while(.not.resflag.and.n.le.naa.and.n.gt.0.and.next(n).gt.0)
               if(qda(n).gt.0.0)then	

!               WE'RE AT A GAUGE AND WE'LL SUBTRACT OUT THE RELEASE
!               RELEASE CAN'T BE GREATER THAN THE GAUGE FLOW
!               nothing is taken out if flow = natural

!     rev. 9.5.66  Oct.  06/09  - NK: fixed bug in flowreset for init flows < 1.0
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

      if(iopt.ge.1)then
         write(55,6105)
         write(55,6100)  ! [in flowreset]
         write(55,6102)
         do n=1,na
           IF(qda(n).ne.0.0.or.qbase(n).ne.0.0)then
!            find the gauge number in this grid           
             Do l=1,no
               if(xxx(n).eq.jx(l).and.yyy(n).eq.iy(l))i=l    
             end do
	       write(55,6101)n,yyy(n),xxx(n),iset(n),da(n),qda(n),qbase(n),i
	     endif
	   end do
	   write(55,*)'In station order:'
	   write(55,*)'Note stations witn no initial flows!!'
!     rev. 9.9.50  Jan.  07/14  - NK: Added zero - initial flow warning
         do l=1,no
           if(iy(l).ge.1.and.jx(l).ge.1.and.iy(l).le.ycount.and.jx(l).le.xcount)then    ! added Jan 23/15  NK   
             n=s(iy(l),jx(l))
             if(n.gt.0.and.n.le.naa)then  ! added Jan 23/15  NK
	         write(55,6101)n,jx(l),iy(l),iset(n),da(n),qda(n),qbase(n),l
               if(qda(n).eq.0.0.and.iopt.ge.1)then
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

      if(debug_output)print*, 'in flowreset at 500'

!     *****************************************************************
!     REV. 8.25 - May.  22/97 -  FIXED ALLOCATING THE BASIN # IN flowreset
!     rev. 9.9.22  Jul.  29/14  - NK: Fixed basin no assignment in flowreset.f
!     ASSIGN THE GRID TO A DRAINAGE BASIN FOR CALC OF SUM PRECIP 
!     Step 3 WORK UPSTREAM:
!     Step 3 WORK UPSTREAM:
!     Step 3 WORK UPSTREAM:
!     Assign basin numbers to grids u/s of gauging stations
      Write(55,*)
      write(55,*)'STEP 3:'
      write(55,*)'Assign u/s grids with basin numbers - d/s gauge no.'
      write(55,*)'Assign init flow to u/s stations - prorated with d.a.'
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

      do n=1,no
         nxtbasin(n)=0      !added (n) nk
      end do

      if(debug_output)print*, 'in flowreset at 600'

      if(iopt.ge.1)write(55,6044)
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

!     rev. 9.04    Jan    16/01 - fixed grid diagnosis in flowreset
!!1111111111111                     STOP 'program stopped in flowreset @65'
                  endif
                  if(nhyd(i,j).gt.0)then
                    nxtbasin(nhyd(i,j))=nhyd(inx,jnx)
                  endif
                  if(iopt.ge.1)then
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
      if(iopt.ge.1)write(55,6045)(l,nxtbasin(l),l=1,no)

      if(debug_output)print*, 'in flowreset at 700'

!     output added Jul. 29/14  NK
      if(iopt.ge.1)then
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
                     qda(n)=qda(nnx)*da(n)/da(nnx)      !  <<<<<<<<<<<<<<<<<<<<<<<
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

      if(debug_output)write(55,6108)
      if(debug_output)write(55,6100)
      if(debug_output)write(55,6102)(n,yyy(n),xxx(n),iset(n),da(n),qda(n),qbase(n),n=1,naa)

      if(iopt.ge.1)then
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

      if(debug_output)print*, 'in flowreset at 800'


!     *****************************************************************
!     STEP 4
!     FILL IN REMAINING GRIDS
!     WORK DOWNSTREAM
      Write(55,*)
      write(55,*)'STEP 4:'
      write(55,*)'Assign d/s grids with basin numbers '
      write(55,*)'Assign init flow to d/s stations - prorated with d.a.'

!     rev. 9.9.22  Jul.  29/14  - NK: Fixed basin no assignment in flowreset.f
      do n=1,naa
!        this added Nov. 21/07  nk
!        to ensure that all grids belong to a basin even if below 
!        a stream flow gauge.
!        CURRENT LOCATION:
!        LOCATION OF DOWNSTREAM ELEMENT:
         nnx=next(n)
!     REV. 10.1.22 Jan.  25/16  - NK: Fixed flowreset for partioa basins
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
      if(iopt.ge.1)then
	   write(55,*)' basin number allocations after working downstream'
         do i=ycount,1,-1
            write(55,6104)(nhyd(i,j),j=1,xcount)
         end do
      endif
!     output added Jul. 29/14  NK
      
!     rev. 9.9.22  Jul.  29/14  - NK: Fixed basin no assignment in flowreset.f
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
      if(iopt.ge.1)then
	   write(55,*)' basin number allocations after working upstream'
	   write(55,*)' to fill in remaining holes'
         do i=ycount,1,-1
            write(55,6104)(nhyd(i,j),j=1,xcount)
         end do
      endif
!     output added Jul. 29/14  NK

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

      if(debug_output)write(55,55099)&
               '      l      i      j     n  next(n)   iset(n)& 
       datemp        qda      datemp        qda'
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
           stop 'Program aborted in flowreset @ 1153'           
         else
           msg1='a'
           if(debug_output)then
                 if(debug_output)write(55,*)&
              '      l    i    j    n next(n) iset(n) datemp',&
              '         qda      datemp        qda'
                  if(debug_output)write(55,5556)msg1,l,i,j,n,next(n),iset(n),datemp(n),&
              qda(n),datemp(nnx),qda(nnx)
           endif
         endif
!        this loop revised Dec. 22/05  nk

!        check added Apr. 26/12 nk
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
           stop 'Program aborted in flowreset @ 1127'
         endif
         do while(n.lt.na.and.iset(nnx).ne.1&
     		        .and.n.gt.0.and.next(n).gt.0)    !??????????????????????
!          DOWNSTREAM FLOW HAS NOT BEEN SET BY A GAUGE
!          IT COULD HAVE BEEN SET BY A GAUGE UP ANOTHER TRIBUTARY 
!          IF ISET=1 THEN FLOW HAS BEEN SET BY A DOWNSTREAM GAUGE
!          AND WE DON'T TOUCH IT.
           nnx=next(n)    ! leave in - we're in an n loop!

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
           if(iopt.ge.1)then
            write(55,5556)msg1,nhyd(i,j),xxx(n),yyy(n),n,next(n),&
                  iset(n),datemp(n),qda(n),datemp(nnx),qda(nnx)
           endif
           n=next(n)
         end do
	 endif   ! inbsn
      end do
      

!     output added Jul. 29/14  NK
      if(iopt99)then
        open(unit=99,file='debug\basin_no.xyz',status='unknown')
        do i=1,ycount
          do j=1,xcount
            write(99,*)float(j)*xdelta+xorigin-0.5*xdelta,&
                    float(i)*ydelta+yorigin-0.5*ydelta,&
                    nhyd(i,j)             
          end do
        end do
        close(unit=99,status='keep')
      
!     output added Jul. 29/14  NK
        do i=1,ycount
          do j=1,xcount
            outarray(i,j)=nhyd(i,j)             
          end do
        end do
!     rev. 9.9.75  Sep.  11/15  - NK: Added basin_no.r2c output to flowreset.f
        author='flowreset'     
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

      if(debug_output)print*, 'in flowreset at 801'

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
      IF(IOPT.GE.1)THEN
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
      if(iopt.ge.1)write(55,*)'Prorating the flows downstream:'
      do n=1,naa
        if(iset(n).eq.2)then
!         IF ISET=2 flow can be modified by tributary
          if(iopt.ge.1)write(55,*)n,qda(n),da(n),datemp(n)
!     REV. 10.1.22 Jan.  25/16  - NK: Fixed flowreset for partial basins
          if(next(n).gt.0)then
            if(iopt99)then        ! added Feb. 2/18 NK
              fatal=.false.
              if(da(next(n)).le.0.0)then
                write(98,98001)'Warning: Downstream drainage area ='&
                ,da(next(n)),&
                ' in grid #, x,y:',n,xxx(n),yyy(n),' ok if outlet grid'
98001            format(a35,i10,a16,3i8,a18)  
                fatal=.true.
              endif
            endif
!     rev. 9.9.31  Oct.  13/14  - NK: Changed flow initialization RE: zero init flows
            if(da(next(n)).gt.0.0)then
!             fixed Jan. 4/15  nk
              qda(next(n))=qda(n)*da(n)/da(next(n))  !<<<<<<<<<<<<<<<<check check check check 
            else
!             this should only happen at the outlets where da is set = 0.0            
             qda(next(n))=qda(n)
            endif
	      if(iopt.ge.1)write(55,*)n,qda(n),da(n),da(next(n))
	      if(iopt.ge.1)write(55,*)
            iset(n)=1
          endif
        endif    !     REV. 10.1.22
      end do
      if(fatal.and.iopt99)then
        print*,'OK if these locations are receiving grids'
        print*,'OK if these locations are receiving grids'
        print*
      endif
      
      if(iopt.ge.1)then
         iasdf=5
         msg=' base flows routed downstream'
         write(55,5551)iasdf,msg
      endif
!     FIXED A BUG HERE ON APRIL 16/96
!     THIS DO EXECTUED ONLY FOR IOPT.GE.1  -  NO GOOD I THINK

      if(iopt.ge.1)write(55,5557)
      do n=1,naa
         i=yyy(n)
         j=xxx(n)
         qdagrd(i,j)=qda(n)
         qbase(n)=qda(n)
         if(iopt.ge.1)then
            write(55,5556)msg1,nhyd(i,j),xxx(n),yyy(n),n,next(n),&
                         iset(n),da(n),qda(n),datemp(n)
         endif
      end do
      if(debug_output)print*, 'in flowreset at 802'
      if(iopt.ge.1)then
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
      if(debug_output)print*, 'in flowreset at 802'

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
!     rev. 9.2.18  Oct.  27/05  - NK: Fixed bug in flowreset (init spike)
!           qda(n)=qda(nnx)*da(n)/datemp(nnx)
            qda(n)=qda(nnx)*da(n)/da(nnx)
          else
            qda(n)=0.1     ! ok Jul. 11/02
          endif
          qbase(n)=qda(n)
        endif
      end do

      if(iopt.ge.1)then
         iasdf=6
         msg=' distributed base flows'
         write(55,5551)iasdf,msg
         msg1='f'
        write(55,5557) 
        do n=1,naa
            i=yyy(n)
            j=xxx(n)
            write(55,5556)msg1,nhyd(i,j),xxx(n),yyy(n),n,next(n),&
                         iset(n),da(n),qda(n),datemp(n)
         end do
      endif
      if(debug_output)print*, 'in flowreset at 803'

!     FIX ON APRIL 17/97    
      do n=1,naa
         i=yyy(n)
         j=xxx(n)
         qdagrd(i,j)=qda(n)
         qbase(n)=qda(n)
      end do

      if(iopt.ge.1)then 
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
      if(debug_output)print*, 'in flowreset at 804'

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
         do while(.not.resflag.and.n.le.naa.and.n.gt.0.and.next(n).gt.0)     ! ???????????????
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

      if(iopt.ge.1) write(55,6109)

      if(debug_output)write(55,6100)
      if(debug_output)write(55,6102)(n,yyy(n),xxx(n),iset(n),da(n),qda(n),qbase(n),n=1,naa)
      if(debug_output)print*, 'in flowreset at 805'


!     check it out!!!!!!!!!!!!!!!!!1
!     iopt in the wrong place??????  changed jan 29/09 nk



      if(iopt.ge.1)then
         iasdf=7
         msg=' reservoir flows added back in'
         write(55,5551)iasdf,msg
	endif

      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        qdagrd(i,j)=qda(n)
      end do

      if(iopt.ge.1)then
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
      if(debug_output)print*, 'in flowreset at 806'

!     REPLACE THE PROPER QDA VALUES FOR LZS INITIATION IN S/R RUNOF5
      do n=1,naa
         qbase(n)=abs(qbase(n))
      end do

      if(iopt.ge.1)then
         write(55,6100)
         write(55,6102)(n,yyy(n),xxx(n),iset(n),da(n),qda(n),qbase(n),n=1,naa)
      endif

!     FIX ON June 19/06  nk    
      do n=1,naa
         i=yyy(n)
         j=xxx(n)
         qdagrd(i,j)=qda(n)
      end do

      if(iopt.ge.1)then 
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

      if(iopt.ge.1)then 
        iasdf=9
        msg=' write nhyd(i,j) '
	  write(55,6110)
        write(55,5551)iasdf,msg
        do i=ycount,1,-1
          write(55,5553)(nhyd(i,j),j=1,xcount)
        end do
      endif
      if(debug_output)print*, 'in flowreset at 808'

!     THIS SECTION CAME FROM RESET.FOR

!     rev  9.1.02  July  12/01  - put in dacheck in flowreset for flag
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

!     rev. 9.1.36  Jan.  28/03  - Fixed wetland init condition in flowreset

!        check for r2 value made in param.for
! * * * * * * * * TS - ADDED WATFLOOD ROUTING INITIALIZATION * * * * * *

         if(manningflg.eq.'y')then
           qch=(cap(n)/rl(n))**1.67*slope(n)/chawid(n)**0.667/r2n(n)
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
!               WETLAND FLOW INITIALIZATION
                over(n)=0.0 
                if(manningflg.eq.'y')then
                   store1(n)=rl(n)*(qo2(n)*chawid(n)**0.667*r2n(n)/slope(n))**.60
	            else
                   store1(n)=rl(n)*(qo2(n)*r2(n)/slope(n))**.75
                endif

                store2(n)=store1(n)
                flowxa(n)=store1(n)/rl(n)
	            hcha1(n)=amax1(0.0,flowxa(n)/chawid(n))!  mar 07/09 mk
	            hcha2(n)=hcha1(n)
	            hwet1(n)=hcha1(n)
	            hwet2(n)=hwet1(n)
!               assumes inbank flow only
	            wstore1(n)=wetwid(n)*rl(n)*hwet1(n)*abs(theta(n))
	            wstore2(n)=wstore1(n)
                satxa(n)=wstore1(n)/rl(n)/abs(theta(n))
	          else
!                CHANNEL FLOW INITIALIZATION
!                 if(aclass(n,classcount-2).gt.0.0) 
                over(n)=0.0
                if(manningflg.eq.'y')then
                   store1(n)=rl(n)*(qo2(n)*chawid(n)**0.667*r2n(n)/slope(n))**.60
	            else
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
                   over(n)=((qo2(n)-qch)*r1n(n)*6.0/slope(n))**.75
!                  the factor of 6.0 is incorportated in r1
                 else
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
                   over(n)=((qo2(n)-qch)*r1n(n)*6.0/slope(n))**.75
!                  the factor of 6.0 is incorportated in r1
                 else
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
         do ii=1,classcount
            r(n,ii)=0.0
         end do
         qr(n)=0.

      end do

      if(debug_output)print*, 'in flowreset at 1123'

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
          write(98,*)'Warning: Zero or negative drainage area found in grid # ',n
          write(98,*)'Warning:To accept this, hit return - but check outcome'
          write(98,*)'Warning:Program may crash'

          qbase(n)=0.0
! * * * * * * *  TS  - added initializations * * * * * * * * 
          if(wetflg.eq.'y')then
	      qiwet1(n)=qbase(n)
	      qiwet2(n)=qbase(n)
	      qowet1(n)=qbase(n)
	      qowet2(n)=qbase(n)
          endif
        endif
          flz2(n)=1.0-(1.0-flz(n))
          pwr2(n)=1.0/pwr(n)       ! used in soilinit 
        if(qbase(n).gt.0.0)then
!         WITHOUT THIS, N IS SET TO 0 !!!!!  FAR OUT & VERY WIERD
          lzs(n)=(qbase(n)/flz2(n)/tdum/frac(n))**pwr2(n)

!     rev. 9.4.13  Jul.  09/07  - NK: modified lzs to account for lake area (flowreset) 
!         There is no lzs under laks & reservoirs so adjust:
          lzs(n)=lzs(n)*(1.0-aclass(n,classcount-1))

!         print*,n,qbase(n),flz2(jj),frac(n),pwr2(jj),lzs(n)
        else
!         print*,n,qbase(n),flz2(jj),frac(n),pwr2(jj)
          lzs(n)=1.0
!         DEFAULT VALUE
        endif
      if(debug_output.and.n.eq.1)write(55,*)'           n        flz         flz2       pwr         pwr2'  
      if(debug_output)write(55,*)n,flz(n),flz2(n),pwr(n),pwr2(n)
      end do
      if(debug_output)print*, 'in flowreset at 1173'

      if(iopt.ge.1)then
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

      if(debug_output)then
        if(debug_output)write(55,6010)
        if(debug_output)write(55,6011)(n,i,j,qda(n),qi1(n),qi2(n),qo1(n),qo2(n),&
              store1(n),store2(n),cap(n),over(n),n=1,na)
      endif

!     initialize reservoir storage
      if(debug_output)then
        write(53,*)'In flowreset'
        write(53,*)'Note: set iopt = 0 to not write this.'
        write(53,*)'Initial reservoir storage & outflow:'
        write(53,*)'     grid#        res#        try1       storage',&
                  '       outflow'
	  endif


      if(numa.eq.0)write(51,*)'Flow initializing completed'

      if(debug_output)print*, 'in flowreset at 1223'


!     rev. 9.5.25  Mar.  20/08  - NK: fixed lake initiation - moved code route -> flowreset
!       this bunch of code was moved from route where it was repeated if timestep < 1 hour
!       This one kept me up until 1 am!

        do n=1,naa
          qowetsum(n)=0.0
          qiwetsum(n)=0.0
          wstoreinit(n)=wstore2(n)  
        end do

        convert=al*al*1000.0  ! converts vol in m^3 to mm on the unit grid

        do n=1,naa
          res(n)=0
        end do
        
        
!       Make a table for the next downstream lake        
!           LOOK TO SEE IF THERE IS A RESERVOIR IN THIS GRID:
        

        do n=1,naa
          do l=1,noresv
            lll=next(n)
!           LOOK TO SEE IF THERE IS A RESERVOIR IN THIS GRID:
            if(yyy(n).eq.ires(l).and.xxx(n).eq.jres(l))then
                
              res(n)=l
              
!             IF THERE IS, SET INITIAL RESERVOIR STORAGE:
!              if(b1(l).gt.0.0.and.b2(l).gt.0.0)then
              if(b1(l).ne.0.0.and.b2(l).ne.0.0)then
!               WE HAVE A NATURAL CONTROL & WE NEED INITIAL
!               RESERVOIR STORAGE
!               INITIAL FLOWS ARE CALC IN SUB SO LEAVE THEM ALONE 
!     rev. 9.1.07  Jan.   3/02  - check that outlet is in a lake
!                 Jan. 03/02  NK
                if(ireach(n).eq.0)then
                  print*,' in flowreset in grid number',n
                  print*,' reach number =',ireach(n)
                  print*,' i.e. lake or reservoir #',l,' outlet' 
                  print*,' is not in a lake'
                  print*,' please make grid part of a lake or  '
                  print*,' move outlet to a grid in a lake'
!                 note that dam locations with releases do not have to be in a lake!! 
                  print*
                  stop 'Program aborted in flowreset @ 170'
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
                  endif
                endif

                write(53,*)
                write(53,*)'Reservoir or lake no:',l
                write(53,53999)'b`s',b1(l),b2(l),b3(l),b4(l),b5(l),b6(l),b7(l)
53999           format(a5,5g12.3,2f12.3) 
                write(53,*)'qo1(n) =',qo1(n),' qo2(n) =',qo2(n)

                if(b1(l).gt.0.0)then
!                 we have a rating curve and calculate initial storage                
                  if(b3(l).eq.0.0)then
!                   use 2 coefficients  
                          write(53,*)l,'power function'
                          store1(n)=(qo2(n)/b1(l))**(1/b2(l))
                          store2(n)=store1(n)
                          try1=qo1(n)
                          write(53,*)n,l,try1,qo1(n),store1(n),'done'
                          
                  else
                    write(53,*)l,'bisection method'
!                   use bisection to get init flow
!                   actually, I made the int. storage a little larger
!                   needs to be fixed
                      store1(n)=1000.
                      try1=0.0
!                     use 2-5 coefficients
                      write(53,*)'         n           l          try1',&    
                          '         qo1         store1'
                      do while(try1.lt.qo1(lll))
!                       keep doubling the res. storage until the corresponding
!                       flow is larger than the initialized streamflow.
                        store1(n)=2.0*store1(n)                      
!     rev. 9.8.89  Oct.  27/13  - NK: Fixed undefined (NAN) problem in flowint
                        if(b5(l).ne.0.00000)then
                          try1=b1(l)*store1(n)+&
                              b2(l)*store1(n)**2+&
                              b3(l)*store1(n)**3+&
                              b4(l)*store1(n)**4+&
                              b5(l)*store1(n)**5
                        elseif(b5(l).eq.0.0000.and.b4(l).ne.0.00000)then
                          try1=b1(l)*store1(n)+&
                              b2(l)*store1(n)**2+&
                              b3(l)*store1(n)**3+&
                              b4(l)*store1(n)**4
                        else
                          try1=b1(l)*store1(n)+&
                              b2(l)*store1(n)**2+&
                              b3(l)*store1(n)**3
                        endif
      if(debug_output)write(53,*)b1(l)*store1(n),b2(l)*store1(n)**2,b3(l)*store1(n)**3
      if(debug_output) write(53,*),b4(l)*store1(n)**4,b5(l)*store1(n)**5
                        write(53,*)n,l,try1,qo1(n),store1(n),'poli'
                        if(try1.lt.0.0)then
                          print*,' trial value for flow =',try1
                          print*,' trial reservoir outflow < 0.0'
                          print*,' polinimial functins not'
                          print*,' monotonically increasing'
                          print*,' for reservoir',l
                          print*,' Please fix the function'
                          print*,' Program may excecute but results' 
                          print*,' are approximate'
                          print*,'                   @ 167/flowreset '
                          print*
                        endif
                      end do
                  endif

                  write(53,*)n,l,try1,qo1(n),store1(n),'done'
                  write(53,*)
                else  ! b1(l)=0.0
!     rev. 9.8.50  Feb.  27/13  - NK: Initialize store1&2() for zero lake outflow
!                 use a default storage - no outflow anyways             
                  store1(n)=100.0
                  store2(n)=100.0

                endif   !if(b1(l).ne.0.0.and.b2(l).ne.0.0)then
      if(debug_output)print*, 'in flowreset at 1224'

!                if(resumflg.ne.'y')then
!                  store1(n)=max(100.0,store1(n))
!                  store2(n)=store1(n)
!                  if(iopt.ge.1)then
!                    write(51,8492)
 !                   write(51,8493)n,l,b1(l),b2(l),qda(n),store1(n)
!                    write(53,8492)
!                    write(53,8493)n,l,b1(l),b2(l),qda(n),store1(n)
!                  endif
!                endif
              endif     !if(b1(l).ne.0.0.and.b2(l).ne.0.0)then

!     rev. 9.9.34  Oct.  17/14  - NK: Added re-compute of lake storage re: new lake levels
!             this section duplicated in sub with this revision #
	        if(b6(l).ne.0.0.or.b7(l).ne.0.0)then
!               initial lake levels can be lower than the datum	        
                write(53,*)
                write(53,*)'Override flow based initial storage:'
                write(53,*)'Override flow based initial storage:'
                write(53,*)'Override flow based initial storage:'
!               Note: these are over written again in rerout for special cases
!               Note: these are over written again in rerout for special cases
!               Note: these are over written again in rerout for special cases
                

      go to 99999
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
                   
                if(store2(n).ge.0.0)then
                  if(b3(l).eq.0.0)then
!     REV. 10.1.21 Jan.  23/16  - NK: Fixed lake init flow bug in flowreset
                    qo1(n)=b1(l)*store2(n)**b2(l)
                  else
                    qo1(n)=store2(n)*(b1(l)+store2(n)*(b2(l)+store2(n)*&
                     (b3(l)+store2(n)*(b4(l)+b5(l)*store2(n)))))
                  endif   
                  qo2(n)=qo1(n)
                else
                  qo2(n)=0.0
                endif
              endif
      99999 continue
!     rev. 9.9.46  Dec.  10/14  - NK: Added check on initial lake outflow
              if(qo2(n)/da(n).gt.0.10)then
                write(98,98002)&
                'Warning: Initial calculated outflow for lake',l,&
                ' is ',qo2(n),' CMS which seems high ',&
                'Check  initial lake level or rating curve'
98002         format(a44,i5,a4,f10.0,a22,a41)          
              endif
      if(debug_output)print*, 'in flowreset at 1225'
                 
!     rev. 9.9.20  Jul.  24/14  - NK: Added dead storage for lakes "store_dead"
!     rev. 9.9.38  Nov.  12/14  - NK: Added LKdepth to ill file
!               Add dead storage:
!             -----------------                
!               here dead storage is added to the live storage     
!               Use only for reservoir with releases
!               b6-b7 = current live storage
!               lake depth is the total lake depth used in the lake evap routine
!               store1 & store2 are TOTAL storage of the lake
!                store_dead(l)=((b7(l)-LKinvert(l)))*lake_area(l)
!                store_dead(l)=amax1(0.0,store_dead(l))
!                LKdepth(l)=b6(l)-LKinvert(l)
                          store1(n)=(qo2(n)/b1(l))**(1/b2(l))
                          store2(n)=store1(n)
                store1(n)=store1(n)+store_dead(l)
                store2(n)=store1(n)
!c              endif
              
                
                
              write(53,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
              write(53,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
              write(53,*)'Resetting flows at end of id= ',id
              write(53,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
              write(53,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
              write(53,*)'Grid #               ',n
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
              write(53,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
            endif
         end do
      end do  
        
      do n=1,naa  
        if(res(n).ne.0)then
        if(rlake(n).lt.1.0.and.rlake(n).gt.0.1)then
            write(98,98003)&
              'Warning: rlake(',ibn(n),') has the value of ',rlake(n),&
                     'Warning: Manning`s n will be decreased'
98003       format(a20,i5,a15,f10.3,a30)            
        endif
        endif
      end do
                        
!     rev. 10.2.20 Apr.  06/18  - NK: Added res_next for NBS
!     This section needs to be here because res(n) needs to be defines as above
         write(53,*)
         write(53,*)'      Derived in flowreset @ 2515'
         write(53,*)
         write(53,*)'      Lake #           Next Lake #'
!        initialize values         
         do l=1,noresv
             res_next(l)=0
         end do
         do l=1,noresv
              i=ires(l)
              j=jres(l)
!             Check that we are in the basin - needed when dealing with
!             sub-watershed models              
              if(i.le.ycount.and.j.le.xcount.and.i.gt.0.and.j.gt.0)then
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
                          if(nresv.gt.0)write(53,53000)&
                              resname(l)(1:12),l,& 
                              res_next(l),resname(res_next(l))(1:12)
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
 8492 format(' ','initialize resvr flow & storage - in flowreset'/&
                 'n,l,b1,b2,qda,store1')
 8493 format(' in route/flowreset:',2i5,4e12.3)
     
      if(debug_output)print*, 'in flowreset at 1226'
      
!     rev. 10.2.07 Nov.  03/17  - NK: New rt_pond subroutine for channel pond routing    
      if(iopt99)write(51,4951)&
         '          #          n          water          gridArea   ',&
         ' channelArea  channel/gridArea iak ireach wetflg pondflg  '      
4951      format(2a58)
      do n=1,naa
          if(pondflg(n))then
             store1(n)=(qo2(n)/rlake(n))**(1.0/1.75)
             store2(n)=store1(n)
          endif
      end do
      if(iopt99)then
          do n=1,naa
              write(51,4952)n,r2n(n),aclass(n,classcount-1),&
                        grid_area(n),chaarea(n),&
                        chaarea(n)/grid_area(n),&
                        ibn(n),ireach(n),&
                        wetland_flag(n),pondflg(n)
4952          format('Grid # :',i5,5f15.3,2i5,5x,l1,5x,l1)
          end do
      endif
!     rev. 9.7.01  Jun.  09/10  - NK: fixed error.xyz & error.r2s

!     create flow station heirarchy
      if(debug_output)write(63,*)'**********************'
      if(debug_output)write(63,*)'**********************'
      if(debug_output)write(63,*)'In flowreset:'     
      if(debug_output)write(63,*)'     flow station    downstream flow station'

      do l=1,no
        if(inbsnflg(l).eq.1)then    ! added nk Nov. 29/10
	    ds_sta(l)=-1
!        what grid no?
          i=iy(l)
          j=jx(l)
	    n=s(i,j)
!         find ds station - start in the next grid
          n2=n
          found_sta=.false.
          do while(.not.found_sta.and.n2.lt.naa&
     		        .and.n2.gt.0.and.next(n).gt.0)         ! ???????????????
	      n2=next(n2)    ! cell # we are in
!           see if there's a staton in this grid
!           loop thru stations
	      ll=0
	      do while(ll.lt.no.and.n2.gt.0)  ! ??????????????????
              ll=ll+1
	        if(l.ne.ll)then
                if(inbsnflg(ll).eq.1)then    ! added nk Jan. 05/11
!                 skip current station
!                 find sta in this grid & get its number
	            ii=iy(ll)
                  jj=jx(ll)
	            nnn1=s(ii,jj)  
	            if(nnn1.eq.n2)then  
!                   found a station
	              ds_sta(l)=ll
	              us_sta(ds_sta(l))=l
	              found_sta=.true.
	            endif
	          endif
	        endif
	      end do
	    end do
	  endif  ! inbsnflg = 1
	end do
!     end
!     rev. 9.7.01  Jun.  09/10  - NK: fixed error.xyz & error.r2s

      if(debug_output)print*, 'in flowreset at 1227'

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
55999   format(999i3)
	  end do
	endif
      
!     rev. 9.8.39  Nov.  26/12  - NK: added check for flow stations in lakes
!     Check that no flow stations are in lakse or reservoirs
!     and only in outlets if they are.
      fatal=.false.
      do l=1,no
!       rev. 9.9.05  Jan.  02/14  - NK: Add check if in-basin in flowreset
        if(inbsnflg(l).ge.1)then
          n=s(iy(l),jx(l))
          if(ireach(n).gt.0.and.res(n).eq.0)then
            if(nopt(l).eq.2)then
              fatal=.true.
              print*
              print*,'Fatal ERROR'
            endif
            write(98,98000)'Warning:Flow station',l,& 
         'is located in a lake or reservoir that is NOT an outlet'
98000 format(a20,i5,a56)           
          endif
        endif
      end do
      print*
      if(routeflg.ne.'q')then
        if(fatal)then
            write(98,*)'Program aborted in flowreset @ 1987'
            stop 'Program aborted in flowreset @ 1987'
        endif
      ENDIF

      firstpass=.false.

      if(iopt.ge.1)then
       write(98,*)'Info: Flow initialization completed'      
       write(63,*)'Flow initialization completed'      
        write(55,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'      
        write(55,*)'Flow initialization completed'      
        write(55,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'  
      endif
            
      RETURN

! FORMATS

  500 format(8f10.3)
  501 format(3i5,4x,a1)
 1011 format(' ',3x,'  i  ires(i) jres(i)    b1(i)     b2(i)    b3(i)     b4(i)')
 1013 format(' ',3x,i3,2i8,5f12.5,a12/)
 1761 format(' Error: prob. cause - gauge not in watershed'/&
               ' gauge no. ',i3,' colum no.',i3,' row no.',i3/) 
 1777 format(' ','in flowreset: l,i,j,n,nhyd/',5i5)
 1778 format(' stream gauge ',i5,' is outside the watershed')
 1779 format('Check the gauge locations in basin/xxxx.str and'/& 
     'in strfw/xxxx.str  Station is ignored!!!'/&
     'This is OK if you are using the sub-watershed option'/&
     'run with iopt = 1 and check the spl.err file for info'/&
     '              ')
 5003 format(2i5,4g10.3,5x,a12)
 5004 format(2i5,5g10.3,5x,a12)
 5008 format(' sub'/'   id    l    i    j    n nlow    qda       area - 1st&low'/)
 5009 format(6i5,2f12.3)
 5551 format('  write #',i2,a30)
 5553 format(999i5)
 5554 format(999f8.0)
 5555 format(999f5.1)
 5556 format(' ',a1,6i7,4f12.3)
5557  format(' msg   bsn    col   row     n  next  n    iset     da',&      
         '       qda         datemp')
 5558 format('    n  row   col       iset     da',&
         '     qda       datemp')
 5559 format(999f6.0)
 6000 format(' ','flowreset called to initialize flows and lzs')
 6001 format(5x,' no of storms =',i5)
 6002 format(35x,' * * * time =',f7.2,' hours * * *')
 6004 format('+',i5,' hrs rainfall available, ',i5,' hrs used')
 6005 format(' Warning: no streamflow stations selected for '/&
     ' error calculation. Enter data in 1st line of .str file.')
 6007 format(' no,nl,mhtot,kt/',4i5)
 6008 format(' id,nr,tj1,mo,conv/',2i5,f5.2,i5,f5.2)
 6010 format(' reset:'/&
       '    n    i    j         qda         qi1          qi2',&             
       '        qo1        qo2      store1      store2         cap',&
       '         over')
 6011 format(3i5,9f12.3)
 6043 format(5I12)
 6044 format('          n          i           j       nhyd   nxtbasin')     
 6045 format(' l,nxtbasin/',2i5)
 6046 format('n=',i5,' nhyd(',2i5,')=',i10,'  check if in watershed')
 6100 format(' [in flowreset]'/&
      '     n    i    j iset         da         qda       qbase   gage')
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
51002 format('     n          da         cap       chaxa      chadep',&
      '      chawid      wetwid        wcap       wetxa  %wet class',&
      '         sl2')     
51003 format(' Channel properties:' )   
99001 format(2f12.3,i5,2x,a12,a50,3f10.0,' %')
99002 format(2f12.3,i5,2x,a12,3f10.0,' %')

      END SUBROUTINE flowreset





