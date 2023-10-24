      SUBROUTINE rdshed(unitNum,flnNum)

!***********************************************************************
!
!  THIS SUBROUTINE READS ALL THE DATA CALCULATED BY SEPERATE
!  PROGRAM CALLED "BASIN"
!          effeciency not to important- called only once
!          a total of 10 resevoirs are allowed.
!          b1,b2,b3,b4, ....... are parameters for  each resevoir,
!
!     REV. 7.2    Sept. 19/94 -  added ireach(n) for dwoper input 
!     REV. 8.97   July  12/99 -  demonstration copy addition
!     REV. 8.99a  Jul.     99 -  lat-long watershed data
!     REV. 8.99c  Oct.   5/99 -  irough -> sl1 input in shed
!
!     REV  9.00   Mar.   2000 -  TS: CONVERSION TO FORTRAN 90
!     REV. 9.02   Oct.   5/00 -  added option to debug on one grid
!     REV. 9.03   Nov.   2000 -  TS: ADDED WATFLOOD SWAMP ROUTING
!     rev. 9.1.06  Oct.  16/01  - nrvr added to area3 to set # river types
!     rev. 9.1.46  Jul.  17/03  - WATFLOOD LITE incorporated 
!     rev. 9.1.57  Jul.  06/04  - NK: Fixed major bug in shed.for max instead of min
!     rev. 9.1.60  Jul.  27/04  - NK: reversed definitions for sl1 & sl2 Int. Slope
!     rev. 9.2.28  Jan.  30/06  - NK: Added low slope a4 for grids with water
!     rev. 9.3.03  Sep.  09/06  - NK: read s(i,j) from table instead of grid
!
!
!  classcount    - max number of permeability classes in an element
!  al       - grid length in m.
!  step     - grid length in km but if al<1000, it is equal to 1
!  na       - total # of squares in the watershed including extra one in
!             which water runs. e.g.lake
!  xxx(n)    - the north-south coordinate of the n'th element
!  yyy(n)   - the east-west coordinate of the n'th element
!  order    - coordinate file,highest to lowest element.
!  da(n)    - drainage element (n).
!  ibn(n)   - basin number 1-5
!  irough(n)- over land flow conveyance class 1-5.
!  ichnl(n) - #channels/square class1-59
!  aimp(n)  - q/age impermeable area.from.01 to .9
!  next(n)  - the element number into which the n'th element drains
!  s(i,j)   - contains the order number for each square (used to be 
!             drainage direction)
!           - impevious area should never be 1.0 or 0.0.

!***********************************************************************

      use area_watflood

!// Added by Dave
      USE EF_Module
      implicit none
	TYPE(ShedParam) :: header
!// End Dave addition

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

        CHARACTER(128):: qstr
        character(20) :: junk,junk20
	  character(1)  :: junk1
        CHARACTER(80) :: comment
!        CHARACTER(79) :: hdrcomment(100)
        INTEGER :: chksum(160),thms(25),gr10k(25),saug(25),hmbr(25),
     *             dffn(25),simpson(25),colum(50),canag(25),noire8k(25)
        INTEGER :: chsm,ndam1,iallcnt5,dummy1,ntest,nchr,iallocate,ios,
     *             latdegmin,latminmin,latdegmax,latminmax,londegmin,
     *             lonminmin,londegmax,lonminmax,i,j,n,ii,igridflg,
     *             nrvr1,ntmp,l,newformat
        real*4  :: sumclass,cintv,conv,classmax
        INTEGER   :: result1,n_hdr_lines
	  logical :: exists

!      DATA thms/0,211,176,148,108,70,42,14,11,0,0,0,0,0,0,0,0,0,0,0,0,
!     *          0,0,0,0/
      DATA gr10k/47,161,222,229,189,142,69,39,20,7,3,0,0,0,0,0,0,0,0,0,
     *            0,0,0,0,0/
!      DATA canag/45,155,204,210,172,143,100,48,4,0,0,0,0,0,0,0,0,0,0,0,
!     *            0,0,0,0,0/
!      DATA saug/0,191,216,259,278,297,303,109,58,0,0,0,0,0,0,0,0,0,0,0,
!     *          0,0,0,0,0/
!      DATA hmbr/199,288,619,643,823,726,565,458,305,302,122,0,0,0,0,0,
!     *            0,0,0,0,0,0,0,0,0/
!      DATA dffn/149,498,623,610,431,396,223,73,0,0,0,0,0,0,0,0,0,0,0,0,
!     *            0,0,0,0,0/
!      DATA simpson/0,26,114,336,431,467,451,588,578,440,597,527,410,
!     *             292,203,0,0,0,0,0,0,0,0,0,0/
!      DATA colum/0,1100, 837, 723, 453, 711, 914,1072,1683,1621,1944,
!     *           2102,2528,2385,2175,1834,2049,2384,2611,2936,3165,2987,
!     *           2468,2709,2506,3315,3425,3132,2762,3366,
!     *           3059,2217,2021,1371,1336,1291, 904, 576,1118,1272, 359,
!     *              0,   0,   0,   0,   0,   0,   0,   0,   0/
      DATA noire8k/0,0,82,179,299,213,117,288,263,222,
     * 206,169,129,44,0,0,0,0,0,0,0,0,0,0,0/
      DATA ntest/-10588/qstr/'copyright n.kouwen'/nchr/18/
	DATA iallcnt5/0/
      
      
!// Added by Dave

! parameter type definitions
	integer   unitNum, flnNum, iStat

! Local variables
	character(4096) line, subString, tmpString
	character(128) keyword, value
	integer lineLen, keyLen, wordCount, attCount
	logical rStat, lineType, foundEndHeader

	character(64) attribName
	integer ai, vi, xi, yi, attLen, error, rank
	real(4) val

! initialize default values
	CALL InitShedParam(header)	

! Set unit and fln number   ! now in argument list
!	unitNum = 31
!	flnNum = 1

	foundEndHeader = .false.

!// End Dave addition
      
      IF(.not.ID.GT.1)THEN      !  changed Jul. 12/04  nk
!	else
!     TS: NNOTE=100 -> DEFINED AS A PARAMETER IN AREA10A
!     TS - ALLOCATION OF AREA10A ARRAY (REMAINDER)
        if(iallcnt5.eq.0)then
          allocate(note(nnote),stat=iAllocate)
          if(iAllocate.ne.0) STOP
     *     'Error with allocation of area10a array in sheda'
          iallcnt5=1
	    endif
      endif

!     AL MUST BE IN INTEGER KILOMETERS
!     COORDINATES OF THE OUTSIDE ELEMENTS OF THE GRID.
!     ALL THE BASIN DATA IS READ ON FILE 9

!     rev. 9.9.35  Oct.  20/14  - NK: Added keyword & file checks
        inquire(FILE=fln(flnnum),EXIST=exists)
        if(.not.exists)then
          print*
          print*,'Fatal error'
          print*,'File not found ',fln(flnnum)(1:50)
          print*,'not found. Please check the event\event.evt file'
          print*,'for correct keyword and/or file name'
          print*
          stop 'Program aborted in read_shd_ef @ 1301'
        endif

!     basin/bsnm_shd.r2c
      open(unit=unitNum ,file=fln(flnNum) ,status='old',iostat=ios)
      if(ios.ne.0)then
        print*
        print*,'Problems in file',flnNum,fln(flnNum)(1:40)
        write(*,99162)fln(flnNum)(1:40)
        write(51,99162)fln(flnNum)(1:40)
99162   format(' Warning: Error opening or reading fln:',a30)
        print*,'iostat code =',ios
        STOP 'program aborted in read_shed_ef.for @ 130'
      else
          write(*,*)'Reading ',fln(flnNum)(1:50)
      endif

!// Added by Dave
	line(1:1) = '#'

	do WHILE((.NOT.foundEndHeader) .AND.
     &	    ((line(1:1) .eq. '#') .OR.
     &		(line(1:1) .eq. ':') .OR.
     &		(LEN_TRIM(line) .eq. 0))) 	

		read(UNIT=unitNum, FMT='((A))', iostat=ios) line	! read a line
		if(ios .eq. -1)then
			write(6,'((A))') 'ERROR: Premature EndOfFile encountered'
			STOP ' Stopped in read_shed_ef'
		end if

		rStat = Detab(line)				! replace tabs with spaces
		line = ADJUSTL(line)		! Get rid of leading white space
		lineLen = LEN_TRIM(line)		! Find the length excluding trailing spaces
          print*,line(1:lineLen)
		if(line(1:1) .eq. ':')then
			wordCount = SplitLine(line, keyword, subString)	! find the keyword
			rStat = ToLowerCase(keyword)
			KeyLen = LEN_TRIM(keyword)

			if(keyword(1:KeyLen) .eq. ':endheader')then
				foundEndHeader = .TRUE.

			else
				iStat = ParseShedParam(header,keyword,keyLen,
     &													subString)
				if(iStat .lt. 0) then
					write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
					write(*,'(2(A))') '   in line: ',line					
					STOP ' Stopped in read_shed_ef'
					return
				else if(iStat .eq. 0) then
!					write(*,'((A), (A))')  'Unrecognized keyword line: ',
!     &										line
				endif
			end if
		end if
	end do
	print*,'Found endheader'
      
      coordsys1=header%r2cp%csp%projection
      datum1=header%r2cp%csp%ellipsoid


	xcount=header%r2cp%xCount
	ycount=header%r2cp%yCount
	xorigin = header%r2cp%xOrigin
	yorigin = header%r2cp%yOrigin
	xdelta = header%r2cp%xDelta
	ydelta = header%r2cp%yDelta	



	al = header%nominalGridSize_AL
	cintv = header%contourInterval
	ntype = header%classCount-1
	classcount = header%classCount

!        every where loops are classcount
!        this needs to be changed  fix fix
!        to reading in the actual # land classes &
!        not no classes -1 for inpervious

	nrvr = header%numRiverClasses
	conv = header%r2cp%unitConv
	na   = header%totalNumOfGrids
	naa = header%numGridsInBasin
	nnprint = header%debugGridNo
      
      print*,'classcount =',classcount
      
	If(nnprint.gt.naa)then
	  print*,'debug grid is outside the watershed'
	  print*,'This can happen if:'
	  print*,'     You chose the wrong number or'
	  print*,'     You are using a subwatershed'
	  print*,'Please edit the shd file so'
	  print*, 'DebugGridNumber < NumGridsInBasin'
	  print*,nnprint,naa
	  print*
	  stop 'Program aborted in read_shd_ef @ 244'
	endif

d       print*,'al,cintv,classcount,nrvr,conv'
d       print*,al,cintv,classcount,nrvr,conv
d       print*,'na,naa,nnprint',na,naa,nnprint
d       print*

!     rev. 10.2.12 Dec.  30/17  - NK: Added frame headers to static r2c files incl. shd file   

	attCount =header%r2cp%ep%attCount 
	CALL LoadAttributeData(header%r2cp%ep, xCount,
     &						yCount, unitNum)

!// End Dave addition
      
      
!       fix fix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!       eventually, ycount & jnmax should be replaced everywhere in the code
c        ycount=ycount
c	  xcount=xcount
      
      print*,IsLatLong(header%r2cp%csp)

!       added ll separation Jul. 27/04  nk
	if(IsLatLong(header%r2cp%csp).eq.(.true.))then
!         lat-long coordinates
          iymin=int(yorigin*60.0)
          iymax=int((yorigin+ycount*ydelta)*60.0) 
          jxmin=int(xorigin*60.0)  
          jxmax=int((xorigin+xcount*xdelta)*60.0)  
	    llflg='y'            ! added Mar. 15/06 nk
          grde=xdelta*60.      ! grid size long in minutes
          grdn=ydelta*60.      ! grid size lat  in minutes
      else
!         utm or cartesian coordinates
          jxmin=int(xorigin/1000.)
          jxmax=jxmin+grde*(xcount-1)
          iymin=int(yorigin/1000.)
	    iymax=iymin+grdn*(ycount-1)
          grde=xdelta/1000.
          grdn=ydelta/1000.  
	    llflg='n'          ! added Mar. 15/06 nk
      endif

	astep=al/1000.
      istep=int(astep)         ! not used in computations
      step2=astep*astep        ! area of centre grid in km
      if(istep.lt.1)istep=1

      if(calling_program_flg.eq.'snw       ')return
      if(calling_program_flg.eq.'moist     ')return

      nastart=1
      naend=naa
      ib=1
      it=ycount-1
      
      
!     rl() and ch_length() are the same thing. ch_length used in bsn

!     rev. 10.2.07 Nov.  03/17  - NK: New rt_pond subroutine for channel pond routing 
!     add pondflg below      
c      if(iallcnt5.eq.1)then
      if(.not.allocated(s))then
       allocate(s(ycount,xcount),            !dummy(ycount,xcount),
     *    xxx(na),yyy(na),da(na),bnkfll(na),slope(na),elev(na),rl(na),
     *    ibn(na),sl1(na),sl2(na),irough(na),ichnl(na),next(na),
     *    ireach(na),frac(na),aclass(na,classcount),glacier_flag(na),
     *    flz(na),pwr(na),r1n(na),r2n(na),mndr(na),aa2(na),aa3(na),
     *    aa4(na),widep(na),theta(na),kcond(na),pool(na),rlake(na),
     *    pondflg(na),fpFactor(na),  
     *    flz2(na),pwr2(na),grid_area(na),flag1(na),     !!diverflg(na),
     *    pot(classcount),potfs(classcount),   ! moved from read_shed_ef may 15/07 nk
     *    stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *     'Error with allocation of area16a arrays in sheda'
!              glacier_flag(na)      added Mar, 28/06  nk
        iallcnt5=2
      endif

!     rev. 9.8.90  Oct.  30/13  - NK: Added fetch to the shd file 
      if(.not.allocated(fetch))then
        allocate(fetch(na,8),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *    'Error with allocation of Fetch vectors in read_shed @ 336'
        iallcnt5=2
      endif

	if(dds_flag.eq.1.or.numa.ge.1)then
        if(.not.allocated(s))allocate(penalty(na),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *     'Error with allocation of penalty in read_shed @ 332'
	endif


!  FOR READING IN DATA, HAVE TO USE UNFORM DATA BALOCKS
!  FOR COMPUTATIONS, CAN USE VARIABLE ROW LENGTHS
!  JL J(LEFT) AND JR J(RIGHT) DEFINE THE OUTER LIMITS OF A GIVEN ROW

      write(51,*)'The rank of each grid: highest =1'

!     rev. 9.3.03  Sep.  09/06  - NK: read s(i,j) from table instead of grid
!       valuse from this table are no longer used. Grid is for info only
!       and can be dropped for next version

c      do i=ycount,1,-1
c	  read(31,*)      ! pass over the data
c      end do

!     ENVIRONMENTAL IMPACT MODE - UP TO 16 CLASSES WITHIN ELMTS

!     REV. 8.99c  Oct.   5/99 -  IROUGH -> SL1 INPUT IN SHED
!     rev. 9.1.60  Jul.  27/04  - NK: reversed definitions for sl1 & sl2 Int. Slope
!       newformat file                   read data:
        Write(51,*)'shd file stuff:'
        read(31,3002)junk,junk20
        write(51,3002)junk,junk20
        read(31,3002)junk,junk20
        write(51,3002)junk,junk20
        print*,'na =',na
        noresv=0
        s=0   ! initialize all values in s(i,j)
        do n=1,na
          read(31,*,iostat=ios)i,next(n),yyy(n),xxx(n),da(n),
     *      bnkfll(n),slope(n),elev(n),rl(n),
     *      ibn(n),sl1(n),ichnl(n),ireach(n),grid_area(n),
     *      (aclass(n,ii),ii=1,classcount)
 8999 format(5x,4i5,3f10.0,2f7.0,i5,f10.5,2i5,3x,f15.0,<classCount>f6.3)
d            write(51,8999)n,next(n),yyy(n),xxx(n),da(n),
d    *      bnkfll(n),slope(n),elev(n),rl(n),
d    *      ibn(n),sl1(n),ichnl(n),ireach(n),grid_area(n),
d    *      (aclass(n,ii),ii=1,classcount)
          if(ios.ne.0)then
            print*,'problems reading ',fln(1)
            print*,'for grid no ',n,'/',na
            print*
d            write(51,8999)n,next(n),yyy(n),xxx(n),da(n),
d    *      bnkfll(n),slope(n),elev(n),rl(n),
d    *      ibn(n),sl1(n),ichnl(n),ireach(n),grid_area(n),
d    *      (aclass(n,ii),ii=1,classcount)
            if(next(n).eq.0)exit
c            stop 'Program aborted in rdshed @ 675'
          endif
!     rev. 9.3.03  Sep.  09/06  - NK: read s(i,j) from table instead of grid
!         to retrieve grid # if i,j are known
	    s(yyy(n),xxx(n))=n
!          print*,n,yyy(n),xxx(n),s(yyy(n),xxx(n))
          frac(n)=grid_area(n)/al/al
          noresv=max(noresv,ireach(n))
        end do
        print*,'Finished reading the watershed file'

!     set the debug grid and class numbers or location
	if(nnprint.eq.0)then
!       for old format
        nnprint=s(ipr,jpr)
      elseif(nnprint.gt.naa)then
        nnprint=naa/2
        ipr=yyy(nnprint)   ! probably not needed anywhere
        jpr=xxx(nnprint)
      else
!       this can happen when a sub-watershed is used
!       for new format
        ipr=yyy(nnprint)   ! probably not needed anywhere
        jpr=xxx(nnprint)
      endif

      iiprint=1
	classmax=0.0
	do ii=1,classcount
	  if(aclass(nnprint,ii).gt.classmax)iiprint=ii
	  classmax=amax1(aclass(nnprint,ii),classmax)
	end do

      write(51,*)
      write(51,5000)nnprint,ipr,jpr
	write(51,5004)iiclass
      write(51,*)

!     for -ve ipr & jpr, compute only on the debug grid & ignore the rest
      if(iopt.ge.1)then
        if(ipr.lt.0.and.jpr.le.0)then
  	    ipr=abs(ipr)
  	    jpr=abs(jpr)
          do n=1,na
      	    if(yyy(n).eq.ipr.and.xxx(n).eq.jpr)then
              nastart=n
              naend=n
            endif
	    end do
        endif
      endif

!     REV. 8.92 - Dec.  24/89 -  CHECK FOR 100% ACLASS COVERAGE
      igridflg=0
      do n=1,naa
        sumclass=0.0
        do ii=1,classcount
          sumclass=sumclass+aclass(n,ii)
        end do
        if(sumclass.ne.1.0)then
          igridflg=1
          write(98,9023)n,yyy(n),xxx(n),sumclass
          do ii=1,classcount
            if(sumclass.gt.0.0)then
              aclass(n,ii)=aclass(n,ii)/sumclass
            else
              write(98,9024)n,yyy(n),xxx(n)
            endif
          end do
        end if
      end do

!     FOR SPLD ONLY:
!     * * * * * * * * * * * * * * * * *
!     include 'demosize.for'

!     Checking data:

      nrvr1=0
      do n=1,naa   
!       moved this to flowinit at one point but then it got 
!       recalculated with each iteration when optimizing.
!       A serious snafu resulting in the convergence problem on opt.

        if(slope(n).lt.0.0)then	
          print*,'In rdshed reading the file :',fln(1)
          print*,'The slope in grid no ',n,' is ',slope(n)
          print*,'Please check the elevations in the map file '
          print*,'or change the slope value in the shd file'
          print*,'The former is recommended as the permanent solution'
          print*
          stop 'Program aborted in rdshed @ 756'
        endif
        slope(n)=sqrt(slope(n))
        elev(n)=0.01*elev(n)
!       check to see how many basins/river classes there are:
        nrvr1=max(nrvr1,ibn(n))
      end do

      if( nrvr.ne.nrvr1)then
!      if(nrvr.lt.nrvr1)then
        print*,'no of river classes specifies in .shd header =',nrvr
        print*,'no found in .shd grid data =',nrvr1
        print*,'the header should = no found in gridded data'
        print*,'the no in the header should be ',nrvr1 
        print*,'Ignore this warning if you are running a sub-watershed'
        print*,' program paused in shed @ 473'
        read(*,*)
      endif
          
      if(nrvr.le.0.or.nrvr.gt.16)then
        print*,' In the .shd file: ',fln(1)
        print*,' Column 16-20 should have the value for nrvr'
        print*,' the number of river classes '
        print*,' It must have a value 1-16'
        print*,' The number of river/basin types should match'
        print*,' in the .par and .shd files'
        print*,' Please enter the correct number in the'
        print*,' .par and .shd files. To get the correct number,'
        print*,' create a new .shd file using bsn.exe and then match'
        print*,' it in the .par file - making sure there are an '
        print*,' equal number of parameter sets for rivers'
        print*
        stop ' Program aborted in rdshd @ 729'
      endif


!  taken out because of problems in snw.exe  nk oct 11/05
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     remove for unix
!      call keychk(qstr,nchr,result1)
!      call userchk(ntest,result1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

      close(unit=31)
!     WE WILL USE THIS UNIT NUMBER AGAIN FOR THE DAMAGE SITE FILE

!      if(iopt.gt.0)then
	   write(51,*)' Note: order not the same as the .shd file YET'
         write(51,6006)
         write(51,6007)

	   write(51,*)' Note: order not the same as the .shd file YET'

         if(classcount.le.0)then
            do n=1,naa
               write(51,6004)n,yyy(n),xxx(n),da(n),bnkfll(n),
     *         slope(n)**2,
     *         elev(n),ibn(n),sl2(n),ichnl(n),next(n),ireach(n),
     *	       grid_area(n),(aclass(n,ii),ii=1,classcount)
            end do
         else
            do n=1,naa
c               write(51,6004)n,yyy(n),xxx(n),da(n),bnkfll(n),
c     *         slope(n)**2,elev(n),ibn(n),sl1(n),
c     *         ichnl(n),next(n),ireach(n),grid_area(n),
c     *         (aclass(n,ii),ii=1,classcount)
            end do
         endif
!      endif

!     section moved from spl
!     REV. 8.99c  Oct.   5/99 -  IROUGH -> SL2 INPUT IN SHED
      if(newformat.ne.1)then
        do n=1,naa

!		  changed sl1(n) to the internal slope. Was sl2(n)  nk jul 27/04
!         rev. 9.1.60  Jul.  27/04  - NK: reversed definitions for sl1 & sl2 Int. Slope

          sl1(n)=cintv*(float(irough(n))+.0001)/al


!         REVISED 99/09/23 TO SET MAX ON SL1()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         rev. 9.1.57  Jul.  06/04  - NK: Fixed major bug in shed.for max instead of min
!         rev. 9.1.57  Jul.  06/04  - NK: Fixed major bug in shed.for max instead of min
!         sl1(n)=amin1(1.0,sl1(n))    << don't use this RR & GR  July/04
!!!           sl1(n)=amax1(1.0,sl1(n))     !  <<< so wrong
!         with this, the flows are way peakier. Way too much interflow
!         and not nearly enough drng.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        end do
      endif

!     WRITE THE MAP INFORMATION TO THE /SPL/SIMOUT/PIC.LST FILE:
      write(56,9005)na,xcount,ycount
      do n=1,na
         write(56,9005)n,yyy(n),xxx(n),next(n)
      end do





!     TS - ALLOCATIONS OF AREA4A ARRAYS
!     classcount for the number of land cover classes
!     nrivertype for the number of channel or basin types
!     moved here from spl9  nk 06/07/00
!     then moved from rdpar nk 28/12/04
!     moved back to rdpar   nk 27/07/06 because needed for bsn.for
!     parameter allocation moved to rdpar  27/07/06 nk  - needed by bsn
!     no good. Moved back here. fix bsn some other way


!      allocate(mndr(nrvr),
!     *r1(nrvr),r2(nrvr),r2low(nrvr),r2hgh(nrvr),r2dlt(nrvr),
!     *aa2(nrvr),aa3(nrvr),aa4(nrvr),rivtype(nrvr),
!     *theta(nrvr),thetadlt(nrvr),thetalow(nrvr),thetahgh(nrvr),
!     *widep(nrvr),widepdlt(nrvr),wideplow(nrvr),widephgh(nrvr),
!     *kcond(nrvr),kconddlt(nrvr),kcondlow(nrvr),kcondhgh(nrvr),
!     *flz(nrvr),flz2(nrvr),flzlow(nrvr),flzhgh(nrvr),flzdlt(nrvr),
!     *pwr(nrvr),pwr2(nrvr),pwrlow(nrvr),pwrhgh(nrvr),pwrdlt(nrvr),
!     *stat=iAllocate)
!      if(iAllocate.ne.0) STOP
!     *   'Error with allocation of area4a arrays in rdshed @ 968'
      
      
      allocate(r1(na),r2(na),stat=iAllocate)
      if(iAllocate.ne.0) STOP
     *   'Error with allocation of r1 & r2 arrays in rdshed @ 968'


!          allocate(r1n(nrvr),r2n(nrvr),
!     *       r2nlow(nrvr),r2nhgh(nrvr),r2ndlt(nrvr),
!     *       stat=iAllocate)
!          if(iAllocate.ne.0)STOP
!     *   'Error with allocation of area4 arrays in rdshed @ 975'



! * * * *  TS  * * * * *
!     TS - ALLOCATIONS OF AREAWET ARRAYS
      if(.not.allocated(wetwid))then
        allocate(wetwid(na),chawid(na),chadep(na),
     *  wstore1(na),wstore2(na),
     *  wcap(na),flowxa(na),chaxa(na),chflxa(na),obflxa(na),totwid(na),
     *  chflfrac(na),obflfrac(na),satxa(na),wetxa(na),hcha1(na),
     *  hcha2(na),hfp(na),hwet1(na),hwet2(na),qin(na),qswevp(na),
     *  qswrain(na),qiwet1(na),qiwet2(na),qowet1(na),qowet2(na),
     *  wetarea(na),chaarea(na),bin_precip(na),wsat(na),wetfrac(na),
     *  qiob1(na),qiob2(na),qoob1(na),qoob2(na),   !rev 9.5
     *  qUS1(na),qUS2(na),  !     REV. 10.1.24 Jan.  30/16  - NK: Added qUS1 & qUS2 for watbal
     *  qiwetsum(na),qowetsum(na),wstoreinit(na),resvstore1(na),  !rev 9.5.1
     *  resvstore2(na),store_error_flag(na),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *   'Warning: error with allocation of  arrays in read_shd @ 676'
	endif
! * * * * * * * * * * *



!     TS - ALLOCATIONS OF AREA1A ARRAYS

!     *xxx(na),yyy(na),da(na),bnkfll(na),slope(na),elev(na),
!     *ibn(na),irough(na),ichnl(na),next(na),ireach(na),frac(na),
!     *aclass(na,classcount),ch_length(na),sl1(na),rl(na),
!     *pot(classcount),potfs(classcount),qlz(na),psmear(na),punused(na),
      if(.not.allocated(qi1))then
        allocate(
     *  qi1(na),qi2(na),qo1(na),qo2(na),qr(na),
!     rev. 10.2.44 Jan.  21/19  - NK: Added qOld() to allow previous weighted outflow on the resume file
     *  qOld(na),        
     *  d2(na),qda(na),cap(na),over(na),
     *  qmax(na),res(na),netinflow(na),netoutflow(na),totwetl(na),
     *  sump(na),store1(na),store2(na),att(na),
     *  qbase(na),   !nreach(na), dim for nreach changed in flowinit
     *  rf(na,classcount),rffs(na,classcount),
!     rev. 10.4.21 Apr.  21/20  = NK Add UZS deficit to wfo file = UZS(class=classcount)
     *  r(na,classcount),effpor(na,classcount),totuzsdef(na),
     *  v(na,classcount),totd1(na),totuzs(na),totsnw(na),qstream(na),
     *  totchnl(na),totgrid(na),storinit(na),d1(na,classcount),
     *  d1fs(na,classcount),uzs(na,classcount),uzsfs(na,classcount),
     *  lzs(na),sumf(na,classcount),sumrechrg(na),
     *  sumffs(na,classcount),qlz(na),sr(classcount),
     *  x4(classcount),x5(classcount),q1(na,classcount),
     *  q1fs(na,classcount),qint(na,classcount),qintfs(na,classcount),
     *  fake(classcount),fakefs(classcount),qdrng(na),qdrngfs(na),
     *  drng(na,classcount),drngfs(na,classcount),sq1(classcount),
     *  sq1fs(classcount),sqint(classcount),sqintfs(classcount),
     *  sdrng(classcount),sdrngfs(classcount),
     *  sexcess(classcount),qstrm(na),
     *  sumq1(na),sumqint(na),sumq1fs(na),sumqintfs(na),
     *  strloss(na),rechrg(na),sumrff(na),
     *  ice_fctr(na),ice_fctr_min(na),ice_fctr_max(n),
     *  stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *  'Warning: error with allocation of  arrays in read_shd @ 710'
      endif
!     rev. 10.1.98 Oct   04/17  - NK: Deal with -ve flows in route
      if(.not.allocated(qqlow))
     *                    allocate(qqlow(na),stat=iAllocate)
      if(iAllocate.ne.0) STOP
     *  'Warning: error with allocation of  qqlow in read_shd @ 710'

c      allocate(rechrg(na),sumrff(na),stat=iAllocate)

!     TS - ALLOCATIONS OF AREAETA ARRAYS (PARTIAL)
!     RAD ALLOCATED IN SHEDA.FOR
!     rev. 9.1.80  Mar.  31/05  - NK: added sublimation   (sublim)
c       allocate(strloss(na),stat=iAllocate)
c      if(iAllocate.ne.0) STOP
c     *   'Warning: error with allocation of areaeta arrays in spl9'

!     OCT30/03 TS:  ADDED SUMQI,SUMQINT,SUMQ1FS AND SUMQINTFS ALLOCATIONS

!     JUN28/06 TS: ADDED DF, DFFS FOR ISOTOPE ROUTINES
      if(.not.allocated(df))then
        allocate(df(na,classcount),dffs(na,classcount),
     *       qdrng2(na,classcount),qdrngfs2(na,classcount),
     *       qdf(na,classcount),qdffs(na,classcount),
     *       dprecip(na,classcount),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *   'Warning: error with allocation of arrays in read_shd @ 731'
	endif

!     rev. 9.6.05  Apr.  06/10  - NK: added store_error_flag for -ve storage grids
	do n=1,naa
	  store_error_flag(n)=.false.
	end do

      return


 9998 ntmp=0
      ndam=0
      write(98,9025)
      write(98,9026)
      close(unit=34)
      RETURN

 9999 ndam=0
      write(98,9026)
      close(unit=34)
      RETURN



99901 write(*,99902)fln(4)
99902 format(' file',a30,' not found for unit 34 - check event file')
      STOP 'program stopped in shed.for at 99902'

99910 write(*,99911)fln(4)
99911 format(' no data found or problems with data in ',a30)
      STOP 'program stopped in shed at 99911'

! FORMATS

 1000 format(i5)
 1002 format(' ',i5,'stream gage locations have been passed over')
 1003 format(' ',i5,'reservoir locations have been passed over')
 1004 format(' ',i5,'damage sites:')
 1098 format(a80)
 1099 format(2i5,1x,a12,7x,4e10.3,f10.3)
 1100 format(' ',2i5,1x,a12,7x,4e10.3,f10.3)
 1101 format(' reading the stream gauge location file: ',a30)
 1102 format(' ',2i5,1x,a12,7x,4e10.3,f10.3/)
 1776 format(' ','l,iys(l),jxs(l)',5i5)
        write(51,*)
 3001 format(a20,i12)
 3002 format(2a20)
 3003 format(a20,f12.3)
 3004 format(a20,a10)

 5000 format(' Debug grid reset to grid number n,row col',3i7)
 5001 format(' ',i5,10f12.3,f12.7)
 5002 format('     n          da         cap       chaxa      chadep',
     *'      chawid      wetwid        wcap       wetxa  %wet class')     
 5003 format(' Channel properties:' )   
 5004 format(' Debug class set to class number ',i7)
5007  format(a1,a79)
 6004 format(1x,3i5,2f10.1,f10.5,f7.0,i5,f10.5,3i5,f12.0,99f5.2)
 6006 format(2x,'basin file:')
 6007 format(4x,'n   yy   xx      da       cap       slope    elv',
     *'     ibn   sl2    ich  next reach gridarea imp area & fractions')
 6010 format(12i5)
 6014 format(i5)
 6020 format(/' Number of river type or basins found =',i5/)


 9000 format(5x,2i5,3f10.5,f7.0,5i5,17f5.2)
 9001 format(999i4)
 9002 format(22x,f10.0)
 9003 format(3i5,f7.1,f10.0)
 9004 format(12i5,2f5.1)
 9005 format(12i5,2f5.0)
 9006 format(8i10)
 9007 format('   iymin     iymax     jxmin     jxmax in minutes')
 9008 format(999i5)
 9021 format(f10.0,f5.0,2i5,15x,17i1)
 9022 format(' al,astep,cintv,nrvr,classcount/ ',f10.0,f8.3,f5.0,
     * 2i5,2x,17i1)
 9023 format(' Warning: area correction in grid(n,i,j)',3i5,f9.5)
 9024 format(' Warning: total area = 0.0 for grid(n,i,j)',3i5)
 9025 format(' Warning: no reservoirs or lakes in bsnm.str file')
 9026 format(' Warning: no damage sites in bsnm.str file')

      RETURN

      END SUBROUTINE rdshed

