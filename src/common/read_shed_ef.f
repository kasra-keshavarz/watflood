      SUBROUTINE read_shed_ef(unitNum,flnNum)

!***********************************************************************
!    Copyright (C) 1987-2018 by Nicholas Kouwen and Dave  
        
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
! ****NOTE****This subroutine works but still requires extensive cleaning up
! **** Dave Watson is still working on it
!***********************************************************************!

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ~~~~NOTE:  This s/r was modified to read the gridded parameter file  
! ~~~~as well as the shed file for watroute by NK  Jul. 5/07
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!***********************************************************************
!     written Mar/06 by Dave Watson
!     derived from rdshed written by Nick Kouwen
!     This subroutine reads the ensim compatible shd file (r2c format).
!	This file is generated by a separate program called "BASIN" (by Nick Kouwen)
!***********************************************************************!
!f
!  ntype    - max number of permeability classes in an element
!  classcount - number of land cover classes  = ntype + 1
!  al       - grid length in m.
!  step     - grid length in km but if al<1000, it is equal to 1
!  na       - total # of squares in the watershed including extra one in
!             which water runs. e.g.lake
!  xxx(n)    - the north-south coordinate of the n'th element
!  yyy(n)   - the east-west coordinate of the n'th element
!  order    - coordinate file,highest to lowest element.
!  da(n)    - drainage element (n).
!  ibn(n)   - river class number (no longer basin number)
!  irough(n)- over land flow conveyance class 1-5.
!  ichnl(n) - #channels/square class1-59
!  aimp(n)  - q/age impermeable area.from.01 to .9
!  next(n)  - the element number into which the n'th element drains
!  s(i,j)   - contains the order number for each square (used to be 
!             drainage direction)
!           - impevious area should never be 1.0 or 0.0.

!***********************************************************************

c     rev. 9.5.20  Mar.  06/08  - NK: added resvstore for iso model
!     rev. 9.7.24  Apr.  20/11  - NK: Added diverflg to indicate if a diversion is in grid
!     rev. 10.1.07 Dec.  02/15  - NK: Added ice_fctr(n) to route 
!     rev. 10.1.11 Dec.  11/15  - NK: Revised ice factor initialization and calculation   

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
c        CHARACTER(79) :: hdrcomment(100)
        INTEGER :: chksum(160),thms(25),gr10k(25),saug(25),hmbr(25),
     *             dffn(25),simpson(25),colum(50),canag(25),noire8k(25)
        INTEGER :: chsm,ndam1,iallcnt5,dummy1,ntest,nchr,iallocate,ios,
     *             latdegmin,latminmin,latdegmax,latminmax,londegmin,
     *             lonminmin,londegmax,lonminmax,i,j,n,ii,igridflg,
     *             nrvr1,ntmp,l,newformat
        real*4  :: sumclass,cintv,conv
        INTEGER(kind=2) :: result1,n_hdr_lines,zone
	  logical :: exists


      real*4, dimension(:,:),   allocatable :: dummy  ! now in area17

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
	integer*4 unitNum, flnNum, iStat

! Local variables
	character*4096 line, subString, tmpString
	character*128 keyword, value
	integer lineLen, keyLen, wordCount, attCount
	logical rStat, lineType, foundEndHeader

	character *64 attribName
	integer ai, vi, xi, yi, attLen, error, rank
	real*4 val

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

!      if(iopt.eq.2)print*,'opening unitnum',unitnum
!	if(iopt.eq.2)print*,'filename =',fln(flnnum)
d      print*,'opening unitnum',unitnum
d	print*,'filename =',fln(flnnum)(1:50)



!     rev. 9.9.35  Oct.  20/14  - NK: Added keyword & file checks
        inquire(FILE=fln(flnnum),EXIST=exists)
        if(.not.exists)then
          print*
          print*,'Fatal error -  file not found'
          print*,fln(flnnum)(1:60)
          print*,'not found. Please check the event\event.evt file'
          print*,'for correct keyword and/or file name'
          print*
          stop 'Program aborted in read_shd_ef @ 1301'
        endif
!       check to see that it is an r2c file        
		if(.not.IsFileTypeR2C(fln(flnnum))) then
                 print*,'Old format shd files not accepted'
                 print*,'Please create EF ????_shd.r2c files & rerun'
!	           stop 'Program aborted in read_shd_ef @ 994'
		endif

!     basin/bsnm_shd.r2c
      open(unit=unitNum ,file=fln(flnNum) ,status='old',iostat=ios)
      if(ios.ne.0)then
        print*
        print*,'Problems in file',flnNum,fln(flnNum)(1:40)
        write(*,99162,iostat=ios)fln(flnNum)(1:40)
        write(51,99162,iostat=ios)fln(flnNum)(1:40)
99162   format(' Warning: Error opening or reading fln:',a30)
        print*,'iostat code =',ios
        STOP 'program aborted in read_shed_ef.for @ 130'
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
	
      coordsys1=header%r2cp%csp%projection
      datum1=header%r2cp%csp%ellipsoid

!     rev. 9.9.62  Mar.  21/15  - NK: Change zone fom character to integer
c!     ugly but effective
c	open(unit=99,file='junk',status='unknown')
c	write(99,99000)header%r2cp%csp%zone
c99000 format(i10)
c	rewind 99
c	read(99,99001)zone1
c99001 format(a10)
c      close(unit=99,status='delete')

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
	na = header%totalNumOfGrids
	naa = header%numGridsInBasin
	nnprint = header%debugGridNo

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

d     if(iopt.eq.2)then
d       print*,'al,cintv,classcount,nrvr,conv'
d       print*,al,cintv,classcount,nrvr,conv
d       print*,'na,naa,nnprint',na,naa,nnprint
d       print*
d     endif

!     rev. 10.2.12 Dec.  30/17  - NK: Added frame headers to static r2c files incl. shd file   

	attCount =header%r2cp%ep%attCount 
	CALL LoadAttributeData(header%r2cp%ep, xCount,
     &						yCount, unitNum)

!// End Dave addition

!       fix fix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!       eventually, ycount & jnmax should be replaced everywhere in the code
c        ycount=ycount
c	  xcount=xcount

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

!     rev. 9.3.04  Oct.  24/06  - NK: routing parameters dim to na in rte
!     rev. 9.5.06  Feb.  05/08  - NK: added pool and pool_o in rdpar & route

!     NOTE:  depending on whether this s/r is for rte or spl, dimension 
!            routing pars differently.  na for rte & nrvr for spl
!            Now we dim all for na for spl & rte so mem is wasted in spl
!            But it is so read_shed_ef is the same for all.
!            Probably not a good idea to change for spl because of opt.

!     rev. 10.2.07 Nov.  03/17  - NK: New rt_pond subroutine for channel pond routing 
!     add pondflg below      
c      if(iallcnt5.eq.1)then
      if(.not.allocated(s))then
       allocate(s(ycount,xcount),dummy(ycount,xcount),
     *    xxx(na),yyy(na),da(na),bnkfll(na),slope(na),elev(na),rl(na),
     *    ibn(na),sl1(na),sl2(na),irough(na),ichnl(na),next(na),
     *    ireach(na),frac(na),aclass(na,classcount),glacier_flag(na),
     *    flz(na),pwr(na),r1n(na),r2n(na),mndr(na),aa2(na),aa3(na),
     *    aa4(na),widep(na),theta(na),kcond(na),pool(na),rlake(na),
     *    pondflg(na),fpFactor(na),  
     *    flz2(na),pwr2(na),grid_area(na),nca_1d(na),flag1(na),     !!diverflg(na),
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

c!     rev. 9.7.24  Apr.  20/11  - NK: Added diverflg to indicate if a diversion is in grid
c      do n=1,na
c        diverflg(n)=.false.
c      end do
      
!     initialize fetch
      do n=1,na
        do ii=1,8
          fetch(n,ii)=0.0
        end do
      end do

!// Added by Dave	
!// First find and copy rank attribute first data over to global array
	do ai=1,attCount-classcount             !<<<<<<<<<<<<<<<<<<<fix fix
		vi = 0
		attribName = header%r2cp%ep%attList(ai)%name
		rStat = ToLowerCase(attribName)
		attLen = LEN_TRIM(attribName)
		if(attribName(1:attLen) .eq. 'rank')then 
			do yi=1,yCount
				do xi=1,xCount
					vi = vi+1
					val = header%r2cp%ep%attList(ai)%val(vi)
					s(yi,xi) = val
				end do
			end do
		end if
	end do


!// Copy attribute data (not classes yet) over to global attributes
	vi = 0
	do yi=1,yCount
		do xi=1,xCount
			vi = vi+1
			rank = s(yi,xi)
			if(rank.gt.0) then
			 do ai=1,attCount-classcount
			 
			  attribName = header%r2cp%ep%attList(ai)%name
			  rStat = ToLowerCase(attribName)
			  attLen = LEN_TRIM(attribName)
			  val = header%r2cp%ep%attList(ai)%val(vi)
			  if(attribName(1:attLen) .eq. 'next')then
				next(rank) = int(val)
			  else if(attribName(1:attLen) .eq. 'da')then
				da(rank) = val
			  else if(attribName(1:attLen).eq.'bankfull')then
				bnkfll(rank) = val
			  else if(attribName(1:attLen).eq.'chnlslope')then
				slope(rank) = val
			  else if(attribName(1:attLen) .eq. 'elev')then
				elev(rank) = val
			  else if(attribName(1:attLen) .eq. 'chnllength')then
				rl(rank) = val
			  else if(attribName(1:attLen) .eq. 'iak')then
				ibn(rank) = val
			  else if(attribName(1:attLen) .eq. 'intslope')then 
				sl1(rank) = val
			  else if(attribName(1:attLen) .eq. 'chnl')then
				ichnl(rank) = val
			  else if(attribName(1:attLen) .eq. 'reach')then
				ireach(rank) = val

!      attributes added by nk  Oct. 1/06  for watroute par file
			  else if(attribName(1:attLen) .eq. 'flz')then
				flz(rank) = val
			  else if(attribName(1:attLen) .eq. 'pwr')then
				pwr(rank) = val
			  else if(attribName(1:attLen) .eq. 'r1n')then
				r1n(rank) = val
			  else if(attribName(1:attLen) .eq. 'r2n')then
				r2n(rank) = val
			  else if(attribName(1:attLen) .eq. 'mndr')then
				mndr(rank) = val
			  else if(attribName(1:attLen) .eq. 'aa2')then
				aa2(rank) = val
			  else if(attribName(1:attLen) .eq. 'aa3')then
				aa3(rank) = val
			  else if(attribName(1:attLen) .eq. 'aa4')then
				aa4(rank) = val
			  else if(attribName(1:attLen) .eq. 'widep')then
				widep(rank) = val
!               end attributes added by nk

			  else if(attribName(1:attLen) .eq. 'gridarea')then
				grid_area(rank) = val
	                  frac(rank)=grid_area(rank)/al/al
!               frac is still used in the code but no longer in the shed file   
			  else if(attribName(1:attLen) .eq. 'frac')then
!				frac(rank) = val
                print*
	          print*,'Error: old format shd file found'
	          print*,'Please create a new bsnm_shd.r2c file using the'
	          print*,'current version of bsn.exe'
	          print*
	          stop 'Program aborted in read_shed_ef.for @ 447'
!               impervious area no longer used  nk
!   		  	  else if(attribName(1:attLen) .eq. 'imperv')then
!					aclass(rank, classcount) = val
!             end no longer used

!             Attributes added by NK. Oct. 2013 for Fetch
!     rev. 9.8.90  Oct.  30/13  - NK: Added fetch to the shd file 
			  else if(attribName(1:attLen) .eq. 'fetchne')then
				fetch(rank,1) = val
			  else if(attribName(1:attLen) .eq. 'fetche')then
				fetch(rank,2) = val
			  else if(attribName(1:attLen) .eq. 'fetchse')then
				fetch(rank,3) = val
			  else if(attribName(1:attLen) .eq. 'fetchs')then
				fetch(rank,4) = val
			  else if(attribName(1:attLen) .eq. 'fetchsw')then
				fetch(rank,5) = val
			  else if(attribName(1:attLen) .eq. 'fetchw')then
				fetch(rank,6) = val
			  else if(attribName(1:attLen) .eq. 'fetchnw')then
				fetch(rank,7) = val
			  else if(attribName(1:attLen) .eq. 'fetchn')then
				fetch(rank,8) = val
				
!             End Attributes added by NK. Oct. 2013 for Fetch

			  end if
			 end do
			endif
		end do
	end do	

      manningflg='y'

d     print*,'attribute count =',attcount
d     print*,'class count =    ',classcount

!// Copy class attribute data over to global attributes
	vi = 0
	do yi=1,yCount
		do xi=1,xCount
			vi = vi+1
			rank = s(yi,xi)
			if(rank.gt.0) then
				do ai=attCount-classcount+1,attCount  
				  val = header%r2cp%ep%attList(ai)%val(vi)
					aclass(rank,ai-(attCount-classcount)) = val
				enddo
			endif
		end do
	end do
      	
!     pick the dominant class in this grid for debug printing
      iiprint=1
	do ii=2,classcount
	  if(aclass(nnprint,ii).gt.aclass(nnprint,iiprint))iiprint=ii
      end do
        
!// Copy rows and col over to global attributes
!// I'm not sure if we need this...check with Nick
	vi = 0
	do yi=1,yCount
		do xi=1,xCount
			vi = vi+1
			rank = s(yi,xi)
			if(rank.gt.0) then
				xxx(rank) = xi
				yyy(rank) = yi
			endif
		end do
	end do	

      do i=1,ycount
	  do j=1,xcount
	    dummy(i,j)=0.0
	  end do
	end do

!     Write the land cover class fractions to spl.txt
      if(iopt99)then
        write(51,50002,iostat=ios)
        do ii=1,classcount
	    write(51,50001,iostat=ios)ii
          do i=1,ycount
	      do j=1,xcount
              n=s(i,j)
              if(n.gt.0)dummy(i,j)=aclass(n,ii) 
            end do
            write(51,50000,iostat=ios)(dummy(i,j),j=1,xcount)
	    end do
	  end do
	endif
50000 format(999f6.2)
50001 format('class no=',i5)
50002 format('Land Cover class fractions:')
 
!// Deallocate the attribute data now that global attributes have been set
	do ai=1,attCount
		deallocate ( header%r2cp%ep%attList(ai)%val, STAT = error )
		if (error.ne.0) STOP 'deallocation error in read_gsm_ef()' 
	end do

!// End Dave addition

!     GRAND RIVER:
      if(iymin.eq.4790.and.jxmin.eq.500)then
         do i=1,25
          chksum(i)=gr10k(i)
        end do
      end if

      ichsm=1     
!     if one is no good, kill run
      do i=1,min(25,ycount)
        chsm=0
        do j=1,xcount
          chsm=chsm+s(i,j)
        end do
!       print*,' i,shksun(i),chsm /',i,chksum(i),chsm
        if(chsm.ne.chksum(i))then
          ichsm=0        ! program will abort
        end if
      end do

!	close this loophole - Oct. 15/03 NK
!     when the top left hand corner of the grid are all zeros
!     the checksums are all zero and the program would run.
      chsm=0
	do i=1,min(25,ycount)
        do j=1,xcount
          chsm=chsm+s(i,j)
        end do
      end do

      if(chsm.eq.0)then
        ichsm=0        ! program will abort
      end if
     
!     rev. 9.1.46  Jul.  17/03  - WATFLOOD LITE incorporated 
!     WATFLOOD LT
!     WATFLOOD LT
      if(ycount.le.7.and.xcount.le.7.and.na.le.15.and.al.le.2000.0)then
!       when ever these conditions are met, the program will run
!       Since the gr10k files exceed these bounds, the messages will
!       appropriate.
        ichsm=3
      endif

	if(nnprint.gt.naa)then
        nnprint=naa/2
        ipr=yyy(nnprint)   ! probably not needed anywhere
        jpr=xxx(nnprint)
      else
!       this can happen when a sub-watershed is used
!       for new format
        ipr=yyy(nnprint)   ! probably not needed anywhere
        jpr=xxx(nnprint)
      endif
      if(iopt.gt.1)write(51,*,iostat=ios)
      if(iopt.gt.1)write(51,5000,iostat=ios)nnprint,ipr,jpr

!     REV. 8.92 - Dec.  24/89 -  CHECK FOR 100% ACLASS COVERAGE
      igridflg=0
      do n=1,naa
        sumclass=0.0
        do ii=1,classcount
          sumclass=sumclass+aclass(n,ii)
        end do
        if(sumclass.ne.1.0)then
          igridflg=1
          if(iopt.gt.1)write(51,9023,iostat=ios)n,yyy(n),xxx(n),sumclass
          do ii=1,classcount
            if(sumclass.gt.0.0)then
              aclass(n,ii)=aclass(n,ii)/sumclass
            else
              if(iopt.gt.1)write(51,9024,iostat=ios)n,yyy(n),xxx(n)
            endif
          end do
        end if
      end do

!     Checking data:
      nrvr1=0
      do n=1,na   
!       moved this to flowinit at one point but then it got 
!       recalculated with each iteration when optimizing.
!       A serious snafu resulting in the convergence problem on opt.

        if(slope(n).lt.0.0)then	
          print*,'In read_shed_ef reading the file :',fln(flnNum)
          print*,'The slope in grid no ',n,' is ',slope(n)
          print*,'row',yyy(n),' col',xxx(n)
          print*,'Please check the elevations in the map file '
          print*,'or change the slope value in the shd file'
          print*,'The former is recommended as the permanent solution'
          print*
          stop 'Program aborted in read_shed_ef @ 756'
        endif
        slope(n)=sqrt(slope(n))
	  sl2(n)=sqrt(sl1(n))

!     rev. 9.5.65  Sep.  26/09  - NK: lapse rate changed from dC per 100 m to dC per m
!        originally lapse rate was given in dC  per 100 m 
!        in the par file. 
c        elev(n)=0.01*elev(n)

!       check to see how many basins/river classes there are:
        nrvr1=max(nrvr1,ibn(n))
      end do

      close(unit=unitNum)
!     WE WILL USE THIS UNIT NUMBER AGAIN FOR THE DAMAGE SITE FILE

      if(iopt.gt.1)then
      write(51,*,iostat=ios)'Grid properties:'
      write(51,6006,iostat=ios)
      write(51,6007,iostat=ios)
      do n=1,naa
        write(51,6004,iostat=ios)n,yyy(n),xxx(n),da(n),bnkfll(n),
     *        slope(n)**2,elev(n),ibn(n),sl1(n),
     *        ichnl(n),next(n),ireach(n),frac(n),
     *        (aclass(n,ii),ii=1,classcount)
      end do
      do n=1,naa
        write(51,*,iostat=ios)n,yyy(n),xxx(n),(fetch(n,i)/1000.0,i=1,8)
      end do
      end if


!     WRITE THE MAP INFORMATION TO THE /SPL/SIMOUT/PIC.LST FILE:
!      write(56,9005)na,xcount,ycount
!      do n=1,na
!         write(56,9005)n,yyy(n),xxx(n),next(n)
!      end do
      
      if(.not.allocated(r1))then
        allocate(r1(na),r2(na),stat=iAllocate)
      endif
      if(iAllocate.ne.0) STOP
     *   'Error with allocation of r1 & r2 arrays in rdshed @ 968'


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

 5000 format(' Debug grid reset to grid number n,row col',3i7)
 6004 format(1x,3i5,f12.1,f9.2,f8.4,f7.0,i5,f7.4,3i5,2x,f5.2,2x,
     *                 <classcount>f5.2)
 6006 format(2x,'basin file:')
 6007 format(4x,'n    yy   xx          da      cap   slope   elv',
     *'   ibn    sl2  ich next reach  frac   classes 1 -> classcount')
 6014 format(i5)
 9005 format(12i5,2f5.0)
 9023 format(' Warning: area correction in grid(n,i,j)',3i5,f9.5)
 9024 format(' Warning: total area = 0.0 for grid(n,i,j)',3i5)
 9025 format(' Warning: no reservoirs or lakes in bsnm.str file')
 9026 format(' Warning: no damage sites in bsnm.str file')

      RETURN

      END SUBROUTINE read_shed_ef

