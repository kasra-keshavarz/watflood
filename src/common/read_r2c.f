      SUBROUTINE read_r2c(unitNum,flnNum,hdrflg,jz,newDataFlag)

      
!***********************************************************************
!     Copyright (c) 1987-2018 Nicholas Kouwen and Dave Watson (NRC) 2006-2018      fjz
        
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
     
! **********************************************************************
!  read_r2c - written May/06 by Dave Watson
!     - Adapdet for general r2c read by NK  Jul. 24/06
!	- Derived from rdtemp written by Nick Kouwen
!	- This subroutine reads the ensim compatible gridded r2c format
!     - First time through it reads the header & the first frame
!       and looks for time stamp on next frame
! **********************************************************************

!     rev. 9.5.13  Feb.  25/08  - NK: changed tolerance for coordinate check to .gt.0.001
!     rev. 9.8.13  Jan.  17/11  - NK: modifications to read_r2c for single frame data
!     rev. 9.9.07  Jan.  10/14  - NK: Overhaul of the frame numbers to EnSim specs
!     rev. 9.9.09  Feb.  24/14  - NK: Fixed reading the time stamp in r2c frame headers
!     rev. 10.4.43 Oct.  10/21  = NK Fixed missing last frames in read_r2c

      use area_watflood

! R2C data module
	use EF_Module
      
      implicit none
	type(TempParam) :: header
	type(FrameRecord) :: frameRec

      CHARACTER(128) :: qstr

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      DIMENSION :: ntoc(9000)

      real*4  :: x1,deltatmp,teltatmp
      integer :: idlast,mhlast,jan,ntocorrect,ios,ios1,
     *           ntoc,i,iend,j,jz,n,
     *           ndeltatt,tem_hr,iostat,iAllocate,iDeallocate,
     *           nh_count,un,fn,xcount_local,ycount_local,
     *           nolines,xcount_max,ycount_max,ntest,result1,nchr
      character(10) ::  fileformat

      character*20  :: junk
	character*1   :: junk2,firstpass(999),hdrflg
	character*6   :: junk0
      DATA idlast/0/
      data firstpass/999*'y'/
      LOGICAL      :: exists,newDataFlag,newFrameFlag
      data newFrameFlag/.false./

     
! Local variables
	integer*4 unitNum, flnNum, iStat
	character*4096 line, subString
	character*128 keyword, value 
	integer lineLen, keyLen, wordCount
	logical rStat, foundEndHeader(999),foundFirstFrame(999),foundFrame(999)
	integer frameCount

	integer frameHour(999) !used to be k
	integer frameHourLast(999) !used to be k
	integer dataHour(999) !used to be nr_count

      DATA ntest/64161/qstr/'timer'/nchr/5/

! Debug line
d     if(iopt.ge.2)then
d     	print*, 'Inside read_r2c at time = ',jz, 'hrs ',fln(flnNum)(1:30)
d         print*,'unit #',unitNum
d     endif

! Initialize default values within frame module
	CALL InitFrameRecord(frameRec)

	foundEndHeader(unitNum) = .false.

! If hdrflg==1 ,then it's the first time through for this event
	if(hdrflg.eq.'1')then

!     rev. 9.5.44  Oct.  27/08  - NK: removed code & obj modules for hasp & rainbow
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     remove for unix
c         if(ichsm.eq.0)then
c           call keychk(qstr,nchr,result1)
c           call userchk(ntest,result1)
c           result1=0
c         endif

		inquire(FILE=fln(flnNum),EXIST=exists)
		if(exists) then
		  open(unit=unitNum,file=fln(flnNum),status='old',iostat=ios)
	    if(ios.ne.0)then
		    print*, 'Problems opening file'
              print*,fln(flnNum)(1:72)
		    STOP ' Stopped in read_r2c @ 89'
		  endif
	      if(iopt.ge.1)print*,'Opened unit'
     *                      ,unitNum,' filename  ',fln(flnNum)(1:40)
		else
		  print*, 'Attempting to open file name '
            print*,fln(flnNum)(1:72)
            print*, 'Fln number  ',flnNum
	      print*, 'Unit number ',unitNum
		  print*, 'but it is not found (in this location)'
		  STOP 'Program aborted in read_r2c @ 95'
          endif
d         print*,'Opened file: ',fln(flnNum)(1:40)          

!     Dave:  this needs to be fixed for this s/r  NK

! Initialize default values within XXXXXXXX module
		call InitTempParam(header)	

! Search for and read r2c file header
		line(1:1) = '#'
		do while((.NOT.foundEndHeader(unitNum)) .AND.
     &	    ((line(1:1) .eq. '#') .OR.
     &		(line(1:1) .eq. ':') .OR.
     &		(len_trim(line) .eq. 0))) 	

			read(unit=unitNum, FMT='((A))', iostat=ios) line	! read a line
c           print*,line(1:70)			
			if(ios .eq. -1)then
				print*, 'ERROR: Premature EndOfFile encountered'
	                  print*, 'FILE=',fln(flnNum)
				STOP ' Stopped in read_r2c @115'
			end if
d             if(iopt.eq.3)print*,line(1:40)
			rStat = Detab(line)				! replace tabs with spaces
			line = adjustl(line)		! Get rid of leading white space
			lineLen = len_trim(line)		! Find the length excluding trailing spaces

			if(line(1:1) .eq. ':')then
				wordCount = SplitLine(line, keyword, subString)	! find the keyword
				rStat = ToLowerCase(keyword)
				KeyLen = len_trim(keyword)

				if(keyword(1:KeyLen) .eq. ':endheader')then
					foundEndHeader(unitNum) = .TRUE.
				else
					! parse the header
					iStat = ParseTempParam(header,keyword,keyLen,
     &													subString)
					if(iStat .lt. 0) then
						write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
						write(*,'(2(A))') '   in line: ',line					
	                              print*, 'FILE=',fln(flnNum)
						STOP ' Stopped in read_r2c @ 137'
						return
					else if(iStat .eq. 0) then
!						write(*,'((A), (A))')
!     &								'Unrecognized keyword line: ',line
					endif
				end if
			end if
		end do

! Assign the parsed parameters to the model variables		
		xcount_temp = header%r2cp%xCount
		ycount_temp = header%r2cp%yCount
		xorigin_temp = header%r2cp%xOrigin
		yorigin_temp = header%r2cp%yOrigin
		xdelta_temp = header%r2cp%xDelta
		ydelta_temp = header%r2cp%yDelta
          deltat_temp = header%r2cp%deltaT
		
!       This section copied from read_swe_ef by NK  Jul. 5/07
!       Validate parameters
!     rev. 9.5.13  Feb.  25/08  - NK: changed tolerance for coordinate check to .gt.0.001
        if(fln(flnNum)(1:7).ne.'..\cmc\')then   
!         Assumes processed CMC Grib2 files are in *watflood\cmc\*          
          if(abs(header%r2cp%xOrigin-xorigin).gt.0.001)
     *        print*,'xorigin_temp.ne.xorigin'
	    if(abs(header%r2cp%yOrigin-yorigin).gt.0.001)
     *        print*,'yorigin_temp.ne.yorigin'
          if(abs(header%r2cp%xDelta-xdelta).gt.0.001)
     *     	print*,'xdelta_temp.ne.xdelta'
	    if(abs(header%r2cp%yDelta-ydelta).gt.0.001)
     *  	    print*,'ydelta_temp.ne.ydelta'
          if(header%r2cp%xCount.ne.xcount)
     *    	print*,'xcount_temp.ne.xcount'
          if(header%r2cp%yCount.ne.ycount)
     *  	    print*,'ycount_temp.ne.ycount'
          if(abs(header%r2cp%xOrigin-xorigin).gt.0.001.or.
     &  	  abs(header%r2cp%yOrigin-yorigin).gt.0.001.or.
     &	    abs(header%r2cp%xDelta-xdelta).gt.0.001.or.
     &	    abs(header%r2cp%yDelta-ydelta).gt.0.001.or.
     &   	  header%r2cp%xCount.ne.xcount.or.
     &  	  header%r2cp%yCount.ne.ycount) then
            print*,'WARNING'
	      PRINT*,'Mismatch between ',fln(flnNum)(1:60)
            print*,'    and SHD files'
            print*,'Check files for origins, deltas and counts'
            print*,'Could be due to # significant digits in header' 
            print*,'   file      shd file'
	      print*,header%r2cp%xOrigin,xorigin
	      print*,header%r2cp%yOrigin,yorigin
	      print*,header%r2cp%xDelta,xdelta
	      print*,header%r2cp%yDelta,ydelta
	      print*,header%r2cp%xCount,xcount
	      print*,header%r2cp%yCount,ycount
	      print*,'Possible problem: pdl file does not match grid'
            print*,'No problem if intended - eg for Grib2 files'
c	  	  pause 'Program paused in read_r2c_ef @ 169'
          endif
        endif


!     rev. 9.9.07  Jan.  10/14  - NK: Overhaul of the frame numbers to EnSim specs
c      go to 99999
c     temporary fix to bypass scanning the file for apparently no reason  
C     24/01/2003 well, apparently there is areason - sometimes you need to know how mane frames there will be
c     as in diff.f    
! Scan data frames
		frameCount = 0

		do while((.NOT.EOF(unitNum)))  ! EOF not used in unix  nk
!		do while((.NOT.found_data_end))   ! this does not work on PC

			read(unit=unitNum, FMT='((A))', iostat=ios) line	! read a line
c	      print*,line(1:72)
			if(ios .eq. -1)then
			 write(6,'((A))') 'ERROR: Premature EndOfFile encountered'
	                  print*, 'FILE=',fln(flnNum)
                        print*,'frame count =',frameCount
	                  print*,'last line read:'
	                  print*,line
				STOP ' Stopped in read_r2c @ 186'
			end if
d             if(iopt.eq.3)print*,line(1:40)

			rStat = Detab(line)				! replace tabs with spaces
			line = adjustl(line)		! Get rid of leading white space
			lineLen = len_trim(line)		! Find the length excluding trailing spaces

			if(line(1:1) .eq. ':')then
				wordCount = SplitLine(line, keyword, subString)	! find the keyword
				rStat = ToLowerCase(keyword)
				KeyLen = len_trim(keyword)
				
				if(keyword(1:KeyLen).eq.':frame')then
					iStat = ParseFrameLine(frameRec,keyword,keyLen,
     &													subString)
					
!	Identify the first and last frame's timestamp 
					if(iStat .lt. 0) then
					   write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
						write(*,'(2(A))') '   in line: ',line					
	                              print*, 'FILE=',fln(flnNum)
						STOP ' Stopped in read_r2c @ 208'
						return
					else if(iStat .eq. 0) then
!						write(*,'((A), (A))')
!     &								'Unrecognized keyword line: ', line
					else if(frameRec%frame.EQ.1) then

						header%startJulianDay =
     &						JDATE(frameRec%tStamp%year,
     &							frameRec%tStamp%month,
     &							frameRec%tStamp%day)
						header%startHour = frameRec%tStamp%hour
					else
						header%endJulianDay =
     &						JDATE(frameRec%tStamp%year,
     &							frameRec%tStamp%month,
     &							frameRec%tStamp%day)
						header%endHour = frameRec%tStamp%hour
					end if
					header%r2cp%frameCount = header%r2cp%frameCount+1
				end if
			end if
		enddo	
	
c      print*,'header%r2cp%frameCount',header%r2cp%frameCount

c99999 continue
c      header%r2cp%frameCount=1
      
      if(header%r2cp%frameCount.ne.0)then
!    	   deltat2 = model timestep in hours 
		deltat2 = 1	
!	  njtemp = number of hours spanned by this tem file
!	  convert to julian days to properly calculate hours spanned
        nhtemp = (header%endJulianDay - header%startJulianDay)*24
     &				 + (header%endHour - header%startHour) + 1

        
!	 Position to start of data (immediately after first frame record but before first frame data)
		REWIND (unitNum)
		foundFirstFrame(unitNum) = .false.
		do WHILE(.NOT.foundFirstFrame(unitNum))
			read(unit=unitNum, FMT='((A))', iostat=ios) line	! read a line
c           print*,line(1:60)
			if(line(1:1) .eq. ':')then
				wordCount = SplitLine(line, keyword, subString)	! find the keyword
				rStat = ToLowerCase(keyword)
				KeyLen = len_trim(keyword)
				if(keyword(1:KeyLen).eq.':frame')then
					iStat = ParseFrameLine(frameRec,keyword,
     &							keyLen,subString)
					foundFirstFrame(unitNum) = .true.
				end if
			end if		      
		enddo	
	endif
cd     print*,'foundfirstframe(unitNum)',foundfirstframe(unitNum)		
      
      	frameHour(unitNum)=
     *  	(JDATE(frameRec%tStamp%year,frameRec%tStamp%month,
     &		frameRec%tStamp%day) - header%startJulianDay) * 24
     &			+ (frameRec%tStamp%hour - header%startHour) + 1
c          print*,frameRec%tStamp%year,frameRec%tStamp%month,
c     *      frameRec%tStamp%day,frameRec%tStamp%hour
!         This is the time stamp of the currently read frame          
         
          
c          write(600,*)'jz       ',jz
c          write(600,*)'frameHour',frameHour(unitNum)
c          write(600,*)'year     ',frameRec%tStamp%year
c          write(600,*)'month    ',frameRec%tStamp%month
c          write(600,*)'day      ',frameRec%tStamp%day
c          write(600,*)'juDay    ',header%startJulianDay
c          write(600,*)'hour     ',frameRec%tStamp%hour
c          write(600,*)'startHour',header%startHour
c          write(600,*)'-------------------------------------------'
          
          
!     rev. 10.1.70 Mar.  03/17  - NK: Added year_now2 etc. for converting Grib2 files 
      
          
              year_now2=frameRec%tStamp%year
              month_now2=frameRec%tStamp%month
              day_now2=frameRec%tStamp%day
              hour_now2=frameRec%tStamp%hour
c      write(800,*)'First frame r2c:',
c     *      year_now2,month_now2,day_now2,hour_now2  !  from read_r2c frame

!         This is the time stamp before read frame          
          year_now1=year_now2
          month_now1=month_now2
          day_now1=day_now2
          hour_now1=hour_now2
          
c!     Added 2019/02/16  NK to skip past the end of the year.
c      if(month_now2.eq.1.and.day_now2.eq.1.and.hour_now2.eq.0)then
c          if(frameHour(unitNum).gt.10)then
c              year_now2=year_now2-1
c              month_now2=12
c              day_now2=31
c              hour_now2=24
c          endif
c      endif


!     this makes the frameHour = 1 as the hour of the first time stamp is 
!     subtracted from each of the following frame. This allows the application 
!     of the data to the time period up to the next frame which is a hell of a good idea
!     
      if(frameHour(unitNum).ne.1)then
        print*,'finished reading 1st frame header in ',fln(flnNum)(1:30)
        print*,'frameHour=',frameHour(unitNum)
        print*,'WARNING:'
        print*,'first frame not labeled as frame # 1'
        print*,'in file ',fln(flnNum)(1:50)
        print*      
      endif

!      set model hour to first hour so we have data as we start     <<<<<<<< nasty!!!
		frameHour(unitNum)=1 
		dataHour(unitNum) = 0 ! set data hour to just before first hour

! firstpass means this is the first tem file of the first event...which is the only time it needs to be done
		if(firstpass(unitNum).eq.'y')then
		  firstpass(unitNum)='n'        
!		  allocate(ttemp(ycount3,xcount3),tempv(na),tempvmin(na),
!     *						rh(na),stat=iAllocate)
          if(.NOT.allocated(inarray))then    
!               inarray not previously allocated
		    allocate(inarray(ycount_temp,xcount_temp),
     *                                         stat=iAllocate)
              print*,1,'Inarray allocated with ',
     *                   ycount_temp,xcount_temp,' unit=',unitNum
			if(iAllocate.ne.0) STOP
     *	   'Error with allocation in read_r2c'
			xcount_max=xcount_temp
			ycount_max=ycount_temp
			xcount_local=xcount_temp
			ycount_local=ycount_temp
                  
          else
!               inarray previously allocated but check to see this 
!               data fits in the allocated memory
 		    if(xcount_temp.gt.xcount_max.or.ycount_temp.gt.ycount_max)then
		      deallocate(inarray,stat=iDeallocate)
                print*,'inarray deallocated'
			if(iDeallocate.ne.0)then
			  print*,'Warning: error with deallocation of inarray'
			endif

			xcount_max=xcount_temp
			ycount_max=ycount_temp
			xcount_local=xcount_temp
			ycount_local=ycount_temp
!           outarray is in areawfo
!           this has to be allocated before calling write_r2c
		    allocate(inarray(ycount_max,xcount_max),
     *                             stat=iAllocate)
              print*,2,'Inarray allocated with ',
     *                   ycount_temp,xcount_temp,' unit=',unitNum
			if(iAllocate.ne.0)then
			  STOP 'Error with allocation of inarray in read_qlz'      
			end if
	    endif
         endif
		endif
! Check for change of grid size (could happen for next event)
        if(allocated(inarray))then    ! check that inarray is allocated
		  if(xcount_temp.gt.xcount_max.
     *                          or.ycount_temp.gt.ycount_max)then
			deallocate(inarray,stat=iDeallocate)
			if(iDeallocate.ne.0) then
			  print*,'Warning: error with deallocation of intarray'
			endif

			xcount_max=xcount_temp
			ycount_max=ycount_temp
!       outarray is in areawfo
!       this has to be allocated before calling write_r2c
		  allocate(inarray(ycount_max,xcount_max),stat=iAllocate)
              print*,3,'Inarray allocated with ',
     *                   ycount_temp,xcount_temp,' unit=',unitNum
			if(iAllocate.ne.0) then
			STOP 'Error with allocation of inarray in read_r2c'      
			end if
		  endif
	      endif

		xcount_temp=xcount_local
		ycount_temp=ycount_local
          
		return
      endif          ! (hdrflg.eq.'1')

          
! **********************************************************
! **********************************************************
	
!     Finished reading the header and first frame   

! **********************************************************
! **********************************************************

! Check for change of grid size (could happen for next event)
        if(.not.allocated(inarray))then    ! check that inarray is allocated
c		  if(xcount_temp.gt.xcount_max.
c     *                          or.ycount_temp.gt.ycount_max)then
c			deallocate(inarray,stat=iDeallocate)
c			if(iDeallocate.ne.0) then
c			  print*,'Warning: error with deallocation of intarray'
c			endif

			xcount_max=max(xcount_max,xcount_temp)
			ycount_max=max(ycount_max,ycount_temp)
!       outarray is in areawfo
!       this has to be allocated before calling write_r2c
		    allocate(inarray(ycount_max,xcount_max),stat=iAllocate)
              print*,4,'Inarray allocated with ',
     *                   ycount_temp,xcount_temp,' unit=',unitNum
		    if(iAllocate.ne.0) then
		      STOP 'Error with allocation of inarray in read_r2c'      
			end if
d			print*,'inarray reallocated with :',ycount_max,xcount_max
d			print*
d			print*
c		  endif
	  endif

		xcount_temp=xcount_local
		ycount_temp=ycount_local
          
cd      print*,'frameHour = ',frameHour(unitNum),', dataHour = ',
cd    *   dataHour(unitNum)
          
          

	 dataHour(unitNum) = jz
c	 dataHour(unitNum) = dataHour(unitNum) + 1
c      datahour(unitNum)=1
c      print*,'jz = ',jz

c      if(foundFirstFrame(unitNum))then
c      print*,jz,frameHour(unitNum),datahour(unitNum)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(frameHour(unitNum).eq.datahour(unitNum))then
          
          
!      needed for DA.f          
c      write(800,*)':::',line(1:50) ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c     write(800,*)'In read_r2c    ',jz,month_now1,day_now1,hour_now1     
          
           
        do i=1,ycount_local
          read(unitNum,*,iostat=ios)(inarray(i,j),j=1,xcount_local)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        end do
!      write(800,*)'from read_r2c :',
!     *          jz,month_now1,day_now1,hour_now1,inarray(1,1) !  from read_r2c frame

      
c        if(iopt.ge.2)then
c          print*,'CG value',inarray(ycount/2,xcount/2)   ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
cd          print*
c          print*,'fln ',fln(flnNum)(1:40),' SW corner data'
c          print*,frameHour(unitNum),datahour(unitNum)
c          do i=1,4
c              write(*,*)(inarray(i,j),j=1,4)
c          end do
c        endif

        if(ios.ne.0)then
          print*,' Error:'
          print*,' At frame # ',frameRec%frame
          print*,' frameHour',frameHour(unitNum)
          print*,' last data read:'
          print*,' row   = ',i,' column= ',j
          print*,' unit number ',unitNum
          print*,' file name   ',fln(flnNum)(1:40)
          print*
          print*,'Likely problem, frame repeated'
          print*,'OR with regl_conv not enough gap filling passes'
          print*,'OR r2s format file saved as single frame r2c'
          print*,'To continue - hit return'
          print*,'To abort, hit "ctrl C"'
!          pause

!     rev. 10.5.04 Mar.  21/23  = NK Changed pause in read_r2c to continue task
          inquire(file='runReport.txt',EXIST=exists)
          open(unit=99,file='runReport.txt',status='unknown') 
          If(exists)then
              do while(.not.eof(99))
                  read(99,*,iostat=ios)line
              end do
          endif
          write(99,*)'Error:~~~~~~~~~~~~~~~~~~~~'
          write(99,*)'file name   ',fln(flnNum)(1:40)
          write(99,*)'Unit number ',unitNum
          write(99,*)'At frame #  ',frameRec%frame
          write(99,*)'frameHour   ',frameHour(unitNum)
          write(99,*)'last data read:'
          write(99,*)'row   = ',i,' column= ',j
          write(99,*)
          write(99,*)'Likely problem, frame repeated'
          write(99,*)'OR with regl_conv not enough gap filling passes'
          write(99,*)'OR r2s format file saved as single frame r2c'
          write(99,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~'
          close(unit=99,status='keep')
          STOP ' program aborted - read_r2c @ 364'
          
c       return             ! this had to be added to read the diff files 20140127 NK
        else
d         if(iopt.ge.2)then
d           print*,'got next frame - in read_r2c on unit=',unitNum
d         endif
          newDataFlag=.true.
        endif
c       endif

!	Read the next frame. If we are at the end of the file...close
        foundFrame(unitNum) = .false.
        do WHILE(.NOT.foundFrame(unitNum))
            
            read(unit=unitNum, FMT='((A))', iostat=ios) line	! read a line

c      write(800,*)'nxt',line(1:50) ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

c           write(800,*)line(1:72) ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!           print*,line(1:72) ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!			if(EOF(unitNum)) then 
            if(ios.lt.0) then   ! EOF not used in unix  nkc

!               NK  Oct. 29/06
!               Added global flag to end event
!               end of event is when eof file found in first file
!               other files may have more data but we can only run shortest
                found_data_end=.true.
                close(unit=unitNum,status='keep',iostat=ios)
              if(dds_flag.eq.0)then
                 write(51,*,iostat=ios)'read_r2c: Closed unit ',
     *                        unitNum,' Filename=  ',fln(flnNum)(1:40)
!                 write(800,*)'read_r2c: Closed unit ',
!     *                        unitNum,' Filename=  ',fln(flnNum)(1:40)
                 write(*,*) 'read_r2c: Closed unit ',
     *          	          unitNum,' Filename=  ',fln(flnNum)(1:40)
                 write(*,*)'Last frame: ',
     *                        year_now2,month_now2,day_now2,hour_now2
              endif
!     rev. 10.4.43 Oct.  10/21  = NK Fixed missing last frames in read_r2c
!             This is the time stamp of the last currently read frame          
              year_now1=year_now2
              month_now1=month_now2
              day_now1=day_now2
              hour_now1=hour_now2
              RETURN
            endif

            if(line(1:1) .eq. ':')then
                wordCount = SplitLine(line, keyword, subString)	! find the keyword
                rStat = ToLowerCase(keyword)
                KeyLen = len_trim(keyword)
                if(keyword(1:KeyLen).eq.':frame')then
                    iStat = ParseFrameLine(frameRec,keyword,
     &												keyLen,subString)
                    foundFrame(unitNum) = .true.
                end if
            end if		      
        enddo	

!       so the grid size can be checked in sub at any time
		xcount_temp=xcount_local
		ycount_temp=ycount_local
          
!         This is the time stamp of the previously read frame          
          

         frameHourLast(unitNum)=frameHour(unitNum) 
         frameDayLast=frameRec%tStamp%day
         
         
         
!	  Determine the next frameHour		
		frameHour(unitNum)=(JDATE(frameRec%tStamp%year,frameRec%tStamp%month,
     &		frameRec%tStamp%day) - header%startJulianDay) * 24
     &			+ (frameRec%tStamp%hour - header%startHour) + 1
          
c          print*,frameRec%tStamp%year,frameRec%tStamp%month,
c     *      frameRec%tStamp%day,frameRec%tStamp%hour

!         This is the time stamp of the currently read frame          
          year_now1=year_now2
          month_now1=month_now2
          day_now1=day_now2
          hour_now1=hour_now2
          
!         This is the time stamp of the next read frame          
          year_now2=frameRec%tStamp%year
          month_now2=frameRec%tStamp%month
          day_now2=frameRec%tStamp%day
          hour_now2=frameRec%tStamp%hour
          
          
!      write(800,*)'read_data_r2c :',
!     *          jz,month_now1,day_now1,hour_now1,inarray(1,1) !  from read_r2c frame

d          write(63,*)jz,' next frame header =',frameHour(unitNum)
          
d        if(iopt.ge.2)then
d          print*,'finished reading frame header in ',fln(flnNum)(1:30)
d          print*,'frameHour  next =',frameHour(unitNum)   
c          write(456,*)year_now2,month_now2,day_now2,hour_now2   
d        endif

c	else
c		do nothing...keep last temperature values...
!     rev. 10.2.24 May   21/18  - NK: Added error message in Read_rain & read_tmp
        if(frameHour(unitNum)-frameHourLast(unitNum).le.0)then
          print*
          print*,'Error:'
          print*,'file name ',fln(flnNum)(1:50)
          print*,'frameHour(unitNum) =',frameHour(unitNum)
          print*,'frameHourLast(unitNum) =',frameHourLast(unitNum)
          print*,'Note: Frame numbers must be monotomically increasing'
          print*,'in an r2c time series file.'
          print*,'year ',frameRec%tStamp%year
          print*,'month',frameRec%tStamp%month
          print*,'day  ',frameRec%tStamp%day
          print*,'hour ',frameRec%tStamp%hour
          print*
          print*,'Possible cause: '
          print*,'month, day or hour repeated in the r2c file'
          print*
        call BEEPQQ(3000,300)
        call SLEEPQQ(100) 
        call BEEPQQ(2000,500)
!       take out for unix - non standard
          print*,'Unfinished conversion - hit enter to continue'
          pause 'or Ctrl C to abort & fix input data'
c          stop 'Program aborted in read_temp_ef @ 411'
        endif

	endif

	RETURN

!     FORMATS

 9993 format(' error reading inarray at hour/line ',2i10/)
	

      END SUBROUTINE read_r2c

