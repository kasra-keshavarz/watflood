      SUBROUTINE write_r2c(un,fn,
     *            no_frames,no_classes,frame_no,class_no,
     *            no_signf)

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
      
!     This file is part of WATFLOOD(R).

!     WATFLOOD is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, Version 3

!     WATFLOOD is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have recieved a copy of the GNU General Public License
!     along with WATFLOOD.  If not, see <https://www.gnu.org/licenses/>.
      
!     rev. 9.1.64  Oct.  03/04  - NK: Coded up new header in ragmet.for
!     rev. 9.5.70  Oct.  11/09  - NK: fixed timer for r2c frames (use year_now)
!     rev. 9.9.23  Aug.  10/14  - NK: Uniform arguments for write_r2c

!***********************************************************************
! - THIS SUBROUTINE OUTPUTS and r2c file.

! - List of arguments:

!   I - un      int        unit number
!   I - fn      int        file number
!   I - itogo   int        no. of hours until next rainfall
!   R - unit_conversion    REAL*4     conversion factor (area2)
!   I - FLN     CHAR*12    file names
!   I - DIR     CHAR*12    directory name
!   I - ocflg   int        open_close flag 1 open file -1 close file
!   I - frmflg  int        frame flag  1 new frame -1 end frame  
!                          0 each call = 1 frame with frame marks
!***********************************************************************

!Usage - argument list:
!
!                'time series'               'static'
!un                    n                        n
!fn                   aaa                      aaa
!no_frames            m|m                      0|0
!no_classes           0|0                      m|m
!frame_no             0|j                      0|0
!class_no             0|0                      0|n
!no_signf             k|k                      k|k

! write header | write data frames
! m, n & k = interger values
! aaa = chracted value
! for time series r2c:
! first call: set no_frames = 0 & frame_no = 0 - this will write the header
! subsequent calls to write the date, set
!             frame_no = # of trhe consecutive frames
!             no_frames = # of frames to be written
!             one file for all events: set no_frames = frame_no + 8784
!                                      the file will be closed when the program stops normally
!             one file for each event: set no_frames = # frames to be writted this event
!                                      so either mhtot OR mhtot/delta t of the write interval
!for static r2c file:
! first call: set the no_classes = the # of attributes to be written and class_no = 0
! subsequent calls to write the data, set
!             class_no = # of the consecutive data sets to be written
!             no_classes - as above

      use area_watflood
	implicit none
	save

      INTEGER       :: un,fn,ii,no_signf,hour_no,hours_togo,
     *                 no_frames,no_classes,frame_no,class_no,
     *                 i,j,ios
      integer       :: nh  ! hour number
	character(20) :: junk
      CHARACTER(10) :: time
      CHARACTER(8)  :: cday

c      write(*,9)un,fn,no_frames,no_classes,frame_no,class_no,no_signf
c9     format(8i8)     

!     FIRST TIME THROUGH THIS SUBROUTINE ONLY
!     OPEN OUTPUT FILE AND WRITE HEADER

c      hour_no=frame_no
c      hour_no=nh

c     if(iopt.eq.2)then
c       write(63,*)'In write_r2c @ L42 - hour_no=',frame_no
c     endif
        
!     this s/r will only allow:
!     multiple frames for 1 class
!                or
!     1 frame for multiple classes
      if(no_frames.gt.1.and.no_classes.gt.1)then
        write(98,*)'Error: Programming error'
        write(98,*)'Error: no_frames > 1  and  no_classes > 1'
        write(98,*)'Error: This is not allowed'
        write(98,*)'Error: Program aborted due to programming error'
        print*,'Programming error'
        print*,'no_frames > 1  and  no_classes > 1'
        print*,'This is not allowed'
        print*
        stop 'Program aborted due to programming error'
!       This can only be cause by misuse of this s/r 
!         in the calling program
      endif

c      if(frame_no.eq.0)then
!     rev. 9.9.23  Aug.  10/14  - NK: Uniform arguments for write_r2c
      if(frame_no.eq.0.and.class_no.eq.0)then

!       write the header ONLY

!       FILE NAMES AND UNIT NUMBERS DIFFER BY 30
!        write(*,1400)fn,fln(fn)
! 1400   format(' opening fln(',i3,'):',a30,'---')
!        write(*,*)

        open(unit=un,file=fln(fn),status='unknown',iostat=ios)
	  if(ios.ne.0)then   ! added Nov. 10/14  nk
	    write(98,*)'Error: Unable to open file ',fln(fn)(1:40)
          write(98,*)'Error: On unit ',un
	    write(98,*)'Error: ios = ',ios
          write(98,*)'Error: Possible cause(s):'
          write(98,*)'Error: file in use by another application'
          write(98,*)'Error: or target directory does not exist'
	    write(98,*)'Error: or New event files req`d for the new'
	    write(98,*)'Error: wind and temp. difference files'
	    write(98,*)'Error: Program aborted in write_r2c @ L83'
          print*
	    print*,'Unable to open file ',fln(fn)(1:40)
          print*,'On unit ',un
	    print*,'ios = ',ios
          print*,'Possible cause(s):'
          print*,'file in use by another application'
          print*,'or target directory does not exist'
	    print*,'or New event files req`d for the new'
	    print*,'wind and temp. difference files'
	    stop 'in write_r2c @ L83'
	  endif 
	  if(iopt.ge.1)then
          PRINT*
          write(98,*)'Info: In write_r2c opened unit=',un
		write(98,*)'Info: filename=',fln(fn)(1:50)
	  endif
        write(un,3005)'########################################'
        write(un,3005)':FileType r2c  ASCII  EnSim 1.0         '
        write(un,3005)'#                                       '
	  write(un,3005)'# DataType               2D Rect Cell   '
        write(un,3005)'#                                       '
        write(un,3005)':Application             WATFLOOD       '
	  write(un,3005)':Version                 2.1.23         '
	  write(un,3020)':WrittenBy          ',author
	  write(un,3001)':Weight used        ',nw
        call date_and_time(cday,time)
	  write(un,3010)':CreationDate       ',
     *       cday(1:4),cday(5:6),cday(7:8),time(1:2),time(3:4)
3010  format(a20,a4,'-',a2,'-',a2,2x,a2,':',a2)
        write(un,3005)'#                                       '
	  write(un,3005)'#---------------------------------------'

        write(un,3005)'#                                       '
	  write(un,3020)':Name               ',name
        write(un,3005)'#                                       '
	  write(un,3004)':Projection         ',coordsys_temp
	  if(coordsys_temp(1:7).eq.'LatLong')then
	    write(un,3004)':Ellipsoid          ',datum_temp
!     rev. 9.9.62  Mar.  21/15  - NK: Change zone from character to integer
	  elseif(coordsys_temp(1:3).eq.'UTM')then
	    write(un,3004)':Ellipsoid          ',datum_temp
          write(un,3001)':Zone               ',zone_temp
	  endif
        write(un,3005)'#                                       '
        write(un,3003)':xOrigin            ',xorigin_temp
        write(un,3003)':yOrigin            ',yorigin_temp
        write(un,3005)'#                                       '
        write(un,3020)':SourceFile         ',source_file_name
        write(un,3005)'#                                       '
c        if(no_frames.eq.1.and.no_classes.ge.1)then
        if(no_classes.ge.1)then
          do i=1,no_classes
      	    write(un,3007)':AttributeName',i,' Class',i
	    end do
        endif
c          write(un,3005)'#                                       '

!       Fixed by adding attribute_count to the write statement  NK  Dec. 30/17
        if(attribute_count.gt.0)then
          do i=1,attribute_count
            write(un,3019)':AttributeName      '
     *               ,attribute_count,attribute_name
            write(un,3019)':AttributeUnits     '
     *               ,attribute_count,attribute_units  
          end do
        endif
        attribute_count=0  ! don't want this lose in the code



!         see note below @***       
        write(un,3005)'#                                       '
        write(un,3001)':xCount             ',xcount_temp
        write(un,3001)':yCount             ',ycount_temp
        write(un,3003)':xDelta             ',xdelta_temp
        write(un,3003)':yDelta             ',ydelta_temp
        if(no_frames.le.1)then
          write(un,3005)'#                                       '
          write(un,3004)':SampleTime         ',startdate,starttime
!          write(un,3004)':StartDate          '
        endif
        write(un,3005)'#                                       '
        if(unit_conversion.ne.0.0)then
          write(un,3003)':UnitConversion     ',unit_conversion
        endif
        if(name.eq.'Snow Water Equivalent                   ')then
	    write(un,3003)':InitHeatDeficit    ',init_heat_deficit
        endif
!       added Sept. 26/09 nk to add lapse rates to header info

	  if(no_hdrcomments.gt.0)then
          write(un,3005)'#                                       '
	    do i=1,no_hdrcomments
	      write(un,3008)hdrcomment(i),hdrcomment_value(i)
            if(iopt.ge.1)write(*,3008)hdrcomment(i),hdrcomment_value(i)
	    end do
	  endif
        write(un,3005)'#                                       '
     	  write(un,3005)':endHeader                              '

        return

	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! fix fix

!       this is just to read in the new format
!       still needs to be programmed for multiple classes

!     the header is written only the first time through
!     this part is called with each class

!     :Frame and :EndFrame lines are written only for time series

      if(data_source.eq.'     ')data_source=source  ! source from area2

c      if(hour_no.ge.1)then
      if(frame_no.ge.1.and.class_no.eq.0)then
      
!     rev. 9.5.70  Oct.  11/09  - NK: fixed timer for r2c frames (use year_now)

!       write the frame header:
!       write the frame header:
!       write the frame header:
          
        if(no_frames.gt.1)write(un,3333)':Frame',
     *     abs(frame_no),abs(frame_no),
     *         year_now,month_now,day_now,hour_now
c        if(no_frames.gt.1)write(900,3333)':Frame',
c     *     abs(frame_no),abs(frame_no),
c     *         year_now,month_now,day_now,hour_now
 3333   format(a6,2i10,3x,'"',i4,'/',i2.2,'/',i2.2,1x,i2.2,
     *                   ':00:00.000"',2x,a5,2i5,2x,a5)
          
        
      endif
c          if(no_frames.gt.1)write(*,3024)':Frame',
c     *     abs(frame_no),abs(frame_no),year1,month_now,day_now,hour_now

!       NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE 
!       The r2c grids are written upside down: 
!                 south on top, north on bottom    !!!!!!!!!
c      if(frame_no.gt.0)then
      if(frame_no.gt.0.or.class_no.gt.0)then
        if(no_signf.eq.0)then
          do i=1,ycount_temp
            write(un,1300)(outarray(i,j),j=1,xcount_temp)
1300        format(9999(1x,f5.0))
          end do          
        elseif(no_signf.eq.1)then   ! swe, precip and temp
          do i=1,ycount_temp
            write(un,1301)(outarray(i,j),j=1,xcount_temp)
1301        format(9999(1x,f5.1))
          end do          
        elseif(no_signf.eq.2)then   ! swe, precip and temp
          do i=1,ycount_temp
            write(un,1302)(outarray(i,j),j=1,xcount_temp)
1302        format(9999(1x,f5.2))
          end do          
        elseif(no_signf.eq.3)then   ! swe, precip and temp
          do i=1,ycount_temp
            write(un,1303)(outarray(i,j),j=1,xcount_temp)
1303        format(9999(1x,f6.3))
          end do          
        elseif(no_signf.eq.4)then   ! swe, precip and temp
          do i=1,ycount_temp
            write(un,1304)(outarray(i,j),j=1,xcount_temp)
1304        format(9999(1x,f7.4))
          end do          
        elseif(no_signf.eq.5)then   ! swe, precip and temp
          do i=1,ycount_temp
            write(un,1305)(outarray(i,j),j=1,xcount_temp)
1305        format(9999(1x,f8.5))
          end do          
        elseif(no_signf.eq.6)then   ! swe, precip and temp
          do i=1,ycount_temp
            write(un,1306)(outarray(i,j),j=1,xcount_temp)
1306        format(9999(1x,f9.6))
          end do          
        elseif(no_signf.eq.8)then   ! flow
          do i=1,ycount_temp
            write(un,1308)(outarray(i,j),j=1,xcount_temp)
1308        format(9999(1x,e10.3))
          end do          
        elseif(no_signf.eq.9)then   ! drain, dsnow
          do i=1,ycount_temp
            write(un,1309)(outarray(i,j),j=1,xcount_temp)
1309        format(9999(1x,f6.2))
          end do          
        elseif(no_signf.eq.10)then   ! drain, dsnow
          do i=1,ycount_temp
            write(un,1310)(outarray(i,j),j=1,xcount_temp)
1310        format(9999(1x,f10.0))
          end do          
        elseif(no_signf.eq.11)then   ! drain, dsnow
          do i=1,ycount_temp
            write(un,1311)(outarray(i,j),j=1,xcount_temp)
1311        format(9999(1x,f8.3))
          end do          
        elseif(no_signf.eq.12)then   ! elv_mean,elv_max
          do i=1,ycount_temp
            write(un,1312)(outarray(i,j),j=1,xcount_temp)
1312        format(9999(1x,f8.2))
          end do          
        else                        ! init soil moisture
          do i=1,ycount_temp
            write(un,1307)(outarray(i,j),j=1,xcount_temp)
1307        format(9999(1x,e12.6))
          end do          
        endif
      endif

!       changed by NK  Aug. 8/14
c        if(class_no.eq.0)write(un,3012)':EndFrame'
      if(frame_no.ge.1.and.class_no.eq.0)write(un,3012)':EndFrame'

      if(frame_no.ge.no_frames.and.class_no.ge.no_classes)then
c      if(nh.eq.no_frames.and.class_no.eq.no_classes)then
        close(unit=un,status='keep')
        if(no_frames.gt.0)then
          write(51,*,iostat=ios)'Closed unit ',un,' Filename=  ',
     *              fln(fn)(1:72),'fr #',frame_no,no_frames   
         write(98,*)'Info: Closed unit ',un,' Filename= ',fln(fn)(1:72),
     *              'fr #',frame_no,no_frames   
        else
          write(51,*,iostat=ios)'Closed unit ',un,' Filename=  ',
     *             fln(fn)(1:72), 'fr #',class_no,no_classes   
         write(98,*)'Info: Closed unit ',un,' Filename= ',
     *              fln(fn)(1:72),'fr #',class_no,no_classes   
        endif
	endif

      return

! FORMATS
 3000 format(a10,i5)
 3001 format(a20,i16)
 3002 format(2a20)
 3003 format(a20,f16.6)
 3004 format(a20,a10,2x,a10)
 3005 format(a40)
 3006 format(a3,a10)
 3007 format(a14,i5,a6,i5)
 3008 format(a40,f16.7)
 3012 format(a9)
 3019 format(a20,i5,5x,a40)
 3020 format(a20,10x,a40)
     
       END SUBROUTINE write_r2c









