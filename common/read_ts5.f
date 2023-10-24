      SUBROUTINE read_ts5(unitNum,flnNum)

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
!   THIS S/R  reads daily time series data 
!
!     rev. 9.9.77  Sep.  11/15  - NK: S/r read_ts5 created
!     rev. 10.2.51 Apr.  03/19  - NK: New section to read ts5 format file for swe
!
!***********************************************************************

      use area_watflood
      use areacg
      implicit none
      
      logical         :: foundEndHeader
      character*512   :: line,junk
      character(1)    :: chr(256)
!Tegan: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      character*10    :: yyyymmdd(367)  ! needs one more becuse read overshoots the list
!Tegan: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer         :: flnNum,ll,i,j,rank,iallocate,unitnum
      integer         :: iDeallocate,ios
      
!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE
      
      open(unit=99,file=fln(flnNum),status='old')
d      print*,'Opened ',fln(flnNum)(1:50)
      foundEndHeader=.false.
      i=0
      do while(.not.line(1:10).eq.':BeginLine')
          read(99,99000)line
d          print*,'xx',line(1:50)
          if(line(1:14).eq.':AttributeName')then
             i=i+1
          endif
      end do
!         allocate the x & y coordinates for the no of stations          
          if(allocated(x_temp))then
            deallocate(x_temp,y_temp,gname_temp)
          endif
          allocate(x_temp(i),y_temp(i),
     *                         gname_temp(i),stat=iAllocate)
          if(iAllocate.ne.0)then
              print*,'Attempting to read ',fln(flnNum)(1:40)
              print*,'allocating x_temp & y_temp for ',i,' locations'
              STOP
     *      'Error with allocation of x_temp & y_temp in read_ts5 @ 38'
          endif
      
      rewind(unit=99)
      j=0    ! use j as the no of locations (cols)
      i=0
      read(99,99000)line
d      print*,'yy',line(1:72)
      do while(.not.foundEndHeader)

99000   format(a256)       
        if(line(1:10).eq.':BeginLine')then
!         found the start of the x-y coordinates        
          read(line,*)line,xcount_temp  ! xcount_temp = the no of locations (cols)
d          print*,xcount_temp,'# station(s) found' 
d      print*,'# data points= =',xcount_temp
        endif  

!     rev. 10.2.51 Apr.  03/19  - NK: New section to read ts5 format file for swe
!       For each location there is a keyword :Point with x & y  
!       this is the only info we need from the header      
        if(line(1:14).eq.':AttributeName')then
          i=i+1
          read(line,*)junk,gname_temp(i)
d          print*,junk(1:14),i,gname_temp(i)
        endif 
!     END rev. 10.2.51 Apr.  03/19  - NK: New section to read ts5 format file for swe
        
        if(line(1:6).eq.':Point')then
          j=j+1
          read(line,*)junk,x_temp(j),y_temp(j)
d          print*,junk(1:6),j,x_temp(j),y_temp(j)
        endif 
        
        read(99,99000)line
d        print*,'jj',line(1:40)
c        if(line(1:8).eq.'#StopLog')ruletype='StopLog'  
c        if(line(1:12).eq.'#TargetLevel')ruletype='TargetLevel'  
c        print*,ruletype
c        pause 43434
        if(line(1:10).eq.':endHeader')foundEndHeader=.true.   
        if(line(1:10).eq.':EndHeader')foundEndHeader=.true.   
      end do
      
!     Check that the # of locations = # specified in BeginLine
      if(xcount_temp.ne.j)then
        print*,'Error - no of stations',j 
        print*,'does not match the'
        print*,'the number in "beginLine"',xcount_temp
        stop 'Program terminated in read_ts5 @ 45'
      endif
      
c      xcount_temp=j   ! # of  data points
        
      if(allocated(inarray))then    ! check that inarray is allocated
        deallocate(inarray,stat=iDeallocate)
        if(iDeallocate.ne.0) then
          print*,'Warning: error with deallocation of intarray'
        endif
      endif	
c      allocate(inarray(367,ycount_temp),stat=iAllocate)
      allocate(inarray(367,xcount_temp),stat=iAllocate)
      if(iAllocate.ne.0) then
        STOP 'Error with allocation of inarray in read_r2c'      
      end if
      
!     rev. 10.1.88 May   23/17  - NK: Fixed Julian_day problems for iso R/W
      i=0       
      ios=0
!      do while(.not.eof(99))
      do while(ios.eq.0)
        i=i+1    !  line or row #
        read(99,99001,iostat=ios)(chr(j),j=1,256)
99001   format(256a1)        
        do j=1,256
            if(chr(j).eq.'/')chr(j)='-'
        end do
        write(line,*)(chr(j),j=1,256)
c        read(99,*)yyyymmdd(i),(inarray(i,j),j=1,xcount_temp)
        read(line,*)
     *      yyyymmdd(i),(inarray(i,j),j=1,xcount_temp)
c        print*,i,yyyymmdd(i),(inarray(i,j),j=1,3)
      end do
c      ycount_temp=i    ! use for older compiler
      ycount_temp=i-1
d      print*,'# data lines read in read_ts5 =',ycount_temp
d      print*,'# columns read =',xcount_temp
      
c      if(mod(year_now,4).eq.0)then
c!       check if there are 366 data points for this leap year
        if(ycount_temp.ne.366)then
          print*,'Insufficient # of data points '
          print*,'Found ',ycount_temp,' lines of data'
          print*
          print*,'Attempting to read ',fln(flnNum)(1:40)
c          pause 'Program paused in read ts5 @ 114'
c        endif
c      else
c!       check if there are 365 data points for this leap year
c        if(ycount_temp.ne.365)then
c          print*
c          print*,'Insufficient # of data points in this year'
c          print*,'Found ',ycount_temp,' lines of data'
c          print*,'Attempting to read ',fln(flnNum)(1:40)
c          pause 'Program paused in read ts5 @ 170'
c        endif
      endif

  999 RETURN

      END SUBROUTINE read_ts5

