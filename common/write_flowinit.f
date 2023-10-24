      SUBROUTINE write_flowinit

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
!     rev. 9.3.05  Nov.  13/06  - NK: adder write_flowinit.for to flowinit.for

!     "author" must be specifief before calling this s/r

      use area_watflood
	implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      logical    :: exists
	integer(4)    :: i,j,n,ios,k,flnnum
      CHARACTER(30) :: sys_command
      CHARACTER(10) :: time
      CHARACTER(8)  :: cday
	character(1)  :: answer
      logical       :: Lnan
!     write the init file for watroute in r2c format



!     rev. 10.2.50 Mar.  22/19  - NK: New resume files written to \resume\*.*
	  fln(99)='resume\flow_init.r2c'
        open(99,file=fln(99),status='unknown',iostat=ios)
        if(ios.ne.0)then    ! added Nov. 10/14  nk
          print*
         write(63,*)'Unable to open file',fln(99)(1:40)
         write(63,*)'Possible cause(s):'
         write(63,*)'file in use by another application'
         write(63,*)'or target directory does not exist'
          stop 'Program aborted in write_flowinit.f @ 47'
        endif

      write(99,3005)'########################################'
      write(99,3005)':FileType r2c  ASCII  EnSim 1.0         '
      write(99,3005)'#                                       '
	write(99,3005)'# DataType               2D Rect Cell   '
      write(99,3005)'#                                       '
      write(99,3005)':Application             EnSimHydrologic'
	write(99,3005)':Version                 2.1.23         '
	write(99,3002)':WrittenBy          ',author
      call date_and_time(cday,time)
	write(99,3010)':CreationDate       ',
     *       cday(1:4),cday(5:6),cday(7:8),time(1:2),time(3:4)
3010  format(a20,a4,'/',a2,'/',a2,2x,a2,':',a2)
      write(99,3005)'#                                       '
        if(tbcflg.eq.'y')then
	    write(99,3011)'#start y/m/d/h:',
     *                   year_start,mo_start,day_start,hour_start
!     REV. 10.1.26 Mar.  23/16  - NK: Fixed comment for spinup period 
	    write(99,3011)'#end   y/m/d/h:',
     *                   year_now,month_now,day_now,hour_now
 3011 format(a15,4i10)
        endif
	write(99,3005)'#---------------------------------------'
      write(99,3020)':SourceFileName     ',fln(6)
      write(99,3001)':ClassCount         ',11
      write(99,3005)'#                                       '
	write(99,3004)':Projection         ',coordsys1
!     rev. 9.9.62  Mar.  21/15  - NK: Change zone from character to integer
	if(coordsys1.eq.'UTM       ')then
          write(99,3001)':Zone               ',zone1
	    write(99,3004)':Ellipsoid          ',datum1
	endif
	if(coordsys1.eq.'LATLONG   ')then
	    write(99,3004)':Ellipsoid          ',datum1
      endif
      write(99,3005)'#                                       '
      write(99,3003)':xOrigin            ',xorigin
      write(99,3003)':yOrigin            ',yorigin
      write(99,3005)'#                                       '
      write(99,3008)':AttributeName 1 qi1          ' 
      write(99,3008)':AttributeName 2 qo1          '  
      write(99,3008)':AttributeName 3 store1       '
      write(99,3008)':AttributeName 4 over         ' 
      write(99,3008)':AttributeName 5 lzs          ' 
      write(99,3008)':AttributeName 6 qOld         ' 
      if(wetflg.eq.'y')then
        write(99,3008)':AttributeName 7 qiwet        ' 
        write(99,3008)':AttributeName 8 qowet        ' 
        write(99,3008)':AttributeName 9 wstore       ' 
        write(99,3008)':AttributeName 10 hcha         ' 
        write(99,3008)':AttributeName 11 hwet        ' 
      endif

      write(99,3005)'#                                       '
      write(99,3001)':xCount             ',xcount
      write(99,3001)':yCount             ',ycount
      write(99,3003)':xDelta             ',xdelta
      write(99,3003)':yDelta             ',ydelta
      write(99,3005)'#                                       '
      write(99,3005)':EndHeader                              '



!	  if(wetflg.eq.'y')then
!          do n=1,na  
!            write(99,*)n,qiwet1(n),qiwet2(n),qowet1(n),qowet2(n),
!     *       wstore1(n),wstore2(n),qswrain(n),qswevp(n),
!     *       qin(n),hcha1(n),hcha2(n),hwet1(n),hwet2(n)
!          end do
!	  endif


 3000 format(a10,i5)
 3001 format(a20,i16)
 3002 format(2a20)
 3003 format(a20,f16.7)
 3004 format(a20,a10,2x,a10)
 3005 format(a40)
 3006 format(a3,a10)
 3007 format(a14,i5,a6,i5)
 3008 format(a30)
c 3012 format(a9)
 3020 format(a20,a40)
 	if(iopt.eq.2)print*, 'in write_flowinit at 1302'

!     initialize p() to male sure there is no junk in the unused grids
      do i=1,ycount
	  do j=1,xcount
	    p(i,j)=0.0
	  end do
	end do
	k=0

!     Initial grid inflow
      k=k+1
	p(1,1)=float(k)
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        p(i,j)=qi1(n)
        Lnan=isnan(p(i,j))
        if(Lnan)then
            print*,'Warning: qi1 @ row,col,grid# ',i,j,n
            print*,'Warning: is not a number (NAN)'
            print*,'qi1(n) set to ',0.0
            p(i,j)=0.0
        endif
      end do
	do i=1,ycount
        write(99,4009)(p(i,j),j=1,xcount)
      end do

!     Initial grid outflow in m^3/s
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        p(i,j)=qo1(n)
        Lnan=isnan(p(i,j))
        if(Lnan)then
            print*,'Warning: qo1 @ row,col,grid# ',i,j,n
            print*,'Warning: is not a number (NAN)'
            print*,'qo1(n) set to ',-999.0
            p(i,j)=0.0
        endif
      end do
	do i=1,ycount
        write(99,4009)(p(i,j),j=1,xcount)
      end do

!     Initial grid channel storage in m^3
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        p(i,j)=store1(n)
        Lnan=isnan(p(i,j))
        if(Lnan)then
            print*,'Warning: store1 @ row,col,grid# ',i,j,n
            print*,'Warning: is not a number (NAN)'
            print*,'store1(n) set to ',-999.0
            p(i,j)=0.0
        endif
      end do
	do i=1,ycount
        write(99,4009)(p(i,j),j=1,xcount)
      end do

!     Initial overbank storage in m^3(used as conditional so include here)
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        p(i,j)=over(n)
        Lnan=isnan(p(i,j))
        if(Lnan)then
            print*,'Warning: over @ row,col,grid# ',i,j,n
            print*,'Warning: is not a number (NAN)'
            print*,'over(n) set to ',-999.0
            p(i,j)=0.0
        endif
      end do
	do i=1,ycount
        write(99,4009)(p(i,j),j=1,xcount)
      end do

!     Initial lower zone storage in mm
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        p(i,j)=lzs(n)
        Lnan=isnan(p(i,j))
        if(Lnan)then
            print*,'Warning: lzs @ row,col,grid# ',i,j,n
            print*,'Warning: is not a number (NAN)'
            print*,'lzs(n) set to ',-999.0
            p(i,j)=0.0
        endif
      end do
	do i=1,ycount
        write(99,4009)(p(i,j),j=1,xcount)
      end do

!     rev. 10.2.44 Jan.  21/19  - NK: Added qOld() to allow previous weighted outflow on the resume file
!     Initial old flow     m**3
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        p(i,j)=qOld(n)
        Lnan=isnan(p(i,j))
        if(Lnan)then
            print*,'Warning: qold @ row,col,grid# ',i,j,n
            print*,'Warning: is not a number (NAN)'
            print*,'qold(n) set to ',-999.0
            p(i,j)=0.0
        endif
      end do
	do i=1,ycount
        write(99,4009)(p(i,j),j=1,xcount)
      end do

      if(wetflg.eq.'y')then

!       Initial wetland inflow        cms
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          p(i,j)=qiwet2(n)
        Lnan=isnan(p(i,j))
        if(Lnan)then
            print*,'Warning: qiwet2 @ row,col,grid# ',i,j,n
            print*,'Warning: is not a number (NAN)'
            print*,'qiwer2(n) set to ',-999.0
            p(i,j)=0.0
        endif
        end do
	  do i=1,ycount
          write(99,4009)(p(i,j),j=1,xcount)
        end do

!       Initial wetland outflow
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          p(i,j)=qowet2(n)
        Lnan=isnan(p(i,j))
        if(Lnan)then
            print*,'Warning: qowet2 @ row,col,grid# ',i,j,n
            print*,'Warning: is not a number (NAN)'
            print*,'qowet21(n) set to ',-999.0
            p(i,j)=0.0
        endif
        end do
	  do i=1,ycount
          write(99,4009)(p(i,j),j=1,xcount)
        end do

!       Initial wetland storage  m**3
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          p(i,j)=wstore2(n)
        Lnan=isnan(p(i,j))
        if(Lnan)then
            print*,'Warning: wstore2 @ row,col,grid# ',i,j,n
            print*,'Warning: is not a number (NAN)'
            print*,'wstore2(n) set to ',-999.0
            p(i,j)=0.0
        endif
        end do
	  do i=1,ycount
          write(99,4009)(p(i,j),j=1,xcount)
        end do

!       Initial channel  m
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          p(i,j)=hcha2(n)
        Lnan=isnan(p(i,j))
        if(Lnan)then
            print*,'Warning: hcha2 @ row,col,grid# ',i,j,n
            print*,'Warning: is not a number (NAN)'
            print*,'hcha2(n) set to ',-999.0
            p(i,j)=0.0
        endif
        end do
	  do i=1,ycount
          write(99,4009)(p(i,j),j=1,xcount)
        end do

!       Initial wetland depth    m
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          p(i,j)=hwet2(n)
        Lnan=isnan(p(i,j))
        if(Lnan)then
            print*,'Warning: hwet2 @ row,col,grid# ',i,j,n
            print*,'Warning: is not a number (NAN)'
            print*,'hwet2(n) set to ',-999.0
            p(i,j)=0.0
        endif
        end do
	  do i=1,ycount
          write(99,4009)(p(i,j),j=1,xcount)
        end do

      endif

      close(unit=99,status='keep')
      print*
      write(63,*)fln(99)(1:40),' written '
	print*

 4003 format(999(' ',i5))
!     rev. 10.5.05 Mar.  25/23  = NK Fixed reporting in wfo file & roundoff in the resume files
 4009 format(999(' ',E15.7))

      return

      end SUBROUTINE write_flowinit

