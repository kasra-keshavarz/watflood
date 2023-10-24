      SUBROUTINE write_lzsinit

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
!***********************************************************************
!     created Jan. 16/07 nk
!***********************************************************************

      use area_watflood
	implicit none

      INTEGER       :: un,fn,ii,no_signf,hour_no,hours_togo,
     *                 no_frames,no_classes,frame_no,class_no,
     *                 i,j,n,ios
	character(20) :: junk
c      character(40) :: attribute_name
c      character(5)  :: data_source
c	REAL          :: conv
      CHARACTER(10) :: time
	character(512):: line(1000)
      CHARACTER(8)  :: cday
	logical       :: exists

!     rev. 9.3.08  Jan.  15/07  - NK: added lzs_init_new.r2c output to sub.for

!     write the lzs_init_new.r2c file
!     added Jan. 15/07 nk
!      author='spl.exe                               '
      coordsys_temp=coordsys1
      zone_temp=zone1
	datum_temp=datum1
	xorigin_temp=xorigin
	yorigin_temp=yorigin
	xcount_temp=xcount
	ycount_temp=ycount
	xdelta_temp=xdelta
	ydelta_temp=ydelta
	name='Lower Zone Storage                      '
	attribute_units='mm                                      ' 
!     attribute_type='Lower Zone Storage                      '
      attribute_name='Lower Zone Storage                      '
      unit_conversion=1.0
!      init_heat_deficit=deffactor  
      source_file_name=fln(10)    ! met file
!     convert local variables to write module variables

      do i=1,ycount
        do j=1,xcount
!         set values outside grid
          outarray(i,j)=0.0
        end do
      end do          

      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        outarray(i,j)=lzs(n)
      end do

c      open(unit=99,file='junk1',status='unknown')
      write(line,99005)year1,mo1,day_now-1,hour_now
99005 format(i4,'/',i2,'/',i2,2x,i2,':00:00')
c      rewind 99
	read(line,99006)startdate,starttime
99006 format(2a10)
c      close(unit=99,status='delete')

c      inquire(FILE='lzs_init.r2c',EXIST=exists)
c      if(exists)then
c!       write a new lzs file - keep old one
c        fln(99)='lzs_init_new.r2c'
c      else
        fln(99)='lzs_init.r2c'     
c      endif                                 

c!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c      call write_r2c(99,99,0,1,0,1,1)   
c      call write_r2c(99,99,0,0,1,1,1)   
c!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     rev. 9.9.22  Jul.  29/14  - NK: Fixed basin no assignment in flowinit.f
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call write_r2c(99,99,0,1,0,0,1)   
      call write_r2c(99,99,0,1,0,1,1)   
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      return

      END SUBROUTINE write_lzsinit
