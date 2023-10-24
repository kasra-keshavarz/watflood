      SUBROUTINE write_tb0(un,fn,nfg,ng,no_signf)

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
      
!     rev. 9.1.64  Jul. 11/08  - NK: Coded up for tb0
!     rev. 9.8.53  Mar.  20/13  - NK: Add Lake St. Joseph diversion algorithm to REROUT.f

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

      use area_watflood
	implicit none

      INTEGER       :: un,fn,ii,no_signf,ng,nfg,i,j,l,ios
	character(20) :: junk
      CHARACTER(10) :: time
      CHARACTER(8)  :: cday

!     FIRST TIME THROUGH THIS SUBROUTINE ONLY
!     OPEN OUTPUT FILE AND WRITE HEADER

!     write the header

!     FILE NAMES AND UNIT NUMBERS DIFFER BY 30
      iopt=1

      if(nfg.eq.1)then
!     GreenKenue uses LatLong - code below uses LATLONG
      if(coordsys_temp.eq.'LatLong   ')coordsys_temp='LATLONG   '
      if(coordsys_temp.eq.'Cartesian ')coordsys_temp='CARTESIAN '

      open(unit=un,file=fln(fn),status='unknown',iostat=ios)
      if(ios.ne.0)then    ! added Nov. 10/14  nk
        print*
        print*,'Unable to open file',fln(fn)(1:60)
        print*,'Possible cause(s):'
        print*,'file in use by another application'
        print*,'or target directory does not exist'
        stop 'Program aborted in write_tb0.f @ 52'
      else
         print*,'Opened unit=',un,' filename=',fln(fn)(1:40)
      endif
      write(un,3005)'########################################'
      write(un,3005)':FileType tb0  ASCII  EnSim 1.0         '
      write(un,3005)'#                                       '
	write(un,3005)'# DataType               Time Series    '
      write(un,3005)'#                                       '
      write(un,3005)':Application             WATFLOOD       '
	write(un,3005)':Version                 2.1.23         '
	write(un,3020)':WrittenBy          ',author
      call date_and_time(cday,time)
	write(un,3010)':CreationDate       ',
     *     cday(1:4),cday(5:6),cday(7:8),time(1:2),time(3:4)
3010  format(a20,a4,'-',a2,'-',a2,2x,a2,':',a2)
      write(un,3005)'#                                       '
	write(un,3005)'#---------------------------------------'
      write(un,3020)':SourceFile         ',source_file_name                                    '
      write(un,3005)'#                                       '
	write(un,3020)':Name               ',name
      write(un,3005)'#                                       '
!     REV. 10.1.38 Jul   28/16  - NK: Added noDataValue to WFO & tb0 files
      write(un,3008)':NoDataValue          ',noDataValue
	write(un,3004)':Projection         ',coordsys_temp
	if(coordsys_temp(1:7).eq.'LATLONG')then
	  write(un,3004)':Ellipsoid          ',datum_temp
	endif
	if(coordsys_temp(1:7).eq.'LatLong')then
	  write(un,3004)':Ellipsoid          ',datum_temp
	endif
!     rev. 9.9.62  Mar.  21/15  - NK: Change zone from character to integer
	if(coordsys_temp(1:3).eq.'UTM')then
	  write(un,3004)':Ellipsoid          ',datum_temp
        write(un,3001)':Zone               ',zone_temp
	endif
      write(un,3005)'#                                       '
!     fixed format for the date      
      if(mo_start.le.9.and.hour_start.le.9)then
        write(un,3024)':StartDate          ',
     *     year_start,mo_start,day_start,hour_start
3024    format(a20,i4,'/0',i1,'/0',i1,'  ',i2,':00:00')  
      elseif(mo_start.le.9.and.hour_start.gt.9)then
        write(un,3025)':StartDate          ',
     *     year_start,mo_start,day_start,hour_start
3025    format(a20,i4,'/0',i1,'/',i2,'  ',i2,':00:00')  
      elseif(mo_start.gt.9.and.hour_start.le.9)then
        write(un,3026)':StartDate          ',
     *     year_start,mo_start,day_start,hour_start
3026    format(a20,i4,'/',i2,'/0',i1,'  ',i2,':00:00')  
      else    
        write(un,3027)':StartDate          ',
     *     year_start,mo_start,day_start,hour_start
      endif
3027  format(a20,i4,'/',i2,'/',i2,'  ',i2,':00:00')      
      write(un,3028)':deltaT             ',deltat_temp
3028  format(a20,i10)   
      write(un,3003)':UnitConversion     ',unit_conversion
      write(un,3005)'#                                       '
	write(un,3005)':ColumnMetaData                         '
      write(un,3021)'   :ColumnUnits     ',(column_units(l),l=1,ng)
	write(un,3021)'   :ColumnType      ',(column_type(l),l=1,ng)
	write(un,3021)'   :ColumnName      ',(gname(l),l=1,ng)
	if(coordsys_temp.eq.'LATLONG   ')then
	  write(un,3023)'   :ColumnLocationX ',(xsta(l),l=1,ng)
	  write(un,3023)'   :ColumnLocationY ',(ysta(l),l=1,ng)
	elseif(coordsys_temp(1:7).eq.'LatLong')then
	  write(un,3023)'   :ColumnLocationX ',(xsta(l),l=1,ng)
	  write(un,3023)'   :ColumnLocationY ',(ysta(l),l=1,ng)
	elseif(coordsys_temp(1:7).eq.'Lambert')then
	  write(un,3029)'   :ColumnLocationX ',(xsta(l),l=1,ng)
	  write(un,3029)'   :ColumnLocationY ',(ysta(l),l=1,ng)
	else
	  write(un,3022)'   :ColumnLocationX ',(xsta(l),l=1,ng)
	  write(un,3022)'   :ColumnLocationY ',(ysta(l),l=1,ng)
	endif
	
!     rev. 9.8.53  Mar.  20/13  - NK: Add Lake St. Joseph diversion algorithm to REROUT.f
!       this is here because the diversion file needs two sets of lat-long
!       one for the source and another for the delivery location      
      if(wrtdiverflg)then
	  if(coordsys_temp.eq.'LATLONG   ')then
	    write(un,3023)'   :ColumnLocationX1 ',(xsta1(l),l=1,ng)
	    write(un,3023)'   :ColumnLocationY1 ',(ysta1(l),l=1,ng)
	  elseif(coordsys_temp(1:7).eq.'LatLong')then
	    write(un,3023)'   :ColumnLocationX1 ',(xsta1(l),l=1,ng)
	    write(un,3023)'   :ColumnLocationY1 ',(ysta1(l),l=1,ng)
	  else
	    write(un,3022)'   :ColumnLocationX1 ',(xsta1(l),l=1,ng)
	    write(un,3022)'   :ColumnLocationY1 ',(ysta1(l),l=1,ng)
  	  endif
      endif        

!     rev. 9.9.66  Apr.  `3/15  - NK: Added options to write_tb0 str files
      if(wrtstrflg)then
        write(un,3022)'   :Coeff1         ',(Coeff1(l),l=1,ng)
        write(un,3022)'   :Coeff2         ',(Coeff2(l),l=1,ng)
        write(un,3022)'   :Coeff3         ',(Coeff3(l),l=1,ng)
        write(un,3022)'   :Coeff4         ',(Coeff4(l),l=1,ng)
        write(un,3022)'   :Value1         ',(Value1(l),l=1,ng)
      endif
        
     	write(un,3005)':endHeader                              '

      return
      
      endif

!     Note:  ng=xcount_temp


!       NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE 
!       The r2c grids are written upside down: 
!                 south on top, north on bottom    !!!!!!!!!!
        if(no_signf.eq.0)then
          do i=1,ycount_temp
            write(un,1300)(outarray(i,j),j=1,xcount_temp)
1300        format(20x,9999(1x,f5.0))
          end do          
        elseif(no_signf.eq.1)then   ! swe, precip and temp
          do i=1,ycount_temp
            write(un,1301)(outarray(i,j),j=1,xcount_temp)
1301        format(20x,9999(1x,f5.1))
          end do          
        elseif(no_signf.eq.2)then   ! swe, precip and temp
          do i=1,ycount_temp
            write(un,1302)(outarray(i,j),j=1,xcount_temp)
1302        format(20x,9999(1x,f5.2))
          end do          
        elseif(no_signf.eq.3)then   ! swe, precip and temp
          do i=1,ycount_temp
            write(un,1303)(outarray(i,j),j=1,xcount_temp)
1303        format(20x,9999(1x,f6.3))
          end do          
        elseif(no_signf.eq.4)then   ! swe, precip and temp
          do i=1,ycount_temp
            write(un,1304)(outarray(i,j),j=1,xcount_temp)
1304        format(20x,9999(1x,f7.4))
          end do          
        elseif(no_signf.eq.5)then   ! swe, precip and temp
          do i=1,ycount_temp
            write(un,1305)(outarray(i,j),j=1,xcount_temp)
1305        format(20x,9999(1x,f8.5))
          end do          
        elseif(no_signf.eq.6)then   ! swe, precip and temp
          do i=1,ycount_temp
            write(un,1306)(outarray(i,j),j=1,xcount_temp)
1306        format(20x,9999(1x,f9.6))
          end do          
        elseif(no_signf.eq.8)then   ! flow
          do i=1,ycount_temp
            write(un,1308)(outarray(i,j),j=1,xcount_temp)
1308        format(20x,9999(1x,e10.3))
          end do          
        elseif(no_signf.eq.10)then   ! drain, dsnow
          do i=1,ycount_temp
            write(un,1310)(outarray(i,j),j=1,xcount_temp)
1310        format(20x,9999(1x,f10.0))
          end do          
        elseif(no_signf.eq.11)then   ! drain, dsnow
          do i=1,ycount_temp
            write(un,1311)(outarray(i,j),j=1,xcount_temp)
1311        format(20x,9999(1x,f8.3))
          end do
                    
        elseif(no_signf.eq.12)then   ! elv_mean,elv_max
          if(xcount_temp.eq.1)then       
            do i=1,ycount_temp
              write(un,1312)outarray(i,1)
            end do
          else
            do i=1,ycount_temp
              write(un,1312)(outarray(i,j),j=1,xcount_temp)
1312        format(20x,9999(1x,f9.3))
             end do          
          endif
          
        elseif(no_signf.eq.13)then   ! NWM precip
          do i=1,ycount_temp
            write(un,1313)(outarray(i,j),j=1,xcount_temp)
1313        format(20x,9999(1x,f9.3))
          end do
                    
        elseif(no_signf.eq.14)then   ! str files
          do i=1,ycount_temp
            write(un,1314)(outarray(i,j),j=1,xcount_temp)
1314        format(20x,9999(1x,f11.3))
          end do
                    
        else                        ! init soil moisture
          do i=1,ycount_temp
            write(un,1307)(outarray(i,j),j=1,xcount_temp)
1307        format(20x,9999(1x,e12.6))
          end do          
        endif


c        close(unit=un,status='keep')
c        write(51,*)'Closed unit ',un,' Filename=  ',fln(fn)
c        write(*,*)'Closed unit ',un,' Filename=  ',fln(fn)

      return

! FORMATS
 3000 format(a10,i5)
 3001 format(a20,i16)
 3002 format(2a20)
 3003 format(a20,f16.7)
 3004 format(a20,a10,2x,a10)
 3005 format(a40)
 3006 format(a3,a10)
 3007 format(a14,i5,a6,i5)
 3008 format(a20,f10.1)
 3012 format(a9)
 3020 format(a20,a40)
 3021 format(a23,999(a9,3x))
 3022 format(a20,999f12.2)
 3023 format(a20,999f12.4)
 3029 format(a20,999f10.0)
     
       END SUBROUTINE write_tb0









