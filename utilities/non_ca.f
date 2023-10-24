      subroutine non_ca()
      
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
!     This s/r calculates the nca in each cell
!     Revised Jan. 21, 2014 for the new bsn_response.txt file ver. 4      

      use area_watflood
      USE area17   ! see note in bsn.f

      implicit none
      
      CHARACTER(30) :: fn1,junk,class_name(99),cover_file_name(999)
      INTEGER      :: ii,jj,z,n,j,i,xcount_nca,ycount_nca
      Real(4)      :: xorigin_nca,yorigin_nca,xdelta_nca,ydelta_nca,
     *                northlimit,eastlimit
      logical      :: exists
      real*4       :: p_long,p_lat

      integer, dimension(:,:), allocatable :: nca_raw
      real,    dimension(:,:), allocatable :: hit_count,nul_count
      
      integer  :: ios
      integer  :: IALLOCATESTATUS

      print*,'********************************************************'
      print*,'*                                                      *'
      print*,'*                  WATFLOOD (TM)                       *'
      print*,'*                                                      *'
      print*,'*           Subroutine non_ca   June 9, 2011           *'
      print*,'*               Revised  Jan. 21, 2015                 *'
      print*,'*                                                      *'
      print*,'*              (c) N. Kouwen, 1972-2015                *'
      print*,'*                                                      *'
      print*,'********************************************************'
      print*
 
      fln(200)='nca.r2s'

      INQUIRE(FILE=fln(200),EXIST=exists)
      if(exists)then
c        effective_nca=1.00
        print*
        print*,'nca.r2s file found'
        print*,'non-contributing areas will be subtracted from frac for'
        print*,'for each cell'
        print*,'You can not subtract nca from frac if you want to split'
        print*,'land cover classes into contributing & non-contributing'
        print*

!       CHECK FOR EXIST DONE ABOVE        
        open(unit=200,file=fln(200),status='unknown',iostat=ios)
        if(ios.ne.0)then
          print*,'problems opening file= ',fln(200)(1:40)
          print*
          stop 'Program aborted in non_da @ 49'
        else
          print*,'opened input file:',fln(200)(1:40)
          print*
        endif
c	  non_ca_flg=.true.
	  allocate(nca(ycount,xcount),
     *	 nul_count(ycount,xcount),hit_count(ycount,xcount),
     *     stat=iAllocateStatus)
        if (iAllocateStatus .ne. 0) STOP 
     *       '**Allocation failed for nca in non_da.f @ 57**' 
d        print*,'nca() allocated for',ycount,xcount
d        pause 'in non_ca  @ 84'
	  do i=1,ycount
          do j=1,xcount
            nca(i,j)=0
	      nul_count(i,j)=0
	      hit_count(i,j)=0
          end do
        end do

      else
	  non_ca_flg=.false.
        return
      endif

!     read the non contributing area r2s file
      REWIND 200
	do i=1,12
	  read(200,*)  !skip over first 12 lines
	end do

      read(200,*)junk,coordsys2  !<<<<<<<<<<<
      if(coordsys2.eq.'LatLong')coordsys2='LATLONG'
      print*,junk,coordsys2

      if(coordsys2(1:3).eq.'UTM')then
        read(200,*)junk,datum2  !<<<<<<<<<<<
        print*,junk,datum2
        read(200,*)junk,zone2  !<<<<<<<<<<<
        print*,junk,zone2
      elseif(coordsys2(1:9).eq.'CARTESIAN')then
!       datum and zone are undefined & not needed
        datum2='unknown'
        zone2='unknown'
      elseif(coordsys2(1:7).eq.'LATLONG')then
        read(200,*)junk,datum2  !<<<<<<<<<<<
        print*,junk,datum2
      else
        print*,' Valid format (UTM, CARTESIAN, or LATLONG) not found'
        print*,'coordsys found=',coordsys2
        stop 'Program aborted in non_ca@ 89'
      end if

!     READ NEW FORMAT FOR HEADER
        read(200,*)junk    ! reads a line starting with #
        print*,junk
        read(200,*)junk,xorigin_nca  !<<<<<<<<<<<
        print*,junk,xorigin_nca
        read(200,*)junk,yorigin_nca  !<<<<<<<<<<<
        print*,junk,yorigin_nca
        read(200,*)
        read(200,*)junk,xcount_nca  !<<<<<<<<<<<
        print*,junk,xcount_nca
        read(200,*)junk,ycount_nca  !<<<<<<<<<<<
        print*,junk,ycount_nca
        read(200,*)junk,xdelta_nca  !<<<<<<<<<<<
        print*,junk,xdelta_nca
        read(200,*)junk,ydelta_nca  !<<<<<<<<<<<
        print*,junk,ydelta_nca
      close(unit=5,status='keep')

	allocate(nca_raw(ycount_nca,xcount_nca),stat=iAllocateStatus)
        if (iAllocateStatus .ne. 0) STOP 
     *       '**Allocation failed for nca in non_da.f @ 57**' 
	do i=1,ycount_nca
        do j=1,xcount_nca
            nca_raw(i,j)=0
        end do
      end do

      rewind 200
	read(200,*)junk
c	print*,junk
	do while(junk(1:10).ne.':EndHeader')
  	  read(200,*)junk
c        print*,junk
	end do

      print*,' reading the nca file '
	do i=1,ycount_nca
	  read(200,*)(nca_raw(i,j),j=1,xcount_nca)
c	  print*,i,'row',(nca_raw(i,j),j=1,3),'....'
	end do
      close(unit=200,status='keep')

!     calculate limits

!     Calculate top     
      northlimit=yorigin_nca+ycount_nca*ydelta_nca
      eastlimit=xorigin_nca+xcount_nca*xdelta_nca

	print*,'Grid extents of non-contributing areas:'
      print*,'xorigin_nca',xorigin_nca
      print*,'eastlimit',eastlimit
      print*,'yorigin_nca',yorigin_nca
      print*,'northlimit',northlimit
      print*

      print*,' counting pixels '
!     find the total number of pixels of each class in each grid
      write(51,*)'   class   col  row   total   # in class'
      do ii=1,ycount_nca
        do jj=1,xcount_nca
          p_long=xorigin_nca+jj*xdelta_nca  ! pixel long
	    p_lat=yorigin_nca+ii*ydelta_nca   ! pixel lat 
c	print*,ii,jj,p_lat,p_long
!         find the host grid
          j=int((p_long-xorigin)/xdelta)+1
	    i=int((p_lat-yorigin)/ydelta)+1
!         check to see if oixel is in the map domain
	    if(i.ge.1.and.i.le.ycount.and.j.ge.1.and.j.le.xcount)then
	      if(nca_raw(ii,jj).eq.1)then
!             count hits
	        hit_count(i,j)=hit_count(i,j)+1
	      else
	        nul_count(i,j)=nul_count(i,j)+1
	      endif
	    endif
	  end do
	end do

      print*,' calculating the nca on each cell '
	do i=1,ycount
	  do j=1,xcount
c	print*,i,j,nul_count(i,j),hit_count(i,j)
	    if((nul_count(i,j)+hit_count(i,j)).gt.0.0)then
	      nca(i,j)=nul_count(i,j)/(nul_count(i,j)+hit_count(i,j))
	    else
	      nca(i,j)=0.0
	    endif
	  end do
	end do
              
      print*,' writing the nca.xyz file '
	open(unit=200,file='nca.xyz',status='unknown')
	do i=1,ycount
	  do j=1,xcount
	    write(200,*)xorigin+float(j)*xdelta-xdelta/2.0,
     *           yorigin+float(i)*ydelta-ydelta/2.0,nca(i,j)
	  end do
	end do
      print*,'nca.xyx written'
	print*
      print*, ' done computing non contributing areas'
      
      return

      END SUBROUTINE non_ca
