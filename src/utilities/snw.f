       PROGRAM snw

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
!      include 'debug.for'

!************************************************************************
! SNWA - THIS PROGRAM GENERATES A SNOW-WATER CONTENT DATA FILE FROM 
!        SNOW COURSE DATA.
!
!	subroutines called:
!		rdevt
!		weight
!		distr
!
! Revision for weq by class - FEB 20/94
!
! Modified by Tricia Stadnyk - September 2000
! Converted common blocks to modules and added dynamically allocated
! run-time arrays as part of Fortran 90 conversion.
!
!  VERSION 2      Dec. 2000  - TS: Added dynamic memory allocation
!  VERSION 2      ???. ????  - NK: Added radius of influence
!  VERSION 3      Mar. 2018  - NK: Read rad. of inf. from the pt2 file
!
! - list of arguments:
!
!   o - rain( , )  real*4   rain gauge precipitation data
!   i - ng         int      number of snow courses
!   i - nhr=ncover int      number of land cover classes
!   i - conv       real*4   conversion factor
!   i - fln( )     char*12  file names
!   o - sta( )     cmplx    gauge locations (grid square numbers)
!   o - ewg( )     real*4   east-west utm coordinates of gauges
!   o - sng( )     real*4   north-south utm coordinates of gauges
!   o - xcount        int      number of east-west grid squares
!   o - ycount        int      number of north-south grid squares
!   o - ng         int      number of snow courses
!   o - radw       int      utm coordinate of west side of grid
!   o - rads       int      utm coordinate of south side of grid
!   o - gname( )   char*12  rain gauge names
!   i - fln( )     char*30  file names
!       rain( )    real*4   water equivalent in snow pack (course,cover)
!
!************************************************************************



!     fix fix
!     in the future, the whole thing should be made floating point 
!     for the  x-y coordinates:  distr & weight

!     first fix weight & distr
!     then the rest:  snw, ragmet and tmp
      use area_watflood

C///////////////////////// 
C// Added by Dave
	USE EF_module
C// End Dave addition
C/////////////////////////

      CHARACTER(14) :: date
      DIMENSION     :: smc5(16) 
	INTEGER       :: iAllocateStatus,iDeallocateStatus,ii,value
      INTEGER(2)       status1
      CHARACTER(1)     stopflg,newfmtflg
      character(8)  :: buf
      character(20) :: junk,junk1,junk2
      character(22) :: in_fln,out_fln
      character(72) :: line
      character(10) :: fileformat,surveydate,surveytime
      logical       :: exists,use_buf,RIflag

      real, dimension(:),     allocatable :: ewg_fl,sng_fl
      
      
        write(6,6001)
        write(6,6002)
        write(6,6003)
        write(6,6004) 
        write(6,6005) 
        write(6,6006) 
        write(6,6007) 
        write(6,6008) 
        write(6,6009) 
        write(6,6002)
        write(6,6001)

 6001 format(1x,'*******************************************************
     ***********************')
 6002 format(1x,'*                                                      
     *                     *')
 6003 format(1x,'*  ww  ww         ww  a    tttttttt ffffff ll       ooo
     *      ooo    ddddd   *')
 6004 format(1x,'*   ww  ww       ww  aaa   tttttttt ffffff ll      oooo
     *o    ooooo   dddddd  *')
 6005 format(1x,'*    ww  ww     ww  aa aa     tt    ff     ll     oo   
     *oo  oo   oo  dd   dd *')
 6006 format(1x,'*     ww  ww   ww  aaa aaa    tt    ffff   ll     oo   
     *oo  oo   oo  dd   dd *')
 6007 format(1x,'*      ww  ww ww  aaaaaaaaa   tt    ffff   ll     oo   
     *oo  oo   oo  dd   dd *')
 6008 format(1x,'*       ww  www  aa       aa  tt    ff     llllll  oooo
     *o    ooooo   dddddd  *')
 6009 format(1x,'*        ww  w  aa         aa tt    ff     llllll   ooo
     *      ooo    ddddd   *')

	print*,'********************************************************'
	print*,'*                                                      *'
	print*,'*      SSSSSSS    NN     NN WW    WW            WW     *'
	print*,'*     SSSSSSSSS   NNN    NN  WW    WW          WW      *'
	print*,'*     SS          NNNN   NN   WW    WW        WW       *'
	print*,'*     SSSSSSSS    NN NN  NN    WW    WW      WW        *'
	print*,'*      SSSSSSSS   NN  NN NN     WW    WW    WW         *'
	print*,'*            SS   NN   NNNN      WW  WWWW  WW          *'
	print*,'*     SSSSSSSSS   NN    NNN       WWWW  WWWW           *'
	print*,'*      SSSSSSS    NN     NN        WW    WW            *'
	print*,'*                                                      *'
	print*,'*                  WATFLOOD (TM)                       *'
	print*,'*           Version 3       Nov. 26, 2018              *'
	print*,'*            (c) N. Kouwen, 1972-2018                  *'
	print*,'*                                                      *'
	print*,'********************************************************'
      print*

!     USE DFLIB
      CALL GETARG(1, buf, status1)
c      if(status1.ne.1)buf=' '
c      if(buf.eq.'1')then
c        stopflg='y'
c      else
c        stopflg='n'
c      endif
      print*,buf,status1
      use_buf=.false.
      if(status1.eq.8)then
        in_fln=buf(1:8)
        use_buf=.true.
        write(line,1300)'snow1\',buf(1:8),'_crs.pt2'
1300    format(a6,a8,a8)        
        read(line,1301)in_fln
1301    format(a22)
        write(line,1300)'snow1\',buf(1:8),'_swe.r2c'
        read(line,1301)out_fln
        print*,in_fln,' ',out_fln
      endif
        
! TS - ALLOCATION STATEMENT TO DYNAMICALLY ALLOCATE ARRAY FOR AREA12A
      allocate(fln(601),stat=iAllocateStatus)
      if (iAllocateStatus .ne. 0) STOP 
     *    '**Allocation failed for area12 in snwa l57**'

      iopt=0
	nw=2

      open(unit=51,file='snw_info.txt',status='unknown',iostat=ios)
        if(ios.ne.0)then
          print*,'fln= snw.msg','  ios=',ios
          print*
          STOP 'Program snw aborted @ 99'
        endif 

! PRINT EVENT TITLE:
      write(51,1000)

      calling_program_flg='snw       '

!     INPUT EVENT PARTICULARS
      fln(99)='event\event.evt'

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c      call rdevt(date,conv,scale,smc5,nr,nhf)
      call read_evt(date,conv1,scale,smc5,nhg,nhf)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      print*
      print*,'An ENSIM compatible r2c file will be created'
      print*
      print*,'SNW.EXE requires a bsnm.shd file to read he grid'
      print*,'dimensions. The snow1\yyyymmdd_swe.r2c file'
      print*,'will have the same dimensions as the watershed'
      print*,'file.'
      print*
      
      if(use_buf)then
          fln(35)=in_fln
      endif
      
      call find_filetype(1)

	if(IsFileTypeR2C(fln(1))) then
        print*
	  print*,'Gone to read '
	  print*,fln(1)(1:60)
	  print*,'with read_shd_ef'
	  print*
	  call read_shed_ef(31,1)	
	  print*,'Finished reading ',fln(1)(1:50)
	  print*,'# classes =',classcount
	  print*
	else
!     rev. 10.5.00 Nov.  07/22  = NK Modificatons for HYPE i/o 
!			call read_shed_hype(31,1)
!			call rdshed(31,1)
			call read_shed_ef(31,1)
      endif

      if(.not.allocated(s))then
       allocate(s(ycount,xcount),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *     'Error with allocation of s(,) in moist @ 213'
	endif

!     outarray is in areawfo
      allocate(outarray(ycount,xcount),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *    'Error with allocation of outarray arrays in snw'      

      snwver=6.0

      write(51,*)
      write(51,*)
      print*,'Please look in snw_info.txt for diagmostic information'
      print*
      write(51,*)'Echo data from ',fln(35)
	write(51,*)
!     
!     test if file already exists. if yes, do not overwrite.
!     unit=265  fln(35)- snow course data file          yyyymmdd_crs.pt2 
      INQUIRE(FILE=fln(35),EXIST=exists)
      IF(exists)THEN
      open(unit=265,file=fln(35),status='unknown',iostat=ios)
        if(ios.ne.0)then
          print*,'fln(35)=',fln(35),'  ios=',ios
          print*
          STOP 'Program snw aborted @ 88'
        endif 
	  print*,'Opened file name 35 ',fln(35)(1:40)
      else
        print*,'Attempting to open the file '
	  print*,fln(35)(1:60)
        print*,'but it is not found (in this location)'
        print*
        print*,'NOTICE:'
        print*,'This program no longer reads the .snw or .crs files'
        print*,'Please convert the event file to version 9.3'
        print*,'and convert the snow1\yymmdd.snw file and '
        print*,'                snow1\yyyymmdd.crs  file to the '
        print*,'                snow1\yyymmdd_crs.pt2 file and'
        print*,'run the snw.exe program.  Use the file names'
        print*,'snow1\yyyymmdd_crs.pt2  and'
        print*,'snow1\yyyymmdd_swe.r2c in the event file as per'
        print*,'example in tghe current manual'
        print*,'Sorry for the inconvenience  14/07/06  NK' 
        print*
        stop 'Program aborted in snw @ 94'
      endif


!     programmd for pt2 file type 14/07/06 nk

      call find_filetype(35)

      print*
      print*,'filetype=',filetype
	print*

!     Read the yyyymmdd_crs.pt2 file


!     Parse the header
      classcount_temp=0
      RIflag=.false.
      if(filetype.eq.'pt2')then
c        do while(junk1.ne.':Name     ')
c        junk1='    '  
        do while(.not.line(1:10).eq.':EndHeader')
          read(265,26500,iostat=ios)line
26500     format(a60)          
          print*,line
          if(line(1:11).eq.':Projection')read(line,*)junk,coordsys2
          if(line(1:5).eq.':Zone')read(line,*)junk,zone2
          if(line(1:5).eq.':Ellipsoid')read(line,*)junk,datum2
          if(coordsys2(1:7).eq.'LATLONG'.or.
     *           coordsys2(1:7).eq.'latlong')coordsys2='LatLong'
          if(line(1:15).eq.':UnitConversion')read(line,*)junk,conv
          if(line(1:16).eq.':InitHeatDeficit')read(line,*)junk,deffactor
          if(line(1:14).eq.':AttributeName')then
!           changed 2021-08-13              
c              read(line,*)junk,junk1,junk2
c              if(junk2(1:5).eq.'Class')then
                 classcount_temp=classcount_temp+1
                 write(51,*)junk1,classcount_temp
                 write(*,*)junk1,classcount_temp
c              endif
          endif
!          if(line(1:14).eq.':AttributeName')then
!              if(junk2(1:15).eq.'RadiusInfluence')then
!                  RIflag=.true.
!              endif
!          endif
        end do
        classcount_temp=classcount_temp-1
        print*,'par file rad_influence=', rad_influence
        print*,'classcount_temp =',classcount_temp

        ng=256
c	  classcount_temp=classcount_temp/2-2
	  write(51,*)'# classes found =',classcount_temp
        write(51,*)'# classcount in shd file = ',classcount
	  write(*,*)'# classes found =',classcount_temp
        write(*,*)'# classcount in shd file = ',classcount
! TS - ALLOCATION STATEMENT TO DYNAMICALLY ALLOCATE ARRAYS FOR AREA16A
        allocate(sng(ng),ewg(ng),xsta(ng),ysta(ng),gname(ng),    !dst(ng),
     *  rrain(ng,classcount),stat=iAllocateStatus)
        if (iAllocateStatus .ne. 0) STOP
     *    '**Allocation failed for area12 in snwa 238**'
     	 
        allocate(sng_fl(ng),ewg_fl(ng),
     *      radius_influence(ng),stat=iAllocateStatus)
        if (iAllocateStatus .ne. 0) STOP
     *    '**Allocation failed in snw 242**'	 

!       Initialize radius of influence:
!       This variable is read from the pt2 file
!       not the par file.        
        do ii=1,ng
            radius_influence(ii)=-1.0
        end do

        is=1
        do while(ios.eq.0)
          if(RIflag)then  
            read(265,*,iostat=ios)ewg_fl(is),sng_fl(is),gname(is),
     *         radius_influence(is),       
     *         (rrain(is,nh),nh=1,classcount_temp)
          else
            read(265,*,iostat=ios)ewg_fl(is),sng_fl(is),gname(is),
     *         (rrain(is,nh),nh=1,classcount_temp)
          endif
	    if(ios.eq.0)then
              write(51,*)'classs #',is
              if(classcount.gt.classcount_temp)then
                do nh=classcount_temp+1,classcount
                    rrain(is,nh)=rrain(is,classcount_temp)
                end do
	          write(51,*) ewg_fl(is),sng_fl(is),gname(is),
     *             radius_influence(is),       
     *             (rrain(is,nh),nh=1,classcount)
              endif
          endif
          is=is+1
        end do
        ng=is-2
        write(51,*)'ios=',ios,'       ',ng

      else
	  print*,'Old formats no longer accepted.'
	  print*,'Please create a pt2 format file'
        print*
	  stop 'Program aborted in snw @ 323'
      endif

      open(unit=999,file='snow_course_locations.xyz',
     *           status='unknown')
      do is=1,ng
        write(999,99001)ewg_fl(is),sng_fl(is),is,gname(is)
        write(*,99001)ewg_fl(is),sng_fl(is),is,gname(is)
99001   format(2f20.6,i5,5x,a12)
      end do
      
      
      close (unit=99,status='keep')
!     fix fix
!     for compatibility with old format
      do is=1,ng
        if(coordsys2.eq.'LATLONG   ')then 
          ewg(is)=int(ewg_fl(is)*60.0)
          sng(is)=int(sng_fl(is)*60.0)
        else
          ewg(is)=int(ewg_fl(is))
          sng(is)=int(sng_fl(is))
         endif
      end do
!     for compatibility with old format
      do is=1,ng
        xsta(is)=(ewg_fl(is)-xorigin)/xdelta+1
        ysta(is)=(sng_fl(is)-yorigin)/ydelta+1
        write(51,2102)is,xsta(is),ysta(is)
	  write(51,*)
      end do

! TS - ALLOCATION STATEMENT TO DYNAMICALLY ALLOCATE ARRAY FOR AREA16A
      allocate(p(ycount,xcount),w(ycount,xcount,ng),
     *           stat=iAllocateStatus)     !sn(ycount,xcount,nh),
      if (iAllocateStatus .ne. 0) STOP
     *    '**Allocation failed for area16 in snwa l140**'

	if(cnv.le.0.0)cnv=1.0
	
! OPEN A SCRATCH FILE TO CHECK COMPUTED WEIGHTS:
      open(unit=65,file='weight.chk',status='unknown')

!     Write the header

      author='snw.exe                               '
      coordsys_temp=coordsys1
      zone_temp=zone1
	datum_temp=datum1
	xorigin_temp=xorigin
	yorigin_temp=yorigin
	xcount_temp=xcount
	ycount_temp=ycount
	xdelta_temp=xdelta
	ydelta_temp=ydelta
	name='Snow Water Equivalent                   '
	attribute_units='mm                                      ' 
      attribute_type='Course                                  '
      unit_conversion=conv
      init_heat_deficit=deffactor  
      source_file_name=fln(35)                                    


      open(unit=99,file='junk1',status='unknown')
      write(99,99005)year1,mo1,day1,hour1
99005 format(i4,'/',i2,'/',i2,2x,i2,':00:00')
      rewind 99
	read(99,99006)startdate,starttime
99006 format(2a10)
      close(unit=99,status='delete')
      
      
      if(use_buf)then
          fln(36)=out_fln
      endif

      
c      radinfl=55

!     write the header
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c      call write_r2c(266,36,1,classcount,0,ii,1)   
      call write_r2c(266,36,0,classcount,0,0,1)   
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


      do ii=1,classcount
         write(51,*)
         write(*,*)
         write(51,1110)ii
         write(*,1110)ii

!        DISTRIBUTE SNOW ACCORDING TO WEIGHTS AND REWRITE YYMMDD.SNW FILE
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
         call weight(ii,ng,xcount,ycount,0.0)
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
         do is=1,ng
	      write(51,*)ewg_fl(is),sng_fl(is),gname(is)
	      do i=ycount,1,-1
	         write(51,9001)(w(i,j,is),j=1,xcount)
            end do
	   end do

         do i=1,ycount
            do j=1,xcount
               p(i,j)=0.0
            end do
	   end do

!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
         call distr(ii,xcount,ycount,ng)
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

!       convert local variables to write module variables
        do i=1,ycount
          do j=1,xcount
            outarray(i,j)=p(i,j)
          end do
        end do          

c      SUBROUTINE write_r2c(conv,un,fn,
c     *            no_frames,no_classes,frame_no,class_no,
c     *            no_signf)
!       write the data      

!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call write_r2c(266,36,0,classcount,0,ii,1)   
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      end do


! FORMATS:

 1000 format(' snow course distribution application',/)
 1002	format(12f10.1)
 1003 format(' water equivalent for station no ',i3,' is')
 1005	format(999f5.0)
 1110 format('+',' distributing the snow course data for class',i3)
 1111 format(' class ',i2)
 1203 format(' ycount=',i5,'  xcount=',i5)
 2000 format(5i5,2f5.1,30x,f5.2)
 2003 format(2i5,2f5.2)
 2004 format('  rgrd rads radn radw rade grdn grde')
 2006	format(' snow course inputs')
 2007	format(' ng ='i5,' classcount ='i5/)
 2100 format(2i5,1x,a12,1x,17f8.2)
 2101 format(' ',2i5,1x,a12/)
 2102 format(' ','sta(',i2,') xsta,ysta :',2f7.1/17f8.1)
 2110 format('  xcount=',i5,' ycount=',i5,' grdn=',f5.1,' grde=',f5.1)
 2111 format(/' class ',i2)
 9001	format(999f5.2)

 3001 format(a20,i12)
 3002 format(2a20)
 3003 format(a20,f12.3)
 3004 format(a20,a10)
 3005 format(999a10)
 3006 format(999f10.0)
 3007 format(999f10.4)
      
      print*
      print*,'DISCLAIMER'
	print*,'The WATFLOOD software and other material supplied' 
	print*,'in connection herewith is furnished by N. Kouwen and the' 
	print*,'University of Waterloo and is accepted by the' 
	print*,'user upon the express understanding that N. Kouwen' 
	print*,'or the University of Waterloo make no warranties, either' 
	print*,'express or implied, concerning the accuracy, completeness,' 
	print*,'reliability, usability, performance, or fitness for any' 
	print*,'particular purpose.' 
	print*
	print*,'The material is provided "as is". The entire risk as to' 
	print*,'its quality and performance is with the user.'
      print*
      print*
      print*,'NOTE:'
      print*,'Open ', fln(36)(1:60)
      print*,'in Green Kenue and check that all grids in the waterhsed'
      print*,'have values. If areas are not covered, increase the '
      print*,'radius of influence.'
      print*
      
      stop 'Normal ending'

c      if(stopflg.eq.'y')then
c        print*,' normal ending'
c        pause ' Hit enter to exit window'
c        print*
c        stop 
c      else
c        stop ' normal ending'
c      endif

      END PROGRAM snw
