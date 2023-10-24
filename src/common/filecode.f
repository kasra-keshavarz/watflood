!***********************************************************************
      SUBROUTINE find_filetype(fn)
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
      

!     written Dec. 13/05 by NK
!     This s/r will find the file type of fln(n)

      USE area_watflood
      implicit none
      save

c      use area12

      character(1)  :: dot
      integer       :: fn,i,ii,ios

!     find file type:
!     first find the dot in the string:
      dot=' '
      i=0
      do while(dot.ne.'.'.and.i.lt.225)
        i=i+1
!     rev. 9.8.11  Dec.  06/11  - NK: removed 30 char limit on find filetype 
c        if(i.gt.30)then
c!         old format files with soil moisture in the met file.
c          filetype='gsm'
c          return
c        endif
        dot=fln(fn)(i:i)
        chars(i)=dot
c        write(*,*),dot
      end do
c      write(51,*),'dot is char no ',i
c      write(51,*),(chars(ii),ii=1,i)
      filetype=fln(fn)(i+1:i+3)
      filename_length=i+3            

c      write(51,*),'filetype=',filetype
c      write(51,*),'filename_length=',filename_length
c      write(51,*)
c      write(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c      write(*,*)'file no=',fn
c      write(*,*)'file name=',fln(fn)(1:60)	 
c      write(*,*)'filetype=',filetype
c      write(*,*)'filename_length=',filename_length
c      write(*,*)

      return

      end SUBROUTINE find_filetype


!***********************************************************************
      subroutine create_filename(fn,new_fn,ext)
!***********************************************************************

!     written Dec. 13/05 by NK
!     This s/r will change the file name of fln(n) to an ensim
!     compatible format

c      use area12

      USE area_watflood
      implicit none

      character(30) :: old_filename,dummy
	character(3)  :: ext
      character(1)  :: dot
      integer       :: fn,i,ii,ios,new_fn

!     find file type:
!     first find the dot in the string:

      old_filename=fln(fn)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call find_filetype(fn)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      write(98,*),'filetype=',filetype
      write(98,*),'filename_length=',filename_length
      write(98,*)
      
      open(unit=99,file='abcxyz',status='unknown',iostat=ios)
!      write(*,*)(chars(i),i=1,filename_length-4),'_',filetype,'.',ext
      write(99,*)(chars(i),i=1,filename_length-4),'_',filetype,'.',ext
      rewind 99
      read(99,99000)fln(new_fn)   ! importan to have the format here
99000 format(a30)
      close(unit=99,status='delete')

      print*
      print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      print*,'IMPORTANT NOTE:'
      print*,'A new filename       ',fln(new_fn) 
      print*,'has been created from ',old_filename
      print*,'in accordance with the new ENSIM compatible file formats'
      print*
      print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      print*


      return

      end subroutine create_filename


!     this is old code - changed to create_filename

!*********************************************************************

c      SUBROUTINE change_filename(fn)
c
c!     written Dec. 13/05 by NK
c!     This s/r will change the file name of fln(n) to an ensim
c!     compatible format
c
c      use area12

c      character(30) :: old_filename
c      character(1)  :: dot
c      integer       :: fn,i,ii,ios
c
c!     find file type:
c!     first find the dot in the string:
c
c      old_filename=fln(fn)
c!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c      call find_filetype(fn)
c!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c      write(98,*),'filetype=',filetype
c      write(98,*),'filename_length=',filename_length
c      write(98,*)
      
c      open(unit=99,file='abcxyz',status='unknown',iostat=ios)
c      write(*,*)(chars(i),i=1,filename_length-4),'_',filetype,'.r2c'
c      write(99,*)(chars(i),i=1,filename_length-4),'_',filetype,'.r2c'
c      rewind 99
c      read(99,99000)fln(fn)   ! importan to have the format here
c99000 format(a30)            ! otherwise name gets cutoff at the /
c      close(unit=99,status='delete')
c
c      print*
c      print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c      print*,'IMPORTANT NOTE:'
c      print*,'The file name        ',old_filename
c      print*,'has been changed to ',fln(fn) 
c      print*,'in accordance with the new ENSIM compatible file formats'
c      print*
c      print*,'You are asked to change the name in the event file'
c      print*,'If you do not change this name in the event file, '
c      print*,'the newly created file can not be used. '
c      print*
c      print*,'The new names are recognized by ENSIM '
c      print*,'allowing the files to be viewed in ENSIM '
c      print*,'as time series animations etc.'
c      print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c      print*
c
c      return
c
c      end subroutine change_filename


!*********************************************************************
