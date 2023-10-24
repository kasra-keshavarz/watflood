!     ********************************************************
!     *                                                      *
!     *                  WATFLOOD (TM)                       *
!     *                                                      *
!     *     Program BSN Version 10.10     May 11, 2017       *
!     *                                                      *
!     *           (c) N. Kouwen, 1972-2017                   *
!     *                                                      *
!     ********************************************************

      PROGRAM bsn

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

!**********************************************************************
! BSNA - 
!
! Modified by FS - July 1999 variables for lan/long coordinates
! Modified by NK - June 8/00 for free format
!
! Modified by Tricia Stadnyk - September 2000
! Converted common blocks to modules and added dynamically allocated
! run-time arrays as part of Fortran 90 conversion.
!
! Modified June 09, 2011 for non-contributing areas NK
!
!     subroutines called:
!               arrange
!               grade
!
!  VERSION 9      Dec. 2000  - TS: Added dynamic memory allocation
!  Version 9.07 - Jun. 28/02 - 
!  Version 9.08 - Nov. 22/02 - added profiles for sub-basins
!  Version 9.46 - Jul. 31/03 - fixed bug in converting areas
!  rev. 9.1.50  Jan.  14/04  - NK: version number added to the wfo_spec.txt file
!  March 09/09 - nk - fixed river profiles and made them xyz files
!  Version 9.80 - Dec. 17/09 - set stuff outside watershed = 0
!  Version 10.1
!  Version 10.2 - Apr. 2011 - incorporated frames with 0 land cover fraction
!  Version 10.3 - Jun. 2011 - modify frac for nca
!  Version 10.4 - NOv. 2011 - new class_combine.csv file
!  Version 10.5 - Jan. 2012 - create elv_means.r2s
!  Version 10.6 - Mar. 2012 - create split landcover for nca
!  Version 10.7 - Jul. 2012 - added inlet grids
!  Version 10.8 - Sep. 2012 - added elevation datum to get rid of -ve elevations
!  Version 10.9 - Nov. 2015 - 'al' calculation revised to centre grid from origin 
!  Version 11.0 - Jul. 2016 - Reordering of grids for // processing 
!  rev. 10.4.21 Apr.  21/20  = NK Add UZS deficit to wfo file = UZS(class=classcount)
!
! - list of arguments:
!
!   latsdeg = lattitude of south boundary in degrees
!   latsmin = lattitude of south boundary in minutes
!   latndeg = lattitude of north boundary in degrees
!   latnmin = lattitude of north boundary in minutes
!   lonwdeg = longitude of west boundary in degrees
!   lonwmin = longitude of west boundary in minutes
!   lonedeg = longitude of east boundary in degrees
!   lonemin = longitude of east boundary in minutes
!   yint    = north/south height of gru in minutes
!   xint    = west/east width of gru in minutes
!   al      = grid length in meters 
!   cintvl  = contour interval in meters 
!   impr    = percent of urban area that is impervious 
!           if impervious area is classified as such, imp will be 100% 
!           if urban area is classified as such, imp will 
!           range from about 25% to 50% depending on popn. density etc. 
!           default for inpr = 100% 
!           - this class is impervious area, not urban
!   itype   = channel type. Entered in parameter file.del fort.
!           itype = 0 then channel has flood plains
!           itype = 1 flood plains are infinite width
!   ntype   = number of land cover classes including 
!             urban/impervious   (changed Nov. 6/06  nk)
!   elvconv = conversion factor from imp. units to si
!
!************************************************************************

      USE area17   ! see note

      use area_watflood

      implicit none

!     note:  some variable in area1 and area17 have the same name 
!     but have different dimensinos. 
!     So here they are called   rank_2d instead of say rank

      CHARACTER(10) :: time
      CHARACTER(8)  :: cday

      CHARACTER(80) :: comment,line
c      CHARACTER(79) :: hdrcomment(100)
      CHARACTER(50) :: junk,filenames(50)
      CHARACTER(256) :: fn1
      CHARACTER(10) :: coordsys,zone                   !datum1,
      INTEGER(2)    :: nrad,ssquery,n_hdr_lines
      INTEGER       :: chksum(999),iverflg,ireach_2d_max
      integer       :: new_class_count,old_class,new_class,
     *                 class_combine_version
      INTEGER       :: latsdeg,latsmin,latndeg,latnmin,lonwdeg,
     *                 lonwmin,lonedeg,lonemin,lastn,
     *                 longest,maxlength,tempsum,startn,m,k,
     *                 nextgrid,nogrids,nullgrids,frame   !,nnn
      REAL          :: yint,xint,latrad,fudge,slopemin,split,
     *                 sngrid,ewgrid,cintv,aimpr,antype,elvconv,perv,
     *                 sumclass,e1,aimp,minelv,tmp_value,elv_datum,
     *                 xx_temp,yy_temp,      
     *                 delta_x,delta_y,
     *                 connected_nca,nca_use,threshold, 
     *                 imperv(99),pervious(99) 
      integer       :: no_errors,error_count,impr,newformat,iiiii,
     *                 newversion,i,j,n,ii,jj,l,lg,ix,iallocate,         !,local,nn
     *                 iallocatestatus,ios,sss,iDeallocateStatus,
     *                 ycount_new,xcount_new,version_no,no_outlets,
     *                 npoints,ncount,i_datum,j_datum,
     *                 nca_classcount,nca_choice

      INTEGER(2)   :: status1
      CHARACTER(1) :: buf,update_flg,junk1,reply,parflg,respflg,
     *                answer,answer_dem,adjust_frac,
     *                split_type
      integer(2)   :: nrows,srows,wcols,ecols,
     *                i1flg,i2flg,i3flg,j1flg,j2flg
      logical      :: exists,inletflg,slopeflg,errflg_l
      
      integer      :: no_inlets,inletgrid(99)
      

      integer*4, dimension(:),    allocatable:: lastgrid
      integer*4, dimension(:,:),  allocatable:: elv_count
      integer*4, dimension(:,:),  allocatable:: array
      real*4,    dimension(:),    allocatable:: chainage
      real*4,    dimension(:,:),  allocatable:: elv_dem
      real*4,    dimension(:,:),  allocatable:: elv_sum,elv_max,elv_min
      real*4,    dimension(:,:,:),allocatable:: aclass_comb
      real*4,    dimension(:,:,:),allocatable:: aclass_read
      CHARACTER(20), dimension(:),allocatable:: class_name
      logical,   dimension(:),    allocatable:: keepflg

      program_name='bsn       '

!     USE DFLIB
      allocate(fln(601),stat=iAllocateStatus)
      if (iAllocateStatus.ne.0) STOP 
     *    '**Allocation failed for area12**' 

      CALL GETARG(1,buf,status1)
      if(status1.ne.1)iopt=status1
      if(buf.ne.'n')buf='y' ! default - use 'bsn64x n' for no frame headers
      
      no_errors=0
c      iopt=2
      nrvr=0
      coordsys='unknown   '
      datum1='unknown   '
      zone='unknown   '

! TS - INITIALIZATIONS MOVED TO LINE308 DUE TO ALLOCATIONS OF ARRAYS

      print*,'********************************************************'
      print*,'*                                                      *'
      print*,'*                  WATFLOOD (TM)                       *'
      print*,'*                                                      *'
      print*,'*     Program BSN Version 10.9     Nov. 23, 2015       *'
      print*,'*               Revised  Aug. 29, 2023                 *'
      print*,'*                                                      *'
      print*,'*           (c) N. Kouwen, 1972-2023                   *'
      print*,'*                                                      *'
      print*,'********************************************************'
      print*
      print*,'Please see file bsn_info.txt for information re: this run'
      print*
      print*,'VERY IMPORTANT CHANGES:'
      print*
      print*,'In the bsnm.map file'
      print*,'the impervious area is now the LAST class - not the first'
      Print*,'The no of classes is now the TOTAL number - including the'
      print*,'impervious class'
      print*
      print*,'Please change the .map file accordingly if you have not'
      print*,'yet done so.  Sorry for the inconvenience  NK'
      print*
      print*,'Program = modified to allow for non-contributing areas'
      print*,'The file nca.r2s is required - see WATFLOOD manual'
      print*
      print*,'Program = modified to allow for mean & max grid elev`s'
      print*,'The file dem.r2s is required - see WATFLOOD manual'
      print*
      print*,'Program = modified to insert frame headers as comments'
      print*,'To use olw way with no frame headers, '
      print*,'Execute the program with argument ie. `bsn64x n`  '
      print*,'                                  OR  `bsn64d n`  '
      print*
d      print*,'NOTE:'
d      print*,'You are running the debug version of bsn'
d      print*,'Frame headers will be written in the new_shd.r2c file'
d      print*,'These can not be read by GreenKenue so please use'
d      print*,'BSN64x.exe or BSN32x.exe for the working copy.'
d      print*
      pause 'Hit enter to continue - Ctrl C to abort'
      print*

      allocate(lastgrid(100),stat=iAllocateStatus)
      if (iAllocateStatus.ne.0) STOP 
     *    '**Allocation failed for lastgrid**' 
     
     
      INQUIRE(FILE='class_combine.csv',EXIST=exists)
      IF(.not.exists)THEN
        print*,'class_combine.csv not found'
        print*,'Would you like to continue without this file?'
        print*,'If yes, classes will be transfereed to the shd file'
        print*,'in the same order as the map file. Not a good idea'
        print*,'unless you have reordered the classes in the map file'
        print*,' with glaciers (if present) wetland, water & impervious'
        print*,'as the last classes.'
        print*
        print*,'y/n ?'
        read*,line
        if(line(1:1).eq.'n')stop 'Program aborted in bsn @ 211'
      endif
      
      INQUIRE(FILE='bsn_responses.txt',EXIST=exists)
      IF(exists)THEN
        print*,'Found bsn_responses.txt'
        open(unit=99 ,file='bsn_responses.txt',status='old',iostat=ios)
        read(99,*,iostat=ios)line
        close(unit=99)
c        print*,'line=',line(1:20)
      else
        line='junk'
        print*,'bsn_responses.txt  NOT found'

      endif
     
      
c      if(exists.and.line(1:9).eq.'version_#')then
      if(line(1:9).eq.'version_#')then
        open(unit=99,file='bsn_responses.txt',status='old',iostat=ios)
        rewind 99
        read(99,*)line,version_no
        read(99,*)line,fn1   
99999   format((a),$)        
        INQUIRE(FILE=fn1,EXIST=exists)
d        print*,line(1:25),fn1(1:40)       
        if(.not.exists)then
          print*,'file:',fn1(1:60) 
          print*,'not found - please check name'
          stop 'Program aborted @ 212'
        endif
        read(99,*)line,fln(2)
        if(fln(2)(1:3).ne.'na ')then
          INQUIRE(FILE=fln(2),EXIST=exists)
d         print*,line(1:25),fln(2)(1:40)       
          if(.not.exists)then
            print*,'file not found - please check name'
            stop 'Program aborted @ 220'
          endif
        endif
        read(99,*)line,author          
        read(99,*)line,no_outlets,(lastgrid(lg),lg=1,no_outlets)
        if(no_outlets.eq.0)lastgrid(1)=0 ! must be initialized          
        read(99,*)line,no_inlets,(inletgrid(lg),lg=1,no_inlets)
        read(99,*)line,split          
        read(99,*)line,split_type  
        read(99,*)line,slopemin          
        read(99,*)line,adjust_frac
        read(99,*)line,nca_choice          
        read(99,*)line,nca_use    
        read(99,*)line,nca_classcount          
        read(99,*)line,answer_dem  
        close(unit=99,status='keep')
        
        print*
        write(*,99100)'version_#                     ',4
        write(*,99101)'map_file_name                 ',fn1
        write(*,99101)'par_file_name                 ',fln(2)
        write(*,99101)'initial                       ',author
        write(*,99102)'no_outlets_&_locations        ',no_outlets,
     *                              (lastgrid(l),l=1,no_outlets)
        write(*,99103)'no_inlets_&_locations         ',no_inlets,
     *                              (inletgrid(l),l=1,no_inlets)
        write(*,99105)'wetland_split_%               ',split
        write(*,99106)'split_type_1|2                ',split_type
        write(*,99105)'min_allowed_slope             ',slopemin
        write(*,99106)'adjust_frac_y|n               ',adjust_frac
        write(*,99100)'nca_choice_1|2                ',nca_choice         
        write(*,99105)'%_to_use(choice_1)            ',nca_use
        write(*,99100)'nca_classes(1-3)(choice_2)    ',nca_classcount
        write(*,99106)'create_max|mean.r2c_y~n       ',answer_dem
        print*
        print*,'If these values are ok'
        print*,'please hit enter to continue'
        print*,'If not, hit ctrl_C & edit or delete bsn_responses.txt'
        print*,'and relaunch this program.'
        print*
        pause 'Waiting .............'
!       pause pause pause pause pause pause pause pause pause pause pause           

      
      else
        if(exists)then      
          print*,'Old format bsn_responses.txt file found'
        endif
        print*,'Please create a new file by answering the following:'
c        close(unit=99)
        open(unit=99 ,file='bsn_responses.txt',status='unknown',
     *                               iostat=ios)
        if(ios.ne.0)then
          print*,'Problem opening bsn_responses.txt'
          print*
          stop 'Program aborted in bsn @ 247'
        else
          print*,'Opened bsn_responses.txt'
          print*
        endif
!       ..................................................        
        print*,'Enter the basin (map) file name:'
        read(*,*)fn1
!       ..................................................        
        print*,'Enter the parameter (par) file name ONLY'
        print*,'if you need a bsnm_par.r2c file for WATROUTE'
        print*,'other wise, enter:    na'
        read(*,*)fln(2)
!       ..................................................        
        print*
        print*,'Enter your name or initials'
        read(*,*)author
!       ..................................................
        print*,'Once you have a shd file for the whole domain'
        print*,'you can extract sub-watersheds to run on their own'
        print*,'I.e. you can remove downstream grids from the modelling'
        print*,'Load the shd file for the whole domain'
        print*,'into GreenKenue and note the rank(s)'
        print*,'of the location(s) where you would like an outlet -'
        print*,'normally at streamflow locations but not necessarily so'        
        print*,'Enter the number of sub-watersheds'
        print*,'(to NOT remove downstream gridsenter 0)'
        read(*,*)no_outlets
!       ..................................................        
!       grid size reduction:     9/05/2002   NK
        print*
        print*,'Enter the outlet grid rank(s) you would '
        print*,'like included in the simulation'
        print*,'These should NOT be the receiving grids!!!!'
        print*
        print*,'Please enter the rank of',no_outlets,' outlet(s):'
        print*
        print*,'example:  1482'
        print*,'example:  1043'
        print*,'example:  1899'
        print*
        if(no_outlets.gt.0)then   
          do l=1,no_outlets
            print*,'enter outlet grid #',l
            read(*,5)lastgrid(l)
5           format(i5)
          end do
        endif
!       ..................................................   
        print*,'Similarly,'     
        print*,'Upstream watersheds can be removed from the modelling'
        print*,'Enter the number of inlet grids'
        print*,'To use the all upstream watershes enter 0'
        print*,'OR enter the no of grids where upstream area is NOT'
        print*,'to be modelled:'
        read(*,*)no_inlets
!       ..................................................        
!       grid size reduction:     9/05/2002   NK
        print*
        print*,'Enter the inlet grid rank(s) you would '
        print*,'like to use for the simulation'
        print*,'These would normally be streamgage locations'
        Print*,'where you could add inflow to be routed downstream'
        print*,'using either the str or div files'
        print*
        print*,'Please enter the rank of',no_inlets,' inlet(s) '
        print*
        print*,' example:   482'
        print*,' example:    43'
        print*,' example:    99'
        print*
        if(no_inlets.gt.0)then   
          do l=1,no_inlets
            print*,'enter inlet grid #',l
            read(*,5)inletgrid(l)
          end do
        endif
!       ..................................................        
        print*
        print*,'Enter the split: % of wetland coupled to channel'
        print*,'only if you have two identical sets of wetland '
        print*,'land cover grids as the 2 classes before the '
        print*,'water class in the land use section of the map file'
        print*,'Enter 0 if you have just 1 block of wetland cover'
        print*
        print*,'Split = ?'
        read*,split
!     ..................................................        
        if(split.gt.0.0)then
          print*
          print*,'NEW  June 22/11'
          print*,'Wetland split options:'
          print*,'Option 1: Old way - split % of wetland - 1% minimum'
          print*,'Option 2: New way - split % =  max % wetland coupled'
          print*,'                  - wetland >  split % not coupled'
          print*,'                  - wetland <= split % all coupled'
          print*
          print*,'Enter 1 or 2'
          read*,split_type
          print*,'You have entered ',split_type
          if(split_type.ne.'1'.and.split_type.ne.'2')then
            print*,'invaled answer'
            print*,'default = 1 assumed'
            split_type='1'
            print*,'Please edit bsn_response.txt to correct this'
            pause 'Hit enter to continue with this entry'
          endif
      else
!         need value in the file        
!         Added Jul. 15/18 NK          
          split_type='1'  
        endif
!       ..................................................        
        print*
        print*,' Often DEM have flat spots filled and you end up with'
        print*,' unwanted flat spots in your river profile'
        print*,' It causes severe flattening of the hydrographs'
        print*,' Enter the minimum allowable river slope'
        print*,' that you have in your sustem - e.g. 0.0001'
        print*,' Min accepted value = 0.0000001'
        print*,' Max value accepted is 1.0 (45 degrees!)'
        print*
!       This value is needed in grade.for
        read*,slopemin
!       ..................................................        
        print*,'Do you want to incorporate'
        print*,'non-contributing areas (nca)   y/n?'
        print*,'To incorporate nca`s an nca.r2s file is required'
        print*
        read(*,*)adjust_frac
        fln(200)='nca.r2s'
        INQUIRE(FILE=fln(200),EXIST=exists)
        if(exists.and.adjust_frac.eq.'y')then
          effective_nca=1.00
          print*
          print*,'nca.r2s file found'
          print*,'You have 2 choices for using the non contr. area'
          print*,'1. simply reduce the value of frac for each cell'
          print*,'   by it`s % non-contributing'
          print*,'   No extra classes.'
          print*,'   nca(s) is (are) eliminated from the modelling'
          print*,'OR'
          print*,'2. divide up to 3 land cover classes into contr.'
          print*,'   and non-contr. using the % non-contr.'
          print*,'   These will be additional classes which need to be'
          print*,'   accounted for in the par file. Thus these nca`s'
          print*,'   can be made to contribute during very wet years'
          print*
          print*,'Continue with adjustment of frac?'
          print*,'1 = choice # 1'
          print*,'2 = choice # 2'
          read*,nca_choice
!  Version 10.6 - Mar. 2012 - create split landcover for nca
!  fixed Dec. 04/12 NK
c            split_no=0
            nca_classcount=0
            if(nca_choice.eq.1)then
              print*,'What % of the nca would you like to use?'
              print*,'acceptable range 0%-100%'
              read*,nca_use
              if(nca_use.gt.100)stop 'Program aborted @49'
              if(nca_use.lt.0)stop 'Program aborted @ 50'
              effective_nca=nca_use/100.0
	        non_ca_flg=.true.
            elseif(nca_choice.eq.2)then	    
!             cannot reduce frac & split classes      
!             alternative way of splitting contributing and non contributing
!             areas: splitting up to 3 dominant classes
              print*,'Would you like to split any classes into '
              print*,'contributing and non-contributing?'
              print*,'You can only split the first `n` classes'
              print*,'in the shd file (not the map file)'
              print*,'e.g. if crops & grass are the first 2'
              print*,'you can split these by answering 2'
              print*,'If you want to split only the first one'
              print*,'enter 1  -  for no split, enter 0'
              print*,'3 is the maximum'
              print*,'How many?'
c              read*,split_no
              read*,nca_classcount
c              if(split_no.gt.3)then
              if(nca_classcount.gt.3)then
                print*
                print*,'WHOAH!!!!'
                print*,'3 is the maximum allowed'
                stop 'Program aborted at 920'
              endif
              print*,'You have elected to split',nca_classcount,
     *                 ' land cover classes'
              print*
            endif
c          endif
        else   ! nca.r2s file not found so:  
          if(adjust_frac.eq.'y')then
!           nca.r2s file not found            
            print*
            print*,'WARNING:'
            print*,'nca.2rs file not found'
            print*
          endif
          nca_choice=0
c          split_no=0
          nca_classcount=0
          nca_use=0
          print*,'non-contributing areas will not be incorporated'
          print*,'in the shed file'
          print*
        endif   ! does fln(200)    nca.r2s   exist?  
!       ..................................................        
        print*,'Do you want to create  new'
        print*,'elev_means.r2c & elev_max.r2c files   y/n?'
        print*,'To create these files, a dem.2rs file is required'
        print*,'which can be created in GK by saving the dem'
        print*,'as an r2s file' 
        print*
        read(*,*)answer_dem
        INQUIRE(FILE='dem.r2s',EXIST=exists)
        IF(exists)THEN
!         read the header
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          print*
          print*,'A DEM file has been found'
          print*,'so mean grid elevations will be calculated'
          print*,'and written to a files: elev_means.r2c'
          print*,'These mean elevations will be used for lapse rate'
          print*,'adjustments to temperature & precip if left in'
          print*,'the basin directory'
          print*
 !         ..................................................       
        endif
        print*
        print*,'A new (Ver. 4) bsn_resonse.txt file will be created'
        print*,'Any old file will be overwritted with:'
        print*
        write(*,99100)'version_#                     ',4
        write(*,99101)'map_file_name                 ',fn1
        write(*,99101)'par_file_name                 ',fln(2)
        write(*,99101)'initial                       ',author
        write(*,99102)'no_outlets_&_locations        ',no_outlets,
     *                              (lastgrid(l),l=1,no_outlets)
        write(*,99103)'no_inlets_&_locations         ',no_inlets,
     *                              (inletgrid(l),l=1,no_inlets)
        write(*,99105)'wetland_split_%               ',split
        write(*,99106)'split_type_1~2                ',split_type
        write(*,99105)'min_allowed_slope             ',slopemin
        write(*,99106)'adjust_frac_y|n               ',adjust_frac
        write(*,99100)'nca_choice_1|2                ',nca_choice         
        write(*,99105)'%_to_use(choice_1)            ',nca_use
        write(*,99100)'nca_classes(1-3)(choice_2)    ',nca_classcount
        write(*,99106)'create_max|mean.r2c_y~n       ',answer_dem
        print*
        
        INQUIRE(FILE='bsn_responses.txt',EXIST=exists)
        if(exists)then
          print*
          print*,'If you wish to keep the existing bsn_responses.txt,'
          print*,'file, please move it now'
          pause ' Waiting  .......  hit return to continue'
!         pause pause pause pause pause pause pause pause pause pause pause           
        endif  
        open(unit=99,file='bsn_responses.txt',status='unknown',
     *       iostat=ios)
        if(ios.ne.0)then
          print*,'Problems opening bsn_responses.txt'
          print*,'Please ensure the file is not already open'
          pause 'hit enter to continue - look for fort.99'
        endif   
        write(99,99100)'version_#                     ',4
        write(99,99101)'map_file_name                 ',fn1
        write(99,99101)'par_file_name                 ',fln(2)
        write(99,99101)'initial                       ',author
        write(99,99102)'no_outlets_&_locations        ',no_outlets,
     *                              (lastgrid(l),l=1,no_outlets)
        write(99,99103)'no_inlets_&_locations         ',no_inlets,
     *                              (inletgrid(l),l=1,no_inlets)
        write(99,99105)'wetland_split_%               ',split
        write(99,99106)'split_type_1~2                ',split_type
        write(99,99105)'min_allowed_slope             ',slopemin
        write(99,99106)'adjust_frac_y|n               ',adjust_frac
        write(99,99100)'nca_choice_1|2                ',nca_choice         
        write(99,99105)'%_to_use(choice_1)            ',nca_use
        write(99,99100)'nca_classes(1-3)(choice_2)    ',nca_classcount
        write(99,99106)'create_max|mean.r2c_y~n       ',answer_dem
        close(unit=99,status='keep')
        print*
        print*,'New bsn_response.txt file written'
c        endif

      endif
c      endif        

!     open the map file:
      INQUIRE(FILE=fn1,EXIST=exists)
      IF(exists)THEN
        open(31,file=fn1,status='old',iostat=ios)
        if(ios.ne.0)then
           print*,' Unable to open the file ',fn1
           print*,' Are you in the proper directory?'
           print*,' Did you enter the proper map file name?'
           print*,' Does the file you named exist?'
           print*
           STOP ' Program aborted'
        endif
      else
        print*,'file ',fn1(1:40),' not found'
        print*
        stop 'Program aborted in bsn.exe @ 324'
      endif

c      if(respflg.eq.'n')then
c        print*,'Enter the parameter (par) file name'
c        print*,'if you want a bsnm_par.r2c file for watroute'
c        print*,'other wise, hit return'
c        read(*,1)fln(2)
c        write(99,1)fln(2)
c    1   format(a256)
c      endif

      if(fln(2)(1:1).eq.'n')then
!       no par.r2c file wanted for watroute
        print*,'bsnm_par.r2c not wanted for watroute'
        parflg='n'
      else
!       basin/bsnm.par
        open(unit=32 ,file=fln(2),status='old',iostat=ios)
        if(ios.ne.0)then
          print*,'Problems on unit 32'
          write(*,99172)fln(2)
          write(51,99172)fln(2)
99172     format(' Warning: Error opening or reading fln:',a30/
     *    ' Probable cause: missing basin/bsnm.par input file')
          print*,'iostat code =',ios
          STOP 'program aborted in bsn.for @ 139'
        endif
        parflg='y'
      endif


!     get rid of the existing new.map file  nk. Jul. 29/04
      INQUIRE(FILE='new.map',EXIST=exists)
      IF(exists)THEN
        open(unit=999,file='new.map',status='old',iostat=ios)
        close(unit=999,status='delete')
        print*,'existing new.map file deleted'
        print*,'replaced by new_format.map '
        print*
        pause 'hit enter to continue'
      endif

      open(36,file='new_format.shd',form='formatted',
     *     access='sequential',recl=2048,status='unknown',
     *     iostat=ios)
      if(ios.ne.0)then
         print*,' Problems opening the NEW1.SHD file '
         print*,' Possible cause: file open or read only'
         print*
         STOP ' Program BSN aborted'
      endif

!     get rid of the existing new.shd file  nk. Jul. 8/04
      INQUIRE(FILE='new.shd',EXIST=exists)
      IF(exists)THEN
        open(unit=37,file='new.shd',status='old',iostat=ios)
        close(unit=37,status='delete')
        print*,'existing new.txt file deleted'
        print*,'replaced by new_format.map and old_format.map'
        print*
        pause 'hit enter to continue'
!       pause pause pause pause pause pause pause pause pause pause pause           
      endif

      open(37,file='old_format.shd',form='formatted',
     *     access='sequential',recl=2048,status='unknown',
     *     iostat=ios)
      if(ios.ne.0)then
         print*,' Problems opening the NEW.SHD file '
         print*,' Possible cause: file open or read only'
         print*
         STOP ' Program BSN aborted'
      endif

      open(38,file='ssdata.out',form='formatted',
     *     access='sequential',recl=2048,status='unknown',
     *     iostat=ios)
      if(ios.ne.0)then
         print*,' Problems opening the ssdata.out file '
         print*,' Possible cause: file open or read only'
         print*
         STOP ' Program BSN aborted'
      endif

      open(51,file='bsn_info.txt',form='formatted',
     *     access='sequential',recl=2048,status='unknown',
     *     iostat=ios)
      if(ios.ne.0)then
         print*,' Problems opening the bsn.err file '
         print*,' Possible cause: file open or read only'
         print*
         STOP ' Program BSN aborted'
      endif

!     Undocumented feature
!     Undocumented feature
!     Undocumented feature
!     Undocumented feature
!     Wetlands can be either coupled or not coupled to the river
!     To have both present, there should be two wetland classes. 
!     Only the last one will be coupled to the river.
!     Create a map file with two wetland classes, the last one = coupled
!     If they can not be differentiated from a land cover map, first
!     create two identical wetland % grids. 
!     Next divide them up according to the desired split, which is entered 
!     on the second line of the old format map file.
!     The dual wetland percentages are computed in the write statement
!     

!     Changed to accept any number of comments at the top of the file
!         or  -- none at all   NK  Sept. 11/05

      newversion=0
      i=0

      no_hdrcomments=1000
      allocate(hdrcomment(no_hdrcomments),
     * 		  hdrcomment_value(no_hdrcomments),stat=iAllocateStatus)
      if(iAllocateStatus.ne.0)STOP 
     *    '**Allocation failed for hdrcomments in ragmet @ 458**' 
      read(31,5007,iostat=ios)junk1,hdrcomment(i+1)
      if(junk1.eq.'#'.or.junk1.eq.':')then
        i=i+1
        print*
        print*, 'Ensim compatible free format map file expected'
        print*
        newversion=1
!       read the New ENSIM format
!       read the New ENSIM format
!       read the New ENSIM format
!       read the New ENSIM format
!       read the New ENSIM format
!          file can start with # or : for first data line
        do while(junk1.eq.'#')
          i=i+1
          read(31,5007,iostat=ios)junk1,hdrcomment(i)
          write(51,5007)junk1,hdrcomment(i)
          write(*,5007)junk1,hdrcomment(i)
          if(ios.ne.0)then
            print*,'error on line ',i,' of file',fn1
            print*,'Please update BSN.EXE'
            print*
            stop 'Program aborted in BSN @ 237'
          endif
        end do
        backspace 5
        n_hdr_lines=i-1

      endif

      rewind 31

      do while(line(1:10).ne.':endHeader')
        read(31,50001)line
50001   format(a80)
        if(line(1:10).eq.':EndHeader')line=':endHeader'
d	  print*,line(1:60)
        if(line(1:11).eq.':Projection')read(line,*)junk,coordsys
        if(line(1:9) .eq.':CoordSys')read(line,*)junk,coordsys
        if(line(1:10).eq.':Ellipsoid')read(line,*)junk,datum1
        if(line(1:6) .eq.':Datum')read(line,*)junk,datum1
        if(line(1:5) .eq.':Zone')read(line,*)junk,zone
        if(line(1:8) .eq.':xOrigin')read(line,*)junk,xorigin
        if(line(1:8) .eq.':yOrigin')read(line,*)junk,yorigin
        if(line(1:7) .eq.':xCount')read(line,*)junk,xcount
        if(line(1:7) .eq.':yCount')read(line,*)junk,ycount
        if(line(1:7) .eq.':xDelta')read(line,*)junk,xdelta
        if(line(1:7) .eq.':yDelta')read(line,*)junk,ydelta
        if(line(1:16).eq.':contourInterval')read(line,*)junk,cintv
        if(line(1:15).eq.':imperviousArea')read(line,*)junk,aimpr
        if(line(1:11).eq.':classCount')read(line,*)junk,antype
        if(line(1:15).eq.':elevConversion')read(line,*)junk,elvconv
      end do

!     GreenKenue uses LatLong - code below uses LATLONG
      if(coordsys.eq.'LatLong   ')coordsys='LATLONG   '
      if(coordsys.eq.'Cartesian'(1:9))coordsys='CARTESIAN '

      if(coordsys(1:3).eq.'UTM')then
        continue
      elseif(coordsys(1:9).eq.'CARTESIAN')then
        continue
      elseif(coordsys(1:7).eq.'LATLONG')then
        continue
      elseif(coordsys(1:7).eq.'LatLong')then
        continue
      else
        print*,' Valid format'
        print*,'(UTM, CARTESIAN, LatLong or LATLONG) not found'
        print*,'format found is ***',coordsys,'***'
        print*
        print*,'If wrong, please assign a Projection (coordsys)'
        print*,'in Green Kenue'
        print*
        stop 'Program aborted @ 407'
      endif

      print*,'projection=',coordsys
      print*,'datum1=',datum1
      print*,'zone=',zone
      print*,'xorigin=',xorigin
      print*,'yorigin=',yorigin
      print*,'xcount=',xcount
      print*,'ycount=',ycount
      print*,'xdelta=',xdelta
      print*,'ydelta=',ydelta
      print*,'cintv=',cintv
      print*,'aimpr=',aimpr
      print*,'antype=',antype
      print*,'elvconv=',elvconv

      ntype=int(antype)-1
      
c      if(respflg.eq.'n')then
c        print*
c        print*,'Enter the split: % of wetland coupled to channel'
c        print*,'only if you have two identical sets of wetland '
c        print*,'land cover gridsas the 2 classes before the '
c        print*,'water class in the land use section of the map file'
c        print*,'Enter 0 if you have just 1 block of wetland cover'
c        print*
c        print*,'Split ='
c        read*,split
c        write(99,*)split
c      endif

!     fixed Dec. 24/08
      split=split/100.0
      if(split.gt.1.0.or.split.lt.0.0)then
        print*,'value for split is outside acceptable range'
        print*
        stop 'Program aborted @ 483'
      endif

c      if(respflg.eq.'n')then
        print*,'Number of classes now includes the impervious class'
        print*,'Number of classes stipulated = ',ntype+1
        print*,'Is this correct?' 
        print*
        pause ' Hit enter to continue'
!       pause pause pause pause pause pause pause pause pause pause pause           
        
c        if(reply.ne.'y')then
c          print*,'Please fix the .map file'
c          print*
c          stop 'Program aborted in bsn @ 713'
c        endif

c      endif

c      pause 100
c      print*,'newversion=',newversion

      if(newversion.eq.0)then
!       write the new_version.map file - conversion from old format
!       write the new_version.map file - conversion from old format
        write(40,5012) 
        write(40,5009)coordsys
        if(coordsys.eq.'UTM       ')then
          write(40,5010)datum1
          write(40,5011)zone
        endif
c     pause 1
        if(coordsys.eq.'LATLONG   ')then
          write(40,5010)datum1
        endif
        write(40,5012)
        if(coordsys.eq.'LATLONG   ')then
          write(40,5013)xorigin
          write(40,5014)yorigin
        else
          write(40,50131)xorigin         !   *1000.0
          write(40,50141)yorigin         !   *1000.0
        endif
c     pause 2
        write(40,5012)
        write(40,5015)xcount
        write(40,5016)ycount
        if(coordsys.eq.'LATLONG   ')then
          write(40,5017)xdelta/60.0
          write(40,5018)ydelta/60.0
        else
          write(40,50171)xdelta
          write(40,50181)ydelta
        endif
c     pause 3
        write(40,5012)
        write(40,5019)cintv
        write(40,5020)aimpr

        write(40,5021)int(antype)+1  ! changed nov 6/06
!        write(40,5021)int(antype)

c      pause 4
        write(40,5022)elvconv   !units not converted for updated format
        write(40,5024)'          File converted from ',fn1
        write(40,5032)                !  :endHeader
      endif

      print*,'before allocating area17'

! TS - ALLOCATE AREA17A ARRAYS
      allocate(elv_2d(ycount,xcount),
     *da_2d(ycount,xcount),da_nca_2d(ycount,xcount),
     *slope_2d(ycount,xcount),rank_2d(ycount+1,xcount+1),
     *frac_2d(ycount,xcount),frac_nca_2d(ycount,xcount),
     *s_2d(ycount,xcount),next_2d(ycount,xcount),
     *last_2d(ycount,xcount),order_2d(ycount,xcount),
     *ielv(ycount,xcount),
     *iak(ycount,xcount),irough_2d(ycount,xcount),
     *ichnl_2d(ycount,xcount),bnkfll_2d(ycount,xcount),
     *ireach_2d(ycount,xcount),ireach_2d_temp(ycount,xcount),
     *yyy(ycount*xcount),xxx(ycount*xcount),ibn(ycount*xcount),
     *chainage(ycount*xcount),channel(ycount*xcount),
     *dummy(ycount,xcount),idummy(ycount,xcount),
     *new(ycount*xcount),ch_length_2d(ycount,xcount),
     *newxxx(xcount*ycount),newyyy(ycount*xcount),
     *newnext(ycount*xcount),
     *newrank(ycount*xcount),sl1(ycount*xcount),
     *       areamet(ycount,xcount),stat=iAllocateStatus)
      if (iAllocateStatus .ne. 0) STOP 
     *    '**Allocation failed for area17 in bsna l303**' 

      print*,'area17 allocated'

! INITIALIZATION - MOVED SEPT/2000 DUE TO ALLOCATIONS:
      do i=1,ycount
         do j=1,xcount
            rank_2d(i,j)=0
            next_2d(i,j)=0
            da_2d(i,j)=0.0
            da_nca_2d(i,j)=0.0
            bnkfll_2d(i,j)=0.0
            slope_2d(i,j)=0.0
            elv_2d(i,j)=0.0
            ch_length_2d(i,j)=0.0
            iak(i,j)=0.0
            ichnl_2d(i,j)=0.0
            ireach_2d(i,j)=0
            frac_2d(i,j)=0.0
         end do
      end do

      if(aimpr.lt.0.0)then
          print*,'WARNING'
          print*,'imperviousArea < 0 in the map file'
          print*,'This will be replaced by 100'
          print*,'I.e. impervious area in the shd file will be = '
          print*,'to impervious area in the map file'
          print*,'Please fix map file to desired value'
          aimpr=100
          pause 'Hit enter to continue - ctrl C to stop'
      endif
      if(elvconv.eq.0.0)elvconv=1.0

      ijk=ycount*xcount 

! ELEVATION(i,j) IS THE ELEVATION OF THE STREAMBED AT A POINT M
! TWEENTHE ENTRANCE AND THE EXIT OF THE STREAM AT A GIVEN STORAGE
! ELEMENT

      read(31,5000)comment
      if(newversion.eq.0)write(40,5000)comment
!      write(*,5000)comment
!      print*

      do i=ycount,1,-1
         if(newversion.eq.0)then
!          old formatted files
           read(31,1001,iostat=ios)(elv_2d(i,j),j=1,xcount)
           write(40,4001,iostat=ios)(elv_2d(i,j),j=1,xcount)
c           write(51,4001,iostat=ios)(elv_2d(i,j),j=1,xcount)
         else
!          new free format files
           read(31,*,iostat=ios)(elv_2d(i,j),j=1,xcount)
         endif
         if(ios.ne.0)then
           print*,' Problems reading the elevations on line',ycount-i+1
           print*
           STOP ' Program aborted'
         endif
      end do

! S(I,J) IS THE DIRECTION OF DRAINAGE OUT OF AN ELEMENT  
!        1=ne  3=se  4=s  5=sw  6=w  7=nw  8=n


!     CONVERT TO METRIC 

      elv_datum=1.0E+30
      do i=1,ycount
         do j=1,xcount
            elv_2d(i,j)=elv_2d(i,j)*elvconv
c            elv_datum=amin1(elv_2d(i,j),elv_datum)
            if(elv_2d(i,j).lt.elv_datum)then
              elv_datum=elv_2d(i,j)
              i_datum=i
              j_datum=j
            endif
         end do
      end do
      elvconv=1.0  !added 11/09/04  NK

!     Added Sept. 27/12  nk
!  Version 10.8 - Sep. 2012 - added elevation datum to get rid of -ve elevations
!     make all elevations +ve
      if(elv_datum.lt.0.0)then
        elv_datum=elv_datum-0.001
        do i=1,ycount
           do j=1,xcount
c            if(int(elv_2d(i,j)).ne.0)then
c             if(elv_2d(i,j).ne.0.000000.and.frac_2d(i,j).gt.0.0)then
             if(frac_2d(i,j).gt.0.0)then
               elv_2d(i,j)=elv_2d(i,j)-elv_datum
             endif  
           end do
        end do
      endif
c     for the receiving grid:      
      elv_2d(i_datum,j_datum)=elv_2d(i_datum,j_datum)-elv_datum


! FRAC IS THE FRACTION OF A SQUARE IN THE BASIN 
      read(31,5000)comment
      if(newversion.eq.0)write(40,5000)comment
!      write(*,5000)comment
      do i=ycount,1,-1
         if(newversion.eq.0)then
           read(31,1052,iostat=ios)(frac_2d(i,j),j=1,xcount)
           write(40,4050)(int(frac_2d(i,j)),j=1,xcount)
c           write(51,4050)(int(frac_2d(i,j)),j=1,xcount)
!          units not converted for updated format
         else
!          new free format files
           read(31,*,iostat=ios)(frac_2d(i,j),j=1,xcount)
c           write(51,4050)(int(frac_2d(i,j)),j=1,xcount)
         endif
         if(ios.ne.0)then
           print*,' Problems reading the fractions in row',i
           print*
           STOP ' Program aborted'
         endif
      end do

!      do i=ycount,1,-1
!         write(*,1056)(frac_2d(i,j),j=1,xcount)
!      end do
!      print*

!      pause '4 - hit enter to continue'

!      write(*,1109)

      write(51,1109)
      write(51,1057)j,j,(j,j=1,xcount)
      write(51,1057)j,j,(jxmin+istep*(j-1),j=1,xcount) 
      do i=ycount,1,-1
         write(51,1054)iymin+istep*(i-1),i,(frac_2d(i,j),j=1,xcount)
      end do

!      pause 41

      if(frac_2d(ycount,1).lt.0.0)then
         do i=ycount,1,-1
            do j=1,xcount
               frac_2d(i,j)=frac_2d(i,j)/(sstep*sstep)
            end do
         end do
      else
!        frac is given in percent so divide by 100
         do i=ycount,1,-1
            do j=1,xcount
               frac_2d(i,j)=frac_2d(i,j)/100.
            end do
         end do
      end if

!     DRAINAGE DIRECTIONS:

      read(31,5000)comment
      if(newversion.eq.0)write(40,5000)comment
!      write(*,5000)comment

      do i=ycount,1,-1
         if(newversion.eq.0)then
           read(31,1003,iostat=ios)(s_2d(i,j),j=1,xcount)
           write(40,1003,iostat=ios)(s_2d(i,j),j=1,xcount)
         else
!          new free format files
           read(31,*,iostat=ios)(s_2d(i,j),j=1,xcount)
         endif
         if(ios.ne.0)then
           print*,' Problems reading the directions on line',i
           print*
           STOP ' Program aborted'
         endif
      end do

!      write(*,1110)
!      write(*,1111)
!      do i=ycount,1,-1
!         write(*,1003) (s_2d(i,j),j=1,xcount)
!      end do
!      print*

!      pause '6 - Hit enter to continue' 

      write(51,1110)
      write(51,1111)
      write(51,1006)j,j,(j,j=1,xcount)
      write(51,1006)j,j,(jxmin+istep*(j-1),j=1,xcount) 
      do i=ycount,1,-1
         write(51,1006)iymin+istep*(i-1),i,(s_2d(i,j),j=1,xcount)
      end do

      ib=1+1
      it=ycount-1
      slopemin=amax1(slopemin,0.0000001)
      slopemin=amin1(slopemin,1.0)

! THE SUBROUTINE GRADEA CALCULATES THE CHANNEL SLOPE AT THE
! OUTLET OF EACH ELEMENT 

      if(coordsys.eq.'UTM       '.or.coordsys.eq.'CARTESIAN ')then
        areamet(1,1)=ydelta*xdelta
        do i=1,ycount
          do j=1,xcount
            areamet(i,j)=areamet(1,1)
          end do
        end do
        al=sqrt(areamet(1,1))
        sstep=al/1000.0
        step2=sstep**2
        istep=int(sstep)
        do i=1,ycount
          do j=1,xcount
            frac_2d(i,j)=frac_2d(i,j)*areamet(1,1)
            ch_length_2d(i,j)=sqrt(areamet(1,1))
          end do
        end do
      else    ! for latlong:
!       Jan. 19/08  NK
!       previously area of centre grid was used to get grid area
!       now grid area is a function of the latitude
        do i=1,ycount
          do j=1,xcount
!           calculate the latitude in radians
            latrad=(yorigin+ydelta*i)/180.*2.*acos(0.0)
!           calculate the grid area
           areamet(i,j)=ydelta*
     *       (111.1360-0.5623*cos(2.0*latrad)+0.0011*cos(4.0*latrad))*
     *                  xdelta*
     *       (111.4172*cos(latrad)-0.094*cos(3.0*latrad)+
     *        0.0002*cos(5.0*latrad))*1000.0*1000.0
!           frac_2d is the grid area in m**2
!  Version 10.9 - Nov. 2015 - 'al' calculation revised to centre grid from origin 
!            moved below
c            al=sqrt(areamet(1,1))
            frac_2d(i,j)=frac_2d(i,j)*areamet(i,j)
            ch_length_2d(i,j)=sqrt(areamet(i,j))
          end do
        end do
!  Version 10.9 - Nov. 2015 - 'al' calculation revised to centre grid from origin 
        al=sqrt(areamet(ycount/2,xcount/2))
        latrad=(yorigin+ydelta*ycount/2)/180.*2.*acos(0.0)
        delta_x=xdelta*
     *        (111.4172*cos(latrad)-0.094*cos(3.0*latrad)+
     *        0.0002*cos(5.0*latrad))
        delta_y=ydelta*
     *        (111.1360-0.5623*cos(2.0*latrad)+0.0011*cos(4.0*latrad))
        sstep=al/1000.0
        step2=sstep**2
        istep=int(sstep)
      endif


!  add a check on the grid squareness
!     Correcting for non contributing area
!     This s/r will only calculate the nca in each grid
!     Adjustments to the frac or redistributing the classes is done below
      if(adjust_frac.eq.'y')then
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
        call non_ca() ! added June 9/11 NK
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
      else
        print*
        print*,'frac will NOT be adjusted for nca'
        print*,'but the class areas may be depending on your answer'
        print*
        non_ca_flg=.false.

      endif
d      print*,'back from non_ca'

!     looking at the nca map, a lot of the mapped nca IS connected
!     there are rivers in the NCA!! and we are way underestimating on the
!     Battle river. So don't use 100% if the mapped NCA
!     effective_nca added Mar. 19/12 nk


!     Changed this to have the grid area = 1 m**2 when whole grid is nca
!     Oct. 11, 2012  NK
      if(non_ca_flg.and.effective_nca.gt.0.0)then
        do i=1,ycount
          do j=1,xcount
          frac_nca_2d(i,j)=amax1(1.0,                          !frac_2d(i,j)*0.001,
     *               frac_2d(i,j)*(1.0-nca(i,j)*effective_nca))
          end do
        end do
      endif

!  Version 10.5 - Jan. 2012 - create elv_means.r2s
!     read the DEM for the calculation of the grid elevation
!     for calculating lapse rates.
     
      INQUIRE(FILE='dem.r2s',EXIST=exists)
      IF(exists)THEN
!       read the header
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(answer_dem.eq.'y')then
          print*
          print*,'A DEM file has been found'
          print*,'so mean grid elevations will be calculated'
          print*,'and written to a files: elev_means.r2c'
          print*,'These mean elevations will be used for lapse rate'
          
          print*,'adjustments to temperature & precip if left in'
          print*,'the basin directory'
          print*
c        print*,'Do you want to create  new'
c        print*,' elev_means.r2c & elev_max.r2c files   y/n?'
c        print*
c        read(*,*)answer_dem
c        if(answer_dem.eq.'y')then
          print*,'reading dem.r2s file'
          fln(99)='dem.r2s'     
          call read_r2s_beta(99,fln(99)) !EnSim  r2c file
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
          print*,xdelta,xdelta_temp,ydelta,ydelta_temp
          npoints=(xdelta/xdelta_temp+1)*(ydelta/ydelta_temp+1)
          print*,'no points in the dem.r2s file:',npoints 
          allocate(elv_dem(ycount,xcount),outarray(ycount,xcount),
     *           elv_sum(ycount,xcount),elv_min(ycount,xcount),
     *           elv_max(ycount,xcount),elv_count(ycount,xcount),
     *                                             stat=iAllocate)
          if(iAllocate.ne.0) STOP
     *           'Error with allocation of elv_dem in bsn @ 905'
          do i=1,ycount
            do j=1,xcount
              elv_sum(i,j)=0.0
              elv_min(i,j)=1000.0
              elv_max(i,j)=-1000.0
              elv_count(i,j)=0
            end do
          end do
!         calbulate the mean grid elevations
          do ii=1,ycount_temp
            if(mod(ii,ycount_temp/10).eq.0)then
              print*,'doing row',ii,'/',ycount_temp
            endif
            do jj=1,xcount_temp
             xx_temp=xorigin_temp+(jj-1)*xdelta_temp
             yy_temp=yorigin_temp+(ii-1)*ydelta_temp
             i=(yy_temp-yorigin)/ydelta+1
             j=(xx_temp-xorigin)/xdelta+1
             if(i.ge.1.and.j.ge.1.and.i.le.ycount.and.j.le.xcount)then
               if(inarray(ii,jj).gt.0.0)then
                 elv_sum(i,j)=elv_sum(i,j)+inarray(ii,jj)
                 elv_count(i,j)=elv_count(i,j)+1
                 elv_min(i,j)=amin1(elv_min(i,j),inarray(ii,jj))
                 elv_max(i,j)=amax1(elv_max(i,j),inarray(ii,jj))
               endif
               elv_count(i,j)=max(1,elv_count(i,j))
             endif
cd             write(876,87600)xx_temp,xorigin,j,yy_temp,yorigin,i
87600        format(2f10.5,i5,2f10.5,i5)          
           end do
          end do
          print*,'calculating means'
          do i=1,ycount
            do j=1,xcount
              if(elv_count(i,j).gt.0)then  
                elv_dem(i,j)=elv_sum(i,j)/float(elv_count(i,j))
              else
                  elv_dem(i,j)=-999.000
              endif
            end do
          end do
          print*,'Done calculating mean elevations'
          print*
!         new elv_means.r2c file written after grid is adjusted for size.
        endif
      else
        if(answer_dem.eq.'y')then
          print*
          print*,'WARNING:'
          print*,'The file dem.r2s is NOT found'
          print*,'so the elv_means.r2c file can not be created.'
          print*,'The implication is that lapse rate computations will'
          print*,'be based on the channel elevations and not the mean'
          print*,'grid elevations'
          print*
          answer_dem='n'
          print*,'answer_dem =',answer_dem
          pause 'Program pause in bsn @ 1017 - hit enter to continue' 
        endif
      endif

d      print*,'gone to arrange'
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
      call arrange()    
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
d      print*,'back from arrange'
d      print*,'gone to grade'
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
      call grade(slopemin,no_errors)  
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
d      print*,'back from grade'

c      if(0.eq.0)then
      if(lastgrid(1).eq.0)then
        lastgrid(1)=na    ! <<<<<<<do not use naa here - shd file gets truncated
        no_outlets=1
      else
!       the user picks the last grid(s) in the watershed. 
!       here we assign the receiving grid(s) which is the next one
        ix=0
        write(51,*)'revised outlet grids to receiving grids:'
        do lg=1,no_outlets
          i=yyy(lastgrid(lg))
          j=xxx(lastgrid(lg))
          write(51,*)lg,lastgrid(lg),' > ',next_2d(i,j)
          lastgrid(lg)=next_2d(i,j)
          if(lastgrid(lg).gt.na)then
            print*,'It looks like you have picked a lastgrid outside'
            print*,'the original modelled watershed - in one of the '
            print*,'receiving grids. Please reselect more upstream'
            print*,'na,naa,lastgrid',na,naa,lastgrid(lg),l
            print*
            ix=-1
          endif
        end do
        write(51,*)
        if(ix.lt.0)stop 'Program aborted'
      endif

! RANK WILL BE READ BY THE VARIABLE 'S' IN PROGRAM
! SIMPLE. NOTE THAT IN THIS PROGRAM S IS THE DRAINAGE
! DIRECTION!

! IF NTYPE .GE.1 THEN ALL VALUES ON IAK(I,J) MUST BE ZERO 
! BUT DATA IS STILL READ IN TO PRESERVE DATA STRUCTURE
! ANY VALUE OF IAK MAY BE ENTERED BUT WHEN NTYPE>0,
! ALL IAK VALUES ARE CHANGED TO 0
! RIVER CLASS:

      print*,xcount,ycount
      allocate(array(na,na))

      read(31,5000)comment
      if(newversion.eq.0)write(40,5000)comment
!      write(*,5000)comment
      do i=ycount,1,-1
         if(newversion.eq.0)then
           read(31,1003,iostat=ios)(iak(i,j),j=1,xcount)
           if(ios.ne.0)then
             print*,' Problems reading the river class on line',i
             print*
             STOP ' Program aborted'
           endif
           write(40,1003,iostat=ios)(iak(i,j),j=1,xcount)
           write(*,1003,iostat=ios)(iak(i,j),j=1,xcount)
         else
!          new free format files
           read(31,*,iostat=ios)(iak(i,j),j=1,xcount)
         endif
      end do
!     Find number of river classes in the map file:
      do i=1,ycount
        do j=1,xcount
          nrvr=max0(nrvr,iak(i,j))
        end do
      end do
      print*,'No of river classes found in the map file =',nrvr
      print*,'This should match the number specified in the par file'
      print*
      pause 'Hit enter if ok'
      
!     iak = ibn 
!     ibn needed in rdpar later
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        ibn(n)=iak(i,j)
!       diagnostic added Sept. 16/11 nk        
        if(ibn(n).le.0.or.ibn(n).gt.na)then
          print*,'iak value not defined rank=',n,ibn(n)
          print*,'please correct the map file'
          pause
        endif
      end do

!     This added April 1/07  nk
!     was undefined earlier
      do n=naa+1,na
        ibn(n)=0
      end do

!     TS - ALLOCATIONS OF AREA4A ARRAYS
!     ntype for the number of land cover classes
!     nrivertype for the number of channel or basin types
!     moved here from spl9  nk 06/07/00
!     then moved from rdpar nk 28/12/04
!     moved back from rdshed 27/07/06 because needed for bsn.for  nk

      if(iopt.ge.1)print*,'nrvr=',nrvr

!     in spl these allocations are done in read_shed_ef
!     but in this creates this file so it has to be done here

      allocate(mndr(na),r1(na),r2(na),r1n(na),r2n(na),
     *aa2(na),aa3(na),aa4(na),rivtype(na),theta(na),widep(na),
     *kcond(na),flz(na),pwr(na),pwr2(nrvr),pool(na),rlake(na),
     *stat=iAllocate)
      if(iAllocate.ne.0) STOP
     *   'Error with allocation of area4 arrays in bsn'
      
      write(51,5000)comment
      write(51,1006)j,j,(j,j=1,xcount)
      write(51,1006)j,j,(jxmin+istep*(j-1),j=1,xcount) 
      do i=ycount,1,-1
         write(51,1006)iymin+istep*(i-1),i,(iak(i,j),j=1,xcount)
      end do
! CONTOUR COUNT FOR INTERNAL SLOPE:
      read(31,5000)comment
      if(newversion.eq.0)write(40,5000)comment
!      write(*,5000)comment
      do i=ycount,1,-1
         if(newversion.eq.0)then
           read(31,1003,iostat=ios)(irough_2d(i,j),j=1,xcount)
           write(40,4003,iostat=ios)(irough_2d(i,j),j=1,xcount)
         else
!          new free format files
           read(31,*,iostat=ios)(irough_2d(i,j),j=1,xcount)
         endif
         if(ios.ne.0)then
           print*,' Problems reading the no countours on line',i
           print*
           STOP ' Program aborted'
         endif
      end do

!      write(*,1125)
!      do i=ycount,1,-1
!         write(*,1003)(irough_2d(i,j),j=1,xcount)
!      end do

      write(51,5000)comment
      write(51,1006)j,j,(j,j=1,xcount)
      write(51,1006)j,j,(jxmin+istep*(j-1),j=1,xcount) 
      do i=ycount,1,-1
        write(51,1006)iymin+istep*(i-1),i,(irough_2d(i,j),j=1,xcount)
      end do

! NO OF CHANNELS
      read(31,5000)comment
      if(newversion.eq.0)write(40,5000)comment
!      write(*,5000)comment
      do i=ycount,1,-1
         if(newversion.eq.0)then
           read(31,1003,iostat=ios)(ichnl_2d(i,j),j=1,xcount)
           write(40,1003,iostat=ios)(ichnl_2d(i,j),j=1,xcount)
         else
!          new free format files
           read(31,*,iostat=ios)( ichnl_2d(i,j),j=1,xcount)
         endif
      end do
         if(ios.ne.0)then
           print*,' Problems reading the no channels on line',i
           print*
           STOP ' Program aborted'
         endif

!      write(*,1126)
!      do i=ycount,1,-1
!         write(*,1003)(ichnl_2d(i,j),j=1,xcount)
!      end do
      write(51,5000)comment
      write(51,1006)j,j,(j,j=1,xcount)
      write(51,1006)j,j,(jxmin+istep*(j-1),j=1,xcount) 
      do i=ycount,1,-1
         write(51,1006)iymin+istep*(i-1),i,(ichnl_2d(i,j),j=1,xcount)
      end do
      
! DWOPER INPUT NODE
      read(31,5000)comment
      if(newversion.eq.0)write(40,5000)comment
!      write(*,5000)comment

      do i=ycount,1,-1
         if(newversion.eq.0)then
           read(31,1003,iostat=ios)(ireach_2d(i,j),j=1,xcount)
           write(40,4003,iostat=ios)(ireach_2d(i,j),j=1,xcount)
         else
!          new free format files
           read(31,*,iostat=ios)(ireach_2d(i,j),j=1,xcount)
         endif
         if(ios.ne.0)then
           print*,' Problems reading the reach no on line',i
           print*
           STOP ' Program aborted'
         endif
      end do

!      write(*,1126)
!      do i=ycount,1,-1
!         write(*,1003)(ireach_2d(i,j),j=1,xcount)
!      end do

      write(51,5000)comment
      write(51,1006)j,j,(j,j=1,xcount)
      write(51,1006)j,j,(jxmin+istep*(j-1),j=1,xcount) 
      do i=ycount,1,-1
        write(51,1006)iymin+istep*(i-1),i,(ireach_2d(i,j),j=1,xcount)
      end do
      print*,'ntype=',ntype
!     pause 7

d     print*,'gone to ftch to calculate fetch for lake evaporation' 
      print*,'Gone to fetch'
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
      call ftch()    
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
      print*,'Back from fetch'
      print*


      allocate(aclass_3d(ycount,xcount,ntype+2),   !one more to accomodate split
     *        aclass_comb(ycount,xcount,ntype+1),
     *        aclass_read(ycount,xcount,ntype+1),
     *        ijtemp(ycount,xcount),
     *        class_name(ntype+2),
     *        keepflg(na),
     *        stat=iAllocateStatus)
      if (iAllocateStatus .ne. 0) STOP 
     *    '**Allocation failed  in bsn @1197**' 



! PERCENTAGES OF CLASSES WHEN NTYPE > 0
      print*,'Reading the class names as listed in the attribute list'
      print*,'of the map file:'
      if(ntype.ge.1)then
!       read the land cover percentages in the map file
!       do ii=1,ntype
        do ii=1,ntype+1
          read(31,*,iostat=ios)class_name(ii)
          comment=class_name(ii)
          if(ios.ne.0)then
            print*,'ios=',ios
            if(ios.eq.-1)then
              print*,'premature end of file found in map file'
              print*,'looking for class no. ii'
              print*,'no classes found =',ii-1
              print*,'classcount in the map file header =',ntype+1
              stop 'Program aborted in bsn @ 1293'
            endif
            print*,'unknown error - program continues'
          endif
          if(newversion.eq.0)write(40,*)class_name(ii)
          write(*,*)'map file class name: ',ii,class_name(ii)
c          write(51,*)class_name(ii)
          do i=ycount,1,-1
            if(newversion.eq.0)then
              read(31,4050,iostat=ios)(ijtemp(i,j),j=1,xcount)
              write(40,4050)(ijtemp(i,j),j=1,xcount)
              do j=1,xcount
                aclass_read(i,j,ii)=float(ijtemp(i,j))
              end do
            else
!             new free format files
              read(31,*,iostat=ios)(ijtemp(i,j),j=1,xcount)
              do j=1,xcount
                aclass_read(i,j,ii)=float(ijtemp(i,j))
              end do
            endif
            if(ios.ne.0)then
              print*,
     *          ' Problems reading the class % on line',i
              print*,' in class= ',comment
              print*
              STOP ' Program aborted'
            endif
          end do

c          write(51,5000)comment
          write(51,*)comment(1:20),' class # read:  ',ii
          write(51,1057)j,j,(j,j=1,xcount)
          write(51,1057)j,j,(jxmin+istep*(j-1),j=1,xcount) 
          do i=ycount,1,-1  
            write(51,1054)iymin+istep*(i-1),i,
     *       (aclass_read(i,j,ii),j=1,xcount)
          end do

        end do
        print*,'Finished reading the class list in the map file'
        
      else
!       NTYPE=0 ... NO LAND COVER AMOUNTS READ IN
!       SET ALL VALUES FOR ACLASS 1 TO 100%
        ntype=1
        do i=1,ycount
          do j=1,xcount
            aclass_read(i,j,1)=100.0
          end do
       end do
      endif

      write(51,*)'map file read in'
      write(51,*)'map file read in'
      write(51,*)'map file read in'

!     this part is to combine multiple classes as given by the 
!     land cover data. For instance, combining savanna, grassland & pasture 
!     into one class.
      INQUIRE(FILE='class_combine.csv',EXIST=exists)
      IF(exists)THEN
        open(unit=99,file='class_combine.csv',status='old')
!       check for file version number
        read(99,*,iostat=ios)junk,class_combine_version
        print*,'class_combine_version ',class_combine_version
        if(ios.ne.0)then
          print*,'Error on first line of the class_combine.csv file'
          print*,'values found',junk,class_combine_version
          stop 'Progran aborted in bsn @ 1410'
        endif
        read(99,*,iostat=ios)junk1
        if(ios.ne.0)then
          print*,'Error on 2nd line of the class_combine.csv file'
          print*,'iostat = ',ios
          print*,'values found:  ',junk
          stop 'Progran aborted in bsn @ 1416'
        endif
        if(class_combine_version.eq.2)then
          print*
          print*,'Version 2 of the class_combine file found'
!         get the new class count:  
          n=0
          new_class_count=0
          write(51,*)'old_class  new_class (assigned)'
          
          do while(.not.eof(99))
            read(99,*,iostat=ios)junk,junk,old_class,new_class
            if(ios.ne.0)then
              print*,'Error in version 2 class_combine.csv file'
              print*,'Check there are no blanks in the names'
              print*,'Last looked at data line #',n+1
              print*
c              stop 'BSN.exe aborted at 1169'
            endif  
            write(51,*)old_class,new_class
            write(*,*)old_class,new_class
            if(ios.eq.0)then
              n=n+1
              new_class_count=max(new_class_count,new_class)
            endif
          end do
          print*
          print*,'In class_combine.csv:'
          print*,'classes found =',n, 'new_class_count=',new_class_count
          print*
        endif
      else    
        INQUIRE(FILE='class_combine.txt',EXIST=exists)
        IF(exists)THEN
          open(unit=99,file='class_combine.txt',status='old')
          print*
          print*,'Version 1 of the class_combine file found'
!         get the new class count:  
          n=0
          new_class_count=0
          write(51,*)'old_class  new_class (assigned)'
          do while(.not.eof(99))
            read(99,*,iostat=ios)old_class,new_class
            write(51,*)old_class,new_class
            if(ios.eq.0)then
              n=n+1
              new_class_count=max(new_class_count,new_class)
            endif
          end do
          print*
          print*,'In class_combine.txt:'
          print*,'classes found =',n, 'new_class_count=',new_class_count
          print*
        endif
      endif
      
!     Added Jul. 16/18  NK
      if(n.ne.ntype+1)then
          print*,'# of entries in the class_combine.csv file '
          print*,'greater than the # of entries in the shd file'
          print*,'     =',ntype+1
          Print*,'# of entries in the shd file is determined by'
          print*,'the # of classes present in the land cover map'
          print*,'# of entries in the class_combine & shd file must'
          print*,'be the same. Please edit the class_combine file'
          print*
          stop 'Program aborted in bsn.f @ 1708'
      endif

      print*,'please check this table carefully'
      print*,'Make sure the numbers march up with the table in the '
      print*,'class_combine.csv file'
      print*,'Warning: / are read as commas so change to _'
      print*,'********************************************'
      pause 'waiting.....'
      
      if(exists)then  ! either one is ok
        deallocate(aclass_comb)
        allocate(aclass_comb(ycount,xcount,n+1),
     *           stat=iAllocateStatus)
        if (iAllocateStatus .ne. 0) STOP 
     *    '**Allocation failed  in bsn @1158**' 

!       initialize
        do i=1,ycount
          do j=1,xcount
c            do jj=1,new_class_count
            do jj=1,n+1
              aclass_comb(i,j,jj)=0.0
            end do
          end do
        end do
        write(51,*)'   Classes read:'
        do jj=1,ntype+1
          write(51,*)'                            class no',jj
          do i=1,ycount
            write(51,876)(aclass_read(i,j,jj),j=1,xcount)
          end do
        end do
        write(51,*)'There will be',ntype+1,' classes in the shd file'
        write(*,*)'There will be',ntype+1,' classes in the shd file'
!       combine the classes  total read = n
        rewind 99
        if(class_combine_version.eq.2)then        
          read(99,*)junk
          print*,junk
          read(99,*)junk
          print*,junk
        endif
        do l=1,n
          if(class_combine_version.eq.2)then        
            read(99,*)junk,junk,jj,ii,class_name(ii)
            print*,jj,ii,class_name(ii)
          else
!           original version          
            read(99,*)jj,ii,class_name(ii)
          endif
c          if(ii.ne.0)then
!           if ii=0, no data was found for this class ( all 0%)
            do i=1,ycount
              do j=1,xcount
!               jj = class number in the map file
                if(jj.gt.0)then
                  aclass_comb(i,j,ii)
     *                =aclass_comb(i,j,ii)+aclass_read(i,j,jj)
                endif
              end do
            end do
c          endif
        end do
        close(unit=99)

!       reassign for printing to file:
        do i=1,ycount
          do j=1,xcount
            do jj=1,new_class_count
              aclass_read(i,j,jj)=aclass_comb(i,j,jj)
            end do
          end do
        end do
!       obtain the new value for ntype
        ntype=new_class_count-1


        write(51,*)
        write(51,*)'New class count classes'
        do jj=1,new_class_count
         write(51,*)'                         newclass no',jj
          do i=ycount,1,-1
            write(51,876)(aclass_read(i,j,jj),j=1,xcount)
          end do
        end do
876     format(999f5.0)
      endif
      print*
      print*,'End of map file reached and classes combined'  
      print*

      INQUIRE(FILE='SubBasinImpPc.txt',EXIST=exists)
      if(exists)then
!       new method - Feb. 02/18   NK
!       variable % imp. area can be used on a subwatershed basis
!       The sub-watersheds are delineated as channel classes          
!       This should give enough flexibility to do this approach
        open(unit=99,file='SubBasinImpPc.txt',status='old',iostat=ios)
!       count the number of entries
        read(99,*)  !read the header
        n=0
        do while(.not.eof(99))
            n=n+1
            read(99,*,iostat=ios)i,imperv(n)
            if(i.eq.n)then
              imperv(n)=imperv(n)*0.01
              pervious(n)=1.0-imperv(n)
              print*,n,imperv(n),pervious(n)
            endif
        end do
        close(unit=99,status='keep')
        print*,'Found ',n-1,' impervious area entries'
        pause 99999
!       This moved impervious area to grassed areas in % given in SubBasinImpPC.txt
!       Badly needed in areas like Toronto where Urban area is not differentiated as pervious & impervious        
        do i=1,ycount
          do j=1,xcount
            if(iak(i,j).gt.0)then
              aclass_read(i,j,1)=aclass_read(i,j,1)+
     *           aclass_read(i,j,ntype+1)*pervious(iak(i,j))
              aclass_read(i,j,ntype+1)=
     *           aclass_read(i,j,ntype+1)*imperv((iak(i,j)))
            endif  
          end do
        end do
      elseif(aimpr.gt.0.0)then
!       Original method:            
!       MODIFY THE BARREN AREAS BY ADDING 100-IMP % OF THE URBAN AREA
!       THAT IS PERVIOUS AND TAKE IMP % OF THE URBAN AREA AS IMPERVIOUS:
!       IMPERVIOUS PERCENTAGE:    
!       move part of the URBAN area to class 1
!       PERVIOUS PERCENTAGE:
        aimpr=aimpr*0.01
        perv=1.-aimpr
        do i=1,ycount
          do j=1,xcount
!           BUG FIXED HERE JUL 19/00 NK
            aclass_read(i,j,1)=
     *           aclass_read(i,j,1)+aclass_read(i,j,ntype+1)*perv
            aclass_read(i,j,ntype+1)=aclass_read(i,j,ntype+1)*aimpr
          end do
        end do
        print*,'Note:    impervious area > 0 in the header'
        print*,int(perv*100),'% of the impervious class (urban)' 
        print*,'has been subtracted from class',ntype+1
        print*,'and added to class 1'
        print*,'Class 1 should be a land cover compatible with'
        print*,'the pervious areas in urban areas (eg. grass)'
        print*
      endif
        
!     remap the land covers:
      do n=1,ntype-2
c        do i=1+srows,ycount-nrows
c          do j=1+wcols,xcount-ecols
        do i=1,ycount
          do j=1,xcount
            aclass_3d(i,j,n)=aclass_read(i,j,n)    
          end do
        end do
      end do

!     split the wetland class in to 2
!     one bog (uncoupled) and 2 fen (coupled)
      if(split.gt.0.0)then
        print*
        print*,'NEW  June 22/11'
        print*,'Wetland split options:'
        print*,'Option 1: Old way - split % of wetland - 1% minimum'
        print*,'Option 2: New way - split % =  max % wetland coupled'
        print*,'                  - wetland >  split % not coupled'
        print*,'                  - wetland <= split % all coupled'
        print*
c        print*,'Enter 1 or 2'
c        read*,reply
        reply=split_type
        print*,'You have entered ',reply
        if(reply.ne.'1'.and.reply.ne.'2')then
          print*,'invaled answer'
          print*,'default =2 assumed'
          reply='2'
          print*
          pause 'Hit enter to continue with this entry'
        endif

!       there will now be one more class
        ntype=ntype+1  
!       reorder class names           
        class_name(ntype+1)=class_name(ntype)
        class_name(ntype)=class_name(ntype-1)
        class_name(ntype-1)=class_name(ntype-2)

      write(51,*)'Wetland stuff:'

c        do i=1+srows,ycount-nrows
c          do j=1+wcols,xcount-ecols
        print*,'        row      column     class%     frac    criteria'
        do i=1,ycount
          do j=1,xcount
!           remap the last 2 classes
            aclass_3d(i,j,ntype)=aclass_read(i,j,ntype-1)    !water
            aclass_3d(i,j,ntype+1)=aclass_read(i,j,ntype)    !imperveous

c            if(reply.eq.'1')then
            if(split_type.eq.'1')then
c	print*,split,aclass_read(i,j,ntype-2)
!     Revision
!     As you can see - changed my mind a couple of times
!     Now if coupled wetland is less than .1 % of the grid area, no coupling
!     to prevent instability                
c              if(aclass_read(i,j,ntype-2)*split.lt.1.0)then
c              if(aclass_read(i,j,ntype-2)*split.lt.0.01)then
c              if(aclass_read(i,j,ntype-2)*split.lt.0.1)then   ! note: frac = interger
                
!     Canged Aug. 010/23 to make it more restrictive.                
             if(aclass_read(i,j,ntype-2)*split*frac_2d(i,j)/
     *                       areamet(i,j).lt.0.1)then   ! note: frac = interger
!               small wetlands are NOT coupled!!
             if(aclass_read(i,j,ntype-2).gt.0.0)then
                print*,i,j,aclass_read(i,j,ntype-2),
     *          frac_2d(i,j)/areamet(i,j),
     *          aclass_read(i,j,ntype-2)*split*frac_2d(i,j)/areamet(i,j)
             endif   
!               I.e. if wetlands area small, they are assumed to be bogs
!               aclass_read(i,j,ntype-2) is left unchanged
!               this way there will be no coupled wetlands lt 1%
                aclass_3d(i,j,ntype-2)=aclass_read(i,j,ntype-2)
                aclass_3d(i,j,ntype-1)=0.0
                
               else
!               wetlands are split into coupled & uncoupled
!               and a class is added 
                aclass_3d(i,j,ntype-2)=
     *                        aclass_read(i,j,ntype-2)*(1.0-split)
                aclass_3d(i,j,ntype-1)=aclass_read(i,j,ntype-2)*split
              endif
            else
              tmp_value=aclass_read(i,j,ntype-2)
              aclass_3d(i,j,ntype-1)=
     *		amin1(split*100.0,aclass_read(i,j,ntype-2))

              aclass_3d(i,j,ntype-2)=
     *        amax1(0.0,aclass_read(i,j,ntype-2)-aclass_3d(i,j,ntype-1))
            endif
              
            write(51,*)rank_2d(i,j),i,j,split,
     *     	   aclass_3d(i,j,ntype-2),aclass_3d(i,j,ntype-1)
     
            if(ireach_2d(i,j).gt.0)then
!             if in a lake or reservoir
!             add the coupled wetland to the lake 
!             as we can't couple wetland to lakes.
!             Changed Mar. 07/08 -nk-
              aclass_3d(i,j,ntype)=
     *                  aclass_3d(i,j,ntype)+aclass_3d(i,j,ntype-1)
              aclass_3d(i,j,ntype-1)=0.0
            endif
          end do
        end do
        print*
        print*,'If values appear above, coupled wetlands removed'
        print*
        print*,'The value entered for "split" =',int(split*100.)
        print*,'Another land cover class was be added to the'
        print*,'new_shd.r2c file'
        Print*,'Please ensure that the par & sdc files are modified'
        print*,'to accomodate the extra wetland class'
        print*
        pause 33333

      else
!       remap the land covers:
        do n=ntype-1,ntype+1
c          do i=1+srows,ycount-nrows
c            do j=1+wcols,xcount-ecols
           do i=1,ycount
            do j=1,xcount
              aclass_3d(i,j,n)=aclass_read(i,j,n)    
            end do
          end do
        end do
      endif
      
      write(51,*)'Wetland split now incorporated'
      frame=0
      do ii=1,ntype+1
      frame=frame+1
          write(51,*)'frame',frame
c        aclass_3d(1,1,ii)=frame
c        do i=1+srows,ycount-nrows
c          write(51,4008)(aclass_3d(i,j,ii),j=1+wcols,xcount-ecols)
        do i=1,ycount
          write(51,4008)(aclass_3d(i,j,ii),j=1,xcount)
        end do
        print*,' frame= ',frame,' written to bsn_info'
      end do

! MODIFY THE PERCENTAGES TO MAKE SURE THEY ADD UP TO 100%
! This is where percentages are converted to fractions !!!!!!!!!!!!!!!!  % -> fraction

      do i=1,ycount
         do j=1,xcount
            sumclass=0.0
            do ii=1,ntype+1
               sumclass=sumclass+aclass_3d(i,j,ii)
            end do
!            if(sumclass.gt.0.0)then 
            if(sumclass.ne.0.0)then   ! changed Mar. 23/06 nk
              do ii=1,ntype+1
                 aclass_3d(i,j,ii)=aclass_3d(i,j,ii)/sumclass
              end do
            endif
         end do
      end do

      write(51,*)'passed line 1375'

!     check that grids with 100% water have been assigned a reach number
!     this check added 24/10/05  nk
      error_count=0      
      do n=1,na
        i=yyy(n)
        j=xxx(n)
        if(aclass_3d(i,j,ntype).ge.0.99999)then
          if(ireach_2d(i,j).eq.0)then
            error_count=error_count+1
            if(error_count.eq.1)then
              print*,'A grid with 100% water has not been assigned '
              print*,'a reach number. Program will crash if you try'
              print*,'to use a resume file'
              write(51,*)'A grid with 100% water has not been assigned'
              write(51,*)'(a) reach number. Program will crash '
              write(51,*)'if you try to use a resume file'
              pause 'Hit enter to continue but you have been warned'
            endif
            write(51,*)'grid,row,col',n,i,j
            print*,'grid,row,col',n,i,j
c            print*
            aclass_3d(i,j,1)=0.01
            aclass_3d(i,j,ntype)=0.99
          endif
        endif
      end do
      if(error_count.ge.1)then
            print*,error_count,' grid(s) with 100% water has(ve)'
        print*,' not been assigned '
            print*,'a reach number(s). The water class has been changed'
        print*,'99% and class 1 has been changed to 1%'
        pause 'Hit enter to continue but you have been warned AGAIN!'
      endif

      write(51,*)'passed line 1411'

     
! READ THE CHANNEL BANKFULL CAPACITIES IF SET UP IN THE MAP FILE:
      read(31,5000,iostat=ios)comment
!      write(40,5000)comment
      write(*,5000)comment
      do i=ycount,1,-1
         read(31,*,iostat=ios)(bnkfll_2d(i,j),j=1,xcount)
!         write(40,*,iostat=ios)(bnkfll_2d(i,j),j=1,xcount)
!         write(*,1053)(bnkfll_2d(i,j),j=1,xcount)
      end do

      write(51,*)'passed line 1424'

      if(ios.ne.0)then
         print*,'ios=',ios
         print*
         print*,' No bankfull values found'
         print*,' Default assumed'
!        WE HAVE NOT FOUND BANKFULL DISCHARGES IN THE MAP FILE
!        SO WE CALCULATE THE DEFAULT BANKFULL CAPACITIES INSTEAD
!        FOR NOW, JUST THE DRAINAGE AREA DIVIDED BY 6
!         write(*,6012)
         do i=ycount,1,-1
            do j=1,xcount
              bnkfll_2d(i,j)=da_2d(i,j)/6.0+0.1
            end do
         end do

       else
         print*
         print*,'Bankfull values found & will be used'
         print*,'Most people do not produce these values'
         print*,'Did you?   Just checking.'
         print*
         pause 'Hit enter to continue'
      endif

      write(51,*)'passed line 1447'



! THESE VALUES ARE NO LONGER USED IN SPL BUT ARE EMBEDDED 
! SHD FILE FORMAT SO ARE KEPT FOR COMPATIBILITY WITH OLD FILES.
      ls=0
      ks=0
      js=0
      ih=0
      local=0

      comment=' Grid ranking - highest =1'
      write(51,1006)
      write(51,5000)comment
      write(51,1006)j,j,(j,j=1,xcount)
      write(51,1006)j,j,(jxmin+istep*(j-1),j=1,xcount) 
      
      do i=ycount,1,-1
        write(51,1058)(rank_2d(i,j),j=1,xcount)
      end do

!      pause 1
      write(51,*)'passed line 1479'


!     error checking
      do n=1,na
         errflg=0 
         i=yyy(n)
         j=xxx(n)
         if(frac_2d(i,j).le.0.0)then
             write(*,1129)i,j,frac_2d(i,j)
             write(51,1129)i,j,frac_2d(i,j)
             errflg=1
         endif

         if(iak(i,j).le.0.or.iak(i,j).gt.16)then
           write(*,1014)n,i,j,elv_2d(i,j)
           write(51,1014)n,i,j,elv_2d(i,j)
           errflg=1
         endif

         if(irough_2d(i,j).eq.0)then
          write(*,1015)n,i,j,elv_2d(i,j)
           write(51,1015)n,i,j,elv_2d(i,j)
           errflg=1
         endif

         if(ichnl_2d(i,j).le.0.or.ichnl_2d(i,j).gt.5)then
           write(*,1016)n,i,j,elv_2d(i,j)
           write(51,1016)n,i,j,elv_2d(i,j)
                  errflg=1
        endif

        if(next_2d(i,j).eq.0)then
          write(*,1017)n,i,j,elv_2d(i,j)
          write(51,1017)n,i,j,elv_2d(i,j)
          write(*,1018)
          write(51,1018)
                  write(*,1019)
                  print*
                  write(51,1019)
                  write(*,10019)
                  write(51,10019)
          errflg=1
        endif

        if(na.le.naa)then
          write(*,*)'fatal error - outlet elements have not been'
          write(*,*)'properly designated'
          write(*,*)'set frac of the recieving element = 0'
          write(*,*)'frac =',frac_2d(xxx(na),yyy(na))
          write(*,*)'elevation of the receiving grid must be > 0.0'
          write(*,*)'location of outlet:'
          write(*,*)'na, x,y  =',na,xxx(na),yyy(na),
     *     xorigin+(xdelta*xxx(na))-xdelta/2,
     *     yorigin+(ydelta*yyy(na))-ydelta/2
          
          write(*,*)'naa,x,y  =',naa,xxx(naa),yyy(naa)
          
          write(51,*)'fatal error - outlet elements have not been'
          write(51,*)'properly designated'
          write(51,*)'set frac of the recieving element = 0'
          write(51,*)'elevation of the receiving grid must be > 0.0'
          write(51,*)'location of outlet:'
          write(51,*)'na  =',na
          write(51,*)'naa =',naa
          STOP ' program aborted in bsn @ 2195'
        endif

        if(next_2d(i,j).le.rank_2d(i,j))then
                  write(51,6015)i,j,rank_2d(i,j),elv_2d(i,j)
                  write(*,6015)i,j,rank_2d(i,j),elv_2d(i,j)
            errflg=1          
        endif
      end do

      if(errflg.gt.0)then
        print*
        print*,'These can be fatal errors  <<<<<<<<<<<<<<<<<<<<<<<<<<<'
        print*,'Please correct the map file to elliminate these errors'
        print*
        pause 'Hit enter to continue'
      endif

      nullgrids=na-naa

      if(iopt.eq.1)pause 2

!     set all next_2d() values =
      do n=na,naa+1,-1                   ! nk nov 20/02
        next_2d(yyy(n),xxx(n))=0
      end do

      write(51,*)'passed line 1549'




!     Correcting for non contributing area
c	call non_ca() ! added June 9/11 NK

c	if(non_ca_flg)then
c	  do i=1,ycount
c	    do j=1,xcount
c       	  frac_2d(i,j)=amax1(frac_2d(i,j)*(1.0-nca(i,j)),
c     *                         frac_2d(i,j)*0.001)
c          end do
c	  end do
c	endif




!     MAPPING THE SUB-BASIN!     MAPPING THE SUB-BASIN!     MAPPING THE SUB-BASIN
!     MAPPING THE SUB-BASIN!     MAPPING THE SUB-BASIN!     MAPPING THE SUB-BASIN
!     MAPPING THE SUB-BASIN!     MAPPING THE SUB-BASIN!     MAPPING THE SUB-BASIN
!     MAPPING THE SUB-BASIN!     MAPPING THE SUB-BASIN!     MAPPING THE SUB-BASIN
!     a sub-basin can only have one outlet!!!!

      write(51,*)
      write(51,*)'list of outlet grids:'
      do lg=1,no_outlets
        write(51,*)lg,lastgrid(lg)
c	  write(*,*)lg,lastgrid(lg)
      end do

      write(51,*)
      write(51,*)'n,lg,lstgrid,i,j,rank,next'
c      write(*,*)
c	write(*,*)'n,lg,lstgrid,i,j,rank,next'

!     set all rank_2d(i,j)  -ve
      do n=1,na
        i=yyy(n)
        j=xxx(n)
c        do lg=1,no_outlets
c          if(n.eq.lastgrid(lg))then
            rank_2d(i,j)=rank_2d(i,j)*-1
            next_2d(i,j)=next_2d(i,j)*-1
            keepflg(n)=.false.
            write(51,*)n,lg,i,j,rank_2d(i,j),next_2d(i,j)
c          endif
c	  end do
      end do
      
!     leave all rank_2d() and next_2d() outside basin -ve  
!     work from the lastgrid up the sub-watershed to turn 
!     rank_2d() and next_2d() positive

!     Set outlet grids +ve
      write(51,*)'lg,lstgrid,i,j,rank,next'
      do lg=1,no_outlets
        i=yyy(lastgrid(lg))
        j=xxx(lastgrid(lg))
!       make sure they are not already +ve
        if(rank_2d(i,j).lt.0.and.next_2d(i,j).lt.0)then
          rank_2d(i,j)=rank_2d(i,j)*-1
          next_2d(i,j)=next_2d(i,j)*-1
c!         keep the outlet/inlet grid
          keepflg(rank_2d(i,j))=.true.   
          write(51,*)lg,lastgrid(lg),i,j,rank_2d(i,j),next_2d(i,j)
        endif
      end do

      if(iopt.eq.1)pause 32

!     now set rank_2d() & next_2d() back to +ve for grids in watershed(s)
!     BUT LAST GRID IS ALREADY +VE !!!!!
!     flag all the grids that are above the outlet. 
      do lg=1,no_outlets
        do n=lastgrid(lg),1,-1
          i=yyy(n)
          j=xxx(n)
!         next conditional took a while to figure out Jul. 11/05
          if(lastgrid(lg).le.naa)then   
!           this treats the receiving grids too
            nextgrid=abs(next_2d(i,j))
            k=yyy(nextgrid)
            l=xxx(nextgrid)
!           this puts the upstream grid in the watershed
c            if(rank_2d(k,l).gt.0)then
!           make sure the grid has not already been made positive
!           this can happen with multiple basins
            if(rank_2d(k,l).gt.0.and.rank_2d(i,j).lt.0)then
c     *         		  .and.next_2d(i,j).lt.0
              next_2d(i,j)=next_2d(i,j)*-1    
              rank_2d(i,j)=rank_2d(i,j)*-1  
              keepflg(n)=.true.   
            endif
          else
!           for the case where the whole watershed is included
!           all receiving grids need to be included
!           but the last grid rank is already +ve
            next_2d(i,j)=next_2d(i,j)*-1    
!           added this conditional - took 2 days!  Sept. 13/05 NK
            if(rank_2d(i,j).lt.0)then
              rank_2d(i,j)=rank_2d(i,j)*-1     
              keepflg(n)=.true.   
            endif
          endif
        end do                             
          end do
          
!  Version 10.7 - Jul. 2012 - added inlet grids
!     Remove grids upstream from inlet grids
      if(no_inlets.gt.0)then
!       temporarily set inlet grid rank -ve      
        do lg=1,no_inlets
          n=inletgrid(lg)
          if(n.le.0)then
              print*,'error: inlet grid #',lg,n
              print*,'is not inthe watershed.'
              stop 'Program aborted in bsn @ 2351'
          endif
          i=yyy(n)
          j=xxx(n)
          rank_2d(i,j)=abs(rank_2d(i,j))*-1  
        end do
!       do the rest
        do lg=1,no_inlets
          do n=inletgrid(lg),1,-1
            i=yyy(n)
            j=xxx(n)
!           don't touch inlet grid
            if(n.ne.inletgrid(lg))then   
              nextgrid=abs(next_2d(i,j))
              k=yyy(nextgrid)
              l=xxx(nextgrid)
!             this takes upstream grid out of the watershed
!             make sure the grid has not already been made positive
!             this can happen with multiple basins
              if(rank_2d(k,l).lt.0)then
                next_2d(i,j)=abs(next_2d(i,j))*-1    
                rank_2d(i,j)=abs(rank_2d(i,j))*-1  
                keepflg(n)=.false.   
              endif
            endif
          end do
        end do
!       set inlet grid rank back to +ve      
        do lg=1,no_inlets
          n=inletgrid(lg)
          i=yyy(n)
          j=xxx(n)
          rank_2d(i,j)=abs(rank_2d(i,j)) 
        end do
      endif
!     Doneth      

      write(51,*)'*****************************************************'
      write(51,*)'rank'
      do i=1,ycount
        write(51,10)(max(0,rank_2d(i,j)),j=1,xcount)
c	  write(10,10)(rank_2d(i,j),j=1,70)
10      format(9995i5)
      end do
      write(51,*)'next'
      do i=1,ycount
        write(51,10)(max(0,next_2d(i,j)),j=1,xcount)
c	  write(10,10)(next_2d(i,j),j=1,70)
      end do

      if(iopt.eq.1)pause 33

      do n=1,na
        newnext(n)=0
        newrank(n)=0
      end do

!     count the number of grids in the sub-watershed (could be all)
!     and create a new coordinate matrix
      nogrids=0
      do n=1,na
        i=yyy(n)
        j=xxx(n)
        if(rank_2d(i,j).gt.0)then
          nogrids=nogrids+1
          newyyy(nogrids)=yyy(n)
          newxxx(nogrids)=xxx(n)
          newrank(n)=nogrids
        endif
      end do

      if(iopt.eq.1)pause 4

!     reassign the next grid numbering
      do n=1,na
        i=yyy(n)
        j=xxx(n)
        if(rank_2d(i,j).gt.0)then
c          if(n.lt.lastgrid(lg))then
          if(n.le.naa)then
c          if(n.lt.naa)then
            newnext(n)=newrank(next_2d(i,j))
          else
            newnext(n)=0
          endif
c      print*,newnext(n),next_2d(i,j)
        endif
      end do

c      if(iopt.eq.1)pause 5

!     reassign the rank number to remaning grids
      do n=1,na
        i=yyy(n)
        j=xxx(n)
        if(rank_2d(i,j).gt.0)then
          rank_2d(i,j)=newrank(n)
        else
          rank_2d(i,j)=0
        endif
c      print*,rank_2d(i,j),newrank(n)
      end do
c      pause 'newrank'

      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        if(rank_2d(i,j).gt.0)then
          next_2d(i,j)=newnext(n)
        else
          next_2d(i,j)=0
        endif
c      print*,rank_2d(i,j),newnext(n)
      end do
c      pause 'newnext'

!     set frac for outlet grids = 0 to bypass in splx
      do n=1,na
        if(newnext(n).eq.0)then
          i=yyy(n)
          j=xxx(n)
          frac_2d(i,j)=0.0
        endif
      end do

c      i=yyy(lastgrid(lg))
c      j=xxx(lastgrid(lg))
c	next_2d(i,j)=0

c	print*,lastgrid(lg),i,j,next_2d(i,j)
c      pause '-----------------------'

      if(iopt.eq.1)pause 6

      
      print*,'Would you like to remove blank rows & columns  y/n?'
      print*,'Answer `n` if you want to use precit & tmp domain set up'
      print*,'for the whole map file domain'
      print*,'Answer `y` if you have created or will set up the '
      print*,'the precip & tmp files to match the sub-basin domain'
      print*
      read*,answer
      nrows=0
      srows=0
      wcols=0
      ecols=0
      if(answer.eq.'y')then
          
!     Remove the blank rows and colums
!     Check for number of blank rows and columns:

      i1flg=0   ! change to 1 as soon as we find an entry
!     scan the rows until we find the first entry     
      do i=1,ycount
        do j=1,xcount
          if(rank_2d(i,j).gt.0.and.i1flg.eq.0)then
!           we have found the first line with data
            i1flg=1
            srows=i-1
          endif
        end do
      end do

      if(iopt.eq.1)pause 7

!     do it in reverse from the top (north) down
      i1flg=0
      do i=ycount,1,-1
        do j=1,xcount
          if(rank_2d(i,j).gt.0.and.i1flg.eq.0)then
!           we have found the first line with data
            i1flg=1
            nrows=ycount-i
          endif
        end do
      end do

      if(iopt.eq.1)pause 8

      j1flg=0   ! change to 1 as soon as we find an entry
!     scan the rows from the west until we find the first entry     
      do j=1,xcount
        do i=1,ycount
          if(rank_2d(i,j).gt.0.and.j1flg.eq.0)then
!           we have found the first line with data
            j1flg=1
            wcols=j-1
          endif
        end do
      end do

      if(iopt.eq.1)pause 9

!     do it in reverse from the east
      j1flg=0
      do j=xcount,1,-1
        do i=1,ycount
          if(rank_2d(i,j).gt.0.and.j1flg.eq.0)then
!           we have found the first line with data
            j1flg=1
            ecols=xcount-j
          endif
        end do
      end do

      if(iopt.eq.1)pause 10

!     This will not reduce the border to less than 2 grids
      nrows=max(0,nrows-2)
      srows=max(0,srows-2)
      wcols=max(0,wcols-2)
      ecols=max(0,ecols-2)

      if(nrows.gt.0.or.srows.gt.0.or.wcols.gt.0.or.ecols.gt.0)then
        print*
        print*
        print*,'Grid size for the .shd file has been reduced by'
        print*
        print*,'No. of rows removed north side = ',nrows
        print*,'No. of rows removed south side = ',srows
        print*,'No. of columns removed west side = ',wcols
        print*,'No. of columns removed east side = ',ecols
        print*
        
! New Oct. 04/2019        
        print*,'Do you wish to proceed with these numbers? y/n'
        print*,'Use  n  if setting up a sub-basin and wish to use'
        print*,'the same domain as the whole shd file and use the'
        print*,'same precip & temperature files '
        read*,answer
        if(answer.eq.'n')then
          print*,'Use same numbers as for the whole domain'  
          print*
          print*,'How many rows to remove in the north?'  
          read*,nrows
          print*,'How many rows to remove in the south?'  
          read*,srows
          print*,'How many columns to remove in the west?'  
          read*,wcols
          print*,'How many colums to remove in the east?'  
          read*,ecols
        endif  
        
        ycount_new=ycount-srows-nrows
        xcount_new=xcount-wcols-ecols

        print*,'xcount=',xcount
        print*,'ycount=',ycount
        print*
        print*,'There will now be ',xcount_new,' columns'
        print*,'                  ',ycount_new,' rows and'
        print*,'This is done to reduce the size of the output files'
        print*
      endif

      endif   !'remove rows & cols ?'
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
!     changed to use the read_par_parser  Dec. 2/11 nk
      classcount=ntype+1
      allocate(nsdc(classcount),snocap(classcount),idump(classcount),
     *        stat=iAllocate)
      if(iAllocate.ne.0) STOP
     *   'Error with allocation of areamelta arrays in bsn @ 158'

      if(iopt.eq.1)pause 11
      if(parflg.eq.'y')then
          call find_filetype(2)
          if(filetype.eq.'csv')then
!           **********************************************************************
            call read_par_parser(32,2)
!           **********************************************************************
            print*,'parameter file read'
          else
            print*
            print*,'Old format .par file no longer accepted'
            print*
            stop 'Program aborted in BSN @ 1975'
          endif
      endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
!     
!     set all values outside the watershed =0
      do i=1+srows,ycount-nrows
        do j=1+wcols,xcount-ecols
          if(rank_2d(i,j).eq.0)then
            next_2d(i,j)=0
            da_2d(i,j)=0
            bnkfll_2d(i,j)=0
            slope_2d(i,j)=0.0
            elv_2d(i,j)=0.0
            ch_length_2d(i,j)=0
            iak(i,j)=0
!           sl1(n)=???
            ichnl_2d(i,j)=0
            ireach_2d(i,j)=0
            frac_2d(i,j)=0.0
            do ii=1,ntype+1
              aclass_3d(i,j,ii)=0.0
            end do
          endif
        end do
      end do

      ireach_2d_max=0
      do i=1+srows,ycount-nrows
        do j=1+wcols,xcount-ecols
          if(rank_2d(i,j).eq.0)then
            ireach_2d_temp(i,j)=ireach_2d(i,j)
            ireach_2d_max=max(ireach_2d_max,ireach_2d_temp(i,j))
          endif
        end do
      end do
      print*,'ireach_2d_max=',ireach_2d_max
      
      
!     This is meant to be used for the first run wheen the river classes have not been split
!     i.e the initial default values.  This allows the wetlands to be set up as special class 
!     right away
c      if(nrvr.eq.1)then
c         print*,'Would you like to split fens into a seperate river'
c         print*,'class?      y/n'
c         read(*,*)answer
c          if(answer.eq.'y')then
c              print*,'What threshold amount - e.g.  0.05 ?'
c              read(*,*)threshold
c              do i=1+srows,ycount-nrows
c                  do j=1+wcols,xcount-ecols
c                      if(aclass_3d(i,j,ntype-1).gt.threshold)then
c                          iak(i,j)=2
c                      endif
c                  end do
c              end do
c              nrvr=2
c          endif
c      endif


! WRITE THE new_format.shd FILE TO DISK:
! WRITE THE new_format.shd FILE TO DISK:
! WRITE THE new_format.shd FILE TO DISK:
! WRITE THE new_format.shd FILE TO DISK:
! WRITE THE new_format.shdD FILE TO DISK:
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

!     Write header for new_format.shd file:

!     transfer the comment from the map to the shd file
      do i=1,n_hdr_lines
        write(36,5007)'#',hdrcomment(i)
      end do
      write(36,5012)
      call date_and_time(cday,time)
      write(36,5026)time(1:2),time(3:4),time(5:6),
     *              cday(7:8),cday(5:6),cday(1:4)
      write(36,5025)fn1  
      write(36,5012)
      write(36,5009)coordsys
      if(coordsys.eq.'LATLONG   ')then
        write(36,5010)datum1
      elseif(coordsys.eq.'UTM       ')then
        write(36,5010)datum1
        write(36,5011)zone
      endif
      write(36,5012)
      if(coordsys.eq.'LATLONG   ')then
        write(36,5013)xorigin+wcols*xdelta
        write(36,5014)yorigin+srows*ydelta
      else
        write(36,50131)xorigin+wcols*xdelta
        write(36,50141)yorigin+srows*ydelta
      endif
      write(36,5012)
      write(36,5015)xcount   !_new
      write(36,5016)ycount   !_new
      write(*,5015)xcount   !_new
      write(*,5016)ycount   !_new

      if(coordsys.eq.'LATLONG   ')then
        write(36,5017)xdelta
        write(36,5018)ydelta
      else
        write(36,50171)xdelta
        write(36,50181)ydelta
      endif
      write(36,5012)
      write(36,50182)al
      write(36,5019)cintv
      write(36,5020)aimpr
      write(36,5021)ntype
      write(36,5027)nrvr
      write(36,5022)elvconv
      write(36,5023)

!     fix fix  do these write(37 's need to be moved down?

      if(na.eq.nogrids)then
        write(37,7001)na,ycount,xcount,ls,ks,js,ijk,ih,local,
     *                       ycount/2,xcount/2,naa
        write(36,5028)na
        write(36,5029)naa
      else
!       for this case there can only be one outlet
        write(37,7001)
     *     nogrids,ycount,xcount,ls,ks,js,ijk,ih,local,
     *                         ycount/2,xcount/2,nogrids-1
        write(36,5028)nogrids
        write(36,5029)nogrids-1
      endif

      write(36,5030)min(na/2,nogrids/2)
      write(36,5033)     !  #  Note: -ve slopes not corrected in this file
      write(36,5032)     !  enddheader

      if(iopt.eq.2)pause 12

!     Old format shed file
!     Old format shed file
!     write header for old_ format.shd file

      write(37,7021) al,cintv,nrvr,ntype,fn1
!      write(37,7021) al,cintv,impr,ntype,fn1
!     impr is used only in bsn.for so not needed in spl
!     use space of nrvr - Oct. 16/01 NK

! CHECK FOR COORDINATE SYSTEM
      if(coordsys.eq.'UTM      ') then
!        FOR UTM COORDINATES
         write(37,7005)iymin+srows,iymax-nrows,jxmin+wcols,jxmax-ecols
      elseif(coordsys.eq.'CARTESIAN ')then
!        for CARTESIAN coordinates
         write(37,7005)iymin+srows,iymax-nrows,jxmin+wcols,jxmax-ecols
      else
!        FOR LAT-LONG COORDINATES
         iiiii=-1

!     the jx and iy's are reversed!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if(nrows.gt.0.or.srows.gt.0.or.wcols.gt.0.or.ecols.gt.0)then
          print*,' Old extent of the grid in minutes:'
          print*,iymin,iymax,jxmin,jxmax
        endif

!     reduce the watflood domain by removing blank rows & colums
c        iymin=iymin+int(real(srows)*yint)
c        iymax=iymax-int(real(nrows)*yint)
c        jxmin=jxmin+int(real(wcols)*xint)
c        jxmax=jxmax-int(real(ecols)*xint)

        iymin=iymin+int(real(srows)*ydelta)
        iymax=iymax-int(real(nrows)*ydelta)
        jxmin=jxmin+int(real(wcols)*xdelta)
        jxmax=jxmax-int(real(ecols)*xdelta)

        if(nrows.gt.0.or.srows.gt.0.or.wcols.gt.0.or.ecols.gt.0)then
          print*,' New extent of the grid in minutes:'
          print*,iymin,iymax,jxmin,jxmax
        endif

        latsdeg=int(real(iymin)/60.0)
        latsmin=iymin-int(real(latsdeg)*60.0)
        latndeg=int(real(iymax)/60.0) 
        latnmin=iymax-int(real(latndeg)*60.0)

        lonwdeg=int(real(jxmin)/60.0)
        lonwmin=jxmin-int(real(lonwdeg)*60.0)
        lonedeg=int(real(jxmax)/60.0)
        lonemin=jxmax-int(real(lonedeg)*60.0)

        if(nrows.gt.0.or.srows.gt.0.or.wcols.gt.0.or.ecols.gt.0)then
          print*,' New extent of the grid in degrees & minutes:'
          print*,latsdeg,latsmin,latndeg,latnmin
          print*,lonwdeg,lonwmin,lonedeg,lonemin  
c          print*,' x interval =',xint,'y interval =',yint
          print*,' x interval =',xdelta,'y interval =',ydelta
          print*
        endif

        write(37,1005)iiiii,iiiii,iiiii,iiiii,latsdeg,latsmin,latndeg,
     *                 latnmin,lonwdeg,lonwmin,lonedeg,lonemin,
     *                 ydelta,xdelta   !yint,xint
      endif
      
      if(iopt.eq.2)pause 13

!     this is already printed above
      write(37,7006) ycount,xcount
!      write(37,7006) (jl(i),jr(i),i=ib,it)

!     for both old & new formats (not r2c)

      write(36,*)'The ranks of each grid: highest = 1'
      do i=ycount,1,-1
        if(na.le.9999)then
c          write(36,1006)(rank_2d(i-nrows,j+wcols),j=1,xcount)
c          write(37,1006)(rank_2d(i-nrows,j+wcols),j=1,xcount)
        else
c          write(36,1058)(rank_2d(i-nrows,j+wcols),j=1,xcount)
c          write(37,1058)(rank_2d(i-nrows,j+wcols),j=1,xcount)
        endif
c        write(51,1006)iymin+istep*(i-1),i,(rank_2d(i,j),j=1,xcount)
      end do

!     splice from move of section above.

      if(iopt.eq.2)pause 14

      write(36,*)'   n, next, row, col,   da,    bankfull, cha_slope,',
     *'     elv, ch_lenth,iak,int_slope,chnl,reach,frac,imperv',
     *' classes 1 -',ntype 

!     added Jun. 1/04  NK
!     changed sl1(n) to the internal slope. 
!     Was sqrt if int. slope but that was misleading when 
!     looking at the shd file.  nk jul 27/04

c      write(1000,*)'contour interval =',cintv,'  in bsn @ 1946'
      do n=1,na
        i=yyy(n)
        j=xxx(n)
        sl1(n)=cintv*(float(irough_2d(i,j))+.0001)/al
c       write(1000,*)n,i,j,irough_2d(i,j),sl1(n)
      end do


!      FIX:  needs to be added to the new format file.
       print*,'na,naa/',na,naa

        do i=1+srows,ycount-nrows
          do j=1+wcols,xcount-ecols
            if(ireach_2d(i,j).gt.0)then
!             if in a lake or reservoir
!             add the coupled wetland to the lake 
!             as we can't copule wetland to lakes.
!             Changed Mar. 07/08 -nk-
              aclass_3d(i,j,ntype)=
     *                  aclass_3d(i,j,ntype)+aclass_3d(i,j,ntype-1)
              aclass_3d(i,j,ntype-1)=0.0
            endif
          end do
        end do

        do i=1,na
          do j=1,na
            array(i,j)=0
          end do
        end do
        do n=1,naa
          if(newnext(n).le.0)then
            print*,'For original rank =',n
            print*,'The new rank =',newrank(n)
            print*,'& the new next rank =',newnext(n)
            print*,'the new next rank can not be =< 0'
c            stop 'BSN aborted @ 2727'
          else       
            array(newrank(n),newnext(n))=n
          endif
        end do
!       note:outlet grids stay where they are        

       do n=1,na
        i=yyy(n)
        j=xxx(n)
        m=iak(i,j)
        if(n.gt.naa)newnext(n)=0
        if(n.eq.lastgrid(lg))da_2d(i,j)=0.0
        if(rank_2d(i,j).gt.0)then
          if(da_2d(i,j).gt.10.0)then
!           new_format.shd
            write(36,1025)newrank(n),newnext(n),
     *        (yyy(n)-srows),(xxx(n)-wcols),
     *        da_2d(i,j),bnkfll_2d(i,j),slope_2d(i,j),
     *        elv_2d(i,j),ch_length_2d(i,j),
     *        iak(i,j),sl1(newrank(n)),ichnl_2d(i,j),ireach_2d(i,j),
     *        frac_2d(i,j)/areamet(i,j),aclass_3d(i,j,ntype+1),
     *        (aclass_3d(i,j,ii),ii=1,ntype)
!
!     *        frac_2d(i,j),(aclass_3d(i,j,ii),ii=1,ntype+1)
!           channels included format
!           old_format.shd
            write(37,1010)newrank(n),yyy(n)-srows,xxx(n)-wcols,
     *        da_2d(i,j),bnkfll_2d(i,j),slope_2d(i,j),elv_2d(i,j),
     *        iak(i,j),irough_2d(i,j),ichnl_2d(i,j),
     *        newnext(n),ireach_2d(i,j),
     *        frac_2d(i,j)/areamet(i,j),aclass_3d(i,j,ntype+1),
     *        (aclass_3d(i,j,ii),ii=1,ntype)
!     *        frac_2d(i,j),(aclass_3d(i,j,ii),ii=1,ntype+1)
          else
!           new_format.shd
            write(36,1026)newrank(n),newnext(n),
     *        yyy(n)-srows,xxx(n)-wcols,
     *        da_2d(i,j),bnkfll_2d(i,j),slope_2d(i,j),
     *        elv_2d(i,j),ch_length_2d(i,j),
     *        iak(i,j),sl1(newrank(n)),ichnl_2d(i,j),ireach_2d(i,j),
     *        frac_2d(i,j)/areamet(i,j),aclass_3d(i,j,ntype+1),
     *        (aclass_3d(i,j,ii),ii=1,ntype)
!           channels included format
!           old_format.shd
            write(37,1011)newrank(n),yyy(n)-srows,xxx(n)-wcols,
     *        da_2d(i,j),bnkfll_2d(i,j),slope_2d(i,j),elv_2d(i,j),
     *        iak(i,j),irough_2d(i,j),ichnl_2d(i,j),
     *        newnext(n),ireach_2d(i,j),
     *        frac_2d(i,j)/areamet(i,j),aclass_3d(i,j,ntype+1),
     *        (aclass_3d(i,j,ii),ii=1,ntype)
!     *        frac_2d(i,j),(aclass_3d(i,j,ii),ii=1,ntype+1)
          endif
        endif
       end do


      close(unit=36,status='keep')
      close(unit=37,status='keep')

        print*,'2nd time'
        print*,'No. of rows removed north side = ',nrows
        print*,'No. of rows removed south side = ',srows
        print*,'No. of columns removed west side = ',wcols
        print*,'No. of columns removed east side = ',ecols
        print*

c      do i=1,na
c        write(600,6090)(array(i,j),j=1,na)
c      end do
c6090  format(<na>i3)      

!     pause 9


!     write the new_shd.r2c file in r2c format
!     write the new_shd.r2c file in r2c format
!     write the new_shd.r2c file in r2c format
!     write the new_shd.r2c file in r2c format
!     write the new_shd.r2c file in r2c format

      open(41,file='new_shd.r2c',status='unknown',iostat=ios)
      if(ios.ne.0)then
         print*,' Problems opening the new_shd.r2c file '
         print*,' Possible cause: file open or read only'
         print*
         STOP ' Program BSN aborted @ 1966'
      endif

      write(41,3005)'########################################'
      write(41,3005)':FileType r2c  ASCII  EnSim 1.0         '
      write(41,3005)'#                                       '
      write(41,3005)'# DataType               2D Rect Cell   '
      write(41,3005)'#                                       '
      write(41,3005)':Application             WATFLOOD       '
      write(41,3005)':Version                 10             '
      write(41,3020)':WrittenBy          ',author
      call date_and_time(cday,time)
      write(41,3011)':CreationDate       ',
     *       cday(1:4),cday(5:6),cday(7:8),time(1:2),time(3:4)
3011  format(a20,a4,'-',a2,'-',a2,2x,a2,':',a2)
      write(41,3005)'#                                       '
      write(41,3005)'#---------------------------------------'
      write(41,3020)':SourceFileName     ',fn1
      write(41,50182)al
      write(41,5019)cintv
      write(41,5020)aimpr
c      write(41,5021)ntype+1+split_no
      write(41,5021)ntype+1+nca_classcount
      write(41,5027)nrvr
      write(41,5022)elvconv
      if(na.eq.nogrids)then
        write(41,5028)na
        write(41,5029)naa
        write(41,5030)na/2
      else
!       for this case there can only be one outlet
        write(41,5028)nogrids
        write(41,5029)nogrids-1
        write(41,5030)nogrids/2
      endif
      write(41,3005)'#                                       '
      write(41,3003)'#effective nca %    ',effective_nca*100.0   
!  Version 10.9 - Nov. 2015 - 'al' calculation revised to centre grid from origin 
      write(41,3003)'#center deltaX  km  ',delta_x     
      write(41,3003)'#center deltaY  km  ',delta_y                 
      write(41,3005)'#                                       '
      write(41,3004)':Projection         ',coordsys
      if(coordsys.eq.'UTM       ')then
          write(41,3004)':Zone               ',zone
      endif
      write(41,3004)':Ellipsoid          ',datum1
      write(41,3005)'#                                       '
      write(41,3003)':xOrigin            ',xorigin+wcols*xdelta
      write(41,3003)':yOrigin            ',yorigin+srows*ydelta
      write(41,3005)'#                                       '
      write(41,3008)':AttributeName 1 Rank         ' 
      write(41,3008)':AttributeName 2 Next         '  
      write(41,3008)':AttributeName 3 DA           '
      write(41,3008)':AttributeName 4 Bankfull     ' 
      write(41,3008)':AttributeName 5 ChnlSlope    ' 
      write(41,3008)':AttributeName 6 Elev         '
      write(41,3008)':AttributeName 7 ChnlLength   ' 
      write(41,3008)':AttributeName 8 IAK          '
      write(41,3008)':AttributeName 9 IntSlope     '
      write(41,3008)':AttributeName 10 Chnl        '
      write(41,3008)':AttributeName 11 Reach       '
      write(41,3008)':AttributeName 12 GridArea    '
!     Added Oct. 2013      
      write(41,3008)':AttributeName 13 FetchNE     '
      write(41,3008)':AttributeName 14 FetchE      '
      write(41,3008)':AttributeName 15 FetchSE     '
      write(41,3008)':AttributeName 16 FetchS      '
      write(41,3008)':AttributeName 17 FetchSW     '
      write(41,3008)':AttributeName 18 FetchW      '
      write(41,3008)':AttributeName 19 FetchNW     '
      write(41,3008)':AttributeName 20 FetchN      '



!  Version 10.6 - Mar. 2012 - create split landcover for nca
!  Version 10.6 - Mar. 2012 - create split landcover for nca
!  Version 10.6 - Mar. 2012 - create split landcover for nca
!  Version 10.6 - Mar. 2012 - create split landcover for nca
c      if(split_no.eq.0)then
      if(nca_classcount.eq.0)then
        do i=1,ntype+1
c          write(41,3009)':AttributeName',12+i,class_name(i)
          write(41,3009)':AttributeName',20+i,class_name(i)
        end do
c      elseif(split_no.eq.1)then
      elseif(nca_classcount.eq.1)then
        write(41,3009)':AttributeName',20+1,class_name(1)
        write(41,3010)':AttributeName',20+2,'nca_',class_name(1)
!       write the rest of the classes        
        do i=2,ntype+1
c          write(41,3009)':AttributeName',13+i,class_name(i)
          write(41,3009)':AttributeName',21+i,class_name(i)
        end do
c      elseif(split_no.eq.2)then
!      fixed bug below:   NK Oct. 15/16 
c      elseif(nca_classcount.eq.3)then
      elseif(nca_classcount.eq.2)then
        write(41,3009)':AttributeName',20+1,class_name(1)
        write(41,3010)':AttributeName',20+2,'nca_',class_name(1)
        write(41,3009)':AttributeName',20+3,class_name(2)
        write(41,3010)':AttributeName',20+4,'nca_',class_name(2)
!       write the rest of the classes        
        do i=3,ntype+1
c          write(41,3009)':AttributeName',14+i,class_name(i)
          write(41,3009)':AttributeName',22+i,class_name(i)
        end do
c      elseif(split_no.eq.3)then
      elseif(nca_classcount.eq.3)then
        write(41,3009)':AttributeName',20+1,class_name(1)
        write(41,3010)':AttributeName',20+2,'nca_',class_name(1)
        write(41,3009)':AttributeName',20+3,class_name(2)
        write(41,3010)':AttributeName',20+4,'nca_',class_name(2)
        write(41,3009)':AttributeName',20+5,class_name(3)
        write(41,3010)':AttributeName',20+6,'nca_',class_name(3)
!       write the rest of the classes        
        do i=4,ntype+1
c          write(41,3009)':AttributeName',15+i,class_name(i)
          write(41,3009)':AttributeName',23+i,class_name(i)
        end do
      else
        print*,'Too many class splits - ony first 3 allowed'
        stop 'Program aborted in BSN @ 2504'   
      endif

! Note<<<<<<<<<<< # split classes does not change ntype!!!

      write(41,3005)'#                                       '
      write(41,3001)':xCount             ',xcount-wcols-ecols
      write(41,3001)':yCount             ',ycount-srows-nrows
      write(41,3003)':xDelta             ',xdelta
      write(41,3003)':yDelta             ',ydelta
      write(41,3005)'#                                       '
      write(41,3005)':EndHeader                              '

      do i=1,ycount
        do j=1,xcount
          dummy(i,j)=0.0
        end do
      end do
      frame=1
c      rank_2d(1,1)=frame

      if(buf.eq.'y')write(41,41000)'# ',frame,'Rank                '
      if(nogrids.le.9999)then
        do i=1+srows,ycount-nrows
          write(41,1058)(rank_2d(i,j),j=1+wcols,xcount-ecols)
        end do
      else
        do i=1+srows,ycount-nrows
          write(41,2058)(rank_2d(i,j),j=1+wcols,xcount-ecols)
        end do
      endif
      print*,' frame= ',frame,' written'

      frame=frame+1 !2
c      next_2d(1,1)=frame
      if(buf.eq.'y')write(41,41000)'# ',frame,'NextGrid            '
      if(nogrids.le.9999)then
        do i=1+srows,ycount-nrows
          write(41,1058)(next_2d(i,j),j=1+wcols,xcount-ecols)
        end do
      else
        do i=1+srows,ycount-nrows
          write(41,2058)(next_2d(i,j),j=1+wcols,xcount-ecols)
        end do
      endif
      print*,' frame= ',frame,' written'

      frame=frame+1 !2
      if(.not.non_ca_flg)then
c       da_2d(1,1)=frame
        if(buf.eq.'y')write(41,41000)'# ',frame,'DrainageArea        '
        do i=1+srows,ycount-nrows
          write(41,4012)(da_2d(i,j),j=1+wcols,xcount-ecols)
        end do
      else
c       da_2d(1,1)=frame
        if(buf.eq.'y')write(41,41000)'# ',frame,'DrainageArea        '
        do i=1+srows,ycount-nrows
          write(41,4012)(da_nca_2d(i,j),j=1+wcols,xcount-ecols)
        end do
        print*,' frame= ',frame,' written'
      endif

      frame=frame+1 !3
c      bnkfll_2d(1,1)=frame
      if(buf.eq.'y')write(41,41000)'# ',frame,'BabnkfullArea       '
      do i=1+srows,ycount-nrows
        write(41,4002)(bnkfll_2d(i,j),j=1+wcols,xcount-ecols)
      end do
      print*,' frame= ',frame,' written'
      
      
      reply='n'
      slopeflg=.false.  !no -ve slopes
      do i=1+srows,ycount-nrows
        do j=1+wcols,xcount-ecols
          if(slope_2d(i,j).le.0.and.frac_2d(i,j).gt.0.00000)then
            slopeflg=.true. !-ve slope found
          endif
        end do
      end do
      if(slopeflg)then
        print*
        print*,'-ve slopes found'
        print*,'You must fix the drainage directions &/or elevations'
        print*,'to fix the problem. However, you can allow bsn.exe to'
        print*,'set these as the min. slope'
        Print*,'Would you like to proceed this way'
        print*,'and accept responsibility  y/n'
        read*,reply
      endif
      if(reply.eq.'y')then
        do i=1+srows,ycount-nrows
          do j=1+wcols,xcount-ecols
            if(slope_2d(i,j).le.0.0.and.frac_2d(i,j).gt.0.00000)then
              slope_2d(i,j)=slopemin
c              print*,'slope fixed at ',i,j 
            endif
          end do
        end do
        print*
        print*,'-ve slopes eliminated'
        print*,'min slope = ',slopemin
        pause 'Hit enter to continue'
      endif
      

      frame=frame+1 !4
c      slope_2d(1,1)=frame
      if(buf.eq.'y')write(41,41000)'# ',frame,'ChnlSlope           '
      do i=1+srows,ycount-nrows
        write(41,4005)(slope_2d(i,j),j=1+wcols,xcount-ecols)
      end do
      print*,' frame= ',frame,' written'

      frame=frame+1 !5
c      elv_2d(1,1)=frame
      if(buf.eq.'y')write(41,41000)'# ',frame,'Elevation           '
      do i=1+srows,ycount-nrows
        write(41,4006)(elv_2d(i,j),j=1+wcols,xcount-ecols)
      end do
      print*,' frame= ',frame,' written'

      frame=frame+1 !6
c      ch_length_2d(1,1)=frame
      if(buf.eq.'y')write(41,41000)'# ',frame,'ChnlLength          '
      do i=1+srows,ycount-nrows
        write(41,4007)(ch_length_2d(i,j),j=1+wcols,xcount-ecols)
      end do
      print*,' frame= ',frame,' written'

      frame=frame+1 !7
c      iak(1,1)=frame
      if(buf.eq.'y')write(41,41000)'# ',frame,'RiverClass          '
      do i=1+srows,ycount-nrows
        write(41,1058)(iak(i,j),j=1+wcols,xcount-ecols)
      end do
      print*,' frame= ',frame,' written'

      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        dummy(i,j)=sl1(n)
      end do
      frame=frame+1 !8
c      dummy(1,1)=frame
      if(buf.eq.'y')write(41,41000)'# ',frame,'IntSlope            '
      do i=1+srows,ycount-nrows
        write(41,4005)(dummy(i,j),j=1+wcols,xcount-ecols)
      end do
      print*,' frame= ',frame,' written'

      frame=frame+1 !9
c      ichnl_2d(1,1)=frame
      if(buf.eq.'y')write(41,41000)'# ',frame,'Chnl                '
      do i=1+srows,ycount-nrows
        write(41,1059)(ichnl_2d(i,j),j=1+wcols,xcount-ecols)
      end do
      print*,' frame= ',frame,' written'   

      frame=frame+1 !10
c      ireach_2d(1,1)=frame
      if(buf.eq.'y')write(41,41000)'# ',frame,'ReachNumber         '
      do i=1+srows,ycount-nrows
        write(41,1058)(ireach_2d(i,j),j=1+wcols,xcount-ecols)
      end do
      print*,' frame= ',frame,' written'

      frame=frame+1 !11
!     write 'gridarea'

      if(.not.non_ca_flg)then
      if(buf.eq.'y')write(41,41000)'# ',frame,'GridArea            '
        do i=1+srows,ycount-nrows
          write(41,4012)(frac_2d(i,j),j=1+wcols,xcount-ecols)
        end do
        print*,' frame= ',frame,' written'
      else
      if(buf.eq.'y')write(41,41000)'# ',frame,'GridArea            '
        do i=1+srows,ycount-nrows
          write(41,4012)(frac_nca_2d(i,j),j=1+wcols,xcount-ecols)
        end do
        print*,' frame= ',frame,' written'
      endif

!     Fetch data added Oct. 2013
      do k=1,8
        do i=1+srows,ycount-nrows
          do j=1+wcols,xcount-ecols
            fetch_2d(i,j,k)=fetch_2d(i,j,k)*al
          end do
        end do
      end do
      do k=1,8
        frame=frame+1 !10
        if(buf.eq.'y')write(41,41001)'# ','fetch_2d               ',k
        do i=1+srows,ycount-nrows
          write(41,1060)(fetch_2d(i,j,k),j=1+wcols,xcount-ecols)
        end do
        print*,' frame= ',frame,' written'
      end do

c      do ii=1,ntype+1 !12
c        frame=frame+1
cc       write(41,41000)'# ',frame,class_name(ii)
c        write(*,41000)'# ',class_name(ii)
c        do i=1+srows,ycount-nrows
c          write(41,4008)(aclass_3d(i,j,ii),j=1+wcols,xcount-ecols)
c        end do
c        print*,' frame= ',frame,' written'
c      end do


!  Version 10.6 - Mar. 2012 - create split landcover for nca
c      if(split_no.eq.0)then
c
c        do ii=1,ntype+1 !12
c          frame=frame+1
cc         write(41,41000)'# ',frame,class_name(ii)
c          write(*,41000)'# ',class_name(ii)
c          do i=1+srows,ycount-nrows
c            write(41,4008)(aclass_3d(i,j,ii),j=1+wcols,xcount-ecols)
c          end do
c          print*,' frame= ',frame,' written'
c        end do

c      if(split_no.ge.1)then
      if(nca_classcount.ge.1)then
      
        print*,wcols,ecols,ycount,xcount
      
!       split the 1st class:    
!       contributing part of the class:    
        ii=1
        frame=frame+1
        if(buf.eq.'y')write(41,41000)'# ',frame,class_name(ii)
c        write(*,41000)'# ',class_name(ii)
        do i=1+srows,ycount-nrows
          print*,i,ii,1+wcols,xcount-ecols
          if(.not.allocated(nca))stop 'nca not allocated'
          write(41,4008)(aclass_3d(i,j,ii)*(1.0-nca(i,j)),
     *          j=1+wcols,xcount-ecols)
        end do
        print*,' frame= ',frame,' written'
          
!       non-contributing part of the class:    
        frame=frame+1
      if(buf.eq.'y')write(41,41002)'# ','nca_',class_name(ii)
        write(*,41002)'# ','nca_',class_name(ii)
        do i=1+srows,ycount-nrows
          write(41,4008)(aclass_3d(i,j,ii)*nca(i,j),
     *              j=1+wcols,xcount-ecols)
        end do
        print*,' frame= ',frame,' written'
      endif

c      if(split_no.ge.2)then
      if(nca_classcount.ge.2)then
!       split the 2nd class:        
!       contributing part of the class:    
        ii=2
        frame=frame+1
        if(buf.eq.'y')write(41,41000)'# ',frame,class_name(ii)
c        write(*,41000)'# ',class_name(ii)
        do i=1+srows,ycount-nrows
          write(41,4008)(aclass_3d(i,j,ii)*(1.0-nca(i,j)),
     *          j=1+wcols,xcount-ecols)
        end do
        print*,' frame= ',frame,' written'
        
!       non-contributing part of the class:    
        frame=frame+1
        if(buf.eq.'y')write(41,41002)'# ','nca_',class_name(ii)
        write(*,41002)'# ','nca_',class_name(ii)
        do i=1+srows,ycount-nrows
          write(41,4008)(aclass_3d(i,j,ii)*nca(i,j),
     *              j=1+wcols,xcount-ecols)
        end do
        print*,' frame= ',frame,' written'
      endif
        
c      if(split_no.eq.3)then
      if(nca_classcount.eq.3)then
!       split the 3nd class:        
!       contributing part of the class:    
        ii=3
        frame=frame+1
        if(buf.eq.'y')write(41,41000)'# ',frame,class_name(ii)
c        write(*,41000)'# ',class_name(ii)
        do i=1+srows,ycount-nrows
          write(41,4008)(aclass_3d(i,j,ii)*(1.0-nca(i,j)),
     *          j=1+wcols,xcount-ecols)
        end do
        print*,' frame= ',frame,' written'
        
!       non-contributing part of the class:    
        frame=frame+1
        if(buf.eq.'y')write(41,41002)'# ','nca_',class_name(ii)
        write(*,41002)'# ','nca_',class_name(ii)
        do i=1+srows,ycount-nrows
          write(41,4008)(aclass_3d(i,j,ii)*nca(i,j),
     *              j=1+wcols,xcount-ecols)
        end do
        print*,' frame= ',frame,' written'
      endif
        
c      if(split_no.gt.0)print*,'end of writing split classes'
      if(nca_classcount.gt.0)print*,'end of writing split classes'
       
c      do ii=1+split_no,ntype+1 !12
      do ii=1+nca_classcount,ntype+1 !12
        frame=frame+1
        if(buf.eq.'y')write(41,41000)'# ',frame,class_name(ii)
c        write(*,41000)'# ',class_name(ii)
        do i=1+srows,ycount-nrows
          write(41,4008)(aclass_3d(i,j,ii),j=1+wcols,xcount-ecols)
        end do
        print*,' frame= ',frame,' written'
      end do

c      do ii=1,ntype+1 !12
c        frame=frame+1
cc       write(41,41000)'# ',frame,class_name(ii)
c        write(*,41000)'# ',class_name(ii)
c        do i=1+srows,ycount-nrows
c          write(41,4008)(aclass_3d(i,j,ii),j=1+wcols,xcount-ecols)
c        end do
c        print*,' frame= ',frame,' written'
c      end do


      close(unit=41,status='keep')

      print*,'new_shd.r2c written'
      print*
      

!  Version 10.5 - Jan. 2012 - create elv_means.r2s
!     write the elv_means.r2c file is dem.r2s avaialble     
      if(answer_dem.eq.'y')then
        author='watflood/bsn                            '
        name='mean grid elevations                    '
        xorigin_temp=xorigin+wcols*xdelta
        yorigin_temp=yorigin+srows*ydelta
        xcount_temp=xcount-wcols-ecols
        ycount_temp=ycount-srows-nrows
        xdelta_temp=xdelta
        ydelta_temp=ydelta
        startdate='unknown   '
        starttime='unknown   '
        unit_conversion=0.0
        attribute_units='masl                                    ' 
!        attribute_type='Runoff                                  '  
        source_file_name='dem.r2s     '     
        no_hdrcomments=0
        do j=1+wcols,xcount-ecols
          do i=1+srows,ycount-nrows
            outarray(i-srows,j-wcols)=elv_dem(i,j)
          end do
        end do
        fln(99)='elv_means.r2c'
!       write the header
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call write_r2c(99,99,0,1,0,0,12)   
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       write the data
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c        call write_r2c(99,99,0,0,1,0,12)   
        call write_r2c(99,99,0,1,0,1,12)   
!       no_frames=2 tricks write_r2c to write frame .....
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print*,'elv_means.r2c written'
      endif

!  Version 10.5 - Jan. 2012 - create elv_max.r2s
!     write the elv_max.r2c file is dem.r2s avaialble     
      if(answer_dem.eq.'y')then
        author='watflood/bsn                            '
        name='max grid elevations                    '
        source_file_name='dem.r2s     '     
        do j=1+wcols,xcount-ecols
          do i=1+srows,ycount-nrows
            outarray(i-srows,j-wcols)=elv_max(i,j)
          end do
        end do
        fln(99)='elv_max.r2c'
!       write the header
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call write_r2c(99,99,0,1,0,0,12)   
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       write the data
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call write_r2c(99,99,0,1,0,1,12)   
!       no_frames=2 tricks write_r2c to write frame .....
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print*,'elv_max.r2c written'
      endif
     
      

      if(parflg.eq.'y')then

!       write the par.r2c file for watroute in r2c format
!       write the par.r2c file for watroute in r2c format
!       write the par.r2c file for watroute in r2c format
!       write the par.r2c file for watroute in r2c format

        open(41,file='new_ch_par.r2c',status='unknown',iostat=ios)
        if(ios.ne.0)then
          print*,' Problems opening the new_par.r2c file '
          print*,' Possible cause: file open or read only'
          print*
          STOP ' Program BSN aborted @ 2358'
        endif

        write(41,3005)'########################################'
        write(41,3005)':FileType r2c  ASCII  EnSim 1.0         '
        write(41,3005)'#                                       '
        write(41,3005)'# DataType               2D Rect Cell   '
        write(41,3005)'#                                       '
        write(41,3005)':Application                    WATFLOOD'
        write(41,3005)':Version                 2.1.23         '
        write(41,3020)':WrittenBy          ',author
        call date_and_time(cday,time)
        write(41,3011)':CreationDate       ',
     *       cday(1:4),cday(5:6),cday(7:8),time(1:2),time(3:4)
        write(41,3005)'#                                       '
        write(41,3005)'#---------------------------------------'
        write(41,3020)':SourceFileName     ',fln(2)
        write(41,50182)al
        write(41,5019)cintv
        write(41,5020)aimpr
        write(41,5021)ntype-ntype  !want classcount to be 0 in bsnm_ch_par.r2c
        write(41,5027)nrvr
        write(41,5022)elvconv
        if(na.eq.nogrids)then
          write(41,5028)na
          write(41,5029)naa
        else
!         for this case there can only be one outlet
          write(41,5028)nogrids
          write(41,5029)nogrids-1
        endif
        write(41,5030)na/2
        write(41,3005)'#                                       '
        write(41,3005)'#                                       '
        write(41,3004)':Projection         ',coordsys
        if(coordsys.eq.'UTM       ')then
            write(41,3004)':Zone               ',zone
        endif
        write(41,3004)':Ellipsoid          ',datum1
        write(41,3005)'#                                       '
        write(41,3003)':xOrigin            ',xorigin+wcols*xdelta
        write(41,3003)':yOrigin            ',yorigin+srows*ydelta
        write(41,3005)'#                                       '
        write(41,3008)':AttributeName 1 Flz         '
        write(41,3008)':AttributeName 2 Pwr         '
        write(41,3008)':AttributeName 3 R1n         '
        write(41,3008)':AttributeName 4 R2n         '
        write(41,3008)':AttributeName 5 mndr        '
        write(41,3008)':AttributeName 6 aa2         '
        write(41,3008)':AttributeName 7 aa3         '
        write(41,3008)':AttributeName 8 aa4         '
        write(41,3008)':AttributeName 9 theta       '
        write(41,3008)':AttributeName 10 widep      '
        write(41,3008)':AttributeName 11 kcond      '
       
       
        write(41,3005)'#          5                             '
        write(41,3001)':xCount             ',xcount-wcols-ecols
        write(41,3001)':yCount             ',ycount-srows-nrows
        write(41,3003)':xDelta             ',xdelta
        write(41,3003)':yDelta             ',ydelta
        write(41,3005)'#                                       '
        write(41,3005)':EndHeader                              '
       
        do i=1,ycount
          do j=1,xcount
            dummy(i,j)=0.0
          end do
        end do
        frame=1

!       note:  iak(i,j)=ibn(n)=basin number
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          dummy(i,j)=flz(n)
        end do
        do i=1+srows,ycount-nrows
          write(41,4009)(dummy(i,j),j=1+wcols,xcount-ecols)
        end do
        print*,' frame= ',frame,' written'
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          dummy(i,j)=pwr(n)
        end do
        do i=1+srows,ycount-nrows
          write(41,4010)(dummy(i,j),j=1+wcols,xcount-ecols)
        end do
        frame=frame+1
        print*,' frame= ',frame,' written'
       
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          dummy(i,j)=R1n(n)
        end do
        do i=1+srows,ycount-nrows
          write(41,4010)(dummy(i,j),j=1+wcols,xcount-ecols)
        end do
        frame=frame+1
        print*,' frame= ',frame,' written'
       
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          dummy(i,j)=R2n(n)
        end do
        do i=1+srows,ycount-nrows
          write(41,4010)(dummy(i,j),j=1+wcols,xcount-ecols)
        end do
        frame=frame+1
        print*,' frame= ',frame,' written'
       
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          dummy(i,j)=mndr(n)
        end do
        do i=1+srows,ycount-nrows
          write(41,4011)(dummy(i,j),j=1+wcols,xcount-ecols)
        end do
        frame=frame+1
        print*,' frame= ',frame,' written'
       
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          dummy(i,j)=aa2(n)
        end do
        do i=1+srows,ycount-nrows
          write(41,4004)(dummy(i,j),j=1+wcols,xcount-ecols)
        end do
        frame=frame+1
        print*,' frame= ',frame,' written'
       
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          dummy(i,j)=aa3(n)
        end do
        do i=1+srows,ycount-nrows
          write(41,4004)(dummy(i,j),j=1+wcols,xcount-ecols)
        end do
        frame=frame+1
        print*,' frame= ',frame,' written'
       
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          dummy(i,j)=aa4(n)
        end do
        do i=1+srows,ycount-nrows
          write(41,4004)(dummy(i,j),j=1+wcols,xcount-ecols)
        end do
        frame=frame+1
        print*,' frame= ',frame,' written'
       
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          dummy(i,j)=theta(n)
        end do
        do i=1+srows,ycount-nrows
          write(41,4010)(dummy(i,j),j=1+wcols,xcount-ecols)
        end do
        frame=frame+1
        print*,' frame= ',frame,' written'
       
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          dummy(i,j)=widep(n)
        end do
        do i=1+srows,ycount-nrows
          write(41,4001)(dummy(i,j),j=1+wcols,xcount-ecols)
        end do
        frame=frame+1
        print*,' frame= ',frame,' written'
       
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          dummy(i,j)=kcond(n)
        end do
        do i=1+srows,ycount-nrows
          write(41,4010)(dummy(i,j),j=1+wcols,xcount-ecols)
        end do
        frame=frame+1
        print*,' frame= ',frame,' written'
        print*
        close(unit=41,status='keep')
        print*,'new_ch_par.r2c written'
        print*
       
!       end of writing the par file for watroute in r2c format
      endif

!     Write the wfo_spec.new file
!     Write the wfo_spec.new file
!     Write the wfo_spec.new file
!     Write the wfo_spec.new file
!     Write the wfo_spec.new file
!     Write the wfo_spec.new file

!     This section is set up so new items can be easily added
!     by the auto-numbering sequence.
!     Don't forget to change the attribute count

      open(unit=99,file='wfo_spec.new',status='unknown',iostat=ios)
      if(ios.ne.0)then
        print*,' unable to open a new wfo_spec.new file'
        print*,' program not aborted but file not written'
        print*
      endif
!     rev. 9.1.50  Jan.  14/04  - NK: version number added to the wfo_spec.txt file
      write(99,99000)       
!     Attribute Count
!  Version 10.6 - Mar. 2012 - create split landcover for nca
c      write(99,99001)13+15*(ntype+1+split_no)
      write(99,99001)17+16*(ntype+1+nca_classcount)        !    change <<<<<<<<<<<<<<<<<!!
!     Reporting time step - can be changed by user in wfo_spec.txt
      write(99,99002)            
!     Start Reporting       Added Jan. 14/04 NK
      write(99,99902)
!     End Reporting       Added Jan. 14/04 NK
      write(99,99903)
!     1 Temperature 
      j=1                   !  j is the attribute number
      write(99,99201)j
!     2 Precipitation
      j=j+1
      write(99,99202)j
!     3 CummulativePrecipitation    Added June 11/03 NK
      j=j+1
      write(99,99203)j
!     4 lower zone storage
      j=j+1
      write(99,99204)j
!     5 ground water discharge
      j=j+1
      write(99,99205)j
!     6 grid runoff
      j=j+1
      write(99,99206)j
!     7 observed outflow      
      j=j+1
      write(99,99207)j
!     8 computed outflow      
      j=j+1
      write(99,99208)j
!     9 weighted swe
      j=j+1
      write(99,99209)j
!     10 wetland depth
      j=j+1
      write(99,99210)j
!     channel depth
      j=j+1
      write(99,99211)j
!     wetland storage
      j=j+1
      write(99,99212)j
!     wetland outflow
      j=j+1
      write(99,99213)j
!     cummulative ET
      j=j+1
      write(99,99214)j
!     Bankfull
      j=j+1
      write(99,99215)j
!     totuzs
      j=j+1
      write(99,99216)j
!     API
      j=j+1
      write(99,99217)j

      
      !     depression storage
c      do i=1,(ntype+1+split_no)

      print*,'NCA class count',nca_classcount
      pause 'waiting'

      do i=1,(ntype+1+nca_classcount)
        j=j+1
        write(99,99014)j,i
      end do
!     depression storage (snow)
c      do i=1,(ntype+1+split_no)
      do i=1,(ntype+1+nca_classcount)
        j=j+1
        write(99,99015)j,i
      end do
! snow water equivalent
c      do i=1,(ntype+1+split_no)
      do i=1,(ntype+1+nca_classcount)
        j=j+1
        write(99,99016)j,i
      end do
!     snow covered area
      do i=1,(ntype+1+nca_classcount)
c      do i=1,(ntype+1+split_no)
        j=j+1
        write(99,99017)j,i
      end do
!     upper zone storage
c      do i=1,(ntype+1+split_no)
c      do i=1,(ntype+1+nca_classcount)
      do i=1,(ntype+nca_classcount)
        j=j+1
        write(99,99018)j,i
      end do
!     rev. 10.4.21 Apr.  21/20  = NK Add UZS deficit to wfo file = UZS(class=classcount)
      j=j+1
      write(99,99030)j,i
!     upper zone storage (snow)
c      do i=1,(ntype+1+split_no)
      do i=1,(ntype+1+nca_classcount)
        j=j+1
        write(99,99019)j,i
      end do
!     surface flow from bare area
c      do i=1,(ntype+1+split_no)
      do i=1,(ntype+1+nca_classcount)
        j=j+1
        write(99,99020)j,i
      end do
!     surface flow from snow covered ground
c      do i=1,(ntype+1+split_no)
      do i=1,(ntype+1+nca_classcount)
        j=j+1
        write(99,99021)j,i
      end do
!     interflow from bare ground
c      do i=1,(ntype+1+split_no)
      do i=1,(ntype+1+nca_classcount)
        j=j+1
        write(99,99022)j,i
      end do
!     interflow from snow covered ground
c      do i=1,(ntype+1+split_no)
      do i=1,(ntype+1+nca_classcount)
        j=j+1
        write(99,99023)j,i
      end do
!     drainage from bare ground
c      do i=1,(ntype+1+split_no)
      do i=1,(ntype+1+nca_classcount)
        j=j+1
        write(99,99024)j,i
      end do
!     drainage from snow covered ground
c      do i=1,(ntype+1+split_no)
      do i=1,(ntype+1+nca_classcount)
        j=j+1
        write(99,99025)j,i
      end do
!     PET 
c      do i=1,(ntype+1+split_no)
      do i=1,(ntype+1+nca_classcount)
        j=j+1
        write(99,99026)j,i
      end do
!     ET 
c      do i=1,(ntype+1+split_no)
      do i=1,(ntype+1+nca_classcount)
        j=j+1
        write(99,99027)j,i
      end do
!     Sublimation from snow 
c      do i=1,(ntype+1+split_no)
      do i=1,(ntype+1+nca_classcount)
        j=j+1
        write(99,99028)j,i
      end do
!     Heat deficit snow 
c      do i=1,(ntype+1+split_no)
      do i=1,(ntype+1+nca_classcount)
        j=j+1
        write(99,99029)j,i
      end do
    
      close(unit=99)

      print*,'wfo_spec.new  written'

!     write the new.pdl file
!     write the new.pdl file
!     write the new.pdl file

      open(unit=33,file='new.pdl',recl=flen,status='unknown')
      if(ios.ne.0)then
        print*,' problems opening file new.pdl   ios=',ios
        stop 'program aborted in inpgrd.for @ 44'
      endif

      nnprint=0
      newformat=0
!     rev. 9.1.58  Jul.  12/04  - NK: New header for the .shd file
      write(33,5012)
      write(33,33001)':FileType             bsnm.pdl'
      write(33,5009)coordsys
      write(33,5010)datum1
      write(33,5011)zone
      write(33,5012)

      if(coordsys.eq.'LATLONG   ')then
        write(33,5013)xorigin+wcols*xdelta
        write(33,5014)yorigin+srows*ydelta
      else
        write(33,50131)xorigin+wcols*xdelta
        write(33,50141)yorigin+srows*ydelta
      endif
      write(33,5012)
      write(33,5015)xcount-wcols-ecols
      write(33,5016)ycount-srows-nrows
      if(coordsys.eq.'LATLONG   ')then
        write(33,5017)xdelta
        write(33,5018)ydelta
      else
        write(33,50171)xdelta
        write(33,50181)ydelta
      endif

      ewgrid=xorigin+float(wcols)*xdelta+
     *           (float(xcount-wcols-ecols)/2.)*xdelta
      sngrid=yorigin+float(srows)*ydelta+
     *           (float(ycount-srows-nrows)/2.)*ydelta

      write(33,5012)
      write(33,33000)':NoPrecipStations     1'
      write(33,5012)
      if(llflg.ne.'y')then
!       FOR UTM 
          write(33,*)ewgrid,sngrid,'  centerville'
      else
!       FOR LAT-LONG
          write(33,*)ewgrid,sngrid,'  centerville'
      endif

      write(33,5012)
      write(33,33000)':NoSnowCourses        1'
      write(33,5012)
      if(llflg.ne.'y')then
!       FOR UTM 
          write(33,*)ewgrid,sngrid,'  centerville'
      else
!       FOR LAT-LONG
          write(33,*)ewgrid,sngrid,'  centervill'
      endif

      write(33,5012)
      write(33,33000)':NoTempStations       1'
      write(33,5012)
      if(llflg.ne.'y')then
!       FOR UTM 
          write(33,*)ewgrid,sngrid,'  centerville'
      else
!       FOR LAT-LONG
        write(33,*)ewgrid,sngrid,'  centerville'
      endif

      write(33,5012)
      write(33,33000)':NoFlowStations       1'
      write(33,5012)
      if(llflg.ne.'y')then
!       FOR UTM 
        write(33,33002)ewgrid,sngrid,'  centerville',
     *	' 0.000E+00 0.000E+00 0.000E+00 0.000E+00  0'
33002   format(2f10.0,a13,a43)
      else
!       FOR LAT-LONG
        write(33,33003)ewgrid,sngrid,'  centerville',
     *	' 0.000E+00 0.000E+00 0.000E+00 0.000E+00  0'
33003   format(2f10.4,a13,a43)
      endif

      write(33,5012)
      write(33,33000)':NoReservoirs         1'
      write(33,5012)
      if(llflg.ne.'y')then
!       FOR UTM 
        write(33,33004)ewgrid,sngrid,'  centerville',
     *' 0.000E+00 0.000E+00 0.000E+00 0.000E+00  0.000E+00'
33004   format(2f10.0,a13,a51)
      else
!       FOR LAT-LONG
        write(33,33005)ewgrid,sngrid,'  centerville',
     *' 0.000E+00 0.000E+00 0.000E+00 0.000E+00  0.000E+00'
33005   format(2f10.4,a13,a51)
      endif

      write(33,5012)
      write(33,33000)':NoDamageSites        1'
      write(33,5012)

33000  format(a23,i5)
33001  format(a30)

      close(unit=33, status='keep')
      print*,'new.pdl  written'
      print*

      if(split.gt.0.00)then
        print*
        print*,'The value entered for "split" =',int(split*100.)
        print*,'Another land cover class was be added to the'
        print*,'new_shd.r2c file'
        Print*,'Please ensure that the par & sdc files are modified'
        print*,'to accomodate the extra wetland class'
        print*
      endif

!     pause 10

!     Added May 14, 2016  NK
!     check for grids with frac = 0 that are NOT outlets
      errflg_l=.false.
c      do n=1,na
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        if(frac_2d(i,j).eq.0.0.and.next_2d(i,j).gt.0)then
          print*,'ERROR:'
          print*,'A grid with frac = 0.0 and with a downstream grid'
          print*,'with frac > 0 is found in row ',i,' and col ',j
          print*,'in grid # ',n
          print*,'Next grid row',i,'col',j,'Nex',next_2d(i,j)
          errflg_l=.true.
        endif
      end do
      if(errflg_l)then
        print*
        print*,'ERROR: with fatal consequences!'
        print*,'Please correct by either setting frac to a proper'
        print*,'value if this(these) grid(s) is(are) in a watershed'
        print*,'OR '
        print*,'set its(their) elevation to 0.0 if not in the watershed'
        print*
        print*,'OR'
        print*,'Receiving grids of 2 or more outlets are not '
        print*,'   at the same - lowest - elevation'
        print*,'All receiving grids must be at the same elv.'
        print*
        stop 'Program aborted in BSN @ 3836'
      endif
      
c	if(na-naa.gt.1)then
         print*,'If you have gotten this far, you probably will have'
        print*,'a good shd file - i.e. there will be a shd file'
        print*
        print*,'The rest of the program tends to work only if you have'
        print*,'a single watershed outlet'
        print*
c        print*,'You can try & continue - hit "enter"'
c        print*
c	  pause 'In bsn @ 3067'
        print*
c	endif

!     only for the mrb for flow 1D reach locations
!     sept. 5/09 nk.
      maxr=0
      do i=ycount,1,-1
          do j=1,xcount
            maxr=(max0(maxr,ireach_2d(i,j)))
        end do
      end do
d     print*,'no of reaches found maxr=',maxr

      open(unit=999,file='1D_reach_locations.xyz',status='unknown')

d     print*,'no_outlets,na,naa',no_outlets,na,naa

c      if(no_outlets.ne.na-naa)then
      if(no_outlets.gt.1)then
        print*,'There is more than one outlet'
        print*,'so river profiles are not computed'
        print*
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
        print*,'*******************************************************'
        print*,'*                                                     *'
        print*,'*                  WATFLOOD (TM)                      *'
        print*,'*                                                     *'
        print*,'*     Program BSN Version 10.10     May 11, 2017      *'
        print*,'*                                                     *'
        print*,'*           (c) N. Kouwen, 1972-2017                  *'
        print*,'*                                                     *'
        print*,'*******************************************************'
        print*
        print*,'Please see bsn_info.txt for information re: this run'
	  stop 'Normal ending @ 3807'
      endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do l=1,maxr
        do n=naa,1,-1
!         work upstream
!         look for next reach
          i=yyy(n)
          j=xxx(n)
          if(next_2d(i,j).ge.1.and.next_2d(i,j).ge.1)then
            ii=yyy(next_2d(i,j))
            jj=xxx(next_2d(i,j))
!           if the next downstream reach is not the same, write
!           coordinates for the reach outlet.
            if(ireach_2d(i,j).eq.l)then
              if(ireach_2d(ii,jj).ne.ireach_2d(i,j))
     *	        write(999,*)j*xdelta+xorigin-xdelta/2.0,
     * 		         	i*ydelta+yorigin-ydelta/2.0,ireach_2d(i,j)
            endif
          endif
        end do
      end do

d     do l=1,no_outlets
d       print*,'l,lastgrid,na,naa  ',l,lastgrid(l),na,naa
d     end do


c      if(na-naa.le.1000)then
c      if(lastgrid(1).eq.naa)then
      if(no_outlets.eq.1)then
!     PROFILES
!     this section only works for a single outlet
      filenames(1)='profil01.dat'
      filenames(2)='profil02.dat'
      filenames(3)='profil03.dat'
      filenames(4)='profil04.dat'
      filenames(5)='profil05.dat'
      filenames(6)='profil06.dat'
      filenames(7)='profil07.dat'
      filenames(8)='profil08.dat'
      filenames(9)='profil09.dat'
      filenames(10)='profil10.dat'
      filenames(11)='river01.xyz'
      filenames(12)='river02.xyz'
      filenames(13)='river03.xyz'
      filenames(14)='river04.xyz'
      filenames(15)='river05.xyz'
      filenames(16)='river06.xyz'
      filenames(17)='river07.xyz'
      filenames(18)='river08.xyz'
      filenames(19)='river09.xyz'
      filenames(20)='river10.xyz'

cd     pause 11

!     compute the chainage of each grid to the outlet
      do n=1,na
        chainage(n)=0.0
        maxlength=0.0
        longest=0
        channel(n)=0
      end do

d     pause 12

!     work up the watershed and compute the chainage to outlet for each grid
      do n=naa,1,-1
        i=yyy(n)
        j=xxx(n)
!       for grids outside the sub-watershed, rank = 0
        if(rank_2d(i,j).gt.0)then  !new
        
          if(abs(newnext(n)).gt.0)then       ! new Jan. 20/15 NK
            chainage(n)=chainage(abs(newnext(n)))+al
          endif
          
!         find the grid with the longest path
          if(chainage(n).gt.maxlength)then
            maxlength=chainage(n)
            longest=n          
          endif
        endif                     !new
      end do

d     pause 13

!     mark off the longest channel
      n=longest
c      do while(n.le.naa)
      do while(n.le.naa.and.n.gt.0)  ! changed Jan. 20/15  NK
        channel(n)=1
        i=yyy(n)
        j=xxx(n)
        n=abs(next_2d(i,j))
      end do

d     pause 14

      open(unit=98,file=filenames(1),status='unknown',iostat=ios)
      open(unit=99,file=filenames(11),status='unknown',iostat=ios)
      if(ios.ne.0)then
        stop 'Program aborted @ 992, can`t open chainage.dat file'
      endif
!     write stuff for the main stem
      do n=1,na
        if(channel(n).eq.1)then
          i=yyy(n)
          j=xxx(n)
          nn=nn+1
          write(98,*)xorigin+j*xdelta-xdelta/2.0,
     *	             yorigin+i*ydelta-ydelta/2.0,
     *		         channel(n),chainage(n),elv_2d(i,j)
          write(99,*)xorigin+j*xdelta-xdelta/2.0,
     *	             yorigin+i*ydelta-ydelta/2.0,
     *		         channel(n),chainage(n),elv_2d(i,j)
!          print*,n,chainage(n),elv_2d(i,j),channel(n)
        endif
      end do
      close(unit=98,status='keep')
      close(unit=99,status='keep')
        print*,'finished writing ',filenames(1)
        print*,'finished writing ',filenames(11)

d      pause 15

!     write stuff for the tributaries
      do k=2,10
!     check the length of each tributary and find longest unmarked:
        longest=0
        maxlength=0.0
        do n=1,naa
          tempsum=al
          m=n
          do while(channel(m).eq.0)
            tempsum=tempsum+al
            i=yyy(m)
            j=xxx(m)
            m=abs(next_2d(i,j))
            if(m.eq.0)then
              print*,'There seems to be a problem in'
              print*,'grid no',rank_2d(i,j),' row',i,' column ',j
              print*,'Check for -ve slope or missing data'
              print*,'in this grid - OK if a receiving grid.'
              print*,'OR'
              print*,'More than one grid drains into the very last grid'
              print*,'Which means there is more than one outlet.'
              print*,'This can be a problem. Edit the map file so only'
              print*,'one grid drains into each receiving outlet grid.'
              print*,'OR'
              print*,'Did you leave blank columns & rows around the'
              print*,'edges of the data? Maybe not'
              print*
              print*,'River profile data not printed'
              print*,'Please check new_format.shd for problems'
              print*
              print*,'new_shd.r2c and new_format.shd have been written'
              print*
              stop 'Program ended in bsn @ 2953'
            endif
          end do
          if(tempsum.gt.maxlength)then
            longest=n
            maxlength=tempsum
          endif          
        end do
!       mark the channel until a previously marked channel is found:
        n=longest
        do while(channel(n).eq.0)
          channel(n)=k
          i=yyy(n)
          j=xxx(n)
          n=abs(next_2d(i,j))  ! still done on the next, not newnext????? 21/11/02
        end do

        open(unit=98,file=filenames(k),status='unknown',iostat=ios)
        open(unit=99,file=filenames(k+10),status='unknown',iostat=ios)
        if(ios.ne.0)then
          stop 'Program aborted @ 992, can`t open chainage.dat file'
        endif
        do n=1,na    !new
          if(rank_2d(i,j).gt.0)then
          if(channel(n).eq.k)then
            i=yyy(n)
            j=xxx(n)
            write(51,*)xorigin+j*xdelta-xdelta/2.0,
     *	             yorigin+i*ydelta-ydelta/2.0,
     *		         channel(n),chainage(n),elv_2d(i,j)
            write(99,*)xorigin+j*xdelta-xdelta/2.0,
     *	             yorigin+i*ydelta-ydelta/2.0,
     *		         channel(n),chainage(n),elv_2d(i,j)
!            print*,n,chainage(n),elv_2d(i,j),channel(n)
            lastn=n
          endif
          endif    !new
        end do

!I DON'T KNOW WHAT THE NEXT FEW LINES DO  NK nOV. 21/02
!        n=abs(newnext(i,j))   !changed
!        i=yyy(n)
!        j=xxx(n)
!        if(rank_2d(i,j).gt.0)then   !new
!        write(51,*)n,chainage(n),elv_2d(i,j),channel(n)
!        write(99,*)j,i,elv_2d(i,j)
! 
!       print*,n,chainage(n),elv_2d(i,j),channel(n)
!        endif      !new



        close(unit=98)
        close(unit=99)
        print*,'finished writing ',filenames(k)
        print*,'finished writing ',filenames(k+10)

      end do        

      endif     ! na-naa.le.1




!     pause 16
      print*
      print*,'No. of errors found in the map file = ',no_errors
      print*,'No. of errors found in the map file = ',no_errors
      print*,'No. of errors found in the map file = ',no_errors
      print*
      if(no_errors.ne.0)then
       print*,' ********* please check the bsn_info.txt file **********'
       print*,' ********* please check the bsn_info.txt file **********'
       print*,' ********* please check the bsn_info.txt file **********'
       print*
      else
        print*,'new_shd.r2c has been written'
        print*,'Please rename new_shd.r2c or replace the bsnm_shd.r2c '
        print*
      endif
c      STOP ' Normal ending'

! TS - ADDED DEALLOCATION STATEMENT FOR AREA17 ARRAYS
c      deallocate(elv_2d,da,slope,frac,aclass,xxx,yyy,rank,s,next,ielv,iak,
c     *irough,ichnl,bnkfll,ireach,dummy,stat=iDeallocateStatus)
c      if (iDeallocateStatus.ne.0) STOP   
c     *    '**Deallocation for AREA17A arrays in bsna failed**' 

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
      print*,'********************************************************'
      print*,'*                                                      *'
      print*,'*                  WATFLOOD (TM)                       *'
      print*,'*                                                      *'
      print*,'*     Program BSN Version 10.10     Aug. 29, 2023      *'
      print*,'*                                                      *'
      print*,'*           (c) N. Kouwen, 1972-2020                   *'
	print*,'*                                                      *'
	print*,'*             Open source code under the               *'
      print*,'*          GNU Lesser General Public License           *'
      print*,'*                                                      *'
      print*,'********************************************************'
      print*
      print*,'Please see file bsn_info.txt for information re: this run'


! FORMATS

 1001 format(999f4.0)
 1002 format(10f8.4)
 1003 format(999i2)
 1004 format(999f6.1)
 1005 format(12i5,2f5.1)
 1006 format(999i4)
 1007 format(12i5,2f5.1)
 1008 format(20f5.3/)
 1009 format(2i5,4f10.5)
 1010 format(3i5,f10.0,f10.2,f10.7,f7.1,5i5,99f5.2)
 1011 format(3i5,f10.3,f10.5,f10.7,f7.1,5i5,99f5.2)
 1012 format(' ','istep is set to 1 - local coordinates are used')
 1013 format('  al cintvl impr ntype elvconv')
 1014 format('   Basin # not coded @ grid #',i6,' @ ',2i6,' elv=',f7.3)
 1015 format('# contours not coded @ grid #',i6,' @ ',2i6,' elv=',f7.3)
 1016 format('# channels not coded @ grid #',i6,' @ ',2i6,' elv=',f7.3)
 1017 format('       next grid = 0 @ grid #',i6,' @ ',2i6,' elv=',f7.3)
 1018 format('Possible cause: wrong drainage direction')
 1019 format('Errors OK if last receiving grid !!!!!!!!!!!!!'/)
10019 format('Please see new_format.shd file for -ve slope location')
 1020 format(f10.0,  f5.0,2i5,2f10.0)
 1021 format(f10.0,  f5.1,2i5,f5.3)
 1022 format(f10.0,  f5.1,2i5,2f10.3)
 1024 format(4i5,f10.3)
 1025 format(4i5,f10.0,f10.2,f10.7,f10.3,f7.0,i5,f10.5,2i5,3x,99f5.2)
 1026 format(4i5,f10.3,f10.5,f10.7,f10.3,f7.0,i5,f10.5,2i5,3x,99f5.2)
 1027 format(4i5,f10.0,f10.2,f10.7,f7.1,f7.0,i5,11e10.3,
     *               f10.5,2i5,3x,17f5.2)
 1028 format(4i5,f10.3,f10.5,f10.7,f7.1,f7.0,i5,11e10.3,
     *               f10.5,2i5,3x,17f5.2)
 1050 format(999f4.2)
 1051 format(999f5.2)
 1052 format(999f4.0)
 1053 format(999f6.0)
 1054 format(2i4,999f5.0)
 1055 format(2i4,999f5.2)
 1056 format(999f6.2)
 1057 format(2i4,999i5)
 1058 format(999i5)
 2058 format(999i6)
 1059 format(999i4)
 1060 format(999f10.0)
 1101 format(/'   ls   ks   js   ih  local')
 1102 format(/'   l  istep cintv local')
 1103 format(/' iymin iymax jxmin jxmax')
 1104 format(/' converted to local coordinates')
 1105 format( ' 1, ycount,iymin,iymax, yint')
 1106 format(/' 1, xcount,jxmin,jxmax, xint')
 1107 format(/' elevations of the channel bottoms')
 1108 format( ' half way along the square')
 1109 format(/' fraction of the elemement within the basin')
 1110 format(/' drainage directions')
 1111 format( ' 1=ne,2=e,3=se,4=s,5=sw,6=w,7=nw,8=n')
 1112 format(/' ib & it - the bottom and top rows for computations')
 1113 format( '           in local coordinates')
 1116 format(/' na - the total number of elements within the basin')
 1117 format( '      including the element accepting the last flow')
 1118 format(/' pairs of local coordinates for each element in order')
 1119 format( ' of descending elevations - ie according to rank')
 1120 format(/' rank - the order in which the calculations are carried')
 1121 format( '        out, from highest to lowest element')
 1122 format(/' drainage area above the outlet of each element')
 1123 format(/' slope of the river from this to the next square')
 1124 format(/' permeability & land use class')
 1125 format(/' number of contours in the element')
 1126 format(/' number of channels traversing the element')
 1127 format(/' % impervious area in each element')
 1128 format(/' the next element - ie the element recieving the flow') 
 1129 format(/' frac_2d(',2i5,')=',f10.3,' - please check')

 3000 format(a10,i5)
 3001 format(a20,i12)
 3002 format(2a20)
 3003 format(a20,f15.5) !changed 15.6 > 15.5 Aug/10 nk
 3004 format(a20,a10,2x,a10)
 3005 format(a40)
 3006 format(a3,a10)
 3007 format(a14,i5,a6,i5)
 3008 format(a30)
 3009 format(a14,i3,1x,a20)
 3010 format(a14,i3,1x,a4,a20)
 3012 format(a9)
 3020 format(a20,a40)

 4001 format(999(' ',f5.0))
 4002 format(999(' ',f10.3))
 4003 format(999(' ',i2))
 4004 format(999(' ',f10.5))
 4005 format(999(' ',f10.7))
 4006 format(999(' ',f8.1))
 4007 format(999(' ',f8.0))
 4008 format(999(' ',f8.3))
 4009 format(999(' ',e10.3))
 4010 format(999(' ',f5.3))
 4011 format(999(' ',f5.1))
 4012 format(999(' ',e12.7))

 4050 format(999(i4))
 4052 format(999(' ',i3))

 5000 format(a80)
 5001 format(a2)
 5002 format(a12,i3,a65)
 5003 format(a20,f12.0)
 5004 format(a20,a10)
 5005 format(a1,a80)
 5006 format(a20,i12)
5007  format(a1,a79)
5009  format(':CoordSys           ',a10)
5010  format(':datum1              ',a10)
5011  format(':Zone               ',a10)
5012  format('#',a80)
5013  format(':xOrigin            ',f12.5) !changed 12.7 > 15.5 Aug/10 nk
5014  format(':yOrigin            ',f12.5) !changed 12.7 > 15.5 Aug/10 nk
50131 format(':xOrigin            ',f12.3)
50141 format(':yOrigin            ',f12.3)
5015  format(':xCount             ',i12)
5016  format(':yCount             ',i12)
5017  format(':xDelta             ',f12.5) !changed 12.7 > 15.5 Aug/10 nk
5018  format(':yDelta             ',f12.5) !changed 12.7 > 15.5 Aug/10 nk
50171 format(':xDelta             ',f12.3)
50181 format(':yDelta             ',f12.3)
50182 format(':NominalGridSize_AL ',f12.3)
5019  format(':ContourInterval    ',f12.3)
5020  format(':ImperviousArea     ',f12.3)
5021  format(':ClassCount         ',i12)
5022  format(':ElevConversion     ',f12.3)
5023  format('#          ')
5024  format('#            ',2a30)
5025  format(':InputFileName      ',a30)
5026  format(':Created     :      ',
     *        2(a2,':'),a2,2x,a2,'-',a2,'-',a4)
5027  format(':NumRiverClasses    ',i12)
5028  format(':TotalNumOfGrids    ',i12)
5029  format(':numGridsInBasin    ',i12)
5030  format(':DebugGridNo        ',i12)
5032  format(':endHeader          ')
5033  format('#                   ')

 6000 format(5x,'assembled basin data file')
 6001 format(' fatal error - outlet elements have not been')
 6002 format(' properly designated')
 6003 format(' set frac of the recieving element = 0'//)
6004  format(' location of outlet:')
6005  format(' grid number  ',i5)
6006  format(' row number   ',i5,f10.3)
6007  format(' column number',i5,f10.3)
 6011 format(' area ratios do not add to 1.00 @ i=',i2,' j=',i2)
 6012 format(' bankfull discharges:')
 6013 format(' zero bankfull discharge encountered')
 6014 format(' probable cause: blank lines at the end of the map file')
 6015 format(' receiving grid higher in grid',3i5,' elv=',f10.3)
 7001 format(12i5)
 7005 format(4i5)
 7006 format(2i5,'             # rows & columns resp.')
 7021 format(f10.0,f5.1,2i5,'     map file used =  ',a30)
41000 format(a1,1x,i5,1x,a20)
41002 format(a1,1x,a4,a20)
41001 format(a1,1x,a20,i5)
99000 format('  7.0 Version Number')   ! added a class
99001 format(i5,' AttributeCount')
99002 format('   24 ReportingTimeStep Hours')
99902 format('    0 Start Reporting Time for GreenKenue (hr)')
99903 format('    0 End Reporting Time for GreenKenue (hr)')
99201 format('1',i4,' Temperature')
99202 format('0',i4,' Precipitation')
99203 format('1',i4,' Cumulative Precipitation ')
99204 format('0',i4,' Lower Zone Storage Class ')  
99205 format('0',i4,' Ground Water Discharge m^3/s ',i2)     
99206 format('0',i4,' Grid Runoff')
99207 format('1',i4,' Observed Outflow')
99208 format('1',i4,' Computed Outflow')
99209 format('1',i4,' Weighted SWE')
99210 format('0',i4,' Wetland Depth')
99211 format('0',i4,' Channel Depth')
99212 format('0',i4,' Wetland Storage in m^3')
99213 format('0',i4,' Wetland Outflow in m^3/s')
99214 format('1',i4,' Cumulative ET')
99215 format('1',i4,' Bankfull')
99216 format('1',i4,' totUZS')
99217 format('1',i4,' API')
      
      
      
99014 format('0',i4,' Depression Storage Class ',i2)   
99015 format('0',i4,' Depression Storage (Snow) Class ',i2)
99016 format('0',i4,' Snow Water Equivalent Class ',i2)
99017 format('0',i4,' Snow Covered Area Class ',i2)
99018 format('0',i4,' Upper Zone Storage Class ',i2)
99030 format('0',i4,' Upper Zone Storage Deficit ',i2)
99019 format('0',i4,' Upper Zone Storage (Snow) Class ',i2)
99020 format('0',i4,' Surface Flow m^3/s Class ',i2)
99021 format('0',i4,' Surface Flow (snow) m^3/s Class ',i2)
99022 format('0',i4,' Interflow m^3/s Class ',i2)
99023 format('0',i4,' Interflow (snow) m^3/s Class',i2)
99024 format('0',i4,' Recharge mm Class ',i2)
99025 format('0',i4,' Recharge mm (snow) Class ',i2)
99026 format('0',i4,' PET (average) mm Class ',i2)
99027 format('0',i4,' ET (cummulative) mm Class ',i2)
99028 format('0',i4,' Sublimation Cummulative) mm (snow) Class ',i2)
99029 format('0',i4,' Heat deficit mm (snow) Class ',i2)
!99030 - see above      
99100 format(a26,i5)
99101 format(a30,a40)
99102 format(a23,<no_outlets+1>i8)
99103 format(a23,<no_inlets+1>i8)
99105 format(a30,f12.7)
99106 format(a30,a1)

      
      END PROGRAM bsn
