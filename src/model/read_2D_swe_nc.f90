! This is part of the netCDF package.
! Copyright 2006 University Corporation for Atmospheric Research/Unidata.
! See COPYRIGHT file for conditions of use.

! This subroutine is based on an example which reads some 4D pressure and
! temperatures. The data file read by this program is produced by
! the companion program pres_temp_4D_wr.f90. It is intended to
! illustrate the use of the netCDF Fortran 90 API.

! This program is part of the netCDF tutorial:
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

! Full documentation of the netCDF Fortran 90 API can be found at:
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90

! $Id: pres_temp_4D_rd.f90,v 1.9 2010/04/06 19:32:09 ed Exp $



  subroutine read_swe_2D_nc(m,fn,ensemble)
  
!     rev. 10.3.08 Mar.  06/20  = NK For FEWS, add snow1\swe.nc for swe updating

!  Argument from CHARM
!  m = sequence # 1, 2, 3, 4,.....,nrecs
!  nrecs = -ve = close file
!  fn = file name number
!  ensemble = index of ensemble member (ignored if there is no ensemble)
!  TEMP_NAME = variabel name e.g. flow, precipitation, temperature
!  TEMP_UNITS = units e.g. cms, mm

  use netcdf
  use area_debug
  use area_watflood
  implicit none
  save

  ! This is the name of the data file we will read.
  character (256) :: FILE_NAME
  integer :: ncid,iAllocate

! Also, I'd like to get rid of NLVLS altogether
  integer, parameter :: NDIMSMAX = 10  !<AM: ideally maximum number possible

  integer :: NDIMS, NRECS
  integer :: NLVLS, NLATS, NLONS


  integer, parameter             :: NOENSEMBLE = -999
  character (len = *), parameter :: ENSEMBLE_NAME = "realization" !< AM: note the conflict with levels ...
  character (len = *), parameter :: LVL_NAME = "level"
  character (len = *), parameter :: row_NAME = "y"
  character (len = *), parameter :: col_NAME = "x"
  character (len = *), parameter :: LAT_NAME = "y"
  character (len = *), parameter :: LON_NAME = "x"
  character (len = *), parameter               :: REC_NAME       = "time"
  character (len = *), parameter               :: TIME_NAME      = "time"
  character (len = *), parameter               :: TIME_STEP_NAME = "time_step"
  character (len = *), parameter               :: TIME_ID_NAME   = "time"
! We will read surface temperature and pressure fields etc. In netCDF
! terminology these are called "variables."
  character (len = *), parameter :: TEMP_NAME = "swe"
  character (len = *), parameter :: TEMP_UNITS = "mm"

  ! The start and count arrays will tell the netCDF library where to
  ! read our data.
  integer :: dimids(NDIMSMAX)
  integer, allocatable :: start(:), count(:)

  ! In addition to the latitude and longitude dimensions, we will also
  ! create latitude and longitude variables which will hold the actual
  ! latitudes and longitudes. Since they hold data about the
  ! coordinate system, the netCDF term for these is: "coordinate
  ! variables."
  real, allocatable :: lats(:), lons(:)
  integer, allocatable :: times(:)                              !, lons(:)
  integer :: lon_varid, lat_varid, rec_varid, lvl_varid
  integer :: lon_dimid, lat_dimid, rec_dimid, lvl_dimid
  integer :: temp_varid
  integer :: time_varid

  ! We recommend that each variable carry a "units" attribute.
  character (len = *), parameter :: UNITS = "units"
!  character (len = *), parameter :: PRES_UNITS = "hPa", TEMP_UNITS = "celsius"
!  character (len = *), parameter :: TEMP_UNITS = "celsius"
  character (len = *), parameter :: LAT_UNITS = "degrees_north"
  character (len = *), parameter :: LON_UNITS = "degrees_east"

  ! Program variables to hold the data we will read in. We will only
  ! need enough space to hold one timestep of data; one record.
  ! Allocate memory for data.
!  real, dimension(:,:,:), allocatable :: pres_in
  real, dimension(:,:,:), allocatable :: temp_in
  real, dimension(:,:), allocatable :: ratio2
  real(4)  :: class_sum

  ! Loop indices
  integer :: lvl, i , j,ios, step, fn, m
  integer :: n,ii
  

  integer           :: ensemble
  integer           :: ensemble_index
!  integer, save     :: number_members            !ask AM: do we need this if done for general save???
  integer           :: number_members             !ask AM: do we need this if done for general save???
  integer           :: ensemble_dimid

! Added for Dt  NK
  real   :: time_step, time
  integer :: time_id, time_step_id,start_time
  integer :: isecnds, iatime(9), ierror

  logical :: firstpass,done

  data firstpass/.true./done/.false./

  if(done)return

  if(firstpass)then
      NLVLS=1                   !<AM: this is actually the dimension of ensemble members ...
      FILE_NAME = fln(fn)

      if(debug_output)write(63,*)'*******************************'
      if(debug_output)write(63,*)'*******************************'
      if(debug_output)write(63,*)'*******************************'
      if(debug_output)write(63,*)'In read_swe_2D_nc'
      if(debug_output)write(63,*)m,fn,TEMP_NAME,' ',TEMP_UNITS
      if(debug_output)write(63,*)fln(fn)(1:40)

!     AM: we need to open the file and read the various dimensions
!         before we can enter the loop ...
!     Open the file.
      if(debug_output)write(63,*)FILE_NAME(1:20), nf90_nowrite, ncid
      call check( nf90_open(FILE_NAME, nf90_nowrite, ncid) )
      if(debug_output)write(63,*)'Opened ',FILE_NAME(1:40)
      write(98,*)'Info: Opened ',FILE_NAME(1:40)
!
!     Get the dimensions
      call check( nf90_inq_varid(ncid, LAT_NAME, lat_varid) )
      call check( nf90_inq_varid(ncid, LON_NAME, lon_varid) )
      if(debug_output)write(63,*)'step 200'

!!!!!!!!!!  call check( nf90_inq_dimid(ncid, LVL_NAME, lvl_dimid) ) !<AM: not a separate variable, only a dimension
      call check( nf90_inq_dimid(ncid, REC_NAME, rec_dimid) ) !<AM: ditto


!     Determine the number of ensemble members - if any
      ierror = nf90_inq_dimid(ncid, ENSEMBLE_NAME, ensemble_dimid) !<AM: ditto
      if ( ierror == nf90_noerr ) then
          ensemble_index = ensemble
      else
          if ( ierror == nf90_ebaddim ) then
              ensemble_index = NOENSEMBLE
          else
              call check( nf90_inq_dimid(ncid, ENSEMBLE_NAME, ensemble_dimid) ) ! Proper error message the lazy way
          endif
      endif

      if(debug_output)write(63,*)ncid, REC_NAME, rec_dimid
      if(debug_output)write(63,*)'step 300'

      call check( nf90_inquire_variable( ncid, lat_varid, ndims=ndims, dimids=dimids ) )
      lat_dimid = dimids(1)
      if(debug_output)write(63,*)'ndims\lat_varid =',ndims,lat_varid

      call check( nf90_inquire_variable( ncid, lon_varid, ndims=ndims, dimids=dimids ) )
      lon_dimid = dimids(1)
      if(debug_output)write(63,*)'ndims\lon_varid =',ndims,lon_varid
      if(debug_output)write(63,*)'step 400'

      call check( nf90_inquire_dimension( ncid, lat_dimid, len=nlats ) )
      if(debug_output)write(63,*)'ycount',nlats,ycount
      call check( nf90_inquire_dimension( ncid, lon_dimid, len=nlons ) )
      if(debug_output)write(63,*)'xcount',nlons,xcount
      call check( nf90_inquire_dimension( ncid, rec_dimid, len=nrecs ) ) !<AM: this is an "unlimited" dimension - we merely get the current size
!!!!!!!!!  call check( nf90_inquire_dimension( ncid, lvl_dimid, len=nlvls ) )
      if(debug_output)write(63,*) 'Lat/lon/levels: ', nlats, nlons, nlvls
      if(debug_output)write(63,*) 'Records:        ', nrecs
      if(debug_output)write(63,*)'step 500'

      if ( ensemble_index /= NOENSEMBLE ) then
          call check( nf90_inquire_dimension( ncid, ensemble_dimid, len=number_members ) )
      else
          number_members = NOENSEMBLE
      endif
!
!     Allocate the arrays we need
      allocate( temp_in(NLONS, NLATS, NLVLS) )   ! for arrays
      allocate(lats(nlats),lons(nlons))
      if(debug_output)write(63,*)'step 600'
!
!     We make it fully flexible, but we know that in the netCDF file the variable
!     temperature is a 4D array. As we use it per time step, the dimension in the
!     program needs to be 3 only.
!
      call check( nf90_inq_varid( ncid, TEMP_NAME, temp_varid) )
      if(debug_output)write(63,*)'TEMP_NAME =',TEMP_NAME
      call check( nf90_inquire_variable( ncid, temp_varid, ndims=ndims ) )
      allocate( count(ndims), start(ndims) )
      if(debug_output)write(63,*)'step 700'

!     insert code for finding the array dimensions
      ! Read the latitude and longitude data.
      call check( nf90_get_var(ncid, lat_varid, lats) )
      call check( nf90_get_var(ncid, lon_varid, lons) )
      if(debug_output)write(63,*)'step 1100'

!      Get the varids of the temperature netCDF variables.
       call check( nf90_inq_varid(ncid, TEMP_NAME, temp_varid) )
       if(debug_output)write(63,*)'step 1200'
!
       start = (/ 1, 1, 1, 1 /)
       if(debug_output)write(63,*)'start',start
       if(debug_output)write(63,*)'nrecs',nrecs
       if(debug_output)write(63,*)'step 1300'

       if(debug_output)write(63,*)'ncid, time_id, time, start',ncid, time_id, time, start
       call check( nf90_inq_varid(ncid, TIME_id_NAME, time_id) )
       if(debug_output)write(63,*)'TIME_id_NAME :',TIME_id_NAME
       if(debug_output)write(63,*)'time_id',time_id
      !  get time:
       call check( nf90_get_var(ncid, time_id, time) )
       
       if(debug_output)write(63,*)'time: ',time
       if(debug_output)write(63,*)'step 1400'

       allocate(times(nrecs))
       do I = 1, nrecs
            call check( nf90_get_var(ncid, time_id, times ,start = start ) )
       enddo
       if(debug_output)write(63,*)(times(i),i=1,nrecs)
       if(debug_output)write(63,*)'step 1500'

       deltaT2=(times(2)-times(1))/60
       epoch_min_start=times(1)-360

       write(98,*)'Info: Start time info in precip.nc:'
       write(98,*)'Info: deltaT_swe =',deltaT2,' hours'
       if(debug_output)write(63,*)'deltaT2 =',deltaT2
! time of the first frame in seconds from Jan. 1/1970 -0600
! FEWS time is at the end of 6hr deltaT
       isecnds=time*60                          !-6*60     ! time of the first frame in seconds from Jan. 1/1970 -0600
!      Get the date & time
      !******************************************
       call PXFLOCALTIME (isecnds, iatime, ierror)
      !******************************************
       if(debug_output)write(63,*)isecnds,iatime,ierror
!     Convert epoch time to yyyymmdd hhmmss 
!     isecnds - input integer with # of seconds since jan 1, 1970 
!     iatime - output array returning the following: 
!        1 - seconds (0 - 61, for leap seconds) 
!        2 - minutes (0 - 59) 
!        3 - hours (0 - 23) 
!        4 - day of the month (1 - 31) 
!        5 - month of the year (1 - 12) 
!        6 - Gregorian year (e.g., 2006) 
!        7 - Day of the week (0 = sunday) 
!        8 - Day of the year (1 - 366) 
!        9 - DST flag (0 = standard, 1 = DST) 
!     ierror - Returns 0 if successful, EINVAL if not.       
       if(debug_output)write(63,6000)iatime(6),iatime(5),iatime(4),iatime(3),iatime(2),iatime(8)
       write(98,6000)iatime(6),iatime(5),iatime(4),iatime(3),iatime(2),iatime(8)
       write(*,6000)iatime(6),iatime(5),iatime(4),iatime(3),iatime(2),iatime(8)
6000   format(' Info: In read_2D_swe_nc: first frame time stamp ',i4,'/',i2.2,'/',i2.2,'   ',i2.2,':',i2.2,'  Julian day ',i3)
       
      ! Read 1 record of NLONS*NLATS*NLVLS values, starting at the beginning
      ! of the record (the (1, 1, 1, rec) element in the netCDF file).
!       count(1:4) = (/ NLONS, NLATS, NLVLS, 1 /)            ! for arrays
        count = (/ NLONS, NLATS, NLVLS, 1 /)            ! for arrays
        start(1:4) = (/ 1, 1, 1, 1 /)

!       Only exact grid matching in FEWS!!!
!        iyoffset=0
!        jxoffset=0
!        xcount2=nlons
!        ycount2=nlats
!        nhtemp=nrecs*deltaT2  
      allocate(ratio2(ycount,xcount),stat=iAllocate)
      if(iAllocate.ne.0)then
         print*,'Error with allocation of ratio2 in read_swe_2D_nc.f @ 63'
         STOP 'Program aborted  @ 595'
      endif
        
        firstpass=.false.
        write(63,*)'End first pass read_2D_swe_nc- initialization completed'

      if(debug_output)write(63,*)'step 1600 - end first pass'

  end if  !    end firstpass    end firstpass    end firstpass    end firstpass

  if(done)return
  if(debug_output)write(63,*)'step 1700'

!  Get the time of the swe update frame:  
       isecnds=time*60+m*deltaT2*3600                          
!      Get the date & time
      !******************************************
       call PXFLOCALTIME (isecnds, iatime, ierror)
      !******************************************
       if(debug_output)write(63,*)'swe update time',isecnds,iatime,ierror
       if(debug_output)write(63,6000)iatime(6),iatime(5),iatime(4),iatime(3),iatime(2),iatime(8)
       write(98,6001)iatime(6),iatime(5),iatime(4),iatime(3),iatime(2),iatime(8)
       write(*,6001)iatime(6),iatime(5),iatime(4),iatime(3),iatime(2),iatime(8)
6001   format(' Info: In read_2D_swe_nc:       frame time stamp ',i4,'/',i2.2,'/',i2.2,'   ',i2.2,':',i2.2,'  Julian day ',i3)
  
  if ( ensemble_index > number_members .and. number_members /= NOENSEMBLE ) then
    write(98,*)'Error: Ensemble member does not exist - number of members is ', number_members
    write(98,*)'Error: program aborted in read_2D_swe_nc @ 278' 
    stop 'Error: program aborted in read_2D_swe_nc @ 278' 
  endif

! If there is an ensemble, the member index is the third dimension and the time
! the fourth. If there is NO ensemble, the time is the third.
  if ( ensemble_index /= NOENSEMBLE ) then
      start(3) = ensemble_index
      start(4) = m
  else
      start(3) = m      ! m = frame number
      start(4) = 1
  endif
  if(debug_output)write(63,*)'step 1800'
  
  
  
! Read the swe data from the file, one record at a time.
! Read the swe data from the file, one record at a time.
! Read the swe data from the file, one record at a time.
! Read the swe data from the file, one record at a time.
! Read the swe data from the file, one record at a time.
!  if(mod(k+deltaT2-1,deltaT2).eq.0)then
!   *********************************************************************  
    call check( nf90_get_var(ncid, temp_varid, temp_in, start = start, &
                                        count = count) )
!   *********************************************************************  
    if(debug_output)write(63,99001)temp_in(jpr,ipr,1),m
99001 format(' Info:  SW corner swe ',f8.2,i10)

    
    do i=1,ycount   ! allows a good check
            write(400,*)(temp_in(j,i,1),j=1,5)
    end do
  if(debug_output)write(63,*)'step 2000'
  
  
      
!   Put into local array
    do i = 1, NLATS
        do j=1,NLONS
             snw(i,j)=temp_in(j, i, 1)
        end do
    end do
    if(debug_output)write(63,*)'frame ',m,'  nrecs = ',nrecs,'swe ',snw(ipr,jpr)
    if(debug_output)write(63,*)'*** SUCCESS reading file ', FILE_NAME(1:20),'@ frame ',m
    write(98,*)'Info: swe from nc file in grid #',nnprint,' = ',snw(ipr,jpr)
!  endif
  if(debug_output)write(63,*)'step 2100'

!     rev. 10.3.07 Mar.  04/20  = NK Fixed weighted swe in wfo file for grids with water
!     rev. 10.3.08 Mar.  06/20  = NK For FEWS, add snow1\swe.nc for swe updating
!     Start by finding the weighted swe for each grid   
  do n=1,naa
       totsnw(n)=0.0
       class_sum=0.0
       do ii=1,classcount
          if(aclass(n,ii).gt.0.0)then
             if(ii.ne.ii_water)then
               class_sum=class_sum+aclass(n,ii)
               totsnw(n)=totsnw(n)+snowc(n,ii)*aclass(n,ii)*sca(n,ii)
             endif
          endif
       end do
       totsnw(n)=totsnw(n)/class_sum
  end do     
    write(98,*)'Info: model se in grid #',nnprint,' = ',totsnw(nnprint)
  if(debug_output)write(63,*)'step 1900'
    
!     rev. 10.3.09 Mar.  07/20  = NK Revise swe updating to maintain relative swe in classes
!     Find the ratio of obs swe / weighted model swe
  do n=1,naa
        i=yyy(n)
        j=xxx(n)
!          if(snw(i,j).gt.1.0)then
          if(totsnw(n).gt.a2)then
!     rev. 10.4.23 Apr.  30/20  = NK repurposed a2 for swe threshold"
!              If there is not much snow, ratio is too erratic   
!              and just use the snow course data                  
               ratio2(i,j)=snw(i,j)/totsnw(n)
          else
               ratio2(i,j)=1.0000
          endif
   end do
    write(98,*)'Info: ratio fews/model in grid #',nnprint,' = ',ratio2(ipr,jpr)
  if(debug_output)write(63,*)'step 2200'
      
!   Adjust and put INTO VECTOR FORMAT
  do ii=1,classcount
      do n=1,naa
         if(aclass(n,ii).gt.0.0)then
            i=yyy(n)
            j=xxx(n)
!             Old way: swe in all classes replaced by swe from the r2c file              
!			snowc(n,ii)=max(0.0,snw(i,j))
            if(debug_output.and.n.eq.nnprint)print*,n,ii,snowc(n,ii),aclass(n,ii),ratio2(i,j)              
!     rev. 10.3.09 Mar.  07/20  = NK Revise swe updating to maintain relative swe in classes
!             Multiply by the ratio so relative swe in classes stay the same
             if(snowc(n,ii).gt.25)then
                 snowc(n,ii)=snowc(n,ii)*ratio2(i,j)
             else
                 snowc(n,ii)=snw(i,j)
             endif
            if(snowc(n,ii).gt.0.0)then
                sca(n,ii)=1.0
                oldsca(n,ii)=1.0
!				set the initial heat deficit for the snow:
!				def(n,ii)=deffactor*snowc(n,ii)               ! <<<< needs to be defined
            else
!     May 10, 2002 Added these lines to avoid underflows later! AB
                sca(n,ii)=0.0
                oldsca(n,ii)=0.0
                def(n,ii)=0.0
            endif
            if(debug_output.and.n.eq.nnprint)print*,n,ii,snowc(n,ii)              
            if(debug_output.and.n.eq.nnprint)print*
      if(n.eq.nnprint)write(98,*)'Info: adjusted swe in grid - class #',nnprint,ii,' = ',snowc(n,ii)
         endif
      end do
  end do

  if(debug_output)write(63,*)'step 2300'
  if(debug_output)write(63,*)'******************************************'
  if(debug_output)write(63,*)'******************************************'
  if(debug_output)write(63,*)'******************************************'
  if(debug_output)write(63,*)'******************************************'
    
  if(m.eq.NRECS)then
!       Close the file.
        call check( nf90_close(ncid) )
        write(98,*)'Info: Closed ',fln(fn)(1:40)
        done=.true.
        return                           ! as part of WF    <<<<<<<<<<<<<<<<<<<<<<<<
  endif
  return

contains
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      write(98,*)'Error: ', trim(nf90_strerror(status))
      write(98,*)'Error: Progran aborted in read_swe_2D_nc'
      write(98,*)'Error: with ',fln(fn)
      write(*,*)'Error: ', trim(nf90_strerror(status))
      write(*,*)'Error: Progran aborted in read_swe_2D_nc'
      write(*,*)'Error: with ',fln(fn)
      stop 
    end if
  end subroutine check

end subroutine read_swe_2D_nc

