      MODULE EF_Module

!***********************************************************************
!    Copyright (C) 2006 by Dave Watson (NRC) and Nicholas Kouwen 
      
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
      
      
!     rev. 9.5.58  Apr.  16/09  - NK: added nudgeflg for forcing gauge flows
!     rev. 9.8.12  Dec.  08/11  - NK: recognize kenueflg in the event file 
!     rev. 9.8.22  Jul.  17/12  - NK: Added resetflg to reset cumm. precip Sept.1
!     rev. 9.8.47  Feb.  04/13  - NK: Headers added for spl & resin csv files
!     rev. 9.8.93  Nov.  12/13  - NK: Added the routing initialization with yyyymmdd_fli.r2c
!     rev. 9.9.00  Dec.  08/13  - NK: Added Lake Evaporation model
      
      
      USE EF_ParseUtilities
        
      TYPE Attribute
            CHARACTER*64 name
            CHARACTER*16 units
            INTEGER index, count
            REAL, DIMENSION(:), ALLOCATABLE :: val
      END TYPE Attribute
      
      TYPE EnSimParam
            TYPE(Attribute) :: attList(1000) !changed 200 -> 1000 nk Mar 14/11
            INTEGER attCount
            INTEGER  fileFormat
            CHARACTER*8 fileType
            CHARACTER*32 application, writtenBy, creationDate, dataType
            CHARACTER*16 version
            CHARACTER*64 name
      END TYPE EnSimParam
      
      TYPE CoordSysParam
            CHARACTER*12 projection, ellipsoid
            INTEGER zone 
            CHARACTER*8 units
      END TYPE CoordSysParam

      TYPE R2CParam
            TYPE(EnSimParam) ep
            TYPE(CoordSysParam) csp
            TYPE(Attribute) attp
            INTEGER  xCount, yCount
            integer  deltaT
            REAL xDelta, yDelta, xOrigin, yOrigin
            REAL angle, unitConv
            INTEGER frameCount
            CHARACTER*32 surveyDate, surveyTime
!     added for reading cumm precip by -nk- Jan. 20/08
            CHARACTER*64 Name
            CHARACTER*32 attributeUnits
      END TYPE R2CParam

      TYPE ShedParam
            TYPE(R2CParam) r2cp
            REAL nominalGridSize_AL, contourInterval, elevConversion
            INTEGER classCount, numRiverClasses, totalNumOfGrids
            INTEGER numGridsInBasin, debugGridNo
            CHARACTER*128 sourceFileName
      END TYPE

!     added for read_event  by -nk- Jan. 20/08
!     don't forget - these are also assigned in area_watflood as stand alone variables
      TYPE EventParam
        REAL smc_init,conv,scale,readscale
        REAL scalesnw,scaletem,readscalesnw
        REAL nhg,nhf,readscaletem,disaggregate
        INTEGER year1,mo1,day1,hour1
        INTEGER deltat_report,nch
        INTEGER spinupevents
        CHARACTER*1 snwflg,sedflg,vapflg
        CHARaCTER*1 smrflg,resinflg,tbcflg
        CHARACTER*1 resumflg,contflg,routeflg
        CHARACTER*1 crseflg,ensimflg,kenueflg,picflg,wetflg
        CHARACTER*1 modelflg,shdflg,trcflg,hdrflg0
        CHARACTER*1 frcflg,initflg,grdflg,ntrlflg
        CHARACTER*1 nudgeflg,resetflg,flowresetflg
        CHARACTER*1 lakeflg,iceflg,icerivflg,icelakeflg
        CHARACTER*1 divertflg,pafflg,fliflg,tb0flg
        CHARACTER*1 xmlflg,nbsflg,fewsflg,netcdfoutflg
        CHARACTER*1 writestr
        CHARACTER*60 infln(100)
      END TYPE

      TYPE RainParam
            TYPE(R2CParam) r2cp
            INTEGER startJulianDay, endJulianDay, startHour, endHour
      END TYPE

      TYPE TempParam
            TYPE(R2CParam) r2cp
            INTEGER startJulianDay, endJulianDay, startHour, endHour
      END TYPE
      
      TYPE GSMParam
            TYPE(R2CParam) r2cp
      END TYPE

      TYPE SWEParam
            TYPE(R2CParam) r2cp
            REAL initHeatDeficit
      END TYPE

      TYPE TimeStamp
            INTEGER year,month,day,hour,min,sec,msec
      END TYPE

      Type FrameRecord
            TYPE(TimeStamp) tStamp
            INTEGER frame, step
            INTEGER hour, hoursToNext
            CHARACTER*5 source
      END TYPE

      TYPE TB0Param
            TYPE(EnSimParam) ep
            TYPE(CoordSysParam) csp
            REAL unitConv
            INTEGER deltaT
            CHARACTER*30 sourceKey     ! added Mar. 21/14 NK for ragmet & tmp
      END TYPE TB0Param

      TYPE FlowParam
            TYPE(TB0Param) tb0p
            INTEGER routingDeltaT
            CHARACTER*2 fillFlag
      END TYPE FlowParam

!     Added for read_tb0  nk  Dec. 2013
      TYPE temp2Param
            TYPE(TB0Param) tb0p
      END TYPE temp2Param
     
      
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
!     Added for ragmet  nk  Jul. 2/08
      TYPE RagParam
            TYPE(TB0Param) tb0p
      END TYPE RagParam
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

      TYPE ResvParam
            TYPE(TB0Param) tb0p
      END TYPE ResvParam

!     Added for lake levels Aug. 02/11
      TYPE LvlParam
            TYPE(TB0Param) tb0p
      END TYPE LvlParam

!     Added for swe report June 18/12 nk
      TYPE CourseParam
            TYPE(TB0Param) tb0p
      END TYPE CourseParam

!     Added for diversions  nk  Jan. 20/09
      TYPE DivParam
            TYPE(TB0Param) tb0p
      END TYPE DivParam

      TYPE ResvinParam
            TYPE(TB0Param) tb0p
      END TYPE ResvinParam

      
      TYPE TB0ColumnMetaData
            INTEGER colCount
            CHARACTER(64), DIMENSION(:), ALLOCATABLE :: colName
            CHARACTER(16), DIMENSION(:), ALLOCATABLE :: colUnits
            CHARACTER(8), DIMENSION(:), ALLOCATABLE :: colType
            REAL, DIMENSION(:), ALLOCATABLE :: colLocX
            REAL, DIMENSION(:), ALLOCATABLE :: colLocY
            REAL, DIMENSION(:), ALLOCATABLE :: colLocX1
            REAL, DIMENSION(:), ALLOCATABLE :: colLocY1
c            REAL, DIMENSION(:), ALLOCATABLE :: colValue1
            REAL, DIMENSION(:), ALLOCATABLE :: elevation
       END TYPE TB0ColumnMetaData

      TYPE FlowColumnMetaData
            TYPE(TB0ColumnMetaData) tb0cmd
            REAL, DIMENSION(:), ALLOCATABLE :: colCoeff1
            REAL, DIMENSION(:), ALLOCATABLE :: colCoeff2
            REAL, DIMENSION(:), ALLOCATABLE :: colCoeff3
            REAL, DIMENSION(:), ALLOCATABLE :: colCoeff4
            REAL, DIMENSION(:), ALLOCATABLE :: colValue1
      END TYPE FlowColumnMetaData

      TYPE Temp2ColumnMetaData
            TYPE(TB0ColumnMetaData) tb0cmd
       END TYPE Temp2ColumnMetaData

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
!     Added for ragmet  nk  Jul. 2/08
      TYPE RagColumnMetaData
            TYPE(TB0ColumnMetaData) tb0cmd
            REAL, DIMENSION(:), ALLOCATABLE :: colValue1
      END TYPE RagColumnMetaData
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

!     rev. 10.1.46 Nov.  08/16  - NK: Changed B1 - 5 to real*8
      TYPE ResvColumnMetaData
            TYPE(TB0ColumnMetaData) tb0cmd
            REAL*8, DIMENSION(:), ALLOCATABLE :: colCoeff1
            REAL*8, DIMENSION(:), ALLOCATABLE :: colCoeff2
            REAL*8, DIMENSION(:), ALLOCATABLE :: colCoeff3
            REAL*8, DIMENSION(:), ALLOCATABLE :: colCoeff4
            REAL*8, DIMENSION(:), ALLOCATABLE :: colCoeff5
            REAL, DIMENSION(:), ALLOCATABLE :: colCoeff6
            REAL, DIMENSION(:), ALLOCATABLE :: colCoeff7
      END TYPE ResvColumnMetaData
      
!     Added for lake levels  nk  Aug. 02/11
      TYPE LvlColumnMetaData
            TYPE(TB0ColumnMetaData) tb0cmd
            REAL, DIMENSION(:), ALLOCATABLE :: colValue1
      END TYPE LvlColumnMetaData

!     Added for diversions  nk  Jan. 20/09
      TYPE DivColumnMetaData
            TYPE(TB0ColumnMetaData) tb0cmd
            REAL, DIMENSION(:), ALLOCATABLE :: colValue1
      END TYPE DivColumnMetaData

      TYPE ResvinColumnMetaData
            TYPE(TB0ColumnMetaData) tb0cmd
            REAL, DIMENSION(:), ALLOCATABLE :: colValue1
      END TYPE ResvinColumnMetaData




!     Added for swe report June 18/12 nk
      TYPE CourseColumnMetaData
            TYPE(TB0ColumnMetaData) tb0cmd
      END TYPE CourseColumnMetaData





C     TYPE (TB0Column) columns(10)
C     TYPE(FlowColumn), DIMENSION(:), ALLOCATABLE :: flowCols

      CONTAINS  
      
      SUBROUTINE InitEnSimParam(ep)
            TYPE(EnSimParam), INTENT(INOUT) :: ep
            ep%fileType = "Unknown"
            ep%dataType = "Unknown"
            ep%application = "Unknown"
            ep%version = "Unknown"
            ep%writtenBy = "Unknown"
            ep%creationDate = "Unknown"
            ep%name = "Unknown"     
            ep%attCount = 0         
      END SUBROUTINE InitEnSimParam

      SUBROUTINE InitCoordSysParam(csp)
            TYPE(CoordSysParam), INTENT(INOUT) :: csp
            csp%projection = "Cartesian"
            csp%ellipsoid = "Unknown"
c            csp%  = 0               !wicked bug found Mar. 21/15 NK
            csp%zone  = 0
      END SUBROUTINE InitCoordSysParam

      SUBROUTINE InitAttribute(attp)
            TYPE(Attribute), INTENT(INOUT) :: attp
            attp%name = "Unknown"
            attp%units = "Unknown"
            attp%index = 0
            attp%count = 0
      END SUBROUTINE InitAttribute

      SUBROUTINE InitR2CParam(r2cp)
            TYPE(R2CParam), INTENT(INOUT) :: r2cp
            CALL InitEnSimParam(r2cp%ep)
            CALL InitCoordSysParam(r2cp%csp)
            r2cp%xCount = 0
            r2cp%yCount = 0
            r2cp%xDelta = 0.0
            r2cp%yDelta= 0.0
            r2cp%angle = 0
            r2cp%xOrigin = 0.0
            r2cp%yOrigin = 0.0
            r2cp%frameCount = 0
            r2cp%unitConv = 1.0
            r2cp%deltaT = 0
            r2cp%surveyDate = "Unknown"
            r2cp%surveyTime = "Unknown"
      END SUBROUTINE InitR2CParam

      SUBROUTINE InitShedParam(shdp)
            TYPE(ShedParam), INTENT(INOUT) :: shdp
            CALL InitR2CParam(shdp%r2cp)
            shdp%nominalGridSize_AL = 0.0
            shdp%contourInterval = 0.0
            shdp%classCount = 0
            shdp%numRiverClasses = 0
            shdp%elevConversion = 1.0
            shdp%totalNumOfGrids = 0
            shdp%numGridsInBasin = 0
            shdp%debugGridNo = 0
      END SUBROUTINE InitShedParam

!     added for read_event  by -nk- Jan. 20/08
      SUBROUTINE InitEventParam(event)
            TYPE(EventParam), INTENT(INOUT) :: event
c           CALL InitR2CParam(evtp%evtp)
            event%year1 = 0
            event%mo1 = 0
            event%hour1 = 0
            event%day1 = 0
            event%deltat_report=1
            event%spinupevents=0
!     fix  initialize flags etc.
      END SUBROUTINE InitEventParam

      SUBROUTINE InitRainParam(rainp)
            TYPE(RainParam), INTENT(INOUT) :: rainp
            CALL InitR2CParam(rainp%r2cp)
            rainp%startJulianDay = 0
            rainp%startHour = 0
      END SUBROUTINE InitRainParam

      SUBROUTINE InitTempParam(tempp)
            TYPE(TempParam), INTENT(INOUT) :: tempp
            CALL InitR2CParam(tempp%r2cp)
      END SUBROUTINE InitTempParam

      SUBROUTINE InitGSMParam(gsmp)
            TYPE(GSMParam), INTENT(INOUT) :: gsmp
            CALL InitR2CParam(gsmp%r2cp)
      END SUBROUTINE InitGSMParam

      SUBROUTINE InitSWEParam(swep)
            TYPE(SWEParam), INTENT(INOUT) :: swep
            CALL InitR2CParam(swep%r2cp)
            swep%initHeatDeficit = 0.0
      END SUBROUTINE InitSWEParam

      SUBROUTINE InitTimeStamp(tStamp)
            TYPE(TimeStamp), INTENT(INOUT) :: tStamp
            tStamp%year = 0
            tStamp%month = 0
            tStamp%day = 0
            tStamp%hour = 0
            tStamp%min = 0
            tStamp%sec = 0
            tStamp%msec = 0
      END SUBROUTINE InitTimeStamp

      SUBROUTINE InitFrameRecord(frameRec)
            TYPE(FrameRecord), INTENT(INOUT) :: frameRec
            CALL InitTimeStamp(frameRec%tStamp)
            frameRec%hour = 0
            frameRec%hoursToNext = 0
            frameRec%source = "     "
      END SUBROUTINE InitFrameRecord

      SUBROUTINE InitTB0Param(tb0p)
            TYPE(TB0Param), INTENT(INOUT) :: tb0p
            CALL InitEnSimParam(tb0p%ep)
            CALL InitCoordSysParam(tb0p%csp)
            tb0p%deltaT = 0
            tb0p%unitConv = 0
      END SUBROUTINE InitTB0Param

      SUBROUTINE InitFlowParam(flowp)
            TYPE(FlowParam), INTENT(INOUT) :: flowp
            CALL InitTB0Param(flowp%tb0p)
            flowp%routingDeltaT = 0
            flowp%fillFlag = "n"
      END SUBROUTINE InitFlowParam

!     added Dec. 2013 nk
      SUBROUTINE InitTemp2Param(temp2p)
            TYPE(Temp2Param), INTENT(INOUT) :: temp2p
            CALL InitTB0Param(temp2p%tb0p)
      END SUBROUTINE InitTemp2Param





!     Added for swe report June 18/12 nk
      SUBROUTINE InitCourseParam(coursep)
            TYPE(CourseParam), INTENT(INOUT) :: coursep
            CALL InitTB0Param(Coursep%tb0p)
      END SUBROUTINE InitCourseParam





!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
!     Added for ragmet  nk  Jul. 2/08
      SUBROUTINE InitRagParam(ragp)
            TYPE(RagParam), INTENT(INOUT) :: ragp
            CALL InitTB0Param(ragp%tb0p)
      END SUBROUTINE InitRagParam
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

      SUBROUTINE InitResvParam(resvp)
            TYPE(ResvParam), INTENT(INOUT) :: resvp
            CALL InitTB0Param(resvp%tb0p)
      END SUBROUTINE InitResvParam

!     Added for lake levels  nk  Aug. 02/11
      SUBROUTINE InitLvlParam(lvlp)
            TYPE(LvlParam), INTENT(INOUT) :: lvlp
            CALL InitTB0Param(lvlp%tb0p)
      END SUBROUTINE InitLvlParam

!     Added for diversions  nk  Jan. 20/09
      SUBROUTINE InitDivParam(divp)
            TYPE(DivParam), INTENT(INOUT) :: divp
            CALL InitTB0Param(divp%tb0p)
      END SUBROUTINE InitDivParam

      SUBROUTINE InitResvinParam(resvinp)
            TYPE(ResvinParam), INTENT(INOUT) :: resvinp
            CALL InitTB0Param(resvinp%tb0p)
      END SUBROUTINE InitResvinParam

      SUBROUTINE InitTB0ColumnMetaData(tb0Cols)
            TYPE(TB0ColumnMetaData), INTENT(INOUT) :: tb0Cols
            tb0Cols%colCount = 0
      END SUBROUTINE InitTB0ColumnMetaData
      
      SUBROUTINE InitFlowColumnMetaData(flowcmd)
            TYPE(FlowColumnMetaData), INTENT(INOUT) :: flowcmd
            CALL InitTB0ColumnMetaData(flowcmd%tb0cmd)
      END SUBROUTINE InitFlowColumnMetaData


!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
!     Added for ragmet  nk  Jul. 2/08
      SUBROUTINE InitRagColumnMetaData(ragcmd)
            TYPE(RagColumnMetaData), INTENT(INOUT) :: ragcmd
            CALL InitTB0ColumnMetaData(ragcmd%tb0cmd)
      END SUBROUTINE InitRagColumnMetaData
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

      SUBROUTINE InitResvColumnMetaData(resvcmd)
        TYPE(ResvColumnMetaData), INTENT(INOUT) :: resvcmd
            CALL InitTB0ColumnMetaData(resvcmd%tb0cmd)
      END SUBROUTINE InitResvColumnMetaData

!     Added for lake levels  nk  Aug. 02/11
      SUBROUTINE InitLvlColumnMetaData(lvlcmd)
        TYPE(LvlColumnMetaData), INTENT(INOUT) :: lvlcmd
            CALL InitTB0ColumnMetaData(lvlcmd%tb0cmd)
      END SUBROUTINE InitLvlColumnMetaData

!     Added for diversions  nk  Jan. 20/09
      SUBROUTINE InitDivColumnMetaData(divcmd)
        TYPE(DivColumnMetaData), INTENT(INOUT) :: divcmd
            CALL InitTB0ColumnMetaData(divcmd%tb0cmd)
      END SUBROUTINE InitDivColumnMetaData

      SUBROUTINE InitResvinColumnMetaData(resvincmd)
        TYPE(ResvColumnMetaData), INTENT(INOUT) :: resvincmd
            CALL InitTB0ColumnMetaData(resvincmd%tb0cmd)
      END SUBROUTINE InitResvinColumnMetaData




!     Added for swe report June 18/12 nk
      SUBROUTINE InitCourseColumnMetaData(Coursecmd)
            TYPE(CourseColumnMetaData), INTENT(INOUT) :: coursecmd
            CALL InitTB0ColumnMetaData(Coursecmd%tb0cmd)
      END SUBROUTINE InitCourseColumnMetaData






      INTEGER FUNCTION TimeSpanHours(timeStamp1, timeStamp2)
            TYPE(TimeStamp), INTENT(INOUT) :: timeStamp1, timeStamp2 
            TimeSpanHours = 0
            return
      END FUNCTION TimeSpanHours


      LOGICAL FUNCTION IsLatLong(csp)
            TYPE(CoordSysParam), INTENT(INOUT) :: csp

            character*12 projString
            integer strLength
            logical rstat

            IsLatLong = .false.

            projString = csp%projection
            rstat = ToLowerCase(projString)
            strLength = LEN_TRIM(projString)
             
      if(strLength.eq.7) then
            if(projString(1:7) .eq. 'latlong') then
                  IsLatLong = .true.
            endif
      endif

      return

      END FUNCTION IsLatLong


C*******************************************************************
C
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
C
      INTEGER FUNCTION ParseShedParam(header, keyword,
     &                                           keyLen, subString)
            TYPE(ShedParam), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount

            if(keyword(1:keyLen) .eq. ':sourcefilename')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseShedParam = -1
                  else
                        read(value, *)    header%sourceFileName
                        ParseShedParam = 1
                  end if
                  return
            else if(keyword(1:KeyLen) .eq. ':nominalgridsize_al')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseShedParam = -1
                  else
                        read(value, *)    header%nominalGridSize_AL
                        ParseShedParam = 1
                  end if
                  return
            else if(keyword(1:KeyLen) .eq. ':contourinterval')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseShedParam = -1
                  else
                        read(value, *)    header%contourInterval
                        ParseShedParam = 1
                  end if
                  return
            else if(keyword(1:KeyLen) .eq. ':classcount')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseShedParam = -1
                  else
                        read(value, *)    header%classCount
                        ParseShedParam = 1
                  end if
                  return      
            else if(keyword(1:KeyLen) .eq. ':elevconversion')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseShedParam = -1
                  else
                        read(value, *)    header%elevConversion
                        ParseShedParam = 1
                  end if
                  return
            else if(keyword(1:KeyLen) .eq. ':numriverclasses')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseShedParam = -1
                  else
                        read(value, *)    header%numRiverClasses
                        ParseShedParam = 1
                  end if
                  return
            else if(keyword(1:KeyLen) .eq. ':totalnumofgrids')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseShedParam = -1
                  else
                        read(value, *)    header%totalNumOfGrids
                        ParseShedParam = 1
                  end if
                  return      
            else if(keyword(1:KeyLen) .eq. ':numgridsinbasin')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseShedParam = -1
                  else
                        read(value, *)    header%numGridsInBasin
                        ParseShedParam = 1
                  end if
                  return                                    
            else if(keyword(1:KeyLen) .eq. ':debuggridno')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseShedParam = -1
                  else
                        read(value, *)    header%debugGridNo
                        ParseShedParam = 1
                  end if
                  return                                    
            end if
            ! if we're here, then keyword not assigned yet
            ! let's look to the parent blocks

            ParseShedParam = ParseR2CParam(header%r2cp, keyword,
     &                                           keyLen, subString)
            return

      END FUNCTION ParseShedParam

C*******************************************************************
C
C
C           Return Value
C           -1 = Problem
C           0 = keyword not found in Type
C         1 = Successfully assigned
C
C
      INTEGER FUNCTION ParseRainParam(header, keyword,
     &                                           keyLen, subString)
            TYPE(RainParam), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount


!     added Feb. 05/08 -nk- for reading cummulative precip files

            if(keyword(1:keyLen) .eq. ':name')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseRainParam = -1
                  else
                        read(value, *)    header%r2cp%name
                        ParseRainParam = 1
                  endif
            elseif(keyword(1:keyLen) .eq. ':attributeunits')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseRainParam = -1
                  else
                        read(value, *)    header%r2cp%attributeunits
                        ParseRainParam = 1
                  endif
            endif


            ParseRainParam = ParseR2CParam(header%r2cp, keyword,
     &                                           keyLen, subString)

      return

      END FUNCTION ParseRainParam

C*******************************************************************
C
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
C
      INTEGER FUNCTION ParseTempParam(header, keyword,
     &                                           keyLen, subString)
            TYPE(TempParam), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount

            ParseTempParam = ParseR2CParam(header%r2cp, keyword,
     &                                           keyLen, subString)
      return
      END FUNCTION ParseTempParam



C*******************************************************************
C
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
C
      INTEGER FUNCTION ParseGSMParam(header, keyword,
     &                                           keyLen, subString)
            TYPE(GSMParam), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount

            ParseGSMParam = ParseR2CParam(header%r2cp, keyword,
     &                                           keyLen, subString)

      return

      END FUNCTION ParseGSMParam


C*******************************************************************
C
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
C
      INTEGER FUNCTION ParseSWEParam(header, keyword,
     &                                           keyLen, subString)
            TYPE(SWEParam), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount

            if(keyword(1:KeyLen) .eq. ':initheatdeficit')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseSWEParam = -1
                  else
                        read(value, *)    header%initHeatDeficit
                        ParseSWEParam = 1
                  end if
                  return
            endif

            ! if we're here, then keyword not assigned yet
            ! let's look to the parent blocks

            ParseSWEParam = ParseR2CParam(header%r2cp, keyword,
     &                                           keyLen, subString)
            return

      END FUNCTION ParseSWEParam




C*******************************************************************
C
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
C
            INTEGER FUNCTION ParseR2CParam(header, keyword,
     &                                           keyLen, subString)
            TYPE(R2CParam), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount

            ParseR2CParam = 0 

            if(keyword(1:KeyLen) .eq. ':xcount')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseR2CParam = -1
                  else
                        read(value, *)    header%xCount
                        ParseR2CParam = 1
                  end if
                  return
            
            else if(keyword(1:KeyLen) .eq. ':ycount')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseR2CParam = -1
                  else
                        read(value, *)    header%yCount
                        ParseR2CParam = 1
                  end if
                  return
            else if(keyword(1:KeyLen) .eq. ':xdelta')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseR2CParam = -1
                  else
                        read(value, *)    header%xdelta
                        ParseR2CParam = 1
                  end if
                  return
            
            else if(keyword(1:KeyLen) .eq. ':ydelta')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseR2CParam = -1
                  else
                        read(value, *)    header%ydelta
                        ParseR2CParam = 1
                  end if
                  return
            
            else if(keyword(1:KeyLen) .eq. ':xorigin')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseR2CParam = -1
                  else
                        read(value, *)    header%xorigin
                        ParseR2CParam = 1
                  end if
                  return
            
            else if(keyword(1:KeyLen) .eq. ':yorigin')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseR2CParam = -1
                  else
                        read(value, *)    header%yorigin
                        ParseR2CParam = 1
                  end if
                  return

            else if(keyword(1:KeyLen) .eq. ':angle')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseR2CParam = -1
                  else
                        read(value, *)    header%angle
                        ParseR2CParam = 1
                  end if
                  return

            else if(keyword(1:KeyLen) .eq. ':unitconversion')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseR2CParam = -1
                  else
                        read(value, *)    header%unitConv
                        ParseR2CParam = 1
                  end if
                  return
            
            else if(keyword(1:KeyLen) .eq. ':surveydate')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseR2CParam = -1
                  else
                        read(value, *)    header%surveyDate
                        ParseR2CParam = 1
                  end if
                  return

            else if(keyword(1:KeyLen) .eq. ':surveytime')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseR2CParam = -1
                  else
                        read(value, *)    header%surveyTime
                        ParseR2CParam = 1
                  end if


                  return

                  
            end if

            ! if we're here, then unknown keyword...nothing is assigned yet
            ! at this point ParseR2CParam = 0
            ! check if it is a core EnSim parameter
            ParseR2CParam = ParseEnSimParam(header, keyword,
     &                                           keyLen, subString)

            if(ParseR2CParam .ne. 0) return 
                  

            ! if we're here, then unknown keyword...nothing is assigned yet
            ! check if it is a CoordSys parameter
            ParseR2CParam = ParseCoordSysParam(header, keyword,
     &                                           keyLen, subString)

            return

      END FUNCTION ParseR2CParam


C*******************************************************************
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
            INTEGER FUNCTION ParseCoordSysParam(header, keyword,
     &                                           keyLen, subString)
            TYPE(R2CParam), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount

            ParseCoordSysParam = 0 

            if(keyword(1:KeyLen) .eq. ':projection')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseCoordSysParam = -1
                  else
                        read(value, *)    header%csp%projection
                        ParseCoordSysParam = 1
                  end if
                  return

            else if(keyword(1:KeyLen) .eq. ':ellipsoid')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseCoordSysParam = -1
                  else
                        read(value, *)    header%csp%ellipsoid
                        ParseCoordSysParam = 1
                  end if
                  return

            else if(keyword(1:KeyLen) .eq. ':zone')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseCoordSysParam = -1
                  else
                        read(value, *)    header%csp%zone
                        ParseCoordSysParam = 1
                  end if
                  return
            end if

            return

                  
            END FUNCTION ParseCoordSysParam

C*******************************************************************
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
            INTEGER FUNCTION ParseCoordSysParamTB0(header, keyword,
     &                                           keyLen, subString)
            TYPE(TB0Param), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount

            ParseCoordSysParamTB0 = 0 
            if(keyword(1:KeyLen) .eq. ':projection')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseCoordSysParamTB0 = -1
                  else
                        read(value, *)    header%csp%projection
                        ParseCoordSysParamTB0 = 1
                  end if
                  return

            else if(keyword(1:KeyLen) .eq. ':ellipsoid')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseCoordSysParamTB0 = -1
                  else
                        read(value, *)    header%csp%ellipsoid
                        ParseCoordSysParamTB0 = 1
                  end if
                  return

            else if(keyword(1:KeyLen) .eq. ':zone')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseCoordSysParamTB0 = -1
                  else
                        read(value, *)    header%csp%zone
                        ParseCoordSysParamTB0 = 1
                  end if
                  return
            end if

            return

            END FUNCTION ParseCoordSysParamTB0



C*******************************************************************
C
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
C
            INTEGER FUNCTION ParseEnSimParam(header, keyword,
     &                                           keyLen, subString)
            TYPE(R2CParam), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount
            logical rstat

            ParseEnSimParam = 0 

            if(keyword(1:KeyLen) .eq. ':filetype')then           ! filetype line
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEnSimParam = -1
                        return
                  else
                        header%ep%fileType = value
                        ParseEnSimParam = 1
                        if(value .NE. "r2c") then                                    ! data type
                              write(6,'((A))') 'ERROR in file: '
                              write(6,'(3(A))') 'Expected: r2c', 
     &                                               ' Found: ',value
                              ParseEnSimParam = -1
                              return
                        end if
                  endif 
                        
                  subString = tmpString
                  if(SplitLine(subString, value, tmpString) .eq. 0)then ! file format
                        ParseEnSimParam = -1
                        return
                  end if
                  rStat = ToLowerCase(value)
                  if(value .eq. 'ascii')then                ! ASCII file format
                        header%ep%fileFormat = 0
                        ParseEnSimParam = 1
                  else if(value .eq. 'binary')then    ! BINARY file format
                        header%ep%fileFormat = 1
                        ParseEnSimParam = 1
                  else                                            ! UNKNOWN file format
                        write(6,'(2(A))')' Stopped: ',
     &                                        'invalid file format'
                        ParseEnSimParam = -1
                  end if
                  return
            
            else if(keyword(1:KeyLen) .eq. ':application')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEnSimParam = -1
                  else
                        read(value, *)    header%ep%application
                        ParseEnSimParam = 1
                  end if
                  return

            else if(keyword(1:KeyLen) .eq. ':writtenby')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEnSimParam = -1
                  else
                        read(value, *)    header%ep%writtenBy
                        ParseEnSimParam = 1
                  end if
                  return

            else if(keyword(1:KeyLen) .eq. ':creationdate')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEnSimParam = -1
                  else
                        read(value, *)    header%ep%creationDate
                        ParseEnSimParam = 1
                  end if
                  return

            else if(keyword(1:KeyLen) .eq. ':datatype')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEnSimParam = -1
                  else
                        read(value, *)    header%ep%dataType
                        ParseEnSimParam = 1
                  end if
                  return

            else if(keyword(1:KeyLen) .eq. ':version')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEnSimParam = -1
                  else
                        read(value, *)    header%ep%version
                        ParseEnSimParam = 1
                  end if
                  return


            else if(keyword(1:KeyLen) .eq. ':name')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEnSimParam = -1
                  else
                        read(value, *)    header%ep%name
                        ParseEnSimParam = 1
                  end if
                  return

            else if(keyword(1:KeyLen) .eq. ':attributename')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                    ParseEnSimParam = -1
                  else
                    if(IsNumber(value)) then
                          header%ep%attCount = header%ep%attCount+1
                          read(value, *)
     &                 header%ep%attList(header%ep%attCount)%index
                    else
                          write(*,'(4(A))')
     &                    'Read ''', TRIM(value), '''',
     &                          ' but expected a number'
                          ParseEnSimParam = -1
                          return
                    end if
                    subString = tmpString
                    if(SplitLine(subString,value,tmpString).eq.0)then
                          ParseEnSimParam = -1
                    else
                      read(value, *)
     &                      header%ep%attList(header%ep%attCount)%name
                    ParseEnSimParam = 1
                    end if
                  endif
                  return
            end if

            return
                  
            END FUNCTION ParseEnSimParam

C*******************************************************************
C
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
C
            INTEGER FUNCTION ParseEnSimParamTB0(header, keyword,
     &                                           keyLen, subString)
            TYPE(TB0Param), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount
            logical rstat

            ParseEnSimParamTB0 = 0 

            if(keyword(1:KeyLen) .eq. ':filetype')then           ! filetype line
               if(SplitLine(subString, value, tmpString) .eq. 0)then
                     ParseEnSimParamTB0 = -1
                     return
               else
                     header%ep%fileType = value
                     ParseEnSimParamTB0 = 1
                     if(value .NE. "tb0") then                                    ! data type
                           write(6,'((A))') 'ERROR in file: '
                           write(6,'(3(A))') 'Expected: TB0', 
     &                                         ' Found: ',value
                           ParseEnSimParamTB0 = -1
                           return
                     end if
               endif 
                     
               subString = tmpString
               if(SplitLine(subString, value, tmpString) .eq. 0)then ! file format
                     ParseEnSimParamTB0 = -1
                     return
               end if
               rStat = ToLowerCase(value)
               if(value .eq. 'ascii')then                ! ASCII file format
                     header%ep%fileFormat = 0
                     ParseEnSimParamTB0 = 1
               else if(value .eq. 'binary')then    ! BINARY file format
                     header%ep%fileFormat = 1
                     ParseEnSimParamTB0 = 1
               else                                            ! UNKNOWN file format
                     write(6,'(2(A))')' Stopped: ',
     &                              'invalid file format'
                     ParseEnSimParamTB0 = -1
               end if
               return
            
            else if(keyword(1:KeyLen) .eq. ':application')then
               if(SplitLine(subString, value, tmpString) .eq. 0)then
                     ParseEnSimParamTB0 = -1
               else
                     read(value, *)    header%ep%application
                     ParseEnSimParamTB0 = 1
               end if
               return

            else if(keyword(1:KeyLen) .eq. ':writtenby')then
               if(SplitLine(subString, value, tmpString) .eq. 0)then
                     ParseEnSimParamTB0 = -1
               else
                     read(value, *)    header%ep%writtenBy
                     ParseEnSimParamTB0 = 1
               end if
               return

            else if(keyword(1:KeyLen) .eq. ':creationdate')then
               if(SplitLine(subString, value, tmpString) .eq. 0)then
                     ParseEnSimParamTB0 = -1
               else
                     read(value, *)    header%ep%creationDate
                     ParseEnSimParamTB0 = 1
               end if
               return

            else if(keyword(1:KeyLen) .eq. ':datatype')then
               if(SplitLine(subString, value, tmpString) .eq. 0)then
                     ParseEnSimParamTB0 = -1
               else
                     read(value, *)    header%ep%dataType
                     ParseEnSimParamTB0 = 1
               end if
               return

            else if(keyword(1:KeyLen) .eq. ':version')then
               if(SplitLine(subString, value, tmpString) .eq. 0)then
                     ParseEnSimParamTB0 = -1
               else
                     read(value, *)    header%ep%version
                     ParseEnSimParamTB0 = 1
               end if
               return


            else if(keyword(1:KeyLen) .eq. ':name')then
               if(SplitLine(subString, value, tmpString) .eq. 0)then
                     ParseEnSimParamTB0 = -1
               else
                     read(value, *)    header%ep%name
                     ParseEnSimParamTB0 = 1
               end if
               return

            else if(keyword(1:KeyLen) .eq. ':attributename')then
               if(SplitLine(subString, value, tmpString) .eq. 0)then
                 ParseEnSimParamTB0 = -1
               else
                 if(IsNumber(value)) then
                       header%ep%attCount = header%ep%attCount+1
                       read(value, *)
     &                 header%ep%attList(header%ep%attCount)%index
                 else
                       write(*,'(4(A))')
     &                 'Read ''', TRIM(value), '''',
     &                       ' but expected a number'
                       ParseEnSimParamTB0 = -1
                       return
                 end if

                 subString = tmpString
                 if(SplitLine(subString, value, tmpString) .eq. 0)then
                       ParseEnSimParamTB0 = -1
                 else
                       read(value, *)
     &                       header%ep%attList(header%ep%attCount)%name
                 ParseEnSimParamTB0 = 1
                 end if
               endif
               return
            end if

            return
                  
            END FUNCTION ParseEnSimParamTB0


C*******************************************************************
C
C
C           Return Value
C           -1 = Problem
C           0 = keyword not found in Type
C         1 = Successfully assigned
C
C
      INTEGER FUNCTION ParseFrameLine(frameRec, keyword,
     &                                           keyLen, subString)
            TYPE(FrameRecord), INTENT(INOUT) :: frameRec 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount

            ParseFrameLine = 0
C Frame     
            if(keyword(1:KeyLen) .eq. ':frame')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseFrameLine = -1
                        return
                  else
                        read(value, *) frameRec%frame
C                       print*,     ' '
C                       print*,     'frame:',frameRec%frame
                  endif

C Step
                  subString = tmpString
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseFrameLine = -1
                        return
                  else
                        read(value, *) frameRec%step
C                       print*,     'step:',frameRec%step
                  end if

C TimeStamp
                  subString = tmpString
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                     ParseFrameLine = -1
                     return
                  else
                     if(ParseTimeStamp(frameRec%tStamp, value).eq.0)then 
                           ParseFrameLine = -1
                           return
                     end if
                  end if

CC Hour=
C                 subString = tmpString
C                 if(SplitLine(subString, value, tmpString) .eq. 0)then
C                       ParseFrameLine = -1
C                       return
C                 end if
C
CC hour
C                 subString = tmpString
C                 if(SplitLine(subString, value, tmpString) .eq. 0)then
C                       ParseFrameLine = -1
C                       return
C                 else
C                       read(value, *) frameRec%hour
C                       print*,     'hour:',frameRec%hour
C                 end if
C
CC hour to next
C                 subString = tmpString
C                 if(SplitLine(subString, value, tmpString) .eq. 0)then
C                       ParseFrameLine = -1
C                       return
C                 else
C                       read(value, *) frameRec%hoursToNext
C                       print*,     'hoursToNext:',frameRec%hoursToNext
C                 end if
C
CC source
C                 subString = tmpString
C                 if(SplitLine(subString, value, tmpString) .eq. 0)then
C                       ParseFrameLine = -1
C                       return
C                 else
C                       read(value, *) frameRec%source
C                       print*,     'source:',frameRec%source
C                 end if

                  ParseFrameLine = 1
            end if

            return

      END FUNCTION ParseFrameLine


C*******************************************************************
C
C
C           Return Value
C           0 = Problem
C         1 = Successfully assigned
C
C
      INTEGER FUNCTION ParseTimeStamp(tStamp, string)
            TYPE(TimeStamp), INTENT(INOUT) :: tStamp 
            CHARACTER*(*), INTENT(INOUT) :: string
            character*128  value
            character*32768 tmpString, localLine
            integer lineLen, rStat,ios
            integer slash, space, colon, point 
      

            localLine = string                  ! work on local copy of the string
            rStat = Detab(localLine)                        ! change any tabs to spaces
            localLine = ADJUSTL(localLine)                  ! get rid of leading spaces
            lineLen = LEN_TRIM(localLine) ! find the new length, excluding trailing spaces
      
            ParseTimeStamp = 0

C a Date field will have a '/'
            slash = INDEX(localLine(1:lineLen),'/')   ! Look for a '/' character

            if(slash.GT.0) then           ! we have a date field
C     Extract the year
                  value = localLine(1:slash-1)
                  read(value, *,iostat=ios) tStamp%year                
                  if(ios.ne.0)then
                      write(98,*)'Invalid time stamp in an r2c file'
                      write(98,*)'Last valid:',
     *                        tStamp%year,tStamp%month,tStamp%day 
                      write(98,*)'Program aborted in EF_module @ 1412'
                      print*,'Invalid time stamp in an r2c file'
                      print*,'Last valid:',
     *                        tStamp%year,tStamp%month,tStamp%day 
                      stop 'Program aborted in EF_module @ 1412'
                  endif
                  localLine = localLine(slash+1:lineLen)
                  localLine = ADJUSTL(localLine)                  ! get rid of leading spaces
                  lineLen = LEN_TRIM(localLine) ! find the new length, excluding trailing spaces

C     Extract the month
                  slash = INDEX(localLine(1:lineLen),'/')   ! Look for a '/' character
                  if(slash.EQ.0) then
                        return
                  end if
                  value = localLine(1:slash-1)
                  read(value, *,iostat=ios) tStamp%month               
                  if(ios.ne.0)then
                      write(98,*)'Invalid time stamp in an r2c file'
                      write(98,*)'Last valid:',
     *                        tStamp%year,tStamp%month,tStamp%day 
                      write(98,*)'Program aborted in EF_module @ 1434'
                      print*,'Invalid time stamp in an r2c file'
                      print*,'Last valid:',
     *                        tStamp%year,tStamp%month,tStamp%day 
                      stop 'Program aborted in EF_module @ 1434'
                  endif
                  localLine = localLine(slash+1:lineLen)
                  localLine = ADJUSTL(localLine)                  ! get rid of leading spaces
                  lineLen = LEN_TRIM(localLine) ! find the new length, excluding trailing spaces

C     Extract the day
                  space = INDEX(localLine(1:lineLen),' ')   ! Look for a ' ' character
                  if(space.EQ.0) then
                        return
                  end if
                  value = localLine(1:space-1)
                  read(value, *,iostat=ios) tStamp%day  
                  if(ios.ne.0)then
                      write(98,*)'Invalid time stamp in an r2c file'
                      write(98,*)'Last valid:',
     *                        tStamp%year,tStamp%month,tStamp%day 
                      write(98,*)'Program aborted in EF_module @ 1451'
                       print*,'Invalid time stamp in an r2c file'
                      print*,'Last valid:',
     *                        tStamp%year,tStamp%month,tStamp%day 
                      stop 'Program aborted in EF_module @ 1451'
                  endif
                  localLine = localLine(space+1:lineLen)
                  localLine = ADJUSTL(localLine)                  ! get rid of leading spaces
                  lineLen = LEN_TRIM(localLine) ! find the new length, excluding trailing spaces
            end if


C Now the time field
            
C     Extract the hour
            colon = INDEX(localLine(1:lineLen),':')   ! Look for a ':' character
            if(colon.EQ.0) then
                  return
            end if 
            value = localLine(1:colon-1)
            read(value, *,iostat=ios) tStamp%hour                
                  if(ios.ne.0)then
                      write(98,*)'Invalid time stamp in an r2c file'
                      write(98,*)'Last valid:',
     *                        tStamp%year,tStamp%month,tStamp%day 
                      write(98,*)'Program aborted in EF_module @ 1472'
                      print*,'Invalid time stamp in an r2c file'
                      print*,'Last valid:',
     *                        tStamp%year,tStamp%month,tStamp%day 
                      stop 'Program aborted in EF_module @ 1472'
                  endif
            localLine = localLine(colon+1:lineLen)
            localLine = ADJUSTL(localLine)                  ! get rid of leading spaces
            lineLen = LEN_TRIM(localLine) ! find the new length, excluding trailing spaces

C     Extract the minute
            if(lineLen.GT.0)then 
                  colon = INDEX(localLine(1:lineLen),':')   ! Look for a ':' character
                  if(colon.EQ.0) then
                        value = localLine(1:lineLen)
                        localLine=''
                  else 
                        value = localLine(1:colon-1)
                        localLine = localLine(colon+1:lineLen)
                  endif
                  localLine = localLine(colon+1:lineLen)
                  localLine = ADJUSTL(localLine)      ! get rid of leading spaces
                  lineLen = LEN_TRIM(localLine) ! find the new length, excluding trailing spaces
                  read(value, *) tStamp%min
            end if                        

C     Extract the sec.msec
C                 first check for msec
            if(lineLen.GT.0)then 
                  point = INDEX(localLine(1:lineLen),'.')   ! Look for a '.' character
                  if(point.GT.0) then 
                        value = localLine(point:lineLen)
                        read(value, *) tStamp%msec                
                        localLine = localLine(1:point-1)
                        localLine = ADJUSTL(localLine)      ! get rid of leading spaces
                        lineLen = LEN_TRIM(localLine) ! find the new length, excluding trailing spaces
                  endif
C                 Now extract the second
                  value = localLine(1:lineLen)
                  read(value, *) tStamp%sec
            end if                        

            ParseTimeStamp = 1

            return

      END FUNCTION ParseTimeStamp

C*******************************************************************
C
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
C
      INTEGER FUNCTION ParseTB0Param(header, keyword,
     &                                           keyLen, subString)
            TYPE(TB0Param), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount

            ParseTB0Param = 0 

            if(keyword(1:KeyLen) .eq. ':deltat')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseTB0Param = -1
                  else
                        read(value, *)    header%deltaT
                        ParseTB0Param = 1
                  end if
                  return
            else if(keyword(1:KeyLen) .eq. ':unitconversion')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseTB0Param = -1
                  else
                        read(value, *)    header%unitConv
                        ParseTB0Param = 1
                  end if
                  return
!                 added Mar. 21/14 NK for ragmet & tmp
            else if(keyword(1:KeyLen) .eq. ':name')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseTB0Param = -1
                  else
                        read(value, *)    header%sourceKey
                        ParseTB0Param = 1
                  end if
                  return
            end if

            ! if we're here, then unknown keyword...nothing is assigned yet
            ! at this point ParseTB0Param = 0
            ! check if it is a core EnSim parameter
            ParseTB0Param = ParseEnSimParamTB0(header, keyword,
     &                                           keyLen, subString)

            if(ParseTB0Param .ne. 0) return 
                  

            ! if we're here, then unknown keyword...nothing is assigned yet
            ! check if it is a CoordSys parameter
            ParseTB0Param = ParseCoordSysParamTB0(header, keyword,
     &                                           keyLen, subString)

            return

      END FUNCTION ParseTB0Param


C*******************************************************************
C
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
C
      INTEGER FUNCTION ParseFlowParam(header, keyword,
     &                                           keyLen, subString)
            TYPE(FlowParam), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount

            if(keyword(1:KeyLen) .eq. ':routingdeltat')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseFlowParam = -1
                  else
                        read(value, *)    header%routingDeltaT
                        ParseFlowParam = 1
                  end if
                  return
            else if(keyword(1:KeyLen) .eq. ':fillflag')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseFlowParam = -1
                  else
                        read(value, *)    header%fillFlag
                        ParseFlowParam = 1
                  end if
                  return
            end if

            ! if we're here, then keyword not assigned yet
            ! let's look to the parent blocks

            ParseFlowParam = ParseTB0Param(header%tb0p, keyword,
     &                                           keyLen, subString)

            return

      END FUNCTION ParseFlowParam

      INTEGER FUNCTION ParseTemp2Param(header, keyword,
     &                                           keyLen, subString)
            TYPE(Temp2Param), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount

            ParseTemp2Param = ParseTB0Param(header%tb0p, keyword,
     &                                           keyLen, subString)
      return
      END FUNCTION ParseTemp2Param


C*******************************************************************
C
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
C
      INTEGER FUNCTION ParseRagParam(header, keyword,
     &                                           keyLen, subString)
            TYPE(RagParam), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount


            ParseRagParam = ParseTB0Param(header%tb0p, keyword,
     &                                           keyLen, subString)
            return

      END FUNCTION ParseRagParam


C*******************************************************************
C
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
C
      INTEGER FUNCTION ParseResvParam(header, keyword,
     &                                           keyLen, subString)
            TYPE(ResvParam), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount


            ! if we're here, then keyword not assigned yet
            ! let's look to the parent blocks

            ParseResvParam = ParseTB0Param(header%tb0p, keyword,
     &                                           keyLen, subString)
            return

      END FUNCTION ParseResvParam


C*******************************************************************
C
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
C
      INTEGER FUNCTION ParseResvinParam(header, keyword,
     &                                           keyLen, subString)
            TYPE(ResvinParam), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount


            ! if we're here, then keyword not assigned yet
            ! let's look to the parent blocks

            ParseResvinParam = ParseTB0Param(header%tb0p, keyword,
     &                                           keyLen, subString)
            return

      END FUNCTION ParseResvinParam



C*******************************************************************
C
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
C
!     rev. 9.5.52  Jan.  20/09  - NK: added reading yyyymmdd_div.pt2 for diversions
      INTEGER FUNCTION ParseDivParam(header, keyword,
     &                                           keyLen, subString)
            TYPE(DivParam), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount

            ! if we're here, then keyword not assigned yet
            ! let's look to the parent blocks

            ParseDivParam = ParseTB0Param(header%tb0p, keyword,
     &                                           keyLen, subString)
            return

      END FUNCTION ParseDivParam

!     rev. 9.8.02  JUne  18/06  - NK: added swe report tb0 file
      INTEGER FUNCTION ParseCourseParam(header, keyword,
     &                                           keyLen, subString)
            TYPE(CourseParam), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount

            ! if we're here, then keyword not assigned yet
            ! let's look to the parent blocks

            ParseCourseParam = ParseTB0Param(header%tb0p, keyword,
     &                                           keyLen, subString)
            return

      END FUNCTION ParseCourseParam





!     rev. 9.8.02  Aug.  02/11  - NK: added lake level tb0 file
      INTEGER FUNCTION ParseLvlParam(header, keyword,
     &                                           keyLen, subString)
            TYPE(LvlParam), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount

            ! if we're here, then keyword not assigned yet
            ! let's look to the parent blocks

            ParseLvlParam = ParseTB0Param(header%tb0p, keyword,
     &                                           keyLen, subString)
            return

      END FUNCTION ParseLvlParam






C*******************************************************************
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
      INTEGER FUNCTION ParseTB0ColumnMetaData(header, keyword,
     &                                           keyLen, subString)
            TYPE(TB0ColumnMetaData), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount,n, ideallocate

            tmpString = subString
            ParseTB0ColumnMetaData = 0

            if(keyword(1:KeyLen) .eq. ':columnlocationx')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseTB0ColumnMetaData = -1
                        return
                  end if

                  header%colCount = wordCount
                  if(allocated(header%colLocX)) then
                        deallocate(header%colLocX,stat=ideallocate)
                  end if
                  allocate(header%colLocX(wordCount))

                  read(tmpString,*)(header%colLocX(n),n=1,wordCount)

c      just checking
c                 do n=1,wordCount
c                       print*,'X',header%colLocX(n)
c                 enddo
c                 pause 11111

                  ParseTB0ColumnMetaData = 1
                  return

            else if(keyword(1:KeyLen) .eq. ':columnlocationy')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseTB0ColumnMetaData = -1
                        return
                  end if

                  header%colCount = wordCount
                  if(allocated(header%colLocY)) then
                        deallocate(header%colLocY,stat=ideallocate)
                  end if
                  allocate(header%colLocY(wordCount))

                  read(tmpString,*)(header%colLocY(n),n=1,wordCount)

C just checking
c                 do n=1,wordCount
c                       print*,'Y',header%colLocY(n)
c                 enddo
c                 pause 11112


                  ParseTB0ColumnMetaData = 1
                  
                  return


            elseif(keyword(1:KeyLen) .eq. ':columnlocationx1')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseTB0ColumnMetaData = -1
                        return
                  end if

                  header%colCount = wordCount
                  if(allocated(header%colLocX1)) then
                        deallocate(header%colLocX1,stat=ideallocate)
                  end if
                  allocate(header%colLocX1(wordCount))

                  read(tmpString,*)(header%colLocX1(n),n=1,wordCount)

c      just checking
c                 do n=1,wordCount
c                       print*,'X1',header%colLocX1(n)
c                 enddo
c                 pause 11113

                  ParseTB0ColumnMetaData = 1
                  return

            else if(keyword(1:KeyLen) .eq. ':columnlocationy1')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseTB0ColumnMetaData = -1
                        
                        return
                  end if

                  header%colCount = wordCount
                  if(allocated(header%colLocY1)) then
                        deallocate(header%colLocY1,stat=ideallocate)
                  end if
                  allocate(header%colLocY1(wordCount))

                  read(tmpString,*)(header%colLocY1(n),n=1,wordCount)

c just checking
c                 do n=1,wordCount
c                       print*,'Y1',header%colLocY1(n)
c                 enddo
c                 pause 11114

                  ParseTB0ColumnMetaData = 1
                  
                  return


      else if(keyword(1:KeyLen) .eq. ':columnname')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseTB0ColumnMetaData = -1
                        return
                  end if

                  header%colCount = wordCount
                  if(allocated(header%colName)) then
                        deallocate(header%colName,stat=ideallocate)
                  end if
                  allocate(header%colName(wordCount))

                  read(tmpString,*)(header%colName(n),n=1,wordCount)

C just checking
c                 do n=1,wordCount
c                       print*,header%colName(n)
c                 enddo

                  ParseTB0ColumnMetaData = 1
                  return

      else if(keyword(1:KeyLen) .eq. ':columclasscount')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseTB0ColumnMetaData = -1
                        return
                  end if

                  header%colCount = wordCount
                  if(allocated(header%colType)) then
                        deallocate(header%colType,stat=ideallocate)
                  end if
                  allocate(header%colType(wordCount))

                  read(tmpString,*)(header%colType(n),n=1,wordCount)

C just checking
c                 do n=1,wordCount
c                       print*,header%colType(n)
c                 enddo

                  ParseTB0ColumnMetaData = 1
                  return


      else if(keyword(1:KeyLen) .eq. ':columnunits')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseTB0ColumnMetaData = -1
                        return
                  end if

                  header%colCount = wordCount
                  if(allocated(header%colUnits)) then
                        deallocate(header%colUnits,stat=ideallocate)
                  end if
                  allocate(header%colUnits(wordCount))

                  read(tmpString,*)(header%colUnits(n),n=1,wordCount)

C just checking
c                 do n=1,wordCount
c                       print*,header%colUnits(n)
c                 enddo

                  ParseTB0ColumnMetaData = 1

!                 need to allocate elevation in case the 
!                 :Elevation 
!                 is not in the tag or rag files
!                 see below
!                 Also have to set default value;
                  if(allocated(header%elevation)) then
                        deallocate(header%elevation,stat=ideallocate)  
                  end if
	  	      allocate(header%elevation(wordCount))
!                 Initialize incase the elevations are missing in the tb0 files
!                 for ragmet & tmp
	            do n=1,wordCount
	              header%elevation(n)=-999.0
	            end do
	              
                  return

!                 added Sept. 24/09 nk
      else if(keyword(1:KeyLen) .eq. ':elevation')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseTB0ColumnMetaData = -1
                        return
                  end if

!                 just making sure it is allocated
!                 this time it is used
                  header%colCount = wordCount
                  if(allocated(header%elevation)) then
                        deallocate(header%elevation,stat=ideallocate)  
                  end if
	  	        allocate(header%elevation(wordCount))

		    	read(tmpString,*)(header%elevation(n),n=1,wordCount)

c just checking
c		    	print*
c			    do n=1,wordCount
c				  print*,'elevation station ',n,header%elevation(n)
c				  print*
c			    enddo

		    	ParseTB0ColumnMetaData = 1
                  return
!                 added Nov. 03/22 nk

            end if

            return

      END FUNCTION ParseTB0ColumnMetaData



C*******************************************************************
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
      INTEGER FUNCTION ParseFlowColumnMetaData(header, keyword,
     &                                           keyLen, subString)
            TYPE(FlowColumnMetaData), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768  tmpString
            integer lineLen, wordCount,n, ideallocate

            tmpString = subString

            if(keyword(1:KeyLen) .eq. ':coeff1')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseFlowColumnMetaData = -1
                        return
                  end if

                  header%tb0cmd%colCount = wordCount
                  if(allocated(header%colCoeff1)) then
                        deallocate(header%colCoeff1,stat=ideallocate)
                  end if
                  allocate(header%colCoeff1(wordCount))
      
                  read(tmpString,*)(header%colCoeff1(n),n=1,wordCount)

C just checking
C                 do n=1,wordCount
C                       print*,header%colCoeff1(n)
C                 enddo

                  ParseFlowColumnMetaData = 1
                  return

            else if(keyword(1:KeyLen) .eq. ':coeff2')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseFlowColumnMetaData = -1
                        return
                  end if

                  header%tb0cmd%colCount = wordCount
                  if(allocated(header%colCoeff2)) then
                        deallocate(header%colCoeff2,stat=ideallocate)
                  end if
                  allocate(header%colCoeff2(wordCount))
                  
                  read(tmpString,*)(header%colCoeff2(n),n=1,wordCount)

C just checking
C                 do n=1,wordCount
C                       print*,header%colCoeff2(n)
C                 enddo

                  ParseFlowColumnMetaData = 1
                  return

            else if(keyword(1:KeyLen) .eq. ':coeff3')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseFlowColumnMetaData = -1
                        return
                  end if

                  header%tb0cmd%colCount = wordCount
                  if(allocated(header%colCoeff3)) then
                        deallocate(header%colCoeff3,stat=ideallocate)
                  end if
                  allocate(header%colCoeff3(wordCount))
                  read(tmpString,*)(header%colCoeff3(n),n=1,wordCount)

C just checking
c                 do n=1,wordCount
c                       print*,header%colCoeff3(n)
c                 enddo

                  ParseFlowColumnMetaData = 1
                  return
            else if(keyword(1:KeyLen) .eq. ':coeff4')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseFlowColumnMetaData = -1
                        return
                  end if
            
                  header%tb0cmd%colCount = wordCount
                  if(allocated(header%colCoeff4)) then
                        deallocate(header%colCoeff4,stat=ideallocate)
                  end if
                  allocate(header%colCoeff4(wordCount))

                  read(tmpString,*)(header%colCoeff4(n),n=1,wordCount)

C just checking
C                 do n=1,wordCount
C                       print*,header%colCoeff4(n)
C                 enddo

                  ParseFlowColumnMetaData = 1
                  return

            else if(keyword(1:KeyLen) .eq. ':value1')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseFlowColumnMetaData = -1
                        return
                  end if

                  header%tb0cmd%colCount = wordCount
                  if(allocated(header%colValue1)) then
                        deallocate(header%colValue1,stat=ideallocate)
                  end if
	  	      allocate(header%colValue1(wordCount))

			read(tmpString,*)(header%colValue1(n),n=1,wordCount)

C just checking
C			do n=1,wordCount
C				print*,header%colValue1(n)
C			enddo

			ParseFlowColumnMetaData = 1
			return


            end if

            ! if we're here, then keyword not assigned yet
            ! let's look to the parent blocks

            ParseFlowColumnMetaData = ParseTB0ColumnMetaData 
     &                    (header%tb0cmd,keyword, keyLen, subString)
            return


            return

      END FUNCTION ParseFlowColumnMetaData



!     rev. 9.5.52  Jan.  20/09  - NK: added reading yyyymmdd_div.pt2 for diversions
C*******************************************************************
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
      INTEGER FUNCTION ParseDivColumnMetaData(header, keyword,
     &                                           keyLen, subString)
            TYPE(DivColumnMetaData), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768  tmpString
            integer lineLen, wordCount,n, ideallocate

            tmpString = subString

            if(keyword(1:KeyLen) .eq. ':value1')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseDivColumnMetaData = -1
                        return
                  end if

                  header%tb0cmd%colCount = wordCount
                  if(allocated(header%colValue1)) then
                        deallocate(header%colValue1,stat=ideallocate)
                  end if
	  	      allocate(header%colValue1(wordCount))

			read(tmpString,*)(header%colValue1(n),n=1,wordCount)

cc just checking
c			do n=1,wordCount
c				print*,header%colValue1(n)
c			enddo
c            pause 4321		

			ParseDivColumnMetaData = 1
			return


            end if

            ! if we're here, then keyword not assigned yet
            ! let's look to the parent blocks

            ParseDivColumnMetaData = ParseTB0ColumnMetaData 
     &                    (header%tb0cmd,keyword, keyLen, subString)
            return


            return

      END FUNCTION ParseDivColumnMetaData











C*******************************************************************
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
      INTEGER FUNCTION ParseRagColumnMetaData(header, keyword,
     &                                           keyLen, subString)
            TYPE(RagColumnMetaData), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768  tmpString
            integer lineLen, wordCount,n, ideallocate

            tmpString = subString

            ! if we're here, then keyword not assigned yet
            ! let's look to the parent blocks

            ParseRagColumnMetaData = ParseTB0ColumnMetaData 
     &                     (header%tb0cmd,keyword, keyLen, subString)
            return

      END FUNCTION ParseRagColumnMetaData



C*******************************************************************
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
      INTEGER FUNCTION ParseCourseColumnMetaData(header, keyword,
     &                                           keyLen, subString)
            TYPE(CourseColumnMetaData), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768  tmpString
            integer lineLen, wordCount,n, ideallocate

            tmpString = subString

            ! if we're here, then keyword not assigned yet
            ! let's look to the parent blocks

            ParseCourseColumnMetaData = ParseTB0ColumnMetaData 
     &                     (header%tb0cmd,keyword, keyLen, subString)
            return

      END FUNCTION ParseCourseColumnMetaData



      INTEGER FUNCTION ParseLvlColumnMetaData(header, keyword,
     &                                           keyLen, subString)
            TYPE(LvlColumnMetaData), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768  tmpString
            integer lineLen, wordCount,n, ideallocate

            tmpString = subString

            ! if we're here, then keyword not assigned yet
            ! let's look to the parent blocks

            ParseLvlColumnMetaData = ParseTB0ColumnMetaData 
     &                     (header%tb0cmd,keyword, keyLen, subString)
            return

      END FUNCTION ParseLvlColumnMetaData








C*******************************************************************
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
      INTEGER FUNCTION ParseResvColumnMetaData(header, keyword,
     &                                           keyLen, subString)
            TYPE(ResvColumnMetaData), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount,n, ideallocate

            tmpString = subString

            if(keyword(1:KeyLen) .eq. ':coeff1')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseResvColumnMetaData = -1
                        return
                  end if

                  header%tb0cmd%colCount = wordCount
                  if(.NOT.allocated(header%colCoeff1)) then
                        allocate(header%colCoeff1(wordCount))
                  end if
                  read(tmpString,*)(header%colCoeff1(n),n=1,wordCount)

C just checking
C                 do n=1,wordCount
C                       print*,header%colCoeff1(n)
C                 enddo

                  ParseResvColumnMetaData = 1
                  return

            else if(keyword(1:KeyLen) .eq. ':coeff2')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseResvColumnMetaData = -1
                        return
                  end if

                  header%tb0cmd%colCount = wordCount
                  if(.NOT.allocated(header%colCoeff2)) then
                        allocate(header%colCoeff2(wordCount))
                  end if

                  read(tmpString,*)(header%colCoeff2(n),n=1,wordCount)

C just checking
C                 do n=1,wordCount
C                       print*,header%colCoeff2(n)
C                 enddo

                  ParseResvColumnMetaData = 1
                  return

            else if(keyword(1:KeyLen) .eq. ':coeff3')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseResvColumnMetaData = -1
                        return
                  end if

                  header%tb0cmd%colCount = wordCount
                  deallocate(header%colCoeff3,stat=ideallocate)
                  allocate(header%colCoeff3(wordCount))

                  read(tmpString,*)(header%colCoeff3(n),n=1,wordCount)

C just checking
c                 do n=1,wordCount
c                       print*,header%colCoeff3(n)
c                 enddo

                  ParseResvColumnMetaData = 1
                  return
            else if(keyword(1:KeyLen) .eq. ':coeff4')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseResvColumnMetaData = -1
                        return
                  end if

                  header%tb0cmd%colCount = wordCount
                  deallocate(header%colCoeff4,stat=ideallocate)
                  allocate(header%colCoeff4(wordCount))

                  read(tmpString,*)(header%colCoeff4(n),n=1,wordCount)

C just checking
C                 do n=1,wordCount
C                       print*,header%colCoeff4(n)
C                 enddo

                  ParseResvColumnMetaData = 1
                  return

            else if(keyword(1:KeyLen) .eq. ':coeff5')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseResvColumnMetaData = -1
                        return
                  end if

                  header%tb0cmd%colCount = wordCount
                  deallocate(header%colCoeff5,stat=ideallocate)
                  allocate(header%colCoeff5(wordCount))

                  read(tmpString,*)(header%colCoeff5(n),n=1,wordCount)

C just checking
C                 do n=1,wordCount
C                       print*,header%colCoeff5(n)
C                 enddo

                  ParseResvColumnMetaData = 1
                  return

            else if(keyword(1:KeyLen) .eq. ':coeff6')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseResvColumnMetaData = -1
                        return
                  end if

                  header%tb0cmd%colCount = wordCount
                  deallocate(header%colCoeff6,stat=ideallocate)
                  allocate(header%colCoeff6(wordCount))

                  read(tmpString,*)(header%colCoeff6(n),n=1,wordCount)

C just checking
c                 do n=1,wordCount
c                       print*,header%colCoeff6(n)
c                 enddo

                  ParseResvColumnMetaData = 1
                  return

            else if(keyword(1:KeyLen) .eq. ':coeff7')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseResvColumnMetaData = -1
                        return
                  end if

                  header%tb0cmd%colCount = wordCount
                  deallocate(header%colCoeff7,stat=ideallocate)
                  allocate(header%colCoeff7(wordCount))

                  read(tmpString,*)(header%colCoeff7(n),n=1,wordCount)

C just checking
c                 do n=1,wordCount
c                       print*,header%colCoeff7(n)
c                 enddo

                  ParseResvColumnMetaData = 1
                  return


            end if

            ! if we're here, then keyword not assigned yet
            ! let's look to the parent blocks

            ParseResvColumnMetaData = ParseTB0ColumnMetaData 
     &                   (header%tb0cmd,keyword, keyLen, subString)
            return

      END FUNCTION ParseResvColumnMetaData




C*******************************************************************
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
      INTEGER FUNCTION ParseResvinColumnMetaData(header, keyword,
     &                                           keyLen, subString)
            TYPE(ResvinColumnMetaData), INTENT(INOUT) :: header 
            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value
            character*32768 tmpString
            integer lineLen, wordCount,n, ideallocate

            tmpString = subString

            if(keyword(1:KeyLen) .eq. ':value1')then
                  wordCount = CountWords(tmpString)
                  if(wordCount.le. 0)then
                        ParseResvinColumnMetaData = -1
                        return
                  end if

                  header%tb0cmd%colCount = wordCount
                  deallocate(header%colValue1,stat=ideallocate)
                  allocate(header%colValue1(wordCount))

                  read(tmpString,*)(header%colValue1(n),n=1,wordCount)

C just checking
C                 do n=1,wordCount
C                       print*,header%colValue1(n)
C                 enddo

                  ParseResvinColumnMetaData = 1
                  return

            end if

            ! if we're here, then keyword not assigned yet
            ! let's look to the parent blocks

            ParseResvinColumnMetaData = ParseTB0ColumnMetaData( 
     &                      header%tb0cmd,keyword, keyLen, subString)
            return

      END FUNCTION ParseResvinColumnMetaData



C*******************************************************************
C
C
C
C
      SUBROUTINE LoadAttributeData(epHeader, xCount, yCount, unitNum)
            TYPE(EnSimParam), INTENT(INOUT) :: epHeader
            INTEGER xCount, yCount, unitNum
            INTEGER index, valCount, n, xi, yi, ios
            logical frameHeaderFlag
            character*12  line
 
            valCount = xCount*yCount

!     rev. 10.2.12 Dec.  30/17  - NK: Added frame headers to static r2c files incl. shd file   
            read(unitNum,32000)line
c            print*,1,line
32000       format(a12)
            if(line(1:1).eq.'#')then
                frameHeaderFlag=.true.
                backspace unitNum
            else
                frameHeaderFlag=.false.
                backspace unitNum
            endif
              
            do n=1, epHeader%attCount
                  ALLOCATE( epHeader%attList(n)%val(valCount))
                  
                  if(frameHeaderFlag)read(unitNum,*)line                
                  
                  do yi=1,yCount
                        index = (yi-1)*xCount
                        read(unitNum,*,iostat=ios)
     &                  (epHeader%attList(n)%val(index+xi),xi=1,xCount)
c     print*,(epHeader%attList(n)%val(index+xi),xi=1,xCount)
c      pause
                  end do
            end do
                  
      END SUBROUTINE LoadAttributeData


C*******************************************************************
C
C
      LOGICAL FUNCTION IsFileTypeR2C(fName)
            
            character*(*) :: fName 
            character*256 tmpString
            integer strLength
            logical rstat

            tmpString = fName
      
            IsFileTypeR2C = .false.

            rstat = ToLowerCase(tmpString)
            strLength = LEN_TRIM(tmpString)
      
            if(tmpString(strLength-2:strLength) .eq. 'r2c') then
                  IsFileTypeR2C = .true.
            endif

            return

      END FUNCTION IsFileTypeR2C

C*******************************************************************
C
C
      LOGICAL FUNCTION IsFileTypeTB0(fName)
            
            character*(*) :: fName 
            character*256 tmpString
            integer strLength
            logical rstat

            tmpString = fName
      
            IsFileTypeTB0 = .false.

            rstat = ToLowerCase(tmpString)
            strLength = LEN_TRIM(tmpString)
      
            if(tmpString(strLength-2:strLength) .eq. 'tb0') then
                  IsFileTypeTB0 = .true.
            endif

            return

      END FUNCTION IsFileTypeTB0

C*******************************************************************
C
C
      LOGICAL FUNCTION GetBoolean(string)
            
            character*(*) :: string 
            character*256 tmpString
            integer strLength
            logical rstat

            tmpString = string
      
            GetBoolean = .false.

            rstat = ToLowerCase(tmpString)
            strLength = LEN_TRIM(tmpString)
      
            if((tmpString(1:1) .eq. 'y') .or.
     &                               (tmpString(1:1) .eq. 't'))  then
                  GetBoolean = .true.
            endif

            return

      END FUNCTION GetBoolean



C*******************************************************************
C
C
C
C
      SUBROUTINE GoToStartOfData(unitNum)
            INTEGER unitNum, ios, lineLen, wordCount, keyLen
            logical rstat, foundEndHeader
            character*32768 line, subString
            character*128 keyword

            foundEndHeader = .false.
            line = '#'
            do WHILE((.NOT.foundEndHeader) .AND.
     &          ((line(1:1) .eq. '#') .OR.
     &            (line(1:1) .eq. ':') .OR.
     &            (LEN_TRIM(line) .eq. 0)))     

                  read(UNIT=unitNum, FMT='((A))', iostat=ios) line      ! read a line

                  rStat = Detab(line)                       ! replace tabs with spaces
                  line = ADJUSTL(line)          ! Get rid of leading white space
                  lineLen = LEN_TRIM(line)            ! Find the length excluding trailing spaces
                  if(line(1:1) .eq. ':')then
                        wordCount = SplitLine(line, keyword, subString) ! find the keyword
                        rStat = ToLowerCase(keyword)
                        KeyLen = LEN_TRIM(keyword)

                        if(keyword(1:KeyLen) .eq. ':endheader')then
                              foundEndHeader = .TRUE.
                        endif
                  end if
            end do
                  
      END SUBROUTINE GoToStartOfData




C*******************************************************************
C
C
      INTEGER FUNCTION CountDataLinesAfterHeader(unitNum)
      
      INTEGER unitNum, ios, lineLen
      character*256 line
      logical rstat

      CountDataLinesAfterHeader = 0

      CALL GoToStartOfData(unitNum)

C Now count the data lines

      DO WHILE(.NOT.EOF(unitNum))

         read(unit=unitNum, FMT='((A))', iostat=ios) line   ! read a line
         rStat = Detab(line)                    ! replace tabs with spaces
         line = ADJUSTL(line)       ! Get rid of leading white space
         lineLen = LEN_TRIM(line)         ! Find the length excluding trailing spaces

         if(lineLen.gt.0 .and. line(1:1) .ne. '#')then
            CountDataLinesAfterHeader=CountDataLinesAfterHeader+1
         end if

      END DO

      END FUNCTION CountDataLinesAfterHeader


C*******************************************************************
C
C
C           Return Value
C           -1 = Problem
C           0 = keyword not foung in Type
C         1 = Successfully assigned
C
C
      INTEGER FUNCTION ParseEventParam(event, keyword,
     &                                           keyLen, subString)

!           added for read_event  by -nk- Jan. 20/08
!           based on read_shed_ef

            TYPE(EventParam), INTENT(INOUT) :: event 

c           TYPE(EventParam), INTENT(INOUT) :: header 

            INTEGER, INTENT(IN) :: keyLen
            CHARACTER*(*), INTENT(INOUT) :: keyword, subString
            character*128  value, rootFileName
            character*32768 tmpString
            integer lineLen, wordCount
!     parse for the numerical values
            if(keyword(1:keyLen) .eq. ':year')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%year1
                        ParseEventParam = 1
                  end if
                  return
            else if(keyword(1:KeyLen) .eq. ':month')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%mo1
                        ParseEventParam = 1
                  end if
                  return
            else if(keyword(1:KeyLen) .eq. ':day')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%day1
                        ParseEventParam = 1
                  end if
                  return
            else if(keyword(1:KeyLen) .eq. ':hour')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%hour1
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':intsoilmoisture')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%smc_init
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':rainconvfactor')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%conv
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen).eq.':eventprecipscalefactor')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%scale
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':precipscalefactor')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%readscale
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':eventsnowscalefactor')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%scalesnw
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':snowscalefactor')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%readscalesnw
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':eventtempscalefactor')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%scaletem
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':tempscalefactor')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%readscaletem
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':disaggregate')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%disaggregate
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':hoursraindata')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%nhg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':hoursflowdata')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%nhf
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':deltat_report')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%deltat_report
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':spinupevents')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%spinupevents
                        ParseEventParam = 1
                  end if
!     parse for the flags
            else if(keyword(1:KeyLen) .eq. ':snwflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%snwflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':sedflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%sedflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':vapflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%vapflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':smrflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%smrflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':resinflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%resinflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':tbcflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%tbcflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':resumflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%resumflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':contflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%contflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':routeflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%routeflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':crseflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%crseflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':ensimflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%ensimflg
                        ParseEventParam = 1
                  end if
!     rev. 9.8.12  Dec.  08/11  - NK: recognize kenueflg in the event file 
!                 in the rest of the code ensimflg is used
            else if(keyword(1:KeyLen) .eq. ':kenueflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
c                        read(value, *)    event%ensimflg
                        read(value, *)    event%kenueflg
                        event%ensimflg=event%kenueflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':picflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%picflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':wetflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%wetflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':modelflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%modelflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':shdflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%shdflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':trcflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%trcflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':frcflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%frcflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':initflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%initflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':grdflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%grdflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':ntrlflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%ntrlflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':nudgeflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%nudgeflg
                        ParseEventParam = 1
                  end if
!     rev. 10.4.42 Sep.  29/21  = NK Added flowresetflg to the event file
            else if(keyword(1:KeyLen) .eq. ':flowresetflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%flowresetflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':resetflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%resetflg
                        ParseEventParam = 1
                  end if
!     rev. 9.8.47  Feb.  04/13  - NK: Headers added for spl & resin csv files
!                 Note: keyword & program variable name not the same
            else if(keyword(1:KeyLen) .eq. ':hdrflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%hdrflg0
                        ParseEventParam = 1
                  end if
!     rev. 9.8.53  Mar.  20/13  - NK: Add Lake St. Joseph diversion algorithm to REROUT.f
            else if(keyword(1:KeyLen) .eq. ':divertflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%divertflg
                        ParseEventParam = 1
                  end if
!     rev. 9.8.81  Sep.  03/13  - NK: Add pafflg and update precip adjustment factors PAF!***********************************************************************
            else if(keyword(1:KeyLen) .eq. ':pafflg')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%pafflg
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':fliflg')then
!     rev. 9.8.93  Nov.  12/13  - NK: Added the routing initialization with yyyymmdd_fli.r2c
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%fliflg
                        ParseEventParam = 1
                  end if
           else if(keyword(1:KeyLen) .eq. ':lakeflg')then
!     rev. 9.9.47  Dec.  23/14  - NK: Added lakeflg for lake evaporation option
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%lakeflg
                        ParseEventParam = 1
                  end if
                  
           else if(keyword(1:KeyLen) .eq. ':iceflg')then
!     rev. 10.1.07 Dec.  02/15  - NK: Added ice_fctr(n) to route 
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%iceflg
                        ParseEventParam = 1
                  end if
                  
           else if(keyword(1:KeyLen) .eq. ':icerivflg')then
!     rev. 10.1.07 Dec.  02/15  - NK: Added ice_fctr(n) to route 
!     rev. 10.1.44 Dec.  02/15  - NK: Reworked icerivflg & icelakeflg
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%icerivflg
                        ParseEventParam = 1
                  end if
                  
           else if(keyword(1:KeyLen) .eq. ':icelakeflg')then
!     rev. 10.1.44 Dec.  02/15  - NK: Reworked icerivflg & icelakeflg
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%icelakeflg
                        ParseEventParam = 1
                  end if
                  
           else if(keyword(1:KeyLen) .eq. ':tb0flg')then
!     REV. 10.1.41 Oct   11/16  - NK: Added tb0flg to write lake_*.tb0 files
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%tb0flg
                        ParseEventParam = 1
                  end if
                  
           else if(keyword(1:KeyLen) .eq. ':xmlflg')then
!     rev. 10.1.64 Jan.  26/17  - NK: Added XML output file 
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%xmlflg
                        ParseEventParam = 1
                  end if
           else if(keyword(1:KeyLen) .eq. ':nbsflg')then
!     rev. 10.2.28 Jul.  08/18  - NK: Revised for TSW = NBS without lake precip/evap'n
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%nbsflg
                        ParseEventParam = 1
                  end if
           else if(keyword(1:KeyLen) .eq. ':fewsflg')then
!     rev. 10.2.28 Jul.  08/18  - NK: Revised for TSW = NBS without lake precip/evap'n
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%fewsflg
                        ParseEventParam = 1
                  end if
           else if(keyword(1:KeyLen) .eq. ':netcdfoutflg')then
!     rev. 10.2.28 Jul.  08/18  - NK: Revised for TSW = NBS without lake precip/evap'n
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%netcdfoutflg
                        ParseEventParam = 1
                  end if
           else if(keyword(1:KeyLen) .eq. ':writestr')then
!     rev. 10.2.28 Jul.  08/18  - NK: Revised for TSW = NBS without lake precip/evap'n
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%writestr
                        ParseEventParam = 1
                  end if
                  
!     parse for the input file names
            else if(keyword(1:KeyLen) .eq. ':basinfilename')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                    if(ExtractFileNameRoot(value,rootFileName).eq.0)then 
                              ParseEventParam = -1
                    else
                              read(rootFileName, *)   event%infln(1)
                              ParseEventParam = 1
                    endif
                  end if
            else if(keyword(1:KeyLen) .eq. ':parfilename')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                    if(ExtractFileNameRoot(value,rootFileName).eq.0)then 
                              ParseEventParam = -1
                    else
                              read(rootFileName, *)   event%infln(2)
                              ParseEventParam = 1
                    endif
                  end if
            else if(keyword(1:KeyLen) .eq. ':pointdatalocations')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(3)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':pointprecip')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(5)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':streamflowdatafile')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(6)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':reservoirreleasefile')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(7)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':reservoirinflowfile')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(8)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':radarfile')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(9)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':griddedrainfile')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(10)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':rawradarfile')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(11)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':clutterfile')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(12)
                        ParseEventParam = 1
                  end if
!     rev. 10.2.33 Sep.  14/18  - NK: Changed unit=42  fln(12) from clutter to model\*.r2c file 
            else if(keyword(1:KeyLen) .eq. ':griddedmodelfile')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(12)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen).eq.':snowcoverdepletioncurve')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(13)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':pointtemps')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(14)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen).eq.':griddedtemperaturefile')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(15)
                        ParseEventParam = 1
                  end if

!   16,17,18 not used  max,min temp & daily snow

            else if(keyword(1:KeyLen) .eq. ':griddednetradiation')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(19)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':pointnetradiation')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(20)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':griddedhumidity')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(21)
                        ParseEventParam = 1
                  end if
c            else if(keyword(1:KeyLen) .eq. ':griddedwind')then
            else if(keyword(1:14) .eq. ':griddedwind   ')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
!                        do nothing but keep as old file have this keyword                  
!                        keep this griddedwind
c                        read(value, *)    event%infln(22)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':griddedlongwave')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(23)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':griddedshortwave')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(24)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':griddedatmpressure')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(25)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':pointhumidity')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(26)
                        ParseEventParam = 1
                  end if
c            else if(keyword(1:KeyLen) .eq. ':pointwind')then
            else if(keyword(1:12) .eq. ':pointwind   ')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
!                        do nothing but keep as old file have this keyword                  
!                        keep this pointwind
c                        read(value, *)    event%infln(27)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':pointlongwave')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(28)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':pointshortwave')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(29)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':pointatmpressure')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(30)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':griddedrunoff')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(31)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':griddedrecharge')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(32)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':griddedleakage')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(33)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':griddedsnowfile')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(34)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':snowcoursefile')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(35)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':griddedinitsnowweq')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(36)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen).eq.':griddedinitsoilmoisture')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(37)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':griddedinitlzs')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(38)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':pointsoilmoisture')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(39)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':waterqualitydatafile')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(40)
                        ParseEventParam = 1
                  end if
!     rev. 9.9.04  Dec.  17/13  - NK: Change over to gridded clamate normals fo diff
            else if(keyword(1:KeyLen) .eq. ':climateNormalsDiff')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(40)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':channelparfile')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(41)
                        ParseEventParam = 1
                  end if

!  42 & 43 not used

            else if(keyword(1:KeyLen) .eq. ':pointsnow')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(44)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':pointdrain')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(45)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':pointdsnow')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(46)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':griddedsnow')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(47)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':griddeddrain')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(48)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':griddeddsnow')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(49)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':initlakelevel')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(50)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':griddedevaporation')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(51)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':diversionflowfile')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(52)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':observedlakelevel')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(53)
                        ParseEventParam = 1
                  end if
!     rev. 10.1.81 May   05/17  - NK: Added snowg\yyyymmdd_swe.tb0 obs. swe
            else if(keyword(1:KeyLen) .eq. ':observedswe')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(54)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':noeventstofollow')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%nch

                        ParseEventParam = 1
                  end if
!     rev. 9.8.93  Nov.  12/13  - NK: Added the routing initialization with yyyymmdd_fli.r2c
            else if(keyword(1:KeyLen) .eq. ':flowinit')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(55)
                        ParseEventParam = 1
                  end if
!     rev. 9.9.00  Dec.  08/13  - NK: Added Lake Evaporation model
! unit=286  fln(56)- point wind speed               yyyymmdd_spd.tb0
! unit=287  fln(57)- point wind direction           yyyymmdd_dir.tb0
! unit=288  fln(58)- gridded wind speed             yyyymmdd_spd.r2c
! unit=289  fln(59)- gridded wind direction         yyyymmdd_dir.r2c
            else if(keyword(1:KeyLen) .eq. ':pointwindspd')then
c            else if(keyword(1:13) .eq. ':pointwindspd')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(56)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':pointwinddir')then
c            else if(keyword(1:13) .eq. ':pointwinddir')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(57)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':griddedwindspd')then
c            else if(keyword(1:15) .eq. ':griddedwindspd')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(58)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':griddedwinddir')then
c            else if(keyword(1:15) .eq. ':griddedwinddir')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(59)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen) .eq. ':climatenormalsdiff')then
!     rev. 9.9.04  Dec.  17/13  - NK: Change over to gridded clamate normals to diff
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(60)
                        ParseEventParam = 1
                  end if
!                 Note: fln(61) is in use for swe                  
!     rev. 9.9.06  Jan.  08/14  - NK: Add daily differences to Harfreaves ETHarg.f
            else if(keyword(1:KeyLen).eq.':griddeddailydifference')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(62)
                        ParseEventParam = 1
                  end if

!     rev. 10.2.11 Dec.  18/17  - NK: 4 files added for BLEND.exe    
            else if(keyword(1:KeyLen).eq.':pointhourlyprecip')then
                 if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(63)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen).eq.':pointdailyprecip')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(64)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen).eq.':griddedhourlyprecip')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(65)
                        ParseEventParam = 1
                  end if
            else if(keyword(1:KeyLen).eq.':griddeddailyprecip')then
                  if(SplitLine(subString, value, tmpString) .eq. 0)then
                        ParseEventParam = -1
                  else
                        read(value, *)    event%infln(68)
                        ParseEventParam = 1
                  end if

            else
                  print*
                  print*,keyword(1:KeyLen)
                  print*,'This key word is not recognized'
                  print*,'Please update this program'
                  print*,'to the latest version'
                  print*
                  print*,'If this exec is an older version & you would'
                  print*,'like to use it on a newer event file that has'
                  print*,'keywords that have been added since this exec'
                  print*,'was built, please edit the event files to'
                  print*,'delete these unrecognized lines'
                  print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                  print*
                  stop 'Program aborted in EF_module @ 2763'
                  return      
            end if

c           ! if we're here, then keyword not assigned yet
c           ! let's look to the parent blocks
c
c           ParseShedParam = ParseR2CParam(header%r2cp, keyword,
c     &                                          keyLen, subString)
            return

      END FUNCTION ParseEventParam



      END MODULE EF_Module


