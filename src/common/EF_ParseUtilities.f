	MODULE EF_ParseUtilities
      
!***********************************************************************
!    Copyright (C) 2006 by Dave Watson NCR 
        
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
      
      
      
	CONTAINS

C***********************************************************************
	integer FUNCTION SplitLine(line, keyword, value)
C***********************************************************************
C Extract the first white space delimited substring from line and place
C it into keyword. Put the rest of line into value.
Ccase

C Do NOT change the incoming line
C
C Return:	0 = the line is empty
C		1 = there is ONLY a keyword
C		2 = there is a keyword and a value
C***********************************************************************

	character*(*) line
	character*(*) keyword, value
c	character*4096 localLine
c	character*32768 localLine
	character*98304 localLine
	integer I, DQ, lineLen
	logical rStat

	SplitLine = 0
	localLine = line			! work on local copy of the string
					
	rStat = Detab(localLine)				! change any tabs to spaces
	localLine = ADJUSTL(localLine)			! get rid of leading spaces
	lineLen = LEN_TRIM(localLine)	! find the new length, excluding trailing spaces

	if(lineLen .eq. 0) then ! we have an empty line
		SplitLine = 0
		return
	end if
			

	I = INDEX(localLine(1:lineLen),' ')	! Find the first space character
	
C Deal with possibility of quotes (strings)
	if(localLine(1:1) .eq. '"')then
		DQ = INDEX(localLine(2:lineLen),'"')+1
		localLine(1:1) = CHAR(32) !replace first occurance of quote with a space
		if(DQ.GT.1) then
			localLine(DQ:DQ) = CHAR(32) !replace 2nd occurance of quote with a space
				
C			localLine = ADJUSTL(localLine)			! get rid of leading spaces
C			lineLen = LEN_TRIM(localLine)
			if(DQ.EQ.lineLen) then
				keyword = localLine(2:lineLen-1)
				value = ''
				SplitLine = 1
				return
			else 
				keyword = localLine(2:DQ-1)
				value = localLine(DQ+1:lineLen)
				value = ADJUSTL(value)
				SplitLine = 2
			end if
			return
		end if
	end if


	if((I.EQ.0) .AND. (lineLen.GT.0))then	! Single word only
		keyword = localLine(1:lineLen)
		value = ''
		SplitLine = 1
	else
		keyword = localLine(1:I)
		value = localLine(I+1:lineLen)
		value = ADJUSTL(value)
		SplitLine = 2
	end if
	
	return
	end FUNCTION


C***********************************************************************
	integer FUNCTION ExtractFileNameRoot(value, filename)
C***********************************************************************
C Extract the root file name from a string that may contain a file path
C
C Return:	0 = the line is empty
C		1 = successful
C
C***********************************************************************

	character*(*) value, filename
	character*4096 localValue
	integer I, DQ, lineLen
	logical rStat

	ExtractFileNameRoot = 0
	localValue = value			! work on local copy of the string
					
	rStat = Deslash(localValue)				! change any slashes to spaces
	localValue = ADJUSTL(localValue)			! get rid of leading spaces
	lineLen = LEN_TRIM(localValue)	! find the new length, excluding trailing spaces

	if(lineLen .eq. 0) then ! we have an empty line
		ExtractFileNameRoot = 0
		return
	end if
			

	I = INDEX(localValue(1:lineLen),' ',BACK = .TRUE.)	! Find the last space character
	

	if((I.EQ.0) .AND. (lineLen.GT.0))then	! Single word only
		filename = localValue(1:lineLen)
		ExtractFileNameRoot = 1
	else
		filename = localValue(I:lineLen)
		filename = ADJUSTL(filename)
		ExtractFileNameRoot = 1
	end if
	
	return
	end FUNCTION





C***********************************************************************
	integer FUNCTION CountWords(line)
C***********************************************************************
C
C Return:	the number of words in the line
C***********************************************************************

	character*(*) line
c	character*4096 localLine
c	character*32768 localLine
	character*98304 localLine
	integer I, DQ, lineLen
	logical rStat

	localLine = line			! work on local copy of the string
					
	rStat = Detab(localLine)				! change any tabs to spaces
	localLine = ADJUSTL(localLine)			! get rid of leading spaces
	lineLen = LEN_TRIM(localLine)	! find the new length, excluding trailing spaces

	CountWords = 0

	if(lineLen .eq. 0) 	return

	I = 1
	do while(I.ne.0)
		CountWords = CountWords+1
		I = INDEX(localLine(1:lineLen),' ')	! Find the first space character
		localLine = localLine(I+1:lineLen)
		localLine = ADJUSTL(localLine)
		lineLen = LEN_TRIM(localLine)
	end do
	
	return
	end FUNCTION




C***********************************************************************
	logical FUNCTION Detab(line)
C***********************************************************************
C Replace all tabs in line with space characters.
C
C Return:	= .false. no change occured
C		= .true. the string was changed
C***********************************************************************

	character*(*) line
	character*1 tab, space
	integer I

	tab = CHAR(9)		! What we are looking for
	space = CHAR(32)	! the substitution value

	Detab = .false.
	do I=1,LEN_TRIM(line)
		if(line(I:I) .eq. tab)then
			line(I:I) = space
			Detab = .true.
		end if
	end do

	return
	end FUNCTION

C***********************************************************************
	logical FUNCTION Deslash(line)
C***********************************************************************
C Replace all slashes (back and forward) in line with space characters.
C
C Return:	= .false. no change occured
C		= .true. the string was changed
C***********************************************************************

	character*(*) line
	character*1 space, backslash,forwardslash 
	integer I

	space = CHAR(32)
	backslash = CHAR(28)
	forwardslash = CHAR(47)

	Deslash = .false.
	do I=1,LEN_TRIM(line)
		if(line(I:I) .eq. backslash)then
			line(I:I) = space
			Deslash = .true.
		else if(line(I:I) .eq. forwardslash)then
			line(I:I) = space
			Deslash = .true.
		end if
	end do

	return
	end FUNCTION



C***********************************************************************
	logical FUNCTION ToLowerCase(line)
C***********************************************************************
C Replace all uppercase character in line with lowercase characters.
C
C Return:	= .false. no change occured
C		= .true. the string was changed
C***********************************************************************

	character*(*) line
	integer I

	ToLowerCase = .false.
	do I=1,LEN_TRIM(line)
		if((line(I:I) .GE. 'A') .AND. (line(I:I) .LE. 'Z'))then
			line(I:I) = CHAR(ICHAR(line(I:I)) + 32)
			ToLowerCase = .true.
		end if
	end do

	return
	end FUNCTION


  
C***********************************************************************
	logical FUNCTION IsNumber(str)
C***********************************************************************
C Check if the string argument is a number
C
C Return:	= .false. not a number
C		= .true. is a number
C***********************************************************************

	character *(*) str 
	integer i  

      IsNumber=.true.
	
      do i=1,LEN_TRIM(str)
		if((str(i:i).lt.'0'.or.str(i:i).gt.'9').and.
     $		str(i:i).ne.'.'.and.str(i:i).ne.'-'.and.
     $		str(i:i).ne.'+')then
			IsNumber=.false.
			return
		end if
      end do

      return  
      end FUNCTION




C***********************************************************************
	INTEGER FUNCTION JDATE (YEAR,MONTH,DAY)
C***********************************************************************
C Conversion from a Gregorian calendar date to a Julian date. Valid for
C any Gregorian calendar date producing a Julian date greater than zero
C
C COMPUTES THE JULIAN DATE (JDATE) GIVEN A GREGORIAN CALENDAR DATE (YR,MTH,DAY).
C***********************************************************************
      INTEGER YEAR,MONTH,DAY,I,J,K
C
      I= YEAR
      J= MONTH
      K= DAY
C
      JDATE= K-32075+1461*(I+4800+(J-14)/12)/4+367*(J-2-(J-14)/12*12)
     &    /12-3*((I+4900+(J-14)/12)/100)/4
C
      RETURN 
      END FUNCTION

 

C***********************************************************************
	SUBROUTINE GDATE (JD, YEAR,MONTH,DAY)
C***********************************************************************
C Conversion from a Julian date to a Gregorian calendar date.
C
C COMPUTES THE GREGORIAN CALENDAR DATE (YEAR,MONTH,DAY) GIVEN THE JULIAN DATE (JD).
C***********************************************************************
      INTEGER JD,YEAR,MONTH,DAY,I,J,K,L,N
C
      L= JD+68569
      N= 4*L/146097
      L= L-(146097*N+3)/4
      I= 4000*(L+1)/1461001
      L= L-1461*I/4+31
      J= 80*L/2447
      K= L-2447*J/80
      L= J/11
      J= J+2-12*L
      I= 100*(N-49)+I+L
C
      YEAR= I
      MONTH= J
      DAY= K
C
      RETURN
      END SUBROUTINE


	END MODULE



