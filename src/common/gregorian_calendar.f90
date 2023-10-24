      subroutine GREGORIAN (JD, YEAR,MONTH,DAY)

!  http://ix23.com/converting-between-julian-dates-and-gregorian-calendar-dates-in-fortran-and-javascript/
!     COMPUTES THE GREGORIAN CALENDAR DATE (YEAR,MONTH,DAY)
!     GIVEN THE JULIAN DATE (JD).
!
      INTEGER JD,YEAR,MONTH,DAY,I,J,K

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

      YEAR= I
      MONTH= J
      DAY= K

      RETURN
    END subroutine     

    
    
    
    
      INTEGER FUNCTION JULIAN(YEAR,MONTH,DAY)

!   COMPUTES THE JULIAN DATE (JD) GIVEN A GREGORIAN CALENDAR
!   DATE (YEAR,MONTH,DAY).

      INTEGER YEAR,MONTH,DAY,I,J,K

      I= YEAR
      J= MONTH
      K= DAY

      julian = K-32075+1461*(I+4800+(J-14)/12)/4+367*(J-2-(J-14)/12*12)&
        /12-3*((I+4900+(J-14)/12)/100)/4

      RETURN
      END