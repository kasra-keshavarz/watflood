        !COMPILER-GENERATED INTERFACE MODULE: Tue Sep  5 14:30:28 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_EVENT__genmod
          INTERFACE 
            SUBROUTINE WRITE_EVENT(DATE,CONV,SCALE,SMC5,MO,NHR,NHF,YY1, &
     &MM1,DD1,HH1,EVTNAME,EVTLENGTH,EVTTYPE)
              CHARACTER(LEN=12) :: DATE
              REAL(KIND=4) :: CONV
              REAL(KIND=4) :: SCALE
              REAL(KIND=4) :: SMC5(5)
              INTEGER(KIND=4) :: MO
              INTEGER(KIND=4) :: NHR
              INTEGER(KIND=4) :: NHF
              CHARACTER(LEN=4) :: YY1
              CHARACTER(LEN=2) :: MM1
              CHARACTER(LEN=2) :: DD1
              CHARACTER(LEN=2) :: HH1
              CHARACTER(LEN=8) :: EVTNAME
              INTEGER(KIND=4) :: EVTLENGTH
              CHARACTER(LEN=8) :: EVTTYPE
            END SUBROUTINE WRITE_EVENT
          END INTERFACE 
        END MODULE WRITE_EVENT__genmod
