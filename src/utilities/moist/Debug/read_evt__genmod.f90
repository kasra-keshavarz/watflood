        !COMPILER-GENERATED INTERFACE MODULE: Sun Oct 22 19:18:51 2017
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_EVT__genmod
          INTERFACE 
            SUBROUTINE READ_EVT(DATE,CONV,SCALE,SMC5,NHG,NHF)
              CHARACTER(LEN=14) :: DATE
              REAL(KIND=4) :: CONV
              REAL(KIND=4) :: SCALE
              REAL(KIND=4) :: SMC5(16)
              INTEGER(KIND=4) :: NHG
              INTEGER(KIND=4) :: NHF
            END SUBROUTINE READ_EVT
          END INTERFACE 
        END MODULE READ_EVT__genmod
