        !COMPILER-GENERATED INTERFACE MODULE: Thu Apr 17 09:35:17 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SHAPE__genmod
          INTERFACE 
            SUBROUTINE SHAPE(XI1,XI2,PSI,DPSI,I34)
              REAL(KIND=4) :: XI1
              REAL(KIND=4) :: XI2
              REAL(KIND=4) :: PSI(8)
              REAL(KIND=4) :: DPSI(2,8)
              INTEGER(KIND=4) :: I34
            END SUBROUTINE SHAPE
          END INTERFACE 
        END MODULE SHAPE__genmod
