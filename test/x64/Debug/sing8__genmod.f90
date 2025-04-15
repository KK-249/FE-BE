        !COMPILER-GENERATED INTERFACE MODULE: Wed Mar 26 15:26:42 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SING8__genmod
          INTERFACE 
            SUBROUTINE SING8(XP,YP,ZP,XQ,YQ,ZQ,H,G,CP,ISING)
              REAL(KIND=4) :: XP
              REAL(KIND=4) :: YP
              REAL(KIND=4) :: ZP
              REAL(KIND=4) :: XQ(8)
              REAL(KIND=4) :: YQ(8)
              REAL(KIND=4) :: ZQ(8)
              COMPLEX(KIND=4) :: H(8)
              COMPLEX(KIND=4) :: G(8)
              REAL(KIND=4) :: CP
              INTEGER(KIND=4) :: ISING
            END SUBROUTINE SING8
          END INTERFACE 
        END MODULE SING8__genmod
