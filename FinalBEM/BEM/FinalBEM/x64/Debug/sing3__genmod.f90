        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  8 19:57:50 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SING3__genmod
          INTERFACE 
            SUBROUTINE SING3(XP,YP,ZP,XQ,YQ,ZQ,H,G,CP,ISING)
              REAL(KIND=4) :: XP
              REAL(KIND=4) :: YP
              REAL(KIND=4) :: ZP
              REAL(KIND=4) :: XQ(4)
              REAL(KIND=4) :: YQ(4)
              REAL(KIND=4) :: ZQ(4)
              COMPLEX(KIND=4) :: H(4)
              COMPLEX(KIND=4) :: G(4)
              REAL(KIND=4) :: CP
              INTEGER(KIND=4) :: ISING
            END SUBROUTINE SING3
          END INTERFACE 
        END MODULE SING3__genmod
