        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 28 19:35:43 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SING4__genmod
          INTERFACE 
            SUBROUTINE SING4(XP,YP,ZP,XQ,YQ,ZQ,H,G,CP,ISING,AA)
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
              REAL(KIND=4) :: AA(8)
            END SUBROUTINE SING4
          END INTERFACE 
        END MODULE SING4__genmod
