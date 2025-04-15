        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar 27 10:54:11 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ELEM8__genmod
          INTERFACE 
            SUBROUTINE ELEM8(XP,YP,ZP,XQ,YQ,ZQ,H,G,CP,AA,QN)
              REAL(KIND=4) :: XP
              REAL(KIND=4) :: YP
              REAL(KIND=4) :: ZP
              REAL(KIND=4) :: XQ(8)
              REAL(KIND=4) :: YQ(8)
              REAL(KIND=4) :: ZQ(8)
              COMPLEX(KIND=4) :: H(8)
              COMPLEX(KIND=4) :: G(8)
              REAL(KIND=4) :: CP
              REAL(KIND=4) :: AA(8)
              REAL(KIND=4) :: QN(3)
            END SUBROUTINE ELEM8
          END INTERFACE 
        END MODULE ELEM8__genmod
