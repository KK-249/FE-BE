        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 23 16:40:50 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CGECO__genmod
          INTERFACE 
            SUBROUTINE CGECO(A,LDA,N,IPVT,RCOND,Z)
              INTEGER(KIND=4) :: LDA
              COMPLEX(KIND=4) :: A(LDA,1)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IPVT(1)
              REAL(KIND=4) :: RCOND
              COMPLEX(KIND=4) :: Z(1)
            END SUBROUTINE CGECO
          END INTERFACE 
        END MODULE CGECO__genmod
