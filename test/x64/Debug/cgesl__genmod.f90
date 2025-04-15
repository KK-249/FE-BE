        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 23 16:40:50 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CGESL__genmod
          INTERFACE 
            SUBROUTINE CGESL(A,LDA,N,IPVT,B,JOB)
              INTEGER(KIND=4) :: LDA
              COMPLEX(KIND=4) :: A(LDA,1)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IPVT(1)
              COMPLEX(KIND=4) :: B(1)
              INTEGER(KIND=4) :: JOB
            END SUBROUTINE CGESL
          END INTERFACE 
        END MODULE CGESL__genmod
