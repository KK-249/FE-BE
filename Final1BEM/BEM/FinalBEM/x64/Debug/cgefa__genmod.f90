        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  8 19:57:49 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CGEFA__genmod
          INTERFACE 
            SUBROUTINE CGEFA(A,LDA,N,IPVT,INFO)
              INTEGER(KIND=4) :: LDA
              COMPLEX(KIND=4) :: A(LDA,1)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IPVT(1)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE CGEFA
          END INTERFACE 
        END MODULE CGEFA__genmod
