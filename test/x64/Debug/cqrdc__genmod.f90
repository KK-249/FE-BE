        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 23 16:40:50 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CQRDC__genmod
          INTERFACE 
            SUBROUTINE CQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,JOB)
              INTEGER(KIND=4) :: LDX
              COMPLEX(KIND=4) :: X(LDX,1)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: P
              COMPLEX(KIND=4) :: QRAUX(1)
              INTEGER(KIND=4) :: JPVT(1)
              COMPLEX(KIND=4) :: WORK(1)
              INTEGER(KIND=4) :: JOB
            END SUBROUTINE CQRDC
          END INTERFACE 
        END MODULE CQRDC__genmod
