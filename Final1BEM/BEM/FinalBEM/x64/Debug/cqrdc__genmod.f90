        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  8 19:57:49 2025
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
