        !COMPILER-GENERATED INTERFACE MODULE: Wed Oct 23 16:40:50 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CQRSL__genmod
          INTERFACE 
            SUBROUTINE CQRSL(X,LDX,N,K,QRAUX,Y,QY,QTY,B,RSD,XB,JOB,INFO)
              INTEGER(KIND=4) :: LDX
              COMPLEX(KIND=4) :: X(LDX,1)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              COMPLEX(KIND=4) :: QRAUX(1)
              COMPLEX(KIND=4) :: Y(1)
              COMPLEX(KIND=4) :: QY(1)
              COMPLEX(KIND=4) :: QTY(1)
              COMPLEX(KIND=4) :: B(1)
              COMPLEX(KIND=4) :: RSD(1)
              COMPLEX(KIND=4) :: XB(1)
              INTEGER(KIND=4) :: JOB
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE CQRSL
          END INTERFACE 
        END MODULE CQRSL__genmod
