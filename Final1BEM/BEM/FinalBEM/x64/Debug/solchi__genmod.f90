        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  8 19:57:50 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLCHI__genmod
          INTERFACE 
            SUBROUTINE SOLCHI(A,NR,NC,B,AUX,AUX2,IPVT,X,N1,N2)
              INTEGER(KIND=4) :: N2
              INTEGER(KIND=4) :: N1
              COMPLEX(KIND=4) :: A(N1,N2)
              INTEGER(KIND=4) :: NR
              INTEGER(KIND=4) :: NC
              COMPLEX(KIND=4) :: B(N1)
              COMPLEX(KIND=4) :: AUX(N2)
              COMPLEX(KIND=4) :: AUX2(N1)
              INTEGER(KIND=4) :: IPVT(N2)
              COMPLEX(KIND=4) :: X(N2)
            END SUBROUTINE SOLCHI
          END INTERFACE 
        END MODULE SOLCHI__genmod
