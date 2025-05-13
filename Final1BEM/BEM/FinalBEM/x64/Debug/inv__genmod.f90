        !COMPILER-GENERATED INTERFACE MODULE: Mon May 12 15:23:19 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INV__genmod
          INTERFACE 
            SUBROUTINE INV(A,N,INFO)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(INOUT) :: A(N,N)
              INTEGER(KIND=4), INTENT(OUT) :: INFO
            END SUBROUTINE INV
          END INTERFACE 
        END MODULE INV__genmod
