      SUBROUTINE PROC1

      use iso_c_binding, only: C_PTR, C_NULL_CHAR, C_LOC, C_DOUBLE, C_SIZE_T
      implicit none

     include "engine.h"
     include "matrix.h"

      integer, parameter      :: index = 4000
      integer                 :: f, ierr
      type(C_PTR)             :: ep         ! MATLAB 引擎指针
      type(C_PTR)             :: mxF        ! mxArray 存放频率
      type(C_PTR)             :: mxVn       ! mxArray 存放返回的振速
      double precision        :: femcc( index )
      integer(C_SIZE_T)       :: total

      integer(C_SIZE_T) :: FREQ1,FREQ2,DF
      integer(C_SIZE_T) :: RHO,C
      integer(C_SIZE_T) :: FREQ,WN,IFR
    

100   WN=2.*3.14159*FREQ/C
      write(6,*) ' '
      write(6,*) ' Frequency=',FREQ,' Hz'
      write(6,*) ' Wavenumber=',WN

      CALL FORM
      CALL SOLVE
      CALL SOLOUT
      CALL FIELD
      CALL POWER
      
      ep = engOpen(C_NULL_CHAR)
      if ( c_associated(ep, C_NULL_PTR) ) then
         print *, 'Cannot start MATLAB engine'
         stop 1
      end if
      
      IF(DF.EQ.0.0.OR.FREQ2.EQ.0.0) RETURN

      FREQ=FREQ+DF
      IFR=IFR+1 
          
      IF(FREQ.GT.FREQ2) THEN
        call engClose(ep)
        RETURN
      ELSE
          mxF = mxCreateDoubleMatrix(1,1,mxREAL)                  !在matlab创建m*n的双精度矩阵
         call mxCopyPtrToReal8(C_LOC(f), mxGetPr(mxF), 1)
         call engPutVariable(ep, C_LOC('freq'), mxF)              !mxF矩阵值赋给freq
         
         call engEvalString(ep, C_LOC('femcc = femaddmass(freq);'))
         
             mxVn = engGetVariable(ep, C_LOC('femcc'))
         if ( c_associated(mxVn, C_NULL_PTR) ) then
           print *, 'Failed to get femcc at freq=', freq
           exit
         end if
         call mxCopyPtrToReal8(mxGetPr(mxVn), C_LOC(femcc(1)), total)
      
         call mxDestroyArray(mxF)
         call mxDestroyArray(mxVn)
        CALL BC
        GO TO 100
      END IF
      call engClose(ep)
      END
