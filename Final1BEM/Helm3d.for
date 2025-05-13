      PROGRAM Helm3D
C
C     3-D BEM CODE FOR the Helmholtz  Equation
C     written by T. W. Wu    timwu@engr.uky.edu
C
      CHARACTER*15 input
      write(*,*) 'Specify input file name.'
      read(*,100) input
100   FORMAT(A15)      
      OPEN(5,FILE=input,STATUS='OLD')
      OPEN(6,FILE='output.dat')
      CALL PREP
      CALL PROC
      STOP
      END
