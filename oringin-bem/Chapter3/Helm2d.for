      PROGRAM HELM2D
C
C     2-D Direct BEM Code for the Helmholtz equation
C     written by T. W. Wu   timwu@engr.uky.edu
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
