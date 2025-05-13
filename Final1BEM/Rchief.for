      SUBROUTINE RCHIEF
C
C---- READ IN CHIEF POINTS
C
      COMMON /CHI/NCH,XCH(3,20)
      WRITE(6,*) ' '
      WRITE(6,*) ' CHIEF POINTS:'
      READ(5,*) NCH
      DO 30 I=1,NCH
        READ(5,*) XCH(1,I),XCH(2,I),XCH(3,I)
        WRITE(6,*) '  X=', XCH(1,I),' Y=', XCH(2,I),' Z=',XCH(3,I)
30    CONTINUE

      RETURN
      END