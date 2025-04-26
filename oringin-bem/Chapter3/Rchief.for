      SUBROUTINE RCHIEF
C
C---- READ IN CHIEF POINTS
C
      COMMON /CHI/NCH,XCH(20),YCH(20)
      WRITE(6,*) ' '
      WRITE(6,*) ' CHIEF POINTS:'
      READ(5,*) NCH
      DO 30 I=1,NCH
        READ(5,*) XCH(I),YCH(I)
        WRITE(6,*) '  X=', XCH(I),' Y=', YCH(I)
30    CONTINUE

      RETURN
      END