      SUBROUTINE RNODE
C
C---- READ IN NODAL POINT DATA
C
      COMMON /CONTR/NNODE,NELEM
      COMMON /NODES/X(2000),Y(2000)
      READ(5,*) NNODE,NREC
      DO 30 NR=1,NREC
      READ(5,*) N1,N2,X1,Y1,X2,Y2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     NNODE =NUMBER OF NODES                                           C
C     NREC =NUMBER OF NODAL POINT DATA RECORDS                         C
C     N1,N2 =NODAL POINT NUMBERS FOR THE FIRST AND THE LAST NODES      C
C     X1,Y1 =COORDINATES OF NODE #N1                                   C
C     X2,Y2 =COORDINATES OF NODE #N2                                   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(N1.EQ.N2) THEN
        X(N1)=X1
        Y(N1)=Y1
      ELSE
        DN=N2-N1
        DX=(X2-X1)/DN
        DY=(Y2-Y1)/DN
        XX=X1-DX
        YY=Y1-DY
        DO 20 I=N1,N2
          XX=XX+DX
          YY=YY+DY
          X(I)=XX
          Y(I)=YY
20      CONTINUE
      END IF
30    CONTINUE
C
C---- PRINT OUT NODAL POINT COORDINATES
C
      WRITE(6,100)
100   FORMAT(/2X,'NODAL POINT COORDINATES:'/)
      DO 102 I=1,NNODE
        WRITE(6,101) I,X(I),Y(I)
101     FORMAT(5X,'NODE #',I4,5X,'X=',F11.3,3X,'Y=',F11.3)
102   CONTINUE
      RETURN
      END
