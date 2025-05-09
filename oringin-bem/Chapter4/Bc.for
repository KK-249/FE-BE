      SUBROUTINE BC
      COMMON /BG/BIG
      COMMON /BCC/CA(4,2000),CB(4,2000),CC(4,2000)
      COMMON /CONTR/NNODE,NELEM
      COMMON /ELEM/NODE(4,2000)
      COMMON /SOL/SP(2000),VN(4,2000)
      COMMON /FRE/FREQ,WN,IFR
      COMPLEX BIG,CA,CB,CC,SP,VN
      COMPLEX TM,TP
      LOGICAL ALARM
C
      BIG=CMPLX(1.0E+15,1.0E+15)
      DO 10 I=1,NNODE
        SP(I)=BIG
10    CONTINUE
C
      DO 20 J=1,4
        DO 20 K=1,NELEM
          VN(J,K)=BIG
20    CONTINUE
C
      IF(IFR.EQ.1) CALL RBC
C
      ALARM=.FALSE.
      DO 50 K=1,NELEM
        I34=4
        IF(NODE(3,K).EQ.NODE(4,K)) I34=3
        DO 50 J=1,I34
	  NOD=NODE(J,K)
	  IF(CA(J,K).EQ.BIG.AND.CB(J,K).EQ.BIG) THEN
	    WRITE(6,*) 'B.C. NOT SPECIFIED ON ELEMENT #',K
	    ALARM=.TRUE.
	    GO TO 120
	  END IF
	  IF(CA(J,K).EQ.0.0.AND.CB(J,K).EQ.0.0) THEN
	    WRITE(6,*) 'CA=0 AND CB=0 ON ELEMENT #',K
	    ALARM=.TRUE.
	    GO TO 120
	  END IF
	  IF(CB(J,K).EQ.0.0) THEN
	    TM=SP(NOD)
	    TP=CC(J,K)/CA(J,K)
	    SP(NOD)=TP
	    IF(TM.NE.BIG.AND.TP.NE.TM) THEN
	      WRITE(6,*) 'MORE THAN ONE PRESSURE AT NODE #',NOD
	      ALARM=.TRUE.
	      GO TO 120
	    END IF
          END IF
50    CONTINUE  
C
C     PRINT OUT BC's
C
120   WRITE(6,123)
123   FORMAT(/2X,'BOUNDARY CONDITIONS:'/)
      DO 130 K=1,NELEM
        WRITE(6,1070) K
1070    FORMAT(3X,'ELEMENT NO.',I4)
        I34=4
        IF(NODE(3,K).EQ.NODE(4,K)) I34=3
        DO 130 J=1,I34
          WRITE(6,1080) J,CA(J,K),CB(J,K),CC(J,K)
1080      FORMAT(4X,'LOCAL NODE #',I3,2X,'CA=',2F6.1,2X,
     &           'CB=',2F6.1,2X,'CC=',2F6.1)
130   CONTINUE
      IF(ALARM) STOP
      RETURN
      END