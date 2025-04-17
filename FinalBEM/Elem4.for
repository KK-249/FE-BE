      SUBROUTINE ELEM4(XP,YP,ZP,XQ,YQ,ZQ,H,G,CP,AA)
C
C  Formulate element coefficient matrices
C
      COMMON /FRE/FREQ,WN

      DIMENSION XII(12),WT(12)
      DIMENSION XQ(4),YQ(4),ZQ(4),QN(3),PSI(8),DPSI(2,8)
      DIMENSION H(4),G(4),AA(4)
      REAL NN1,NN2,NN3
      COMPLEX H,G,GREEN,DGDN,IK

      C1=1.0/12.56637061
      IK=CMPLX(0.,1.)*WN

      DO 100 I=1,4
          H(I)=0.0
          G(I)=0.0
          AA(I)=0
100   CONTINUE
      CALL GETINT(NINP,XII,WT)
C
C  Integration loop
C
      DO 450 INP1=1,NINP
	 DO 400 INP2=1,NINP
	  XI1=XII(INP1)
	  XI2=XII(INP2)
          CALL SHAPE(XI1,XI2,PSI,DPSI,4)
          XX=0.0
          YY=0.0
	  ZZ=0.0
          DXDS1=0.0
	  DXDS2=0.0
          DYDS1=0.0
	  DYDS2=0.0
	  DZDS1=0.0
	  DZDS2=0.0
	  DO 200 I=1,4
            XX=XX+XQ(I)*PSI(I)
            YY=YY+YQ(I)*PSI(I)
	    ZZ=ZZ+ZQ(I)*PSI(I)
            DXDS1=DXDS1+XQ(I)*DPSI(1,I)
	    DYDS1=DYDS1+YQ(I)*DPSI(1,I)
	    DZDS1=DZDS1+ZQ(I)*DPSI(1,I)
	    DXDS2=DXDS2+XQ(I)*DPSI(2,I)
            DYDS2=DYDS2+YQ(I)*DPSI(2,I)
	    DZDS2=DZDS2+ZQ(I)*DPSI(2,I)		 
200       CONTINUE
          NN1=DYDS1*DZDS2-DZDS1*DYDS2
	  NN2=-(DXDS1*DZDS2-DZDS1*DXDS2)
	  NN3=DXDS1*DYDS2-DYDS1*DXDS2
          DETJ=SQRT(NN1**2+NN2**2+NN3**2)
	  QN(1)=NN1/DETJ
	  QN(2)=NN2/DETJ
	  QN(3)=NN3/DETJ
          RX=XX-XP
          RY=YY-YP
	  RZ=ZZ-ZP
          R=SQRT(RX*RX+RY*RY+RZ*RZ)
          DRDN=(RX*QN(1)+RY*QN(2)+RZ*QN(3))/R
          FACTOR=DETJ*WT(INP1)*WT(INP2)

          GREEN=C1*CEXP(-IK*R)/R*FACTOR
          DGDN=-(1.0+IK*R)*GREEN/R*DRDN
          DGDN0=-C1/(R*R)*DRDN*FACTOR

          DO 300 I=1,4
            H(I)=H(I)+PSI(I)*DGDN
            G(I)=G(I)+PSI(I)*GREEN
            AA(I)=AA(I)+PSI(I)*FACTOR
300       CONTINUE
        CP=CP-DGDN0 
400     CONTINUE
450   CONTINUE
      RETURN
      END