      SUBROUTINE ELEMT(XP,YP,NL,KINDI,XQ,YQ,H,G,CP)
C
C  Integrate element coefficient vectors
C
      DIMENSION XQ(4),YQ(4),H(4),G(4),QN(2),PSI(4),
     &          DPSI(4),XII(12),WT(12)
      COMMON /FRE/FREQ,WN
C
      COMPLEX H,G,GREEN,DGDN,HANK0,HANK1,C1
      EXTERNAL BESJ0,BESY0,BESJ1,BESY1
C
      C1=-CMPLX(0.,1.)/4.0
      C0=-1.0/6.2831853

      DO 100 I=1,NL
          H(I)=0.0
          G(I)=0.0
100   CONTINUE
C
      CALL GETINT(KINDI,NINP,XII,WT)
C
C  Integration loop
C
      DO 400 INP=1,NINP
        CALL SHAPE(XII(INP),KINDI,PSI,DPSI)
        XX=0.0
        YY=0.0
        DXDS=0.0
        DYDS=0.0
        DO 200 I=1,NL
          XX=XX+XQ(I)*PSI(I)
          YY=YY+YQ(I)*PSI(I)
          DXDS=DXDS+XQ(I)*DPSI(I)
          DYDS=DYDS+YQ(I)*DPSI(I)
200     CONTINUE
        DETJ=SQRT(DXDS**2+DYDS**2)
        QN(1)=DYDS/DETJ
        QN(2)=-DXDS/DETJ
        RX=XX-XP
        RY=YY-YP
        R=SQRT(RX*RX+RY*RY)
        DRDN=(QN(1)*RX+QN(2)*RY)/R
        HANK0=BESJ0(WN*R)-CMPLX(0.,1.)*BESY0(WN*R)
        HANK1=BESJ1(WN*R)-CMPLX(0.,1.)*BESY1(WN*R)
        GREEN=C1*HANK0*DETJ*WT(INP)
        DGDN=-C1*WN*HANK1*DRDN*DETJ*WT(INP)

        DO 300 I=1,NL
          H(I)=H(I)+PSI(I)*DGDN
          G(I)=G(I)+PSI(I)*GREEN
300     CONTINUE

        DGDN0=C0*DRDN/R*DETJ*WT(INP)
        CP=CP-DGDN0

400   CONTINUE

      RETURN
      END