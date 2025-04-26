      SUBROUTINE SING(XP,YP,NL,KINDI,XQ,YQ,H,G,CP,ISING)
      REAL LNG
      DIMENSION H(4),G(4),PSI(4),DPSI(4),ETA(12),QN(2)
      DIMENSION XQ(4),YQ(4),WT(12)
      COMPLEX H,G,GREEN,DGDN,HANK0,HANK1,C1
      COMMON /FRE/FREQ,WN
      EXTERNAL BESJ0,BESY0,BESJ1,BESY1
C
      C1=-CMPLX(0.,1.)/4.0
      C0=-1.0/6.2831853

      DO 100 I=1,NL
        H(I)=0.0
        G(I)=0.0
100   CONTINUE
C
      CALL GETINT(KINDI,NINP,ETA,WT)
C
c***LNG=LENGTH BETWEEN NODES IN THE +1 TO -1 (PSI) COORDINATE SYSTEM
C
      LNG=2./KINDI
      XIP=-1.+(ISING-1)*LNG
C
C***THIS LOOP IS FOR THE PURPOSE OF INTEGRATING TO THE RIGHT OF P
C
      IF(XIP.NE.1.) THEN

      DO 400 INP=1,NINP
          A=SQRT(1.-XIP)/2.
          Z=A+A*ETA(INP)
          XII=Z**2+XIP
         CALL SHAPE(XII,KINDI,PSI,DPSI)
         XX=0.
         YY=0.
         DXDS=0.
         DYDS=0.

         DO 200 I=1,NL
            XX=XX+XQ(I)*PSI(I)
            YY=YY+YQ(I)*PSI(I)
            DXDS=DXDS+XQ(I)*DPSI(I)
            DYDS=DYDS+YQ(I)*DPSI(I)
200      CONTINUE         

         DETJ=SQRT(DXDS**2+DYDS**2)
         QN(1)=DYDS/DETJ
         QN(2)=-DXDS/DETJ
         RX=XX-XP
         RY=YY-YP
         R=SQRT(RX**2+RY**2)
         DZDE=(SQRT(1.-XIP))/2.
         DRDN=(QN(1)*RX+QN(2)*RY)/R

         HANK0=BESJ0(WN*R)-CMPLX(0.,1.)*BESY0(WN*R)
         HANK1=BESJ1(WN*R)-CMPLX(0.,1.)*BESY1(WN*R)
         GREEN=C1*HANK0*DETJ*WT(INP)*2.*Z*DZDE
         DGDN=-C1*WN*HANK1*DRDN*DETJ*WT(INP)*2.*Z*DZDE

         DO 300 I=1,NL
           H(I)=H(I)+PSI(I)*DGDN
           G(I)=G(I)+PSI(I)*GREEN
300      CONTINUE

         DGDN0=C0*DRDN/R*DETJ*WT(INP)*2.*Z*DZDE
         CP=CP-DGDN0

400   CONTINUE
      ENDIF
C
C*** THIS LOOP IS FOR THE PURPOSE OF INTEGRATING TO THE LEFT OF P
C
      IF(XIP.NE.-1.) THEN
      
      DO 450 INP=1,NINP
          A=SQRT(1.+XIP)/2.
          Z=A+A*ETA(INP)
          XII=XIP-Z**2
         CALL SHAPE(XII,KINDI,PSI,DPSI)
         XX=0.
         YY=0.
         DXDS=0.
         DYDS=0.
         
         DO 250 I=1,NL
            XX=XX+XQ(I)*PSI(I)
            YY=YY+YQ(I)*PSI(I)
            DXDS=DXDS+XQ(I)*DPSI(I)
            DYDS=DYDS+YQ(I)*DPSI(I)
250      CONTINUE
         DETJ=SQRT(DXDS**2+DYDS**2)
         QN(1)=DYDS/DETJ
         QN(2)=-DXDS/DETJ

         RX=XX-XP
         RY=YY-YP
         R=SQRT(RX**2+RY**2)
         DZDE=(SQRT(1.+XIP))/2.
         DRDN=(QN(1)*RX+QN(2)*RY)/R

         HANK0=BESJ0(WN*R)-CMPLX(0.,1.)*BESY0(WN*R)
         HANK1=BESJ1(WN*R)-CMPLX(0.,1.)*BESY1(WN*R)
         GREEN=C1*HANK0*DETJ*WT(INP)*2.*Z*DZDE
         DGDN=-C1*WN*HANK1*DRDN*DETJ*WT(INP)*2.*Z*DZDE

         DO 350 I=1,NL
            H(I)=H(I)+PSI(I)*DGDN
            G(I)=G(I)+PSI(I)*GREEN
350      CONTINUE

         DGDN0=C0*DRDN/R*DETJ*WT(INP)*2.*Z*DZDE
         CP=CP-DGDN0

450   CONTINUE
      ENDIF

      RETURN
      END
