      SUBROUTINE FIELD
C
C  Evaluate Field-Point Solution
C
      COMMON /FPOINTS/NF,XF(100),YF(100),FS(100)
      COMMON /CONTR/NNODE,NELEM
      COMMON /ELEM/KIND(2000),NODE(4,2000)
      COMMON /SOL/SP(2000),VN(4,2000)
      COMMON /NODES/X(2000),Y(2000)

      DIMENSION QN(2),PSI(4),DPSI(4),XII(12),WT(12)
      COMMON /FLUID/RHO,C
      COMMON /FRE/FREQ,WN,IFR
      COMMON /SCA/AM,RAD
      COMMON /PROB/INT,ISC
C
      COMPLEX SP,VN,FS,PP,VV,IK,IRHOW
      COMPLEX GREEN,DGDN,HANK0,HANK1,C1
      EXTERNAL BESJ0,BESY0,BESJ1,BESY1
C
      WRITE(6,*) ' '
      WRITE(6,*) '  Field-Point Solution:'
      WRITE(6,*) ' '

      C1=-CMPLX(0.,1.)/4.0
      IK=CMPLX(0.,1.)*WN
      IRHOW=CMPLX(0.,1.)*RHO*6.283185*FREQ
C
      DO 900 IP=1,NF
        FS(IP)=0.0
        XP=XF(IP)
        YP=YF(IP)

        DO 800 K=1,NELEM
          KINDI=KIND(K)
          NL=KINDI+1
C
          CALL GETINT(KINDI,NINP,XII,WT)

          DO 400 INP=1,NINP
            CALL SHAPE(XII(INP),KINDI,PSI,DPSI)
            XX=0.0
            YY=0.0
            DXDS=0.0
            DYDS=0.0
            PP=0.0
            VV=0.0
            DO 200 I=1,NL
              J=NODE(I,K)
              XX=XX+X(J)*PSI(I)
              YY=YY+Y(J)*PSI(I)
              DXDS=DXDS+X(J)*DPSI(I)
              DYDS=DYDS+Y(J)*DPSI(I)
              PP=PP+SP(J)*PSI(I)
              VV=VV+VN(I,K)*PSI(I)
200         CONTINUE

            DETJ=SQRT(DXDS**2+DYDS**2)
            QN(1)=DYDS/DETJ
            QN(2)=-DXDS/DETJ
            RX=XX-XP
            RY=YY-YP
            R=SQRT(RX*RX+RY*RY)
            DRDN=(QN(1)*RX+QN(2)*RY)/R
            HANK0=BESJ0(WN*R)-CMPLX(0.,1.)*BESY0(WN*R)
            HANK1=BESJ1(WN*R)-CMPLX(0.,1.)*BESY1(WN*R)
            GREEN=C1*HANK0
            DGDN=-C1*WN*HANK1*DRDN

            FS(IP)=FS(IP)+(-IRHOW*GREEN*VV-PP*DGDN)*DETJ*WT(INP)

400       CONTINUE
800     CONTINUE
C
        IF(ISC.EQ.1) THEN
          FS(IP)=FS(IP)+AM*CEXP(-IK*(COS(RAD)*XP+SIN(RAD)*YP))
        END IF
C
        WRITE(6,850) XF(IP),YF(IP),FS(IP),CABS(FS(IP))
850     FORMAT(2X,'X=',F6.2,2X,'Y=',F6.2,2X,'P=',2F9.3,2X,'AMPL=',F9.3)
c
900   CONTINUE

      RETURN
      END