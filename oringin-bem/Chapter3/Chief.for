      SUBROUTINE CHIEF
C
C  This subroutine formulates the CHIEF equations
C
      COMMON /PROB/INT,ISC
      COMMON /BG/BIG
      COMMON /CONTR/NNODE,NELEM
      COMMON /NODES/X(2000),Y(2000)
      COMMON /ELEM/KIND(2000),NODE(4,2000)
      COMMON /MATRX/A(2020,2000),B(2020),XS(2000)
      COMMON /SOL/SP(2000),VN(4,2000)
      COMMON /CHI/NCH,XCH(20),YCH(20)
      COMMON /SCA/AM,RAD
      COMMON /FLUID/RHO,C
      COMMON /FRE/FREQ,WN,IFR
      DIMENSION H(4),G(4),XQ(4),YQ(4)
C
      COMPLEX BIG,A,B,XS,SP,VN,H,G,IK,IRHOW
C
      IK=CMPLX(0.,1.)*WN
      IRHOW=CMPLX(0.,1.)*RHO*6.283185*FREQ
C
C
C  P LOOP
C
      DO 500 IP=1,NCH
        XP=XCH(IP)
        YP=YCH(IP)
        CP=0.0
        IROW=NNODE+IP
C
C  ELEMENT LOOP
C
        DO 400 K=1,NELEM
          KINDI=KIND(K)
          NL=KINDI+1
          ISING=0
          DO 300 J=1,NL
            IQ=NODE(J,K)
            XQ(J)=X(IQ)
            YQ(J)=Y(IQ)
300       CONTINUE
          CALL ELEMT(XP,YP,NL,KINDI,XQ,YQ,H,G,CP)

          DO 310 I=1,NL
            G(I)=-IRHOW*G(I)
310       CONTINUE

          CALL ASSMB(K,NL,IROW,H,G)

400     CONTINUE

        IF(ISC.EQ.1) THEN
          B(IROW)=B(IROW)+AM*CEXP(-IK*(COS(RAD)*XP+SIN(RAD)*YP))
        END IF

500   CONTINUE
      RETURN
      END
