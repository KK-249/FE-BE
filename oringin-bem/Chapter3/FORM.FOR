      SUBROUTINE FORM
C
C  This subroutine formulates the coefficient matrix
C
      COMMON /PROB/INT,ISC
      COMMON /BG/BIG
      COMMON /FLUID/RHO,C
      COMMON /FRE/FREQ,WN,IFR
      COMMON /CONTR/NNODE,NELEM
      COMMON /NODES/X(2000),Y(2000)
      COMMON /ELEM/KIND(2000),NODE(4,2000)
      COMMON /MATRX/A(2020,2000),B(2020),XS(2000)
      COMMON /SOL/SP(2000),VN(4,2000)
      COMMON /CHI/NCH,XCH(20),YCH(20)
      COMMON /SCA/AM,RAD
C
      DIMENSION H(4),G(4),XQ(4),YQ(4)
C
      COMPLEX BIG,A,B,XS,SP,VN,H,G,IK,IRHOW
C
      IK=CMPLX(0.,1.)*WN
      IRHOW=CMPLX(0.,1.)*RHO*6.283185*FREQ
C
      DO 100 I=1,NNODE+NCH
        B(I)=0.0
        DO 100 J=1,NNODE
          A(I,J)=0.0
100   CONTINUE
C
C  P LOOP
C
      DO 500 IP=1,NNODE
        XP=X(IP)
        YP=Y(IP)
        CP=0.0
        IF(INT.EQ.0) CP=1.0
C
C  ELEMENT LOOP
C
        DO 400 K=1,NELEM
          KINDI=KIND(K)
          NL=KINDI+1
          ISING=0
          DO 300 J=1,NL
            IQ=NODE(J,K)
            IF(IQ.EQ.IP) ISING=J
            XQ(J)=X(IQ)
            YQ(J)=Y(IQ)
300       CONTINUE

          IF(ISING.EQ.0) THEN
            CALL ELEMT(XP,YP,NL,KINDI,XQ,YQ,H,G,CP)
          ELSE
            CALL SING(XP,YP,NL,KINDI,XQ,YQ,H,G,CP,ISING)
          END IF

          DO 310 I=1,NL
            G(I)=-IRHOW*G(I)
310       CONTINUE

          CALL ASSMB(K,NL,IP,H,G)

400     CONTINUE
C
C   Contribution from C(P)
C
        IF(SP(IP).EQ.BIG) THEN
          A(IP,IP)=A(IP,IP)+CP
        ELSE
          B(IP)=B(IP)-CP*SP(IP)
        END IF
C
        IF(ISC.EQ.1) THEN
          B(IP)=B(IP)+AM*CEXP(-IK*(COS(RAD)*XP+SIN(RAD)*YP))
        END IF
C
500   CONTINUE

      IF(NCH.NE.0) CALL CHIEF

      RETURN
      END