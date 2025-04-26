      SUBROUTINE FORM
C
C  This subroutine formulates the coefficient matrix
C
      COMMON /PROB/INT,ISC,ISYM
      COMMON /FLUID/RHO,C
      COMMON /FRE/FREQ,WN,IFR
      COMMON /BG/BIG
      COMMON /CONTR/NNODE,NELEM
      COMMON /NODES/X(3,2000)
      COMMON /ELEM/NODE(4,2000)
      COMMON /MATRX/A(2020,2000),B(2020),XS(2000)
      COMMON /SOL/SP(2000),VN(4,2000)
      COMMON /CHI/NCH,XCH(3,20)
      COMMON /SCA/AM,AF,BE,GA

      DIMENSION XQ(4),YQ(4),ZQ(4),H(4),G(4),P(3)
      CHARACTER*1 CHECK
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

       NLOOP=1
       IF(ISYM.NE.0) NLOOP=2
C
C  P LOOP
C
      DO 500 IP=1,NNODE
       P(1)=X(1,IP)
       P(2)=X(2,IP)
       P(3)=X(3,IP)
       CP=0.0
       IF(INT.EQ.0) CP=1.0
C
C  Symmetry Loop
C
       DO 450 LOOP=1,NLOOP
        IF(LOOP.EQ.2) P(ISYM)=-X(ISYM,IP)
        XP=P(1)
        YP=P(2)
        ZP=P(3)
        CHECK='N'
        IF(LOOP.EQ.1) CHECK='Y'
        IF(LOOP.EQ.2.AND.P(ISYM).EQ.0.0) CHECK='Y'
C
C  ELEMENT LOOP
C
         DO 400 K=1,NELEM
          ISING=0
          I34=4
          IF(NODE(3,K).EQ.NODE(4,K)) I34=3           
          DO 300 J=1,I34
            IQ=NODE(J,K)
            IF(IQ.EQ.IP.AND.CHECK.EQ.'Y') ISING=J
            XQ(J)=X(1,IQ)
            YQ(J)=X(2,IQ)
            ZQ(J)=X(3,IQ)
300       CONTINUE

          IF(I34.EQ.4) THEN
            IF(ISING.EQ.0) THEN
              CALL ELEM4(XP,YP,ZP,XQ,YQ,ZQ,H,G,CP)
            ELSE
              CALL SING4(XP,YP,ZP,XQ,YQ,ZQ,H,G,CP,ISING)
            END IF
          ELSE
            IF(ISING.EQ.0) THEN
              CALL ELEM3(XP,YP,ZP,XQ,YQ,ZQ,H,G,CP)
            ELSE
              CALL SING3(XP,YP,ZP,XQ,YQ,ZQ,H,G,CP,ISING)
            END IF
          END IF
         
          DO 310 I=1,I34
            G(I)=-IRHOW*G(I)
310       CONTINUE
   
          CALL ASSMB(K,IP,H,G,I34)

400      CONTINUE
450     CONTINUE

C
C   Contribution from C(P)
C
        IF(SP(IP).EQ.BIG) THEN
          A(IP,IP)=A(IP,IP)+CP
        ELSE
          B(IP)=B(IP)-CP*SP(IP)
        END IF


        IF(ISC.EQ.1) THEN
          XP=X(1,IP)
          YP=X(2,IP)
          ZP=X(3,IP)
          B(IP)=B(IP)+AM*CEXP(-IK*(AF*XP+BE*YP+GA*ZP))
        END IF

500   CONTINUE

      IF(NCH.NE.0) CALL CHIEF

      RETURN
      END
