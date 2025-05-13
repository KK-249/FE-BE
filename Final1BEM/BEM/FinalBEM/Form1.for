      SUBROUTINE FORM1
C
C  This subroutine formulates the coefficient matrix
C
      COMMON /PROB/INT,ISC,ISYM
      COMMON /FLUID/RHO,C
      COMMON /FRE/FREQ,WN,IFR
      COMMON /BG/BIG
      COMMON /CONTR/NNODE,NELEM
      COMMON /NODES/X(3,4000)
      COMMON /ELEM/NODE(8,4000)
      COMMON /MATRX/A(3000,3000),B(3000),XS(3000)
      COMMON /SOL/SP(4000),VN(8,4000)
      COMMON /CHI/NCH,XCH(3,20)
      COMMON /SCA/AM,AF,BE,GA

      DIMENSION XQ(8),YQ(8),ZQ(8),H(8),G(8),P(3),QN(3),AA(8)
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
        !IF(LOOP.EQ.2.AND.P(ISYM).EQ.0.0) CHECK='Y'
C
C  ELEMENT LOOP
C
         DO 400 K=1,NELEM
          ISING=0
          I34=8
          IF(NODE(3,K).EQ.NODE(4,K)) I34=3           
          DO 300 J=1,I34
            IQ=NODE(J,K)
            IF(IQ.EQ.IP.AND.CHECK.EQ.'Y')  ISING=J
            XQ(J)=X(1,IQ)
            YQ(J)=X(2,IQ)
            ZQ(J)=X(3,IQ)
300       CONTINUE

         IF(I34.EQ.8) THEN
             IF(ISING.EQ.0) THEN
               CALL ELEM8(XP,YP,ZP,XQ,YQ,ZQ,H,G,CP,AA,QN)
             ELSE
               CALL SING8(XP,YP,ZP,XQ,YQ,ZQ,H,G,CP,ISING,AA)
             END IF
        END IF
         
          DO 310 I=1,I34
              G(I)=-IRHOW*G(I)
310       CONTINUE
   
          CALL ASSMB(K,IP,H,G,I34)

400      CONTINUE
450      CONTINUE

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
