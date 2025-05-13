      SUBROUTINE CHIEF
C
C  This subroutine formulates the CHIEF equations
C
      COMMON /PROB/INT,ISC,ISYM
      COMMON /BG/BIG
      COMMON /FLUID/RHO,C
      COMMON /FRE/FREQ,WN,IFR
      COMMON /CONTR/NNODE,NELEM
      COMMON /NODES/X(3,4000)
      COMMON /ELEM/NODE(8,4000)
      COMMON /MATRX/A(2020,2000),B(2020),XS(2000)
      COMMON /SOL/SP(4000),VN(4,4000)
      COMMON /CHI/NCH,XCH(3,20)
      COMMON /SCA/AM,AF,BE,GA
      DIMENSION H(4),G(4),XQ(4),YQ(4),ZQ(4),P(3)
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
       P(1)=XCH(1,IP)
       P(2)=XCH(2,IP)
       P(3)=XCH(3,IP)
       CP=0.0
       IROW=NNODE+IP

       NLOOP=1
       IF(ISYM.NE.0) NLOOP=2
C
C  Symmetry Loop
C
       DO 450 LOOP=1,NLOOP
        IF(LOOP.EQ.2) P(ISYM)=-XCH(ISYM,IP)
        XP=P(1)
        YP=P(2)
        ZP=P(3)
C
C  ELEMENT LOOP
C
        DO 400 K=1,NELEM
          I34=4
          IF(NODE(3,K).EQ.NODE(4,K)) I34=3
          DO 300 J=1,I34
            IQ=NODE(J,K)
            XQ(J)=X(1,IQ)
            YQ(J)=X(2,IQ)
            ZQ(J)=X(3,IQ)
300       CONTINUE

          IF(I34.EQ.4) THEN
            CALL ELEM4(XP,YP,ZP,XQ,YQ,ZQ,H,G,CP)
          ELSE
            CALL ELEM3(XP,YP,ZP,XQ,YQ,ZQ,H,G,CP)
          END IF

          DO 310 I=1,I34
            G(I)=-IRHOW*G(I)
310       CONTINUE    
       
          CALL ASSMB(K,IROW,H,G,I34)

400     CONTINUE
450    CONTINUE


       IF(ISC.EQ.1) THEN
         XP=XCH(1,IP)
         YP=XCH(2,IP)
         ZP=XCH(3,IP)
         B(IROW)=B(IROW)+AM*CEXP(-IK*(AF*XP+BE*YP+GA*ZP))
       END IF

500   CONTINUE
      RETURN
      END
