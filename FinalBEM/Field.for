      SUBROUTINE FIELD
C
C  Evaluate Field-Point Solution
C
      COMMON /PROB/INT,ISC,ISYM
      COMMON /FLUID/RHO,C
      COMMON /FPOINTS/NF,XF(3,100),FS(100)
      COMMON /CONTR/NNODE,NELEM
      COMMON /ELEM/NODE(8,4000)
      COMMON /SOL/SP(4000),VN(4,4000)
      COMMON /NODES/X(3,4000)
      COMMON /FRE/FREQ,WN,IFR
      COMMON /SCA/AM,AF,BE,GA
      DIMENSION XII(12),WT(12)
      DIMENSION QN(3),PSI(8),DPSI(2,8),P(3)
      REAL NN1,NN2,NN3
      COMPLEX IK,SP,VN,FS,PP,VV,GREEN,DGDN,IRHOW

      IF(NF.EQ.0) RETURN

      WRITE(6,*) ' '
      WRITE(6,*) '  Field-Point Solution:'
      WRITE(6,*) ' '

      C1=1.0/12.56637061
      IK=CMPLX(0.,1.)*WN
      IRHOW=CMPLX(0.,1.)*RHO*6.283185*FREQ

      NLOOP=1
      IF(ISYM.NE.0) NLOOP=2

      DO 900 IP=1,NF
       FS(IP)=0.0
       P(1)=XF(1,IP)
       P(2)=XF(2,IP)
       P(3)=XF(3,IP)
       DO 840 LOOP=1,NLOOP
        IF(LOOP.EQ.2) P(ISYM)=-XF(ISYM,IP)
        XP=P(1)
        YP=P(2)
        ZP=P(3)
        DO 800 K=1,NELEM
        I34=4
        IF(NODE(3,K).EQ.NODE(4,K)) I34=3
        CALL GETINT(NINP,XII,WT)
        DO 450 INP1=1,NINP
	 DO 400 INP2=1,NINP
          IF(I34.EQ.4) THEN
	    XI1=XII(INP1)
	    XI2=XII(INP2)
            CONST=1.0
          ELSE
            ETA1=XII(INP1)
            ETA2=XII(INP2)
            RHO1=0.5+0.5*ETA1
            THE=0.5+0.5*ETA2
            XI1=RHO1*(1.0-THE)
            XI2=RHO1*THE
            CONST=0.25*RHO1
          END IF
          CALL SHAPE(XI1,XI2,PSI,DPSI,I34)
          XX=0.0
          YY=0.0
	  ZZ=0.0
          DXDS1=0.0
	  DXDS2=0.0
          DYDS1=0.0
	  DYDS2=0.0
	  DZDS1=0.0
	  DZDS2=0.0
          PP=0.0
          VV=0.0
	  DO 200 I=1,I34
            J=NODE(I,K)
            XX=XX+X(1,J)*PSI(I)
            YY=YY+X(2,J)*PSI(I)
	    ZZ=ZZ+X(3,J)*PSI(I)
            DXDS1=DXDS1+X(1,J)*DPSI(1,I)
	    DYDS1=DYDS1+X(2,J)*DPSI(1,I)
	    DZDS1=DZDS1+X(3,J)*DPSI(1,I)
	    DXDS2=DXDS2+X(1,J)*DPSI(2,I)
            DYDS2=DYDS2+X(2,J)*DPSI(2,I)
	    DZDS2=DZDS2+X(3,J)*DPSI(2,I)
            PP=PP+SP(J)*PSI(I)
            VV=VV+VN(I,K)*PSI(I)		 
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

          GREEN=C1*CEXP(-IK*R)/R
          DGDN=-(1.0+IK*R)*GREEN/R*DRDN

          FACTOR=DETJ*WT(INP1)*WT(INP2)*CONST

          FS(IP)=FS(IP)+(-IRHOW*GREEN*VV-PP*DGDN)*FACTOR

400      CONTINUE
450     CONTINUE
800     CONTINUE
840    CONTINUE

       XP=XF(1,IP)
       YP=XF(2,IP)
       ZP=XF(3,IP)

       IF(ISC.EQ.1) THEN
         FS(IP)=FS(IP)+AM*CEXP(-IK*(AF*XP+BE*YP+GA*ZP))
       END IF
       
       AMP=CABS(FS(IP))
       WRITE(6,850) XP,YP,ZP,FS(IP),AMP
850    FORMAT(2X,'X=',F6.2,' Y=',F6.2,' Z=',F6.2,' P=',2F9.3,' AMPL=',F9
     1 .3)

900   CONTINUE
      RETURN
      END