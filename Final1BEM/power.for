      SUBROUTINE POWER
C
C  Evaluate sound power and radiation efficiency
C
      COMMON /PROB/INT,ISC,ISYM
      COMMON /FLUID/RHO,C
      COMMON /CONTR/NNODE,NELEM
      COMMON /ELEM/NODE(8,2000)
      COMMON /SOL/SP(2000),VN(4,2000)
      COMMON /NODES/X(3,4000)
      COMMON /FRE/FREQ,WN,IFR
C      COMMON /SCA/AM,RAD
C
      DIMENSION XII(12),WT(12)
      DIMENSION QN(3),PSI(8),DPSI(2,8),P(3)
      REAL NN1,NN2,NN3
      COMPLEX SP,VN,FS,PP,VV
      
      IF(INT.EQ.1.OR.ISC.NE.0) RETURN

      V2=0.0
      POW=0.0

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

          FACTOR=DETJ*WT(INP1)*WT(INP2)*CONST

          V2=V2+CABS(VV)**2*FACTOR
          POW=POW+0.5*REAL(CONJG(PP)*(-VV))*FACTOR

400      CONTINUE
450     CONTINUE
800   CONTINUE

      IF(ISYM.NE.0) THEN
        POW=POW*2.
        V2=V2*2.
      END IF
    
      SIG=POW/(0.5*RHO*C*V2)

      WRITE(6,*) ' '
      WRITE(6,*) 'Total Radiated Sound Power=', POW
      WRITE(6,*) 'Radiation Efficiency=', SIG

      RETURN
      END
