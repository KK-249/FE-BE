      SUBROUTINE POWER
C
C  Evaluate sound power and radiation efficiency
C
      COMMON /CONTR/NNODE,NELEM
      COMMON /ELEM/KIND(2000),NODE(4,2000)
      COMMON /SOL/SP(2000),VN(4,2000)
      COMMON /NODES/X(2000),Y(2000)
      COMMON /FLUID/RHO,C
      COMMON /PROB/INT,ISC

      COMPLEX SP,VN,PP,VV
      REAL LENGTH
      DIMENSION QN(2),PSI(4),DPSI(4),XII(12),WT(12)
C
      IF(INT.EQ.1.OR.ISC.NE.0) RETURN

      V2=0.0
      POW=0.0
      LENGTH=0.0

      DO 800 K=1,NELEM
        KINDI=KIND(K)
        NL=KINDI+1
C
        CALL GETINT(KINDI,NINP,XII,WT)

        DO 400 INP=1,NINP
          CALL SHAPE(XII(INP),KINDI,PSI,DPSI)
          DXDS=0.0
          DYDS=0.0
          PP=0.0
          VV=0.0
          DO 200 I=1,NL
            J=NODE(I,K)
            DXDS=DXDS+X(J)*DPSI(I)
            DYDS=DYDS+Y(J)*DPSI(I)
            PP=PP+SP(J)*PSI(I)
            VV=VV+VN(I,K)*PSI(I)
200       CONTINUE

          DETJ=SQRT(DXDS**2+DYDS**2)
          QN(1)=DYDS/DETJ
          QN(2)=-DXDS/DETJ

          FACTOR=DETJ*WT(INP)

          V2=V2+CABS(VV)**2*FACTOR
          POW=POW+0.5*REAL(CONJG(PP)*(-VV))*FACTOR
          LENGTH=LENGTH+FACTOR

400     CONTINUE
800   CONTINUE
C
      SIG=POW/(0.5*RHO*C*V2)

      WRITE(6,*) ' '
      WRITE(6,*) 'Total Length of Boundary=',LENGTH
      WRITE(6,*) 'Total Radiated Sound Power=', POW
      WRITE(6,*) 'Radiation Efficiency=', SIG

      RETURN
      END