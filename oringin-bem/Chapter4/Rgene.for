      SUBROUTINE RGENE
C
C---- Read in general information
C
      COMMON /PROB/INT,ISC,ISYM
      COMMON /FLUID/RHO,C
      COMMON /FRANGE/FREQ1,FREQ2,DF
      COMMON /FRE/FREQ,WN,IFR
      COMMON /SCA/AM,AF,BE,GA
C
      READ(5,*) INT,ISC,ISYM
      IF(ISC.EQ.1) THEN
        READ(5,*) AM,AF,BE,GA
      END IF
      READ(5,*) RHO,C
      READ(5,*) FREQ1,FREQ2,DF
      IFR=1
      FREQ=FREQ1
C
C---- PRINT OUT general information
C
      IF(INT.EQ.0) WRITE(6,*) '**** Exterior Problem ****'
      IF(INT.EQ.1) WRITE(6,*) '**** Interior Problem ****'
      IF(ISYM.NE.0) THEN
        WRITE(6,50) ISYM
50      FORMAT('Plane of Symmetry Perpendicular to X',I1,' AXIS') 
      END IF
      IF(ISC.EQ.1) THEN
        WRITE(6,*) '**** Scattering Problem ****'
        WRITE(6,*) 'Amplitude of Incident Wave =',AM
        WRITE(6,*) 'Direction Cosines of the Incident Wave',AF,BE,GA
      END IF
      WRITE(6,*) 'Fluid Density=',RHO
      WRITE(6,*) 'Speed of Sound=',C
      WRITE(6,*) 'First Frequency=',FREQ1,' Hz'
      IF(FREQ2.NE.0.0.AND.DF.NE.0.0) THEN
        WRITE(6,*) 'Last Frequency =',FREQ2,' Hz'
        WRITE(6,*) 'Frequency Step =',DF, ' Hz'  
      END IF    
      RETURN
      END
