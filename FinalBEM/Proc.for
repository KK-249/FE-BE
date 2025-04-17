      SUBROUTINE PROC
C
C---- PROCESSOR
C
      COMMON /FRANGE/FREQ1,FREQ2,DF
      COMMON /FLUID/RHO,C
      COMMON /FRE/FREQ,WN,IFR

100   WN=2.*3.14159*FREQ/C
      write(6,*) ' '
      write(6,*) ' Frequency=',FREQ,' Hz'
      write(6,*) ' Wavenumber=',WN

      CALL FORM
      CALL SOLVE
      CALL SOLOUT
      CALL FIELD
      CALL POWER

      IF(DF.EQ.0.0.OR.FREQ2.EQ.0.0) RETURN

      FREQ=FREQ+DF
      IFR=IFR+1
      IF(FREQ.GT.FREQ2) THEN
        RETURN
      ELSE
        CALL BC
        GO TO 100
      END IF

      END
