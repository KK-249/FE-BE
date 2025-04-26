      SUBROUTINE RBC
      COMMON /BG/BIG
      COMMON /BCC/CA(4,2000),CB(4,2000),CC(4,2000)
      COMMON /CONTR/NNODE,NELEM
      COMMON /ELEM/KIND(2000),NODE(4,2000)

      COMPLEX CA,CB,CC,BIG
C
      DO 20 J=1,4
        DO 20 K=1,NELEM
          CA(J,K)=BIG
          CB(J,K)=BIG
          CC(J,K)=BIG
20    CONTINUE
C
      READ(5,*) NREC
      DO 100 I=1,NREC
        READ(5,*) K1,K2,NOD,CA1,CA2,CB1,CB2,CC1,CC2
        JS=NOD
        JE=NOD
        DO 50 K=K1,K2
          KINDI=KIND(K)
          IF(NOD.EQ.0) THEN
            JS=1
            JE=KINDI+1
          END IF
          DO 40 J=JS,JE
            CA(J,K)=CMPLX(CA1,CA2)
            CB(J,K)=CMPLX(CB1,CB2)
            CC(J,K)=CMPLX(CC1,CC2)
40        CONTINUE
50      CONTINUE
100   CONTINUE
      RETURN
      END