      SUBROUTINE JBC
      COMMON /BG/BIG
      COMMON /BCC/CA(8,4000),CB(8,4000),CC(8,4000)
      COMMON /CONTR/NNODE,NELEM
      COMMON /ELEM/NODE(8,4000)

      COMPLEX CA,CB,CC,BIG
      logical :: found
      integer i, j,  K, NOD
C
      DO 20 J=1,8
        DO 20 K=1,NELEM
          CA(J,K)=0
          CB(J,K)=1
          CC(J,K)=0
20    CONTINUE
C
      
       READ(5,*) NREC
      !NREC = 2666
      DO 100 I=1,NREC
        READ(5,*) NOD,CA1,CA2,CB1,CB2,CC1,CC2
        CALL Find_Element(NODE, 8, 888, NOD, J, K, found)
          
           ! CA(J,K)=CMPLX(CA1,CA2)
            !CB(J,K)=CMPLX(CB1,CB2)
            CC(J,K)=CMPLX(CC1,CC2)
 

            
100   CONTINUE
      RETURN
      END