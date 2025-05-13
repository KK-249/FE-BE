      SUBROUTINE ASSMB(K,IP,H,G,I34)
      DIMENSION H(8),G(8)
      COMMON /ELEM/NODE(8,4000)
      COMMON/BG/BIG
      COMMON /MATRX/A(4040,4000),B(4040),XS(4000)
      COMMON /SOL/SP(4000),VN(8,4000)
      COMMON /BCC/CA(8,4000),CB(8,4000),CC(8,4000)
      COMPLEX H,G,A,B,XS,CA,CB,CC,BIG,SP,VN

      DO 80 J=1,I34
        IQ=NODE(J,K)
        IF(SP(IQ).NE.BIG) THEN
          B(IP)=B(IP)-H(J)*SP(IQ)
          IF(CB(J,K).NE.0.0) THEN
            B(IP)=B(IP)+G(J)*(CC(J,K)-CA(J,K)*SP(IQ))/CB(J,K)
          ELSE
            A(IP,IQ)=A(IP,IQ)-G(J)
          END IF
        ELSE
          A(IP,IQ)=A(IP,IQ)+H(J)+G(J)*CA(J,K)/CB(J,K)
          B(IP)=B(IP)+G(J)*CC(J,K)/CB(J,K)
        END IF
80    CONTINUE
      RETURN
      END