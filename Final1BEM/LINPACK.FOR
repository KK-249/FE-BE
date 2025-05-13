C                                                                       
C---- COMBINED COMPLEX MATRIX SOLVER (GAUSS ELIMINATION + QR DECOMP.)   
C---- 11/7/91                                                           
C                                                                       
      SUBROUTINE CGECO(A,LDA,N,IPVT,RCOND,Z)                            
      INTEGER LDA,N,IPVT(1)                                             
      COMPLEX A(LDA,1),Z(1)                                             
      REAL RCOND                                                        
C                                                                       
C     CGECO FACTORS A COMPLEX MATRIX BY GAUSSIAN ELIMINATION            
C     AND ESTIMATES THE CONDITION OF THE MATRIX.                        
C                                                                       
C     IF  RCOND  IS NOT NEEDED, CGEFA IS SLIGHTLY FASTER.               
C     TO SOLVE  A*X = B , FOLLOW CGECO BY CGESL.                        
C     TO COMPUTE  INVERSE(A)*C , FOLLOW CGECO BY CGESL.                 
C     TO COMPUTE  DETERMINANT(A) , FOLLOW CGECO BY CGEDI.               
C     TO COMPUTE  INVERSE(A) , FOLLOW CGECO BY CGEDI.                   
C                                                                       
C     ON ENTRY                                                          
C                                                                       
C        A       COMPLEX(LDA, N)                                        
C                THE MATRIX TO BE FACTORED.                             
C                                                                       
C        LDA     INTEGER                                                
C                THE LEADING DIMENSION OF THE ARRAY  A .                
C                                                                       
C        N       INTEGER                                                
C                THE ORDER OF THE MATRIX  A .                           
C                                                                       
C     ON RETURN                                                         
C                                                                       
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS         
C                WHICH WERE USED TO OBTAIN IT.                          
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE       
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER          
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.       
C                                                                       
C        IPVT    INTEGER(N)                                             
C                AN INTEGER VECTOR OF PIVOT INDICES.                    
C                                                                       
C        RCOND   REAL                                                   
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .        
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS       
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE             
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . 
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION     
C                           1.0 + RCOND .EQ. 1.0                        
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING           
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF         
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE          
C                UNDERFLOWS.                                            
C                                                                       
C        Z       COMPLEX(N)                                             
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.  
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS      
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT           
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .                    
C                                                                       
C     LINPACK. THIS VERSION DATED 08/14/78 .                            
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
C                                                                       
C     SUBROUTINES AND FUNCTIONS                                         
C                                                                       
C     LINPACK CGEFA                                                     
C     BLAS CAXPY,CDOTC,CSSCAL,SCASUM                                    
C     FORTRAN ABS,AIMAG,AMAX1,CMPLX,CONJG,REAL                          
C                                                                       
C     INTERNAL VARIABLES                                                
C                                                                       
      COMPLEX CDOTC,EK,T,WK,WKM                                         
      REAL ANORM,S,SCASUM,SM,YNORM                                      
      INTEGER INFO,J,K,KB,KP1,L                                         
C                                                                       
      COMPLEX ZDUM,ZDUM1,ZDUM2,CSIGN1                                   
      REAL CABS1                                                        
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))                  
      CSIGN1(ZDUM1,ZDUM2) = CABS1(ZDUM1)*(ZDUM2/CABS1(ZDUM2))           
C                                                                       
C     COMPUTE 1-NORM OF A                                               
C                                                                       
      ANORM = 0.0E0                                                     
      DO 10 J = 1, N                                                    
         ANORM = AMAX1(ANORM,SCASUM(N,A(1,J),1))                        
   10 CONTINUE                                                          
C                                                                       
C     FACTOR                                                            
C                                                                       
      CALL CGEFA(A,LDA,N,IPVT,INFO)                                     
C                                                                       
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .              
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  CTRANS(A)*Y = E . 
C     CTRANS(A)  IS THE CONJUGATE TRANSPOSE OF A .                      
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL           
C     GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(U)*W = E .             
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.            
C                                                                       
C     SOLVE CTRANS(U)*W = E                                             
C                                                                       
      EK = (1.0E0,0.0E0)                                                
      DO 20 J = 1, N                                                    
         Z(J) = (0.0E0,0.0E0)                                           
   20 CONTINUE                                                          
      DO 100 K = 1, N                                                   
         IF (CABS1(Z(K)) .NE. 0.0E0) EK = CSIGN1(EK,-Z(K))              
         IF (CABS1(EK-Z(K)) .LE. CABS1(A(K,K))) GO TO 30                
            S = CABS1(A(K,K))/CABS1(EK-Z(K))                            
            CALL CSSCAL(N,S,Z,1)                                        
            EK = CMPLX(S,0.0E0)*EK                                      
   30    CONTINUE                                                       
         WK = EK - Z(K)                                                 
         WKM = -EK - Z(K)                                               
         S = CABS1(WK)                                                  
         SM = CABS1(WKM)                                                
         IF (CABS1(A(K,K)) .EQ. 0.0E0) GO TO 40                         
            WK = WK/CONJG(A(K,K))                                       
            WKM = WKM/CONJG(A(K,K))                                     
         GO TO 50                                                       
   40    CONTINUE                                                       
            WK = (1.0E0,0.0E0)                                          
            WKM = (1.0E0,0.0E0)                                         
   50    CONTINUE                                                       
         KP1 = K + 1                                                    
         IF (KP1 .GT. N) GO TO 90                                       
            DO 60 J = KP1, N                                            
               SM = SM + CABS1(Z(J)+WKM*CONJG(A(K,J)))                  
               Z(J) = Z(J) + WK*CONJG(A(K,J))                           
               S = S + CABS1(Z(J))                                      
   60       CONTINUE                                                    
            IF (S .GE. SM) GO TO 80                                     
               T = WKM - WK                                             
               WK = WKM                                                 
               DO 70 J = KP1, N                                         
                  Z(J) = Z(J) + T*CONJG(A(K,J))                         
   70          CONTINUE                                                 
   80       CONTINUE                                                    
   90    CONTINUE                                                       
         Z(K) = WK                                                      
  100 CONTINUE                                                          
      S = 1.0E0/SCASUM(N,Z,1)                                           
      CALL CSSCAL(N,S,Z,1)                                              
C                                                                       
C     SOLVE CTRANS(L)*Y = W                                             
C                                                                       
      DO 120 KB = 1, N                                                  
         K = N + 1 - KB                                                 
         IF (K .LT. N) Z(K) = Z(K) + CDOTC(N-K,A(K+1,K),1,Z(K+1),1)     
         IF (CABS1(Z(K)) .LE. 1.0E0) GO TO 110                          
            S = 1.0E0/CABS1(Z(K))                                       
            CALL CSSCAL(N,S,Z,1)                                        
  110    CONTINUE                                                       
         L = IPVT(K)                                                    
         T = Z(L)                                                       
         Z(L) = Z(K)                                                    
         Z(K) = T                                                       
  120 CONTINUE                                                          
      S = 1.0E0/SCASUM(N,Z,1)                                           
      CALL CSSCAL(N,S,Z,1)                                              
C                                                                       
      YNORM = 1.0E0                                                     
C                                                                       
C     SOLVE L*V = Y                                                     
C                                                                       
      DO 140 K = 1, N                                                   
         L = IPVT(K)                                                    
         T = Z(L)                                                       
         Z(L) = Z(K)                                                    
         Z(K) = T                                                       
         IF (K .LT. N) CALL CAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)            
         IF (CABS1(Z(K)) .LE. 1.0E0) GO TO 130                          
            S = 1.0E0/CABS1(Z(K))                                       
            CALL CSSCAL(N,S,Z,1)                                        
            YNORM = S*YNORM                                             
  130    CONTINUE                                                       
  140 CONTINUE                                                          
      S = 1.0E0/SCASUM(N,Z,1)                                           
      CALL CSSCAL(N,S,Z,1)                                              
      YNORM = S*YNORM                                                   
C                                                                       
C     SOLVE  U*Z = V                                                    
C                                                                       
      DO 160 KB = 1, N                                                  
         K = N + 1 - KB                                                 
         IF (CABS1(Z(K)) .LE. CABS1(A(K,K))) GO TO 150                  
            S = CABS1(A(K,K))/CABS1(Z(K))                               
            CALL CSSCAL(N,S,Z,1)                                        
            YNORM = S*YNORM                                             
  150    CONTINUE                                                       
         IF (CABS1(A(K,K)) .NE. 0.0E0) Z(K) = Z(K)/A(K,K)               
         IF (CABS1(A(K,K)) .EQ. 0.0E0) Z(K) = (1.0E0,0.0E0)             
         T = -Z(K)                                                      
         CALL CAXPY(K-1,T,A(1,K),1,Z(1),1)                              
  160 CONTINUE                                                          
C     MAKE ZNORM = 1.0                                                  
      S = 1.0E0/SCASUM(N,Z,1)                                           
      CALL CSSCAL(N,S,Z,1)                                              
      YNORM = S*YNORM                                                   
C                                                                       
      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM                         
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0                               
      RETURN                                                            
      END                                                               
      SUBROUTINE CGEFA(A,LDA,N,IPVT,INFO)                               
      INTEGER LDA,N,IPVT(1),INFO                                        
      COMPLEX A(LDA,1)                                                  
C                                                                       
C     CGEFA FACTORS A COMPLEX MATRIX BY GAUSSIAN ELIMINATION.           
C                                                                       
C     CGEFA IS USUALLY CALLED BY CGECO, BUT IT CAN BE CALLED            
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.          
C     (TIME FOR CGECO) = (1 + 9/N)*(TIME FOR CGEFA) .                   
C                                                                       
C     ON ENTRY                                                          
C                                                                       
C        A       COMPLEX(LDA, N)                                        
C                THE MATRIX TO BE FACTORED.                             
C                                                                       
C        LDA     INTEGER                                                
C                THE LEADING DIMENSION OF THE ARRAY  A .                
C                                                                       
C        N       INTEGER                                                
C                THE ORDER OF THE MATRIX  A .                           
C                                                                       
C     ON RETURN                                                         
C                                                                       
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS         
C                WHICH WERE USED TO OBTAIN IT.                          
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE       
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER          
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.       
C                                                                       
C        IPVT    INTEGER(N)                                             
C                AN INTEGER VECTOR OF PIVOT INDICES.                    
C                                                                       
C        INFO    INTEGER                                                
C                = 0  NORMAL VALUE.                                     
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR       
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES        
C                     INDICATE THAT CGESL OR CGEDI WILL DIVIDE BY ZERO  
C                     IF CALLED.  USE  RCOND  IN CGECO FOR A RELIABLE   
C                     INDICATION OF SINGULARITY.                        
C                                                                       
C     LINPACK. THIS VERSION DATED 08/14/78 .                            
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
C                                                                       
C     SUBROUTINES AND FUNCTIONS                                         
C                                                                       
C     BLAS CAXPY,CSCAL,ICAMAX                                           
C     FORTRAN ABS,AIMAG,REAL                                            
C                                                                       
C     INTERNAL VARIABLES                                                
C                                                                       
      COMPLEX T                                                         
      INTEGER ICAMAX,J,K,KP1,L,NM1                                      
C                                                                       
      COMPLEX ZDUM                                                      
      REAL CABS1                                                        
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))                  
C                                                                       
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING                        
C                                                                       
      INFO = 0                                                          
      NM1 = N - 1                                                       
      IF (NM1 .LT. 1) GO TO 70                                          
      DO 60 K = 1, NM1                                                  
         KP1 = K + 1                                                    
C                                                                       
C        FIND L = PIVOT INDEX                                           
C                                                                       
         L = ICAMAX(N-K+1,A(K,K),1) + K - 1                             
         IPVT(K) = L                                                    
C                                                                       
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED          
C                                                                       
         IF (CABS1(A(L,K)) .EQ. 0.0E0) GO TO 40                         
C                                                                       
C           INTERCHANGE IF NECESSARY                                    
C                                                                       
            IF (L .EQ. K) GO TO 10                                      
               T = A(L,K)                                               
               A(L,K) = A(K,K)                                          
               A(K,K) = T                                               
   10       CONTINUE                                                    
C                                                                       
C           COMPUTE MULTIPLIERS                                         
C                                                                       
            T = -(1.0E0,0.0E0)/A(K,K)                                   
            CALL CSCAL(N-K,T,A(K+1,K),1)                                
C                                                                       
C           ROW ELIMINATION WITH COLUMN INDEXING                        
C                                                                       
            DO 30 J = KP1, N                                            
               T = A(L,J)                                               
               IF (L .EQ. K) GO TO 20                                   
                  A(L,J) = A(K,J)                                       
                  A(K,J) = T                                            
   20          CONTINUE                                                 
               CALL CAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)                  
   30       CONTINUE                                                    
         GO TO 50                                                       
   40    CONTINUE                                                       
            INFO = K                                                    
   50    CONTINUE                                                       
   60 CONTINUE                                                          
   70 CONTINUE                                                          
      IPVT(N) = N                                                       
      IF (CABS1(A(N,N)) .EQ. 0.0E0) INFO = N                            
      RETURN                                                            
      END                                                               
      SUBROUTINE CGESL(A,LDA,N,IPVT,B,JOB)                              
      INTEGER LDA,N,IPVT(1),JOB                                         
      COMPLEX A(LDA,1),B(1)                                             
C                                                                       
C     CGESL SOLVES THE COMPLEX SYSTEM                                   
C     A * X = B  OR  CTRANS(A) * X = B                                  
C     USING THE FACTORS COMPUTED BY CGECO OR CGEFA.                     
C                                                                       
C     ON ENTRY                                                          
C                                                                       
C        A       COMPLEX(LDA, N)                                        
C                THE OUTPUT FROM CGECO OR CGEFA.                        
C                                                                       
C        LDA     INTEGER                                                
C                THE LEADING DIMENSION OF THE ARRAY  A .                
C                                                                       
C        N       INTEGER                                                
C                THE ORDER OF THE MATRIX  A .                           
C                                                                       
C        IPVT    INTEGER(N)                                             
C                THE PIVOT VECTOR FROM CGECO OR CGEFA.                  
C                                                                       
C        B       COMPLEX(N)                                             
C                THE RIGHT HAND SIDE VECTOR.                            
C                                                                       
C        JOB     INTEGER                                                
C                = 0         TO SOLVE  A*X = B ,                        
C                = NONZERO   TO SOLVE  CTRANS(A)*X = B  WHERE           
C                            CTRANS(A)  IS THE CONJUGATE TRANSPOSE.     
C                                                                       
C     ON RETURN                                                         
C                                                                       
C        B       THE SOLUTION VECTOR  X .                               
C                                                                       
C     ERROR CONDITION                                                   
C                                                                       
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A   
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY  
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER       
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE     
C        CALLED CORRECTLY AND IF CGECO HAS SET RCOND .GT. 0.0           
C        OR CGEFA HAS SET INFO .EQ. 0 .                                 
C                                                                       
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 
C     WITH  P  COLUMNS                                                  
C           CALL CGECO(A,LDA,N,IPVT,RCOND,Z)                            
C           IF (RCOND IS TOO SMALL) GO TO ...                           
C           DO 10 J = 1, P                                              
C              CALL CGESL(A,LDA,N,IPVT,C(1,J),0)                        
C        10 CONTINUE                                                    
C                                                                       
C     LINPACK. THIS VERSION DATED 08/14/78 .                            
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
C                                                                       
C     SUBROUTINES AND FUNCTIONS                                         
C                                                                       
C     BLAS CAXPY,CDOTC                                                  
C     FORTRAN CONJG                                                     
C                                                                       
C     INTERNAL VARIABLES                                                
C                                                                       
      COMPLEX CDOTC,T                                                   
      INTEGER K,KB,L,NM1                                                
C                                                                       
      NM1 = N - 1                                                       
      IF (JOB .NE. 0) GO TO 50                                          
C                                                                       
C        JOB = 0 , SOLVE  A * X = B                                     
C        FIRST SOLVE  L*Y = B                                           
C                                                                       
         IF (NM1 .LT. 1) GO TO 30                                       
         DO 20 K = 1, NM1                                               
            L = IPVT(K)                                                 
            T = B(L)                                                    
            IF (L .EQ. K) GO TO 10                                      
               B(L) = B(K)                                              
               B(K) = T                                                 
   10       CONTINUE                                                    
            CALL CAXPY(N-K,T,A(K+1,K),1,B(K+1),1)                       
   20    CONTINUE                                                       
   30    CONTINUE                                                       
C                                                                       
C        NOW SOLVE  U*X = Y                                             
C                                                                       
         DO 40 KB = 1, N                                                
            K = N + 1 - KB                                              
            B(K) = B(K)/A(K,K)                                          
            T = -B(K)                                                   
            CALL CAXPY(K-1,T,A(1,K),1,B(1),1)                           
   40    CONTINUE                                                       
      GO TO 100                                                         
   50 CONTINUE                                                          
C                                                                       
C        JOB = NONZERO, SOLVE  CTRANS(A) * X = B                        
C        FIRST SOLVE  CTRANS(U)*Y = B                                   
C                                                                       
         DO 60 K = 1, N                                                 
            T = CDOTC(K-1,A(1,K),1,B(1),1)                              
            B(K) = (B(K) - T)/CONJG(A(K,K))                             
   60    CONTINUE                                                       
C                                                                       
C        NOW SOLVE CTRANS(L)*X = Y                                      
C                                                                       
         IF (NM1 .LT. 1) GO TO 90                                       
         DO 80 KB = 1, NM1                                              
            K = N - KB                                                  
            B(K) = B(K) + CDOTC(N-K,A(K+1,K),1,B(K+1),1)                
            L = IPVT(K)                                                 
            IF (L .EQ. K) GO TO 70                                      
               T = B(L)                                                 
               B(L) = B(K)                                              
               B(K) = T                                                 
   70       CONTINUE                                                    
   80    CONTINUE                                                       
   90    CONTINUE                                                       
  100 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE  CSCAL(N,CA,CX,INCX)                                   
C                                                                       
C     SCALES A VECTOR BY A CONSTANT.                                    
C     JACK DONGARRA, LINPACK,  3/11/78.                                 
C                                                                       
      COMPLEX CA,CX(1)                                                  
      INTEGER I,INCX,N,NINCX                                            
C                                                                       
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1)GO TO 20                                             
C                                                                       
C        CODE FOR INCREMENT NOT EQUAL TO 1                              
C                                                                       
      NINCX = N*INCX                                                    
      DO 10 I = 1,NINCX,INCX                                            
        CX(I) = CA*CX(I)                                                
   10 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C        CODE FOR INCREMENT EQUAL TO 1                                  
C                                                                       
   20 DO 30 I = 1,N                                                     
        CX(I) = CA*CX(I)                                                
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      INTEGER FUNCTION ISAMAX(N,SX,INCX)                                
C                                                                       
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.            
C     JACK DONGARRA, LINPACK, 3/11/78.                                  
C                                                                       
      REAL SX(1),SMAX                                                   
      INTEGER I,INCX,IX,N                                               
C                                                                       
      ISAMAX = 0                                                        
      IF( N .LT. 1 ) RETURN                                             
      ISAMAX = 1                                                        
      IF(N.EQ.1)RETURN                                                  
      IF(INCX.EQ.1)GO TO 20                                             
C                                                                       
C        CODE FOR INCREMENT NOT EQUAL TO 1                              
C                                                                       
      IX = 1                                                            
      SMAX = ABS(SX(1))                                                 
      IX = IX + INCX                                                    
      DO 10 I = 2,N                                                     
         IF(ABS(SX(IX)).LE.SMAX) GO TO 5                                
         ISAMAX = I                                                     
         SMAX = ABS(SX(IX))                                             
    5    IX = IX + INCX                                                 
   10 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C        CODE FOR INCREMENT EQUAL TO 1                                  
C                                                                       
   20 SMAX = ABS(SX(1))                                                 
      DO 30 I = 2,N                                                     
         IF(ABS(SX(I)).LE.SMAX) GO TO 30                                
         ISAMAX = I                                                     
         SMAX = ABS(SX(I))                                              
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      REAL FUNCTION SASUM(N,SX,INCX)                                    
C                                                                       
C     TAKES THE SUM OF THE ABSOLUTE VALUES.                             
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.                   
C     JACK DONGARRA, LINPACK, 3/11/78.                                  
C                                                                       
      REAL SX(1),STEMP                                                  
      INTEGER I,INCX,M,MP1,N,NINCX                                      
C                                                                       
      SASUM = 0.0E0                                                     
      STEMP = 0.0E0                                                     
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1)GO TO 20                                             
C                                                                       
C        CODE FOR INCREMENT NOT EQUAL TO 1                              
C                                                                       
      NINCX = N*INCX                                                    
      DO 10 I = 1,NINCX,INCX                                            
        STEMP = STEMP + ABS(SX(I))                                      
   10 CONTINUE                                                          
      SASUM = STEMP                                                     
      RETURN                                                            
C                                                                       
C        CODE FOR INCREMENT EQUAL TO 1                                  
C                                                                       
C                                                                       
C        CLEAN-UP LOOP                                                  
C                                                                       
   20 M = MOD(N,6)                                                      
      IF( M .EQ. 0 ) GO TO 40                                           
      DO 30 I = 1,M                                                     
        STEMP = STEMP + ABS(SX(I))                                      
   30 CONTINUE                                                          
      IF( N .LT. 6 ) GO TO 60                                           
   40 MP1 = M + 1                                                       
      DO 50 I = MP1,N,6                                                 
        STEMP = STEMP + ABS(SX(I)) + ABS(SX(I + 1)) + ABS(SX(I + 2))    
     *  + ABS(SX(I + 3)) + ABS(SX(I + 4)) + ABS(SX(I + 5))              
   50 CONTINUE                                                          
   60 SASUM = STEMP                                                     
      RETURN                                                            
      END                                                               
      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)                            
C                                                                       
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.                            
C     USES UNROLLED LOOP FOR INCREMENTS EQUAL TO ONE.                   
C     JACK DONGARRA, LINPACK, 3/11/78.                                  
C                                                                       
      REAL SX(1),SY(1),SA                                               
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N                                 
C                                                                       
      IF(N.LE.0)RETURN                                                  
      IF (SA .EQ. 0.0) RETURN                                           
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               
C                                                                       
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                
C          NOT EQUAL TO 1                                               
C                                                                       
      IX = 1                                                            
      IY = 1                                                            
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 
      DO 10 I = 1,N                                                     
        SY(IY) = SY(IY) + SA*SX(IX)                                     
        IX = IX + INCX                                                  
        IY = IY + INCY                                                  
   10 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C        CODE FOR BOTH INCREMENTS EQUAL TO 1                            
C                                                                       
C                                                                       
C        CLEAN-UP LOOP                                                  
C                                                                       
   20 M = MOD(N,4)                                                      
      IF( M .EQ. 0 ) GO TO 40                                           
      DO 30 I = 1,M                                                     
        SY(I) = SY(I) + SA*SX(I)                                        
   30 CONTINUE                                                          
      IF( N .LT. 4 ) RETURN                                             
   40 MP1 = M + 1                                                       
      DO 50 I = MP1,N,4                                                 
        SY(I) = SY(I) + SA*SX(I)                                        
        SY(I + 1) = SY(I + 1) + SA*SX(I + 1)                            
        SY(I + 2) = SY(I + 2) + SA*SX(I + 2)                            
        SY(I + 3) = SY(I + 3) + SA*SX(I + 3)                            
   50 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      REAL FUNCTION SDOT(N,SX,INCX,SY,INCY)                             
C                                                                       
C     FORMS THE DOT PRODUCT OF TWO VECTORS.                             
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.                  
C     JACK DONGARRA, LINPACK, 3/11/78.                                  
C                                                                       
      REAL SX(1),SY(1),STEMP                                            
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N                                 
C                                                                       
      STEMP = 0.0E0                                                     
      SDOT = 0.0E0                                                      
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               
C                                                                       
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                
C          NOT EQUAL TO 1                                               
C                                                                       
      IX = 1                                                            
      IY = 1                                                            
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 
      DO 10 I = 1,N                                                     
        STEMP = STEMP + SX(IX)*SY(IY)                                   
        IX = IX + INCX                                                  
        IY = IY + INCY                                                  
   10 CONTINUE                                                          
      SDOT = STEMP                                                      
      RETURN                                                            
C                                                                       
C        CODE FOR BOTH INCREMENTS EQUAL TO 1                            
C                                                                       
C                                                                       
C        CLEAN-UP LOOP                                                  
C                                                                       
   20 M = MOD(N,5)                                                      
      IF( M .EQ. 0 ) GO TO 40                                           
      DO 30 I = 1,M                                                     
        STEMP = STEMP + SX(I)*SY(I)                                     
   30 CONTINUE                                                          
      IF( N .LT. 5 ) GO TO 60                                           
   40 MP1 = M + 1                                                       
      DO 50 I = MP1,N,5                                                 
        STEMP = STEMP + SX(I)*SY(I) + SX(I + 1)*SY(I + 1) +             
     *   SX(I + 2)*SY(I + 2) + SX(I + 3)*SY(I + 3) + SX(I + 4)*SY(I + 4)
   50 CONTINUE                                                          
   60 SDOT = STEMP                                                      
      RETURN                                                            
      END                                                               
      SUBROUTINE  SSCAL(N,SA,SX,INCX)                                   
C                                                                       
C     SCALES A VECTOR BY A CONSTANT.                                    
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO 1.                     
C     JACK DONGARRA, LINPACK, 3/11/78.                                  
C                                                                       
      REAL SA,SX(1)                                                     
      INTEGER I,INCX,M,MP1,N,NINCX                                      
C                                                                       
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1)GO TO 20                                             
C                                                                       
C        CODE FOR INCREMENT NOT EQUAL TO 1                              
C                                                                       
      NINCX = N*INCX                                                    
      DO 10 I = 1,NINCX,INCX                                            
        SX(I) = SA*SX(I)                                                
   10 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C        CODE FOR INCREMENT EQUAL TO 1                                  
C                                                                       
C                                                                       
C        CLEAN-UP LOOP                                                  
C                                                                       
   20 M = MOD(N,5)                                                      
      IF( M .EQ. 0 ) GO TO 40                                           
      DO 30 I = 1,M                                                     
        SX(I) = SA*SX(I)                                                
   30 CONTINUE                                                          
      IF( N .LT. 5 ) RETURN                                             
   40 MP1 = M + 1                                                       
      DO 50 I = MP1,N,5                                                 
        SX(I) = SA*SX(I)                                                
        SX(I + 1) = SA*SX(I + 1)                                        
        SX(I + 2) = SA*SX(I + 2)                                        
        SX(I + 3) = SA*SX(I + 3)                                        
        SX(I + 4) = SA*SX(I + 4)                                        
   50 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE CAXPY(N,CA,CX,INCX,CY,INCY)                            
C                                                                       
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.                            
C     JACK DONGARRA, LINPACK, 3/11/78.                                  
C                                                                       
      COMPLEX CX(1),CY(1),CA                                            
      INTEGER I,INCX,INCY,IX,IY,N                                       
C                                                                       
      IF(N.LE.0)RETURN                                                  
      IF (ABS(REAL(CA)) + ABS(AIMAG(CA)) .EQ. 0.0 ) RETURN              
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               
C                                                                       
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                
C          NOT EQUAL TO 1                                               
C                                                                       
      IX = 1                                                            
      IY = 1                                                            
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 
      DO 10 I = 1,N                                                     
        CY(IY) = CY(IY) + CA*CX(IX)                                     
        IX = IX + INCX                                                  
        IY = IY + INCY                                                  
   10 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C        CODE FOR BOTH INCREMENTS EQUAL TO 1                            
C                                                                       
   20 DO 30 I = 1,N                                                     
        CY(I) = CY(I) + CA*CX(I)                                        
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      COMPLEX FUNCTION CDOTC(N,CX,INCX,CY,INCY)                         
C                                                                       
C     FORMS THE DOT PRODUCT OF TWO VECTORS, CONJUGATING THE FIRST       
C     VECTOR.                                                           
C     JACK DONGARRA, LINPACK,  3/11/78.                                 
C                                                                       
      COMPLEX CX(1),CY(1),CTEMP                                         
      INTEGER I,INCX,INCY,IX,IY,N                                       
C                                                                       
      CTEMP = (0.0,0.0)                                                 
      CDOTC = (0.0,0.0)                                                 
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                               
C                                                                       
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                
C          NOT EQUAL TO 1                                               
C                                                                       
      IX = 1                                                            
      IY = 1                                                            
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                 
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                 
      DO 10 I = 1,N                                                     
        CTEMP = CTEMP + CONJG(CX(IX))*CY(IY)                            
        IX = IX + INCX                                                  
        IY = IY + INCY                                                  
   10 CONTINUE                                                          
      CDOTC = CTEMP                                                     
      RETURN                                                            
C                                                                       
C        CODE FOR BOTH INCREMENTS EQUAL TO 1                            
C                                                                       
   20 DO 30 I = 1,N                                                     
        CTEMP = CTEMP + CONJG(CX(I))*CY(I)                              
   30 CONTINUE                                                          
      CDOTC = CTEMP                                                     
      RETURN                                                            
      END                                                               
      SUBROUTINE  CSSCAL(N,SA,CX,INCX)                                  
C                                                                       
C     SCALES A COMPLEX VECTOR BY A REAL CONSTANT.                       
C     JACK DONGARRA, LINPACK, 3/11/78.                                  
C                                                                       
      COMPLEX CX(1)                                                     
      REAL SA                                                           
      INTEGER I,INCX,N,NINCX                                            
C                                                                       
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1)GO TO 20                                             
C                                                                       
C        CODE FOR INCREMENT NOT EQUAL TO 1                              
C                                                                       
      NINCX = N*INCX                                                    
      DO 10 I = 1,NINCX,INCX                                            
        CX(I) = CMPLX(SA*REAL(CX(I)),SA*AIMAG(CX(I)))                   
   10 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C        CODE FOR INCREMENT EQUAL TO 1                                  
C                                                                       
   20 DO 30 I = 1,N                                                     
        CX(I) = CMPLX(SA*REAL(CX(I)),SA*AIMAG(CX(I)))                   
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      REAL FUNCTION SCASUM(N,CX,INCX)                                   
C                                                                       
C     TAKES THE SUM OF THE ABSOLUTE VALUES OF A COMPLEX VECTOR AND      
C     RETURNS A SINGLE PRECISION RESULT.                                
C     JACK DONGARRA, LINPACK, 3/11/78.                                  
C                                                                       
      COMPLEX CX(1)                                                     
      REAL STEMP                                                        
      INTEGER I,INCX,N,NINCX                                            
C                                                                       
      SCASUM = 0.0E0                                                    
      STEMP = 0.0E0                                                     
      IF(N.LE.0)RETURN                                                  
      IF(INCX.EQ.1)GO TO 20                                             
C                                                                       
C        CODE FOR INCREMENT NOT EQUAL TO 1                              
C                                                                       
      NINCX = N*INCX                                                    
      DO 10 I = 1,NINCX,INCX                                            
        STEMP = STEMP + ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)))            
   10 CONTINUE                                                          
      SCASUM = STEMP                                                    
      RETURN                                                            
C                                                                       
C        CODE FOR INCREMENT EQUAL TO 1                                  
C                                                                       
   20 DO 30 I = 1,N                                                     
        STEMP = STEMP + ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)))            
   30 CONTINUE                                                          
      SCASUM = STEMP                                                    
      RETURN                                                            
      END                                                               
      INTEGER FUNCTION ICAMAX(N,CX,INCX)                                
C                                                                       
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.            
C     JACK DONGARRA, LINPACK, 3/11/78.                                  
C                                                                       
      COMPLEX CX(1)                                                     
      REAL SMAX                                                         
      INTEGER I,INCX,IX,N                                               
      COMPLEX ZDUM                                                      
      REAL CABS1                                                        
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))                  
C                                                                       
      ICAMAX = 0                                                        
      IF( N .LT. 1 ) RETURN                                             
      ICAMAX = 1                                                        
      IF(N.EQ.1)RETURN                                                  
      IF(INCX.EQ.1)GO TO 20                                             
C                                                                       
C        CODE FOR INCREMENT NOT EQUAL TO 1                              
C                                                                       
      IX = 1                                                            
      SMAX = CABS1(CX(1))                                               
      IX = IX + INCX                                                    
      DO 10 I = 2,N                                                     
         IF(CABS1(CX(IX)).LE.SMAX) GO TO 5                              
         ICAMAX = I                                                     
         SMAX = CABS1(CX(IX))                                           
    5    IX = IX + INCX                                                 
   10 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C        CODE FOR INCREMENT EQUAL TO 1                                  
C                                                                       
   20 SMAX = CABS1(CX(1))                                               
      DO 30 I = 2,N                                                     
         IF(CABS1(CX(I)).LE.SMAX) GO TO 30                              
         ICAMAX = I                                                     
         SMAX = CABS1(CX(I))                                            
   30 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      subroutine cqrdc(x,ldx,n,p,qraux,jpvt,work,job)                   
      integer ldx,n,p,job                                               
      integer jpvt(1)                                                   
      complex x(ldx,1),qraux(1),work(1)                                 
c                                                                       
c     cqrdc uses householder transformations to compute the qr          
c     factorization of an n by p matrix x.  column pivoting             
c     based on the 2-norms of the reduced columns may be                
c     performed at the users option.                                    
c                                                                       
c     on entry                                                          
cc                                                                      
c        x       complex(ldx,p), where ldx .ge. n.                      
c                x contains the matrix whose decomposition is to be     
c                computed.                                              
c                                                                       
c        ldx     integer.                                               
c                ldx is the leading dimension of the array x.           
c                                                                       
c        n       integer.                                               
c                n is the number of rows of the matrix x.               
c                                                                       
c        p       integer.                                               
c                p is the number of columns of the matrix x.            
c                                                                       
c        jpvt    integer(p).                                            
c                jpvt contains integers that control the selection      
c                of the pivot columns.  the k-th column x(k) of x       
c                is placed in one of three classes according to the     
c                value of jpvt(k).                                      
c                                                                       
c                   if jpvt(k) .gt. 0, then x(k) is an initial          
c                                      column.                          
c                                                                       
c                   if jpvt(k) .eq. 0, then x(k) is a free column.      
c                                                                       
c                   if jpvt(k) .lt. 0, then x(k) is a final column.     
c                                                                       
c                before the decomposition is computed, initial columns  
c                are moved to the beginning of the array x and final    
c                columns to the end.  both initial and final columns    
c                are frozen in place during the computation and only    
c                free columns are moved.  at the k-th stage of the      
c                reduction, if x(k) is occupied by a free column        
c                it is interchanged with the free column of largest     
c                reduced norm.  jpvt is not referenced if               
c                job .eq. 0.                                            
c                                                                       
c        work    complex(p).                                            
c                work is a work array.  work is not referenced if       
c                job .eq. 0.                                            
c                                                                       
c        job     integer.                                               
c                job is an integer that initiates column pivoting.      
c                if job .eq. 0, no pivoting is done.                    
c                if job .ne. 0, pivoting is done.                       
c                                                                       
c     on return                                                         
c                                                                       
c        x       x contains in its upper triangle the upper             
c                triangular matrix r of the qr factorization.           
c                below its diagonal x contains information from         
c                which the unitary part of the decomposition            
c                can be recovered.  note that if pivoting has           
c                been requested, the decomposition is not that          
c                of the original matrix x but that of x                 
c                with its columns permuted as described by jpvt.        
c                                                                       
c        qraux   complex(p).                                            
c                qraux contains further information required to recover 
c                the unitary part of the decomposition.                 
c                                                                       
c        jpvt    jpvt(k) contains the index of the column of the        
c                original matrix that has been interchanged into        
c                the k-th column, if pivoting was requested.            
c                                                                       
c     linpack. this version dated 08/14/78 .                            
c     g.w. stewart, university of maryland, argonne national lab.       
c                                                                       
c     cqrdc uses the following functions and subprograms.               
c                                                                       
c     blas caxpy,cdotc,cscal,cswap,scnrm2                               
c     fortran abs,aimag,amax1,cabs,cmplx,csqrt,min0,real                
c                                                                       
c     internal variables                                                
c                                                                       
      integer j,jp,l,lp1,lup,maxj,pl,pu                                 
      real maxnrm,scnrm2,tt                                             
      complex cdotc,nrmxl,t                                             
      logical negj,swapj                                                
c                                                                       
      complex csign,zdum,zdum1,zdum2                                    
      real cabs1                                                        
      csign(zdum1,zdum2) = cabs(zdum1)*(zdum2/cabs(zdum2))              
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))                  
c                                                                       
      pl = 1                                                            
      pu = 0                                                            
      if (job .eq. 0) go to 60                                          
c                                                                       
c        pivoting has been requested.  rearrange the columns            
c        according to jpvt.                                             
c                                                                       
         do 20 j = 1, p                                                 
            swapj = jpvt(j) .gt. 0                                      
            negj = jpvt(j) .lt. 0                                       
            jpvt(j) = j                                                 
            if (negj) jpvt(j) = -j                                      
            if (.not.swapj) go to 10                                    
               if (j .ne. pl) call cswap(n,x(1,pl),1,x(1,j),1)          
               jpvt(j) = jpvt(pl)                                       
               jpvt(pl) = j                                             
               pl = pl + 1                                              
   10       continue                                                    
   20    continue                                                       
         pu = p                                                         
         do 50 jj = 1, p                                                
            j = p - jj + 1                                              
            if (jpvt(j) .ge. 0) go to 40                                
               jpvt(j) = -jpvt(j)                                       
               if (j .eq. pu) go to 30                                  
                  call cswap(n,x(1,pu),1,x(1,j),1)                      
                  jp = jpvt(pu)                                         
                  jpvt(pu) = jpvt(j)                                    
                  jpvt(j) = jp                                          
   30          continue                                                 
               pu = pu - 1                                              
   40       continue                                                    
   50    continue                                                       
   60 continue                                                          
c                                                                       
c     compute the norms of the free columns.                            
c                                                                       
      if (pu .lt. pl) go to 80                                          
      do 70 j = pl, pu                                                  
         qraux(j) = cmplx(scnrm2(n,x(1,j),1),0.0e0)                     
         work(j) = qraux(j)                                             
   70 continue                                                          
   80 continue                                                          
c                                                                       
c     perform the householder reduction of x.                           
c                                                                       
      lup = min0(n,p)                                                   
      do 200 l = 1, lup                                                 
         if (l .lt. pl .or. l .ge. pu) go to 120                        
c                                                                       
c           locate the column of largest norm and bring it              
c           into the pivot position.                                    
c                                                                       
            maxnrm = 0.0e0                                              
            maxj = l                                                    
            do 100 j = l, pu                                            
               if (real(qraux(j)) .le. maxnrm) go to 90                 
                  maxnrm = real(qraux(j))                               
                  maxj = j                                              
   90          continue                                                 
  100       continue                                                    
            if (maxj .eq. l) go to 110                                  
               call cswap(n,x(1,l),1,x(1,maxj),1)                       
               qraux(maxj) = qraux(l)                                   
               work(maxj) = work(l)                                     
               jp = jpvt(maxj)                                          
               jpvt(maxj) = jpvt(l)                                     
               jpvt(l) = jp                                             
  110       continue                                                    
  120    continue                                                       
         qraux(l) = (0.0e0,0.0e0)                                       
         if (l .eq. n) go to 190                                        
c                                                                       
c           compute the householder transformation for column l.        
c                                                                       
            nrmxl = cmplx(scnrm2(n-l+1,x(l,l),1),0.0e0)                 
            if (cabs1(nrmxl) .eq. 0.0e0) go to 180                      
               if (cabs1(x(l,l)) .ne. 0.0e0)                            
     *            nrmxl = csign(nrmxl,x(l,l))                           
               call cscal(n-l+1,(1.0e0,0.0e0)/nrmxl,x(l,l),1)           
               x(l,l) = (1.0e0,0.0e0) + x(l,l)                          
c                                                                       
c              apply the transformation to the remaining columns,       
c              updating the norms.                                      
c                                                                       
               lp1 = l + 1                                              
               if (p .lt. lp1) go to 170                                
               do 160 j = lp1, p                                        
                  t = -cdotc(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)            
                  call caxpy(n-l+1,t,x(l,l),1,x(l,j),1)                 
                  if (j .lt. pl .or. j .gt. pu) go to 150               
                  if (cabs1(qraux(j)) .eq. 0.0e0) go to 150             
                     tt = 1.0e0 - (cabs(x(l,j))/real(qraux(j)))**2      
                     tt = amax1(tt,0.0e0)                               
                     t = cmplx(tt,0.0e0)                                
                     tt = 1.0e0                                         
     *                    + 0.05e0*tt*(real(qraux(j))/real(work(j)))**2 
                     if (tt .eq. 1.0e0) go to 130                       
                        qraux(j) = qraux(j)*csqrt(t)                    
                     go to 140                                          
  130                continue                                           
                        qraux(j) = cmplx(scnrm2(n-l,x(l+1,j),1),0.0e0)  
                        work(j) = qraux(j)                              
  140                continue                                           
  150             continue                                              
  160          continue                                                 
  170          continue                                                 
c                                                                       
c              save the transformation.                                 
c                                                                       
               qraux(l) = x(l,l)                                        
               x(l,l) = -nrmxl                                          
  180       continue                                                    
  190    continue                                                       
  200 continue                                                          
      return                                                            
      end                                                               
c                                                                       
c                                                                       
c                                                                       
      real function scnrm2( n, cx, incx)                                
      logical imag, scale                                               
      integer          next                                             
      real         cutlo, cuthi, hitest, sum, xmax, absx, zero, one     
      complex      cx(1)                                                
      data         zero, one /0.0e0, 1.0e0/                             
c                                                                       
c     unitary norm of the complex n-vector stored in cx() with storage  
c     increment incx .                                                  
c     if    n .le. 0 return with result = 0.                            
c     if n .ge. 1 then incx must be .ge. 1                              
c                                                                       
c           c.l.lawson , 1978 jan 08                                    
c                                                                       
c     four phase method     using two built-in constants that are       
c     hopefully applicable to all machines.                             
c         cutlo = maximum of  sqrt(u/eps)  over all known machines.     
c         cuthi = minimum of  sqrt(v)      over all known machines.     
c     where                                                             
c         eps = smallest no. such that eps + 1. .gt. 1.                 
c         u   = smallest positive no.   (underflow limit)               
c         v   = largest  no.            (overflow  limit)               
c                                                                       
c     brief outline of algorithm..                                      
c                                                                       
c     phase 1    scans zero components.                                 
c     move to phase 2 when a component is nonzero and .le. cutlo        
c     move to phase 3 when a component is .gt. cutlo                    
c     move to phase 4 when a component is .ge. cuthi/m                  
c     where m = n for x() real and m = 2*n for complex.                 
c                                                                       
c     values for cutlo and cuthi..                                      
c     from the environmental parameters listed in the imsl converter    
c     document the limiting values are as follows..                     
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)                         
c                   thus cutlo = 2**(-51) = 4.44089e-16                 
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.          
c                   thus cuthi = 2**(63.5) = 1.30438e19                 
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.             
c                   thus cutlo = 2**(-33.5) = 8.23181d-11               
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19                    
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /                        
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /                        
      data cutlo, cuthi / 4.441e-16,  1.304e19 /                        
c                                                                       
      if(n .gt. 0) go to 10                                             
         scnrm2  = zero                                                 
         go to 300                                                      
c                                                                       
   10 assign 30 to next                                                 
      sum = zero                                                        
      nn = n * incx                                                     
c                                                 begin main loop       
      do 210 i=1,nn,incx                                                
         absx = abs(real(cx(i)))                                        
         imag = .false.                                                 
         go to next,(30, 50, 70, 90, 110)                               
   30 if( absx .gt. cutlo) go to 85                                     
      assign 50 to next                                                 
      scale = .false.                                                   
c                                                                       
c                        phase 1.  sum is zero                          
c                                                                       
   50 if( absx .eq. zero) go to 200                                     
      if( absx .gt. cutlo) go to 85                                     
c                                                                       
c                                prepare for phase 2.                   
      assign 70 to next                                                 
      go to 105                                                         
c                                                                       
c                                prepare for phase 4.                   
c                                                                       
  100 assign 110 to next                                                
      sum = (sum / absx) / absx                                         
  105 scale = .true.                                                    
      xmax = absx                                                       
      go to 115                                                         
c                                                                       
c                   phase 2.  sum is small.                             
c                             scale to avoid destructive underflow.     
c                                                                       
   70 if( absx .gt. cutlo ) go to 75                                    
c                                                                       
c                     common code for phases 2 and 4.                   
c                     in phase 4 sum is large.  scale to avoid overflow.
c                                                                       
  110 if( absx .le. xmax ) go to 115                                    
         sum = one + sum * (xmax / absx)**2                             
         xmax = absx                                                    
         go to 200                                                      
c                                                                       
  115 sum = sum + (absx/xmax)**2                                        
      go to 200                                                         
c                                                                       
c                                                                       
c                  prepare for phase 3.                                 
c                                                                       
   75 sum = (sum * xmax) * xmax                                         
c                                                                       
   85 assign 90 to next                                                 
      scale = .false.                                                   
c                                                                       
c     for real or d.p. set hitest = cuthi/n                             
c     for complex      set hitest = cuthi/(2*n)                         
c                                                                       
      hitest = cuthi/float( n )                                         
c                                                                       
c                   phase 3.  sum is mid-range.  no scaling.            
c                                                                       
   90 if(absx .ge. hitest) go to 100                                    
         sum = sum + absx**2                                            
  200 continue                                                          
c                  control selection of real and imaginary parts.       
c                                                                       
      if(imag) go to 210                                                
         absx = abs(aimag(cx(i)))                                       
         imag = .true.                                                  
      go to next,(  50, 70, 90, 110 )                                   
c                                                                       
  210 continue                                                          
c                                                                       
c              end of main loop.                                        
c              compute square root and adjust for scaling.              
c                                                                       
      scnrm2 = sqrt(sum)                                                
      if(scale) scnrm2 = scnrm2 * xmax                                  
  300 continue                                                          
      return                                                            
      end                                                               
c                                                                       
c                                                                       
      subroutine  cswap (n,cx,incx,cy,incy)                             
c                                                                       
c     interchanges two vectors.                                         
c     jack dongarra, linpack, 3/11/78.                                  
c                                                                       
      complex cx(1),cy(1),ctemp                                         
      integer i,incx,incy,ix,iy,n                                       
c                                                                       
      if(n.le.0)return                                                  
      if(incx.eq.1.and.incy.eq.1)go to 20                               
c                                                                       
c       code for unequal increments or equal increments not equal       
c         to 1                                                          
c                                                                       
      ix = 1                                                            
      iy = 1                                                            
      if(incx.lt.0)ix = (-n+1)*incx + 1                                 
      if(incy.lt.0)iy = (-n+1)*incy + 1                                 
      do 10 i = 1,n                                                     
        ctemp = cx(ix)                                                  
        cx(ix) = cy(iy)                                                 
        cy(iy) = ctemp                                                  
        ix = ix + incx                                                  
        iy = iy + incy                                                  
   10 continue                                                          
      return                                                            
c                                                                       
c       code for both increments equal to 1                             
   20 do 30 i = 1,n                                                     
        ctemp = cx(i)                                                   
        cx(i) = cy(i)                                                   
        cy(i) = ctemp                                                   
   30 continue                                                          
      return                                                            
      end                                                               
      subroutine cqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)      
      integer ldx,n,k,job,info                                          
      complex x(ldx,1),qraux(1),y(1),qy(1),qty(1),b(1),rsd(1),xb(1)     
c                                                                       
c     cqrsl applies the output of cqrdc to compute coordinate           
c     transformations, projections, and least squares solutions.        
c     for k .le. min(n,p), let xk be the matrix                         
c                                                                       
c            xk = (x(jpvt(1)),x(jpvt(2)), ... ,x(jpvt(k)))              
c                                                                       
c     formed from columnns jpvt(1), ... ,jpvt(k) of the original        
c     n x p matrix x that was input to cqrdc (if no pivoting was        
c     done, xk consists of the first k columns of x in their            
c     original order).  cqrdc produces a factored unitary matrix q      
c     and an upper triangular matrix r such that                        
c                                                                       
c              xk = q * (r)                                             
c                       (0)                                             
c                                                                       
c     this information is contained in coded form in the arrays         
c     x and qraux.                                                      
c                                                                       
c     on entry                                                          
c                                                                       
c        x      complex(ldx,p).                                         
c               x contains the output of cqrdc.                         
c                                                                       
c        ldx    integer.                                                
c               ldx is the leading dimension of the array x.            
c                                                                       
c        n      integer.                                                
c               n is the number of rows of the matrix xk.  it must      
c               have the same value as n in cqrdc.                      
c                                                                       
c        k      integer.                                                
c               k is the number of columns of the matrix xk.  k         
c               must nnot be greater than min(n,p), where p is the      
c               same as in the calling sequence to cqrdc.               
c                                                                       
c        qraux  complex(p).                                             
c               qraux contains the auxiliary output from cqrdc.         
c                                                                       
c        y      complex(n)                                              
c               y contains an n-vector that is to be manipulated        
c               by cqrsl.                                               
c                                                                       
c        job    integer.                                                
c               job specifies what is to be computed.  job has          
c               the decimal expansion abcde, with the following         
c               meaning.                                                
c                                                                       
c                    if a.ne.0, compute qy.                             
c                    if b,c,d, or e .ne. 0, compute qty.                
c                    if c.ne.0, compute b.                              
c                    if d.ne.0, compute rsd.                            
c                    if e.ne.0, compute xb.                             
c                                                                       
c               note that a request to compute b, rsd, or xb            
c               automatically triggers the computation of qty, for      
c               which an array must be provided in the calling          
c               sequence.                                               
c                                                                       
c     on return                                                         
c                                                                       
c        qy     complex(n).                                             
c               qy conntains q*y, if its computation has been           
c               requested.                                              
c                                                                       
c        qty    complex(n).                                             
c               qty contains ctrans(q)*y, if its computation has        
c               been requested.  here ctrans(q) is the conjugate        
c               transpose of the matrix q.                              
c                                                                       
c        b      complex(k)                                              
c               b contains the solution of the least squares problem    
c                                                                       
c                    minimize norm2(y - xk*b),                          
c                                                                       
c               if its computation has been requested.  (note that      
c               if pivoting was requested in cqrdc, the j-th            
c               component of b will be associated with column jpvt(j)   
c               of the original matrix x that was input into cqrdc.)    
c                                                                       
c        rsd    complex(n).                                             
c               rsd contains the least squares residual y - xk*b,       
c               if its computation has been requested.  rsd is          
c               also the orthogonal projection of y onto the            
c               orthogonal complement of the column space of xk.        
c                                                                       
c        xb     complex(n).                                             
c               xb contains the least squares approximation xk*b,       
c               if its computation has been requested.  xb is also      
c               the orthogonal projection of y onto the column space    
c               of x.                                                   
c                                                                       
c        info   integer.                                                
c               info is zero unless the computation of b has            
c               been requested and r is exactly singular.  in           
c               this case, info is the index of the first zero          
c               diagonal element of r and b is left unaltered.          
c                                                                       
c     the parameters qy, qty, b, rsd, and xb are not referenced         
c     if their computation is not requested and in this case            
c     can be replaced by dummy variables in the calling program.        
c     to save storage, the user may in some cases use the same          
c     array for different parameters in the calling sequence.  a        
c     frequently occuring example is when one wishes to compute         
c     any of b, rsd, or xb and does not need y or qty.  in this         
c     case one may identify y, qty, and one of b, rsd, or xb, while     
c     providing separate arrays for anything else that is to be         
c     computed.  thus the calling sequence                              
c                                                                       
c          call cqrsl(x,ldx,n,k,qraux,y,dum,y,b,y,dum,110,info)         
c                                                                       
c     will result in the computation of b and rsd, with rsd             
c     overwriting y.  more generally, each item in the following        
c     list contains groups of permissible identifications for           
c     a single callinng sequence.                                       
c                                                                       
c          1. (y,qty,b) (rsd) (xb) (qy)                                 
c                                                                       
c          2. (y,qty,rsd) (b) (xb) (qy)                                 
c                                                                       
c          3. (y,qty,xb) (b) (rsd) (qy)                                 
c                                                                       
c          4. (y,qy) (qty,b) (rsd) (xb)                                 
c                                                                       
c          5. (y,qy) (qty,rsd) (b) (xb)                                 
c                                                                       
c          6. (y,qy) (qty,xb) (b) (rsd)                                 
c                                                                       
c     in any group the value returned in the array allocated to         
c     the group corresponds to the last member of the group.            
c                                                                       
c     linpack. this version dated 08/14/78 .                            
c     g.w. stewart, university of maryland, argonne national lab.       
c                                                                       
c     cqrsl uses the following functions and subprograms.               
c                                                                       
c     blas caxpy,ccopy,cdotc                                            
c     fortran abs,aimag,min0,mod,real                                   
c                                                                       
c     internal variables                                                
c                                                                       
      integer i,j,jj,ju,kp1                                             
      complex cdotc,t,temp                                              
      logical cb,cqy,cqty,cr,cxb                                        
c                                                                       
      complex zdum                                                      
      real cabs1                                                        
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))                  
c                                                                       
c     set info flag.                                                    
c                                                                       
      info = 0                                                          
c                                                                       
c     determine what is to be computed.                                 
c                                                                       
      cqy = job/10000 .ne. 0                                            
      cqty = mod(job,10000) .ne. 0                                      
      cb = mod(job,1000)/100 .ne. 0                                     
      cr = mod(job,100)/10 .ne. 0                                       
      cxb = mod(job,10) .ne. 0                                          
      ju = min0(k,n-1)                                                  
c                                                                       
c     special action when n=1.                                          
c                                                                       
      if (ju .ne. 0) go to 40                                           
         if (cqy) qy(1) = y(1)                                          
         if (cqty) qty(1) = y(1)                                        
         if (cxb) xb(1) = y(1)                                          
         if (.not.cb) go to 30                                          
            if (cabs1(x(1,1)) .ne. 0.0e0) go to 10                      
               info = 1                                                 
            go to 20                                                    
   10       continue                                                    
               b(1) = y(1)/x(1,1)                                       
   20       continue                                                    
   30    continue                                                       
         if (cr) rsd(1) = (0.0e0,0.0e0)                                 
      go to 250                                                         
   40 continue                                                          
c                                                                       
c        set up to compute qy or qty.                                   
c                                                                       
         if (cqy) call ccopy(n,y,1,qy,1)                                
         if (cqty) call ccopy(n,y,1,qty,1)                              
         if (.not.cqy) go to 70                                         
c                                                                       
c           compute qy.                                                 
c                                                                       
            do 60 jj = 1, ju                                            
               j = ju - jj + 1                                          
               if (cabs1(qraux(j)) .eq. 0.0e0) go to 50                 
                  temp = x(j,j)                                         
                  x(j,j) = qraux(j)                                     
                  t = -cdotc(n-j+1,x(j,j),1,qy(j),1)/x(j,j)             
                  call caxpy(n-j+1,t,x(j,j),1,qy(j),1)                  
                  x(j,j) = temp                                         
   50          continue                                                 
   60       continue                                                    
   70    continue                                                       
         if (.not.cqty) go to 100                                       
c                                                                       
c           compute ctrans(q)*y.                                        
c                                                                       
            do 90 j = 1, ju                                             
               if (cabs1(qraux(j)) .eq. 0.0e0) go to 80                 
                  temp = x(j,j)                                         
                  x(j,j) = qraux(j)                                     
                  t = -cdotc(n-j+1,x(j,j),1,qty(j),1)/x(j,j)            
                  call caxpy(n-j+1,t,x(j,j),1,qty(j),1)                 
                  x(j,j) = temp                                         
   80          continue                                                 
   90       continue                                                    
  100    continue                                                       
c                                                                       
c        set up to compute b, rsd, or xb.                               
c                                                                       
         if (cb) call ccopy(k,qty,1,b,1)                                
         kp1 = k + 1                                                    
         if (cxb) call ccopy(k,qty,1,xb,1)                              
         if (cr .and. k .lt. n) call ccopy(n-k,qty(kp1),1,rsd(kp1),1)   
         if (.not.cxb .or. kp1 .gt. n) go to 120                        
            do 110 i = kp1, n                                           
               xb(i) = (0.0e0,0.0e0)                                    
  110       continue                                                    
  120    continue                                                       
         if (.not.cr) go to 140                                         
            do 130 i = 1, k                                             
               rsd(i) = (0.0e0,0.0e0)                                   
  130       continue                                                    
  140    continue                                                       
         if (.not.cb) go to 190                                         
c                                                                       
c           compute b.                                                  
c                                                                       
            do 170 jj = 1, k                                            
               j = k - jj + 1                                           
               if (cabs1(x(j,j)) .ne. 0.0e0) go to 150                  
                  info = j                                              
c           ......exit                                                  
                  go to 180                                             
  150          continue                                                 
               b(j) = b(j)/x(j,j)                                       
               if (j .eq. 1) go to 160                                  
                  t = -b(j)                                             
                  call caxpy(j-1,t,x(1,j),1,b,1)                        
  160          continue                                                 
  170       continue                                                    
  180       continue                                                    
  190    continue                                                       
         if (.not.cr .and. .not.cxb) go to 240                          
c                                                                       
c           compute rsd or xb as required.                              
c                                                                       
            do 230 jj = 1, ju                                           
               j = ju - jj + 1                                          
               if (cabs1(qraux(j)) .eq. 0.0e0) go to 220                
                  temp = x(j,j)                                         
                  x(j,j) = qraux(j)                                     
                  if (.not.cr) go to 200                                
                     t = -cdotc(n-j+1,x(j,j),1,rsd(j),1)/x(j,j)         
                     call caxpy(n-j+1,t,x(j,j),1,rsd(j),1)              
  200             continue                                              
                  if (.not.cxb) go to 210                               
                     t = -cdotc(n-j+1,x(j,j),1,xb(j),1)/x(j,j)          
                     call caxpy(n-j+1,t,x(j,j),1,xb(j),1)               
  210             continue                                              
                  x(j,j) = temp                                         
  220          continue                                                 
  230       continue                                                    
  240    continue                                                       
  250 continue                                                          
      return                                                            
      end                                                               
c                                                                       
c                                                                       
c                                                                       
c                                                                       
      subroutine  ccopy(n,cx,incx,cy,incy)                              
c                                                                       
c     copies a vector, x, to a vector, y.                               
c     jack dongarra, linpack, 3/11/78.                                  
c                                                                       
      complex cx(1),cy(1)                                               
      integer i,incx,incy,ix,iy,n                                       
c                                                                       
      if(n.le.0)return                                                  
      if(incx.eq.1.and.incy.eq.1)go to 20                               
c                                                                       
c        code for unequal increments or equal increments                
c          not equal to 1                                               
c                                                                       
      ix = 1                                                            
      iy = 1                                                            
      if(incx.lt.0)ix = (-n+1)*incx + 1                                 
      if(incy.lt.0)iy = (-n+1)*incy + 1                                 
      do 10 i = 1,n                                                     
        cy(iy) = cx(ix)                                                 
        ix = ix + incx                                                  
        iy = iy + incy                                                  
   10 continue                                                          
      return                                                            
c                                                                       
c        code for both increments equal to 1                            
c                                                                       
   20 do 30 i = 1,n                                                     
        cy(i) = cx(i)                                                   
   30 continue                                                          
      return                                                            
      end                                                               
