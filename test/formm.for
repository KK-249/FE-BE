      SUBROUTINE FORM
C
C  This subroutine formulates the coefficient matrix
C
      COMMON /PROB/INT,ISC,ISYM
      COMMON /FLUID/RHO,C
      COMMON /FRE/FREQ,WN,IFR
      COMMON /BG/BIG
      COMMON /CONTR/NNODE,NELEM
      COMMON /NODES/X(3,3000)
      COMMON /ELEM/NODE(8,3000)
      COMMON /MATRX/A(3000,3000),B(3000),XS(3000)
      COMMON /SOL/SP(3000),VN(4,3000)
      COMMON /CHI/NCH,XCH(3,20)
      COMMON /SCA/AM,AF,BE,GA

      DIMENSION XQ(8),YQ(8),ZQ(8),H(8),G(8),P(3),QN(3)
      DIMENSION A_sys(4000,4000),G_sys(4000,4000),CP_sys(4000,4000)
      COMPLEX H_sys(4000,4000)
      DIMENSION AA(8),T(8000,3000)
      
      CHARACTER*1 CHECK
C
      COMPLEX BIG,A,B,XS,SP,VN,H,G,IK,IRHOW
C     WN Wavenumbr 波数
      IK=CMPLX(0.,1.)*WN
      IRHOW=CMPLX(0.,1.)*RHO*6.283185*FREQ            !IRHOW=i*rho*w=i*rho*2pi*f 振动 波传播过程中常见
C
      H_sys(:,:)=0
      A_sys(:,:)=0
      G_sys(:,:)=0
      CP_sys(:,:)=0
      
      DO 100 I=1,NNODE+NCH        !B(1:26)=0  A(1:26,1:26)=0
        B(I)=0.0
        DO 100 J=1,NNODE
          A(I,J)=0.0
100   CONTINUE

       NLOOP=1
       IF(ISYM.NE.0) NLOOP=2      !无论关于哪个平面对称 NLPOOP=2
C
C  P LOOP
C
       DO 500 IP=1,NNODE           ! P被赋予X的值 即各个节点的物理坐标
           P(1)=X(1,IP)
           P(2)=X(2,IP)
           P(3)=X(3,IP)
           CP=0.0
       IF(INT.EQ.0) CP=1.0        !INT != 0 外场问题
C
C  Symmetry Loop
C
       DO 450 LOOP=1,NLOOP
        IF(LOOP.EQ.2) P(ISYM)=-X(ISYM,IP)
        XP=P(1)
        YP=P(2)
        ZP=P(3)
        CHECK='N'
        IF(LOOP.EQ.1) CHECK='Y'
        IF(LOOP.EQ.2.AND.P(ISYM).EQ.0.0) CHECK='Y'
C  ELEMENT LOOP
C
       DO 400 K=1,NELEM
          ISING=0
          I34=4
       IF(NODE(3,K).EQ.NODE(4,K)) I34=3           
       DO 300 J=1,I34
          IQ=NODE(J,K)
           IF(IQ.EQ.IP.AND.CHECK.EQ.'Y') ISING=J  !IQ为单元中的一个节点编号 IP也为节点编号
             XQ(J)=X(1,IQ)
             YQ(J)=X(2,IQ)
             ZQ(J)=X(3,IQ)
300    CONTINUE
       
      !修改部分  
       
       
            
        IF(I34.EQ.4) THEN
            IF(ISING.EQ.0) THEN
              CALL ELEM4(XP,YP,ZP,XQ,YQ,ZQ,H,G,CP,AA,QN)
              ELSE
             CALL SING4(XP,YP,ZP,XQ,YQ,ZQ,H,G,CP,ISING)
            END IF
          ELSE
            IF(ISING.EQ.0) THEN
              CALL ELEM3(XP,YP,ZP,XQ,YQ,ZQ,H,G,CP)
            ELSE
              CALL SING3(XP,YP,ZP,XQ,YQ,ZQ,H,G,CP,ISING)
            END IF
          END IF
          
         
          DO 310 I=1,I34
          G(I)=-IRHOW*G(I)/(CMPLX(0.,1.)*6.283185*FREQ)        !数值离散积分后的G
310       CONTINUE
   
          CALL ASSMB(K,IP,H,G,I34)
          
          
          
             DO n = 1,4
             
       A_sys(NODE(n,K),NODE(n,K))=A_sys(NODE(n,K),NODE(n, K))+AA(n)
       H_sys(IP,NODE(n,K))=H_sys(IP,NODE(n,K))+H(n)
       G_sys(IP,NODE(n,K))=G_sys(IP,NODE(n,K))+G(n)
       !T(3*NODE(n,K)-2,NODE(n,K))=QN(1)
       !T(3*NODE(n,K)-1,NODE(n,K))=QN(2)
       !T(3*NODE(n,K)-0,NODE(n,K))=QN(3)
       
             END DO
          
             CP_sys(IP,IP)=CP_sys(IP,IP)+CP

          
400      CONTINUE
450    CONTINUE
       
      !H_sys(:,:)=H_sys(:,:)+CP_sys(:,:)
      
C
C   Contribution from C(P)
C
        IF(SP(IP).EQ.BIG) THEN
          A(IP,IP)=A(IP,IP)+CP
        ELSE
          B(IP)=B(IP)-CP*SP(IP)
        END IF


        IF(ISC.EQ.1) THEN
          XP=X(1,IP)
          YP=X(2,IP)
          ZP=X(3,IP)
          B(IP)=B(IP)+AM*CEXP(-IK*(AF*XP+BE*YP+GA*ZP))
        END IF

500    CONTINUE  
       
      unit = 10
      OPEN(unit,FILE='mat_outputa.txt',STATUS='REPLACE',ACTION='WRITE')
          DO i = 1, NNODE
             WRITE(unit, '( *(F10.4,1X) )') (A_sys(i,j), j = 1, NNODE)
          END DO
      CLOSE(unit)
      
      unit1 = 11
      unit4 = 14
      OPEN(unit1,FILE='mat_outputh.txt',STATUS='REPLACE',ACTION='WRITE')
      OPEN(unit4,FILE='mat_outputh1.txt',STATUS='REPLACE',ACTION='WRITE'
     1 )
          
      DO i = 1, NNODE
        ! 写入实部
       WRITE(unit1, '( *(F20.8,1X) )') (REAL(H_sys(i,j)), j = 1, NNODE)
        ! 写入虚部
       WRITE(unit4, '( *(F20.8,1X) )') (AIMAG(H_sys(i,j)), j = 1, NNODE)
      END DO
      
      CLOSE(unit1)
      CLOSE(unit4)
       

      
      unit2 = 12
      OPEN(unit2,FILE='mat_outputb.txt',STATUS='REPLACE',ACTION='WRITE')
          DO i = 1, NNODE
             WRITE(unit2, '( *(F20.4,1X) )') (G_sys(i,j), j = 1, NNODE)
          END DO
      CLOSE(unit2)
      
      unit3 = 13
      OPEN(unit3,FILE='mat_outputT.txt',STATUS='REPLACE',ACTION='WRITE')
          DO i = 1, 3*NNODE
             WRITE(unit3, '( *(F20.4,1X) )') (T(i,j), j = 1, NNODE)
          END DO
      CLOSE(unit3)
      
      unit5 = 15
      OPEN(unit5,FILE='mat_outputCP.txt',STATUS='REPLACE',ACTION='WRITE'
     1 )
          DO i = 1, NNODE
             WRITE(unit5, '( *(F20.4,1X) )') (CP_sys(i,j), j = 1, NNODE)
          END DO
      CLOSE(unit5)
      
      PRINT *, 'Matrices have been written to matrices_output.txt'
      IF(NCH.NE.0) CALL CHIEF

      RETURN
      END
