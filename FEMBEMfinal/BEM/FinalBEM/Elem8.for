      SUBROUTINE ELEM8(XP,YP,ZP,XQ,YQ,ZQ,H,G,CP,AA,QN)
C
C  Formulate element coefficient matrices
C
      COMMON /FRE/FREQ,WN

      DIMENSION XII(12),WT(12)
      DIMENSION XQ(8),YQ(8),ZQ(8),QN(3),PSI(8),DPSI(2,8)   
      DIMENSION H(8),G(8),AA(8)
      REAL NN1,NN2,NN3
      COMPLEX H,G,GREEN,DGDN,IK

      C1=1.0/12.56637061
      IK=CMPLX(0.,1.)*WN

      DO 100 I=1,8
          H(I)=0.0
          G(I)=0.0
         AA(I)=0.0
100   CONTINUE
      CALL GETINT(NINP,XII,WT)
!     调用 GETINT 子程序获取积分点数量     
C
C  Integration loop
C
      DO 450 INP1=1,NINP                  !1~4
	 DO 400 INP2=1,NINP
	  XI1=XII(INP1)
	  XI2=XII(INP2)
          CALL SHAPE(XI1,XI2,PSI,DPSI,8)
          XX=0.0
          YY=0.0
	    ZZ=0.0
          DXDS1=0.0
	    DXDS2=0.0
          DYDS1=0.0
	    DYDS2=0.0
	    DZDS1=0.0
	    DZDS2=0.0
	  DO 200 I=1,8
!计算全局坐标 以及对各个局部坐标的导数
          XX=XX+XQ(I)*PSI(I)
          YY=YY+YQ(I)*PSI(I)
          ZZ=ZZ+ZQ(I)*PSI(I)
          DXDS1=DXDS1+XQ(I)*DPSI(1,I)
	    DYDS1=DYDS1+YQ(I)*DPSI(1,I)
	    DZDS1=DZDS1+ZQ(I)*DPSI(1,I)
	    DXDS2=DXDS2+XQ(I)*DPSI(2,I)
          DYDS2=DYDS2+YQ(I)*DPSI(2,I)
	    DZDS2=DZDS2+ZQ(I)*DPSI(2,I)		 
200       CONTINUE
          NN1=DYDS1*DZDS2-DZDS1*DYDS2
	    NN2=-(DXDS1*DZDS2-DZDS1*DXDS2)
	    NN3=DXDS1*DYDS2-DYDS1*DXDS2
          DETJ=SQRT(NN1**2+NN2**2+NN3**2)         !雅可比的值
	    QN(1)=NN1/DETJ                          !单位法向量
	    QN(2)=NN2/DETJ
	    QN(3)=NN3/DETJ
          RX=XX-XP                                !声源或观测点到积分点的距离
          RY=YY-YP
	    RZ=ZZ-ZP
          R=SQRT(RX*RX+RY*RY+RZ*RZ)               !矢径
          DRDN=(RX*QN(1)+RY*QN(2)+RZ*QN(3))/R     !距离R沿着单元法向方向的方向导数  dR/dn
          FACTOR=DETJ*WT(INP1)*WT(INP2)           
          
          GREEN=C1*CEXP(-IK*R)/R*FACTOR           !cexp  复指数
          DGDN=-(1.0+IK*R)*GREEN/R*DRDN           !法向导数
          DGDN0=-C1/(R*R)*DRDN*FACTOR
          
          CP=CP-DGDN0 
          
          DO 300 I=1,8
            H(I)=H(I)+PSI(I)*DGDN
            G(I)=G(I)+PSI(I)*GREEN
            AA(I)=AA(I)+PSI(I)*FACTOR 
300       CONTINUE
           
      
        
400     CONTINUE
450   CONTINUE
      RETURN
      END