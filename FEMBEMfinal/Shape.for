      SUBROUTINE SHAPE(XI1,XI2,PSI,DPSI,I34)
C
C---- SHAPE FUNCTIONS FOR BOUNDARY ELEMENTS
C
      DIMENSION PSI(8),DPSI(2,8)
C
C     4-NODE QUADRILATERAL ELEMENT
C     
      IF(I34.EQ.8) THEN
        PSI(1)=0.25*(1+XI1)*(1+XI2)*(XI1+XI2-1)
        PSI(2)=0.25*(1-XI1)*(1+XI2)*(-XI1+XI2-1)
        PSI(3)=0.25*(1-XI1)*(1-XI2)*(-XI1-XI2-1)
        PSI(4)=0.25*(1+XI1)*(1-XI2)*(XI1-XI2-1)
        PSI(5)=0.5*(1-XI1*XI1)*(1+XI2)
        PSI(6)=0.5*(1-XI1)*(1-XI2*XI2)
        PSI(7)=0.5*(1-XI1*XI1)*(1-XI2)
        PSI(8)=0.5*(1+XI1)*(1-XI2*XI2)
        DPSI(1,1)=0.25*(1+XI2)*(XI1+XI2-1)+0.25*(1+XI1)*(1+XI2)
        DPSI(1,2)=-0.25*(1+XI2)*(-XI1+XI2-1)-0.25*(1-XI1)*(1+XI2)
        DPSI(1,3)=-0.25*(1-XI2)*(-XI1-XI2-1)-0.25*(1-XI1)*(1-XI2)
        DPSI(1,4)=0.25*(1-XI2)*(XI1-XI2-1)+0.25*(1+XI1)*(1-XI2)
        DPSI(1,5)=-XI1*(1+XI2)
        DPSI(1,6)=-0.5+0.5*XI2*XI2
        DPSI(1,7)=-XI1*(1-XI2)
        DPSI(1,8)=0.5-0.5*XI2*XI2
        DPSI(2,1)=0.25*(1+XI1)*(XI1+XI2-1)+0.25*(1+XI1)*(1+XI2)
        DPSI(2,2)=0.25*(1-XI1)*(-XI1+XI2-1)+0.25*(1-XI1)*(1+XI2)
        DPSI(2,3)=-0.25*(1-XI1)*(-XI1-XI2-1)-0.25*(1-XI1)*(1-XI2)
        DPSI(2,4)=-0.25*(1+XI1)*(XI1-XI2-1)-0.25*(1+XI1)*(1-XI2)
        DPSI(2,5)=0.5-0.5*XI1*XI1
        DPSI(2,6)=-2*(0.5-0.5*XI1)*XI2
        DPSI(2,7)=-0.5+0.5*XI1*XI1
        DPSI(2,8)=-2*(0.5+0.5*XI1)*XI2
      END IF   
          
      IF(I34.EQ.4) THEN
        PSI(1)=0.25*(1.-XI1)*(1.-XI2)
	PSI(2)=0.25*(1.+XI1)*(1.-XI2)
	PSI(3)=0.25*(1.+XI1)*(1.+XI2)
	PSI(4)=0.25*(1.-XI1)*(1.+XI2)
	DPSI(1,1)=0.25*(XI2-1.)
	DPSI(1,2)=0.25*(1.-XI2)
	DPSI(1,3)=0.25*(1.+XI2)
	DPSI(1,4)=0.25*(-1.-XI2)
	DPSI(2,1)=0.25*(XI1-1.)
	DPSI(2,2)=0.25*(-1.-XI1)
	DPSI(2,3)=0.25*(1.+XI1)
	DPSI(2,4)=0.25*(1.-XI1)
      END IF
C
C     3-NODE TRIANGULAR ELEMENT
C
      IF(I34.EQ.3) THEN
        PSI(1)=XI1
	PSI(2)=XI2
	PSI(3)=1.-XI1-XI2
	DPSI(1,1)=1.0
	DPSI(1,2)=0.0
	DPSI(1,3)=-1.0
	DPSI(2,1)=0.0
	DPSI(2,2)=1.0
	DPSI(2,3)=-1.0
      END IF
      RETURN
      END