      SUBROUTINE SHAPE(XI,N,PSI,DPSI)
C
C---- SHAPE FUNCTIONS FOR BOUNDARY ELEMENTS
C
      DIMENSION PSI(4),DPSI(4)
      GO TO (10,20,30) N
C
C     Linear elements
C
 10   PSI(1)=0.5*(1.0-XI)
      PSI(2)=0.5*(1.0+XI)
      DPSI(1)=-0.5
      DPSI(2)=0.5
      RETURN
C
C     Quadratic elements
C
 20   PSI(1)=0.5*XI*(XI-1.0)
      PSI(2)=1.0-XI**2
      PSI(3)=0.5*XI*(XI+1.0)
      DPSI(1)=XI-0.5
      DPSI(2)=-2.0*XI
      DPSI(3)=XI+0.5
      RETURN
C
C     Cubic elements
C
30    PSI(1)=9./16.*(1./9.-XI**2)*(XI-1.0)
      PSI(2)=27./16.*(1.0-XI**2)*(1./3.-XI)
      PSI(3)=27./16.*(1.0-XI**2)*(1./3.+XI)
      PSI(4)=-9./16.*(1./9.-XI**2)*(XI+1.0)
      DPSI(1)=-9./16.*(3.*XI**2-2.*XI-1./9.)
      DPSI(2)=27./16.*(3.*XI**2-2./3.*XI-1.)
      DPSI(3)=27./16.*(-3.*XI**2-2./3.*XI+1.)
      DPSI(4)=-9./16.*(-3.*XI**2-2.*XI+1./9.)
      RETURN
      END