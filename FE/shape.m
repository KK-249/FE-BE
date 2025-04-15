function [PSI,DPSI] = shape(XI1,XI2,I34)

if I34 == 4
    PSI(1)=0.25*(1.-XI1)*(1.-XI2);
    PSI(2)=0.25*(1.+XI1)*(1.-XI2);
    PSI(3)=0.25*(1.+XI1)*(1.+XI2);
    PSI(4)=0.25*(1.-XI1)*(1.+XI2);
    DPSI(1,1)=0.25*(XI2-1.);    
    DPSI(1,2)=0.25*(1.-XI2);
    DPSI(1,3)=0.25*(1.+XI2);
    DPSI(1,4)=0.25*(-1.-XI2);
    DPSI(2,1)=0.25*(XI1-1.);
    DPSI(2,2)=0.25*(-1.-XI1);
    DPSI(2,3)=0.25*(1.+XI1);
    DPSI(2,4)=0.25*(1.-XI1);
end

if I34 == 3
    PSI(1)=XI1;
    PSI(2)=XI2;
    PSI(3)=1.-XI1-XI2;
    DPSI(1,1)=1.0;
    DPSI(1,2)=0.0;
    DPSI(1,3)=-1.0;
    DPSI(2,1)=0.0;
    DPSI(2,2)=1.0;
    DPSI(2,3)=-1.0; 
end