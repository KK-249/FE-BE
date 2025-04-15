function [H,G,CP,AA,QN]=elem4(XP,YP,ZP,XQ,YQ,ZQ)

C1 = 1.0/12.56637061;
WN = 2*3.14*5/1500;         %波数=角频率除以声速
IK = 1i*WN;
No_point1 = 4;              %高斯积分点
No_point2 = 4;


H = zeros(4,1);
G = zeros(4,1);
CP = 1;
AA = zeros(4,1);

[point1,weight1]=GaussPoint1(No_point1);
[point2,weight2]=GaussPoint1(No_point2);

for i = 1:size(point1,2)
    for j = 1:size(point2,2)
        
        [PSI,DPSI] = shape(point1(i),point2(j),size(XQ,1));
          XX=0.0;
          YY=0.0;
          ZZ=0.0;
          DXDS1=0.0;
          DXDS2=0.0;
          DYDS1=0.0;
          DYDS2=0.0;
          DZDS1=0.0;
          DZDS2=0.0;
          for ii = 1:size(XQ,1)
              XX=XX+XQ(ii,1)*PSI(ii);
              YY=YY+YQ(ii,1)*PSI(ii);
              ZZ=ZZ+ZQ(ii,1)*PSI(ii);
              DXDS1=DXDS1+XQ(ii)*DPSI(1,ii);
              DYDS1=DYDS1+YQ(ii)*DPSI(1,ii);
              DZDS1=DZDS1+ZQ(ii)*DPSI(1,ii);
              DXDS2=DXDS2+XQ(ii)*DPSI(2,ii);
              DYDS2=DYDS2+YQ(ii)*DPSI(2,ii);
              DZDS2=DZDS2+ZQ(ii)*DPSI(2,ii);	 
          end
            
            NN1=DYDS1*DZDS2-DZDS1*DYDS2;
            NN2=-(DXDS1*DZDS2-DZDS1*DXDS2);
            NN3=DXDS1*DYDS2-DYDS1*DXDS2;
            DETJ=sqrt(NN1*NN1+NN2*NN2+NN3*NN3);         
            QN(1)=NN1/DETJ;                         
            QN(2)=NN2/DETJ;
            QN(3)=NN3/DETJ;
            RX=XX-XP;                                
            RY=YY-YP;
            RZ=ZZ-ZP;
            R=sqrt(RX*RX+RY*RY+RZ*RZ);               
            DRDN=(RX*QN(1)+RY*QN(2)+RZ*QN(3))/R;     
            FACTOR=DETJ*weight1(i)*weight2(j);           
          
            GREEN=C1*exp(-IK*R)/R*FACTOR;          
            DGDN=-(1.0+IK*R)*GREEN/R*DRDN;          
            DGDN0=-C1/(R*R)*DRDN*FACTOR;
          
            CP=CP-DGDN0; 
            
            for jj = 1:4
                H(jj,1)=H(jj,1)+PSI(jj)*DGDN;
                G(jj,1)=G(jj,1)+PSI(jj)*GREEN;
                AA(jj,1)=AA(jj,1)+PSI(jj)*FACTOR; 
            end
    end
end



