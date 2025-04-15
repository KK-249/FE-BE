clc
clear
No_INTpoint_x=3;        %高斯积分的点数
No_INTpoint_y=3;
No_INTpoint_z=3;
A=0.0025;
a=0.01;
t=0.001;
R0=0.04115;
% R0=0.07;
Theta=9*pi/180;
% b=0.5*R0*Theta;
density=2800;
Niu=0.33;
E=70e9;
h=t;
alpha=30*pi/180;
% alpha=0*pi/180;

densityp=5440;
Niu1p=0.31;
Niu2p=0.162;
E1p=30.336e9;
E2p=15.875e9;
Gp=5.515e9;
hp=0.0003;H=h+1*hp;

% bt=1;at=(0.5*h)/(0.5*h+hp);
% bm=(0.5*h)/(0.5*h+hp);am=-(0.5*h)/(0.5*h+hp);
% bb=-(0.5*h)/(0.5*h+hp);ab=-1;

bt=1;at=(0.5*h)/(0.5*h+hp);
bm1=(0.5*h)/(0.5*h+hp);am1=-1;
bm2=1;am2=-(0.5*h)/(0.5*h+hp);
bb=-(0.5*h)/(0.5*h+hp);ab=-1;
% AA=2*a*b;

% e31=-6.6;e32=-6.6;
d31=-1.7e-10;d32=-1.7e-10;
d33=4e-10;
e24=17;e15=17;

E33=1.3e-8;
npa=2;
% nps=9;
% nn=[0];
% nn=[450,521,550,580,605,625,635,655];
nn=[0];
% nna=[450,521,550,580,610,681,710,740];
% nna=[521,681];
nna=[0];
% nns=[450,475,521,550,580,605,625,635,655];
% nns=[450,521,550,580,581,605,625,635,655];
nns=[0];
% nns=[1:800];
nps=size(nns,2);
jdxp=10;
jdxq=10;
jdx=jdxp+jdxq;               %number of nodes in x direction
R=R0+2*A*jdxp*tan(alpha);
b=R*Theta;
AA=2*a*b;
jdx1=2*jdx+1;  %41
jdy=40;               %number of nodes in y direction
jdy1=2*jdy;    %80
jdy2=jdy;      %40
k(1:5*(jdy1*(jdx+1)+jdy*jdx)-5*jdy1,1:5*(jdy1*(jdx+1)+jdy*jdx)-5*jdy1)=0;     %system stiffness matrix
m(1:5*(jdy1*(jdx+1)+jdy*jdx)-5*jdy1,1:5*(jdy1*(jdx+1)+jdy*jdx)-5*jdy1)=0;    %system mass matrix
ka(1:5*(jdy1*(jdx+1)+jdy*jdx)-5*jdy1,npa)=0; 
ks(nps,1:5*(jdy1*(jdx+1)+jdy*jdx)-5*jdy1)=0; 
en(1:(jdx)*(jdy),1:8)=0;       %element node  
for ni=1:jdy
    for nj=1:jdx
        if ni~=jdy
        en(ni+(nj-1)*(jdy),4)=(2*ni-1)+(nj-1)*jdy1+(nj-1)*jdy2;
        en(ni+(nj-1)*(jdy),7)=2*ni+(nj-1)*jdy1+(nj-1)*jdy2;
        en(ni+(nj-1)*(jdy),3)=2*ni+1+(nj-1)*jdy1+(nj-1)*jdy2;
        en(ni+(nj-1)*(jdy),1)=(2*ni-1)+(nj)*jdy1+(nj)*jdy2;
        en(ni+(nj-1)*(jdy),5)=2*ni+(nj)*jdy1+(nj)*jdy2;
        en(ni+(nj-1)*(jdy),2)=2*ni+1+(nj)*jdy1+(nj)*jdy2;
        
        else if ni==jdy
        en(ni+(nj-1)*(jdy),4)=(2*ni-1)+(nj-1)*jdy1+(nj-1)*jdy2;
        en(ni+(nj-1)*(jdy),7)=2*ni+(nj-1)*jdy1+(nj-1)*jdy2;
        en(ni+(nj-1)*(jdy),3)=1+(nj-1)*jdy1+(nj-1)*jdy2;
        en(ni+(nj-1)*(jdy),1)=(2*ni-1)+(nj)*jdy1+(nj)*jdy2;
        en(ni+(nj-1)*(jdy),5)=2*ni+(nj)*jdy1+(nj)*jdy2;
        en(ni+(nj-1)*(jdy),2)=1+(nj)*jdy1+(nj)*jdy2;
      
            end
        end
    end
end
for ni=1:jdy
    for nj=1:jdx
        if ni~=jdy
        en(ni+(nj-1)*(jdy),8)=ni+(nj)*jdy1+(nj-1)*jdy2;
        en(ni+(nj-1)*(jdy),6)=ni+1+(nj)*jdy1+(nj-1)*jdy2;

        
        else if ni==jdy
        en(ni+(nj-1)*(jdy),8)=ni+(nj)*jdy1+(nj-1)*jdy2;
        en(ni+(nj-1)*(jdy),6)=1+(nj)*jdy1+(nj-1)*jdy2;
        
            end
        end
    end
end

for ni=1:jdy1  %80 圆锥壳的全局坐标上表面
    for nj=1:jdxp+1      %11
    
        coo(ni+(nj-1)*(jdy1+jdy2),1)=2*A*(nj-1)-0.5*h*sin(alpha);       
        coo(ni+(nj-1)*(jdy1+jdy2),2)=((R0+2*A*(nj-1)*tan(alpha))+0.5*h*cos(alpha))*sin(0.5*Theta*(ni-1));   
        coo(ni+(nj-1)*(jdy1+jdy2),3)=((R0+2*A*(nj-1)*tan(alpha))+0.5*h*cos(alpha))*cos(0.5*Theta*(ni-1));
%         ZHU(ni+(nj-1)*(jdy1+jdy2),1)=((R0+2*A*(nj-1)*tan(alpha)))*Theta;
%         ZHU(ni+(nj-1)*(jdy1+jdy2),2)=(2*A*cos(alpha));
        
    end
end
for ni=1:jdy2       %40
    for nj=1:jdxp          %10
             
        coo(ni+nj*(jdy1)+(nj-1)*(jdy2),1)=A+2*A*(nj-1)-0.5*h*sin(alpha);       
        coo(ni+nj*(jdy1)+(nj-1)*(jdy2),2)=((R0+(A+2*A*(nj-1))*tan(alpha))+0.5*h*cos(alpha))*sin(Theta*(ni-1));   
        coo(ni+nj*(jdy1)+(nj-1)*(jdy2),3)=((R0+(A+2*A*(nj-1))*tan(alpha))+0.5*h*cos(alpha))*cos(Theta*(ni-1));
%         ZHU(ni+nj*(jdy1)+(nj-1)*(jdy2),1)=((R0+(A+2*A*(nj-1))*tan(alpha)))*Theta;
%         ZHU(ni+nj*(jdy1)+(nj-1)*(jdy2),2)=(2*A*cos(alpha));
    end
end
    
%圆柱壳的全局坐标   上表面
for ni=1:jdy1
    for nj=jdxp+2:jdx+1
%              2*a*(21-10-1)
        coo(ni+(nj-1)*(jdy1+jdy2),1)=2*a*(nj-jdxp-1)+2*A*(jdxp);       
        coo(ni+(nj-1)*(jdy1+jdy2),2)=(R+0.5*h)*sin(0.5*Theta*(ni-1));   
        coo(ni+(nj-1)*(jdy1+jdy2),3)=(R+0.5*h)*cos(0.5*Theta*(ni-1));
%         ZHU(ni+(nj-1)*(jdy1+jdy2),1)=(R)*Theta;
%         ZHU(ni+(nj-1)*(jdy1+jdy2),2)=2*a;
    end
end

for ni=1:jdy2
    for nj=jdxp+1:jdx
             
        coo(ni+nj*(jdy1)+(nj-1)*(jdy2),1)=a+2*a*(nj-jdxp-1)+2*A*(jdxp);       
        coo(ni+nj*(jdy1)+(nj-1)*(jdy2),2)=(R+0.5*h)*sin(Theta*(ni-1));   
        coo(ni+nj*(jdy1)+(nj-1)*(jdy2),3)=(R+0.5*h)*cos(Theta*(ni-1));
%         ZHU(ni+nj*(jdy1)+(nj-1)*(jdy2),1)=(R)*Theta;
%         ZHU(ni+nj*(jdy1)+(nj-1)*(jdy2),2)=2*a;
    end
end






for ni=1:jdy1
    for nj=1:jdxp+1
             
        co(ni+(nj-1)*(jdy1+jdy2),1)=2*A*(nj-1)+0.5*h*sin(alpha);       
        co(ni+(nj-1)*(jdy1+jdy2),2)=((R0+2*A*(nj-1)*tan(alpha))-0.5*h*cos(alpha))*sin(0.5*Theta*(ni-1));   
        co(ni+(nj-1)*(jdy1+jdy2),3)=((R0+2*A*(nj-1)*tan(alpha))-0.5*h*cos(alpha))*cos(0.5*Theta*(ni-1));


    end
end

for ni=1:jdy2
    for nj=1:jdxp
             
        co(ni+nj*(jdy1)+(nj-1)*(jdy2),1)=A+2*A*(nj-1)+0.5*h*sin(alpha);       
        co(ni+nj*(jdy1)+(nj-1)*(jdy2),2)=((R0+(A+2*A*(nj-1))*tan(alpha))-0.5*h*cos(alpha))*sin(Theta*(ni-1));   
        co(ni+nj*(jdy1)+(nj-1)*(jdy2),3)=((R0+(A+2*A*(nj-1))*tan(alpha))-0.5*h*cos(alpha))*cos(Theta*(ni-1));

    end
end
for ni=1:jdy1
    for nj=jdxp+2:jdx+1
             
        co(ni+(nj-1)*(jdy1+jdy2),1)=2*a*(nj-jdxp-1)+2*A*(jdxp);     
        co(ni+(nj-1)*(jdy1+jdy2),2)=(R-0.5*h)*sin(0.5*Theta*(ni-1));   
        co(ni+(nj-1)*(jdy1+jdy2),3)=(R-0.5*h)*cos(0.5*Theta*(ni-1));

    end
end

for ni=1:jdy2
    for nj=jdxp+1:jdx
             
        co(ni+nj*(jdy1)+(nj-1)*(jdy2),1)=a+2*a*(nj-jdxp-1)+2*A*(jdxp);    
        co(ni+nj*(jdy1)+(nj-1)*(jdy2),2)=(R-0.5*h)*sin(Theta*(ni-1));   
        co(ni+nj*(jdy1)+(nj-1)*(jdy2),3)=(R-0.5*h)*cos(Theta*(ni-1));

    end
end





for ni=1:jdy1
    for nj=1:jdxp+1
             
        coop(ni+(nj-1)*(jdy1+jdy2),1)=2*A*(nj-1)-(0.5*h+hp)*sin(alpha);       
        coop(ni+(nj-1)*(jdy1+jdy2),2)=((R0+2*A*(nj-1)*tan(alpha))+(0.5*h+hp)*cos(alpha))*sin(0.5*Theta*(ni-1));   
        coop(ni+(nj-1)*(jdy1+jdy2),3)=((R0+2*A*(nj-1)*tan(alpha))+(0.5*h+hp)*cos(alpha))*cos(0.5*Theta*(ni-1));

    end
end
for ni=1:jdy2
    for nj=1:jdxp
             
        coop(ni+nj*(jdy1)+(nj-1)*(jdy2),1)=A+2*A*(nj-1)-(0.5*h+hp)*sin(alpha);       
        coop(ni+nj*(jdy1)+(nj-1)*(jdy2),2)=((R0+(A+2*A*(nj-1))*tan(alpha))+(0.5*h+hp)*cos(alpha))*sin(Theta*(ni-1));   
        coop(ni+nj*(jdy1)+(nj-1)*(jdy2),3)=((R0+(A+2*A*(nj-1))*tan(alpha))+(0.5*h+hp)*cos(alpha))*cos(Theta*(ni-1));

    end
end
for ni=1:jdy1
    for nj=jdxp+2:jdx+1
             
        coop(ni+(nj-1)*(jdy1+jdy2),1)=2*a*(nj-jdxp-1)+2*A*(jdxp);       
        coop(ni+(nj-1)*(jdy1+jdy2),2)=(R+(0.5*h+hp))*sin(0.5*Theta*(ni-1));   
        coop(ni+(nj-1)*(jdy1+jdy2),3)=(R+(0.5*h+hp))*cos(0.5*Theta*(ni-1));

    end
end

for ni=1:jdy2
    for nj=jdxp+1:jdx
             
        coop(ni+nj*(jdy1)+(nj-1)*(jdy2),1)=a+2*a*(nj-jdxp-1)+2*A*(jdxp);       
        coop(ni+nj*(jdy1)+(nj-1)*(jdy2),2)=(R+(0.5*h+hp))*sin(Theta*(ni-1));   
        coop(ni+nj*(jdy1)+(nj-1)*(jdy2),3)=(R+(0.5*h+hp))*cos(Theta*(ni-1));

    end
end






for ni=1:jdy1
    for nj=1:jdxp+1
             
        cop(ni+(nj-1)*(jdy1+jdy2),1)=2*A*(nj-1)+(0.5*h+hp)*sin(alpha);       
        cop(ni+(nj-1)*(jdy1+jdy2),2)=((R0+2*A*(nj-1)*tan(alpha))-(0.5*h+hp)*cos(alpha))*sin(0.5*Theta*(ni-1));   
        cop(ni+(nj-1)*(jdy1+jdy2),3)=((R0+2*A*(nj-1)*tan(alpha))-(0.5*h+hp)*cos(alpha))*cos(0.5*Theta*(ni-1));

    end
end

for ni=1:jdy2
    for nj=1:jdxp
             
        cop(ni+nj*(jdy1)+(nj-1)*(jdy2),1)=A+2*A*(nj-1)+(0.5*h+hp)*sin(alpha);       
        cop(ni+nj*(jdy1)+(nj-1)*(jdy2),2)=((R0+(A+2*A*(nj-1))*tan(alpha))-(0.5*h+hp)*cos(alpha))*sin(Theta*(ni-1));   
        cop(ni+nj*(jdy1)+(nj-1)*(jdy2),3)=((R0+(A+2*A*(nj-1))*tan(alpha))-(0.5*h+hp)*cos(alpha))*cos(Theta*(ni-1));

    end
end
for ni=1:jdy1
    for nj=jdxp+2:jdx+1
             
        cop(ni+(nj-1)*(jdy1+jdy2),1)=2*a*(nj-jdxp-1)+2*A*(jdxp);     
        cop(ni+(nj-1)*(jdy1+jdy2),2)=(R-(0.5*h+hp))*sin(0.5*Theta*(ni-1));   
        cop(ni+(nj-1)*(jdy1+jdy2),3)=(R-(0.5*h+hp))*cos(0.5*Theta*(ni-1));

    end
end

for ni=1:jdy2
    for nj=jdxp+1:jdx
             
        cop(ni+nj*(jdy1)+(nj-1)*(jdy2),1)=a+2*a*(nj-jdxp-1)+2*A*(jdxp);    
        cop(ni+nj*(jdy1)+(nj-1)*(jdy2),2)=(R-(0.5*h+hp))*sin(Theta*(ni-1));   
        cop(ni+nj*(jdy1)+(nj-1)*(jdy2),3)=(R-(0.5*h+hp))*cos(Theta*(ni-1));

    end
end



cooo=(coo+co)/2;

COO=zeros(size(co,1),3);
for i=1:jdy*jdx     %1~800
    CC=ismember(i,nns);
    if   CC~=0
        for j=1:8
            COO(en(i,j),:)=1.02*coop(en(i,j),:);
         end
    end
end
    
x2=co(:,1);
y2=co(:,2);
z2=co(:,3);
x3=COO(:,1);
y3=COO(:,2);
z3=COO(:,3);



% 
for i=1:jdx+1
  for j=1:jdy  
X2(2*i-1,j)=x2((2*j-1)+3*jdy*(i-1));
Y2(2*i-1,j)=y2((2*j-1)+3*jdy*(i-1));
Z2(2*i-1,j)=z2((2*j-1)+3*jdy*(i-1));

X3(2*i-1,j)=x3((2*j-1)+3*jdy*(i-1));
Y3(2*i-1,j)=y3((2*j-1)+3*jdy*(i-1));
Z3(2*i-1,j)=z3((2*j-1)+3*jdy*(i-1));
  end
end
for i=1:jdx
    for j=1:jdy 
X2(2*i,j)=x2(2*jdy+j+3*jdy*(i-1));
Y2(2*i,j)=y2(2*jdy+j+3*jdy*(i-1));
Z2(2*i,j)=z2(2*jdy+j+3*jdy*(i-1));

X3(2*i,j)=x3(2*jdy+j+3*jdy*(i-1));
Y3(2*i,j)=y3(2*jdy+j+3*jdy*(i-1));
Z3(2*i,j)=z3(2*jdy+j+3*jdy*(i-1));
    end
end
X2(:,jdy+1)=X2(:,1);
Y2(:,jdy+1)=Y2(:,1);
Z2(:,jdy+1)=Z2(:,1);


figure(1)
surf(X2,Y2,Z2);
hold on
surf(X3,Y3,Z3);
hold on
% x1=coo(:,1);
% y1=coo(:,2);
% z1=coo(:,3);
% x2=co(:,1);
% y2=co(:,2);
% z2=co(:,3);
% figure(1)
% plot3(x1,y1,z1);
% hold on
% plot3(x2,y2,z2);
disp(1:(jdy1*(jdx+1)+jdy*jdx),1:5)=1;     % node displacement
disp((jdy1*(jdx+1)+jdy*jdx)-jdy1+1:(jdy1*(jdx+1)+jdy*jdx),1:5)=0;
dof=0;                   %degree of freedom
k11=zeros(6400);
m11=zeros(6400);

for ni=1:1280 %2480
    for nj=1:5
        if disp(ni,nj)~=0
            dof=dof+1;
            disp(ni,nj)=dof;
        end
    end
end
jdzb=coo;
jdzb1=co;

dybh=en;
index(1:40)=0; % vector sontaining system dofs of nodes in each element. 

for loopi=1:jdy*jdx  %单元编号

    dyhm=loopi;
    [ek,theta,xv2i,xv1i,xv3i,zmtemp,v3i,D,jtemp]=shellek(E,Niu,t,dyhm,jdzb,jdzb1,dybh);   
    [em]=shellem(zmtemp,v3i,density,theta,t,xv2i,xv1i,No_INTpoint_x,No_INTpoint_y,No_INTpoint_z,jtemp);
 
 end

    for zi=1:8
        index((zi-1)*5+1)=disp(en(loopi,zi),1);
        index((zi-1)*5+2)=disp(en(loopi,zi),2);
        index((zi-1)*5+3)=disp(en(loopi,zi),3);
        index((zi-1)*5+4)=disp(en(loopi,zi),4);
        index((zi-1)*5+5)=disp(en(loopi,zi),5);

    end
    for jx=1:40
        for jy=1:40
            if(index(jx)*index(jy)~=0)
                  k11(index(jx),index(jy))=k11(index(jx),index(jy))+ek(jx,jy);
                  m11(index(jx),index(jy))=m11(index(jx),index(jy))+em(jx,jy);
                 
      end
   end
end


[v11,d11]=eigs(k11(1:1200*5,1:1200*5),m11(1:1200*5,1:1200*5));
tempd=diag(d11);
[d11,sortindex]=sort(tempd);
v11=v11(:,sortindex);
v1=v11(1:600*5,1:600*5);
% v1=v1(1:1500,1:1500);

P=zeros(1:680*5,1:680*5);
for i=1:600*5
    for j=1:600*5
        P(i,j)=v1(i,j);
    end
end

for i=600*5:600*5+80*5
    for j=600*5:600*5+80*5
        P(i,j)=eye(80);
    end
end

for i=1:600*5
    for j=600*5:(600+80)*5
        P(i,j)=-k11(i,i)'*k11(i,j);
    end
end

for i=600*5:(600+80)*5
    for j=1:600*5
        P(i,j)=[0];
    end
end

%
%第二个子结构
dof=0;                   %degree of freedom
k22=zeros(1:1200*5,1:1200*5);
m22=zeros(1:1200*5,1:1200*5);

for ni=1200:2400 %2480
    for nj=1:5
        if disp(ni,nj)~=0
            dof=dof+1;
            disp(ni,nj)=dof;
        end
    end
end
jdzb=coo;
jdzb1=co;

dybh=en;
index(1:40)=0; % vector sontaining system dofs of nodes in each element. 

for loopi=1:jdy*jdx  %单元编号

    dyhm=loopi;
    [ek,theta,xv2i,xv1i,xv3i,zmtemp,v3i,D,jtemp]=shellek(E,Niu,t,dyhm,jdzb,jdzb1,dybh);   
    [em]=shellem(zmtemp,v3i,density,theta,t,xv2i,xv1i,No_INTpoint_x,No_INTpoint_y,No_INTpoint_z,dyhm,dybh);
 
 end

    for zi=1:8
        index((zi-1)*5+1)=disp(en(loopi,zi),1);
        index((zi-1)*5+2)=disp(en(loopi,zi),2);
        index((zi-1)*5+3)=disp(en(loopi,zi),3);
        index((zi-1)*5+4)=disp(en(loopi,zi),4);
        index((zi-1)*5+5)=disp(en(loopi,zi),5);

    end
    for jx=1:40
        for jy=1:40
            if(index(jx)*index(jy)~=0)
                  k22(index(jx),index(jy))=k22(index(jx),index(jy))+ek(jx,jy);
                  m22(index(jx),index(jy))=m22(index(jx),index(jy))+em(jx,jy);
                 
      end
   end
end

%第二个子结构的第一次坐标变换矩阵
[v22,d22]=eigs(k22(81*5:1200*5,81*5:1200*5),m22(81*5:1200*5,81*5:1200*5));
tempd=diag(d22);
[d22,sortindex]=sort(tempd);
v22=v22(:,sortindex);
v2=v22(1:600*5,1:600*5);
% v1=v1(1:1500,1:1500);

P1=zeros(1:680*5,1:680*5); %变换矩阵 2-1
for i=1:600*5
    for j=1:600*5
        P1(i,j)=v2(80+i,80+j);
    end
end

for i=600*5:680*5
    for j=600*5:680*5
        P1(i,j)=eye(80);
    end
end

for i=1:600*5
    for j=600*5:680*5
        P1(i,j)=-k22(1200*5+1-i,j)'*k22(1200*5+1-i,j);
    end
end

for i=600*5:680*5
    for j=1:600*5
        P1(i,j)=zeros;
    end
end

nm=10;
[v,d] =eigs(k,m,nm,'SM');  %v特征向量 d特征值
tempd=diag(d);
[d,sortindex]=sort(tempd);
omega=sqrt(d);
v=v(:,sortindex);
mode_number=1:nm;
frequency(mode_number)=sqrt(d(mode_number))/(2*pi);
Frequency=real(frequency);
V=real(v);

xvec=v;
mn=3;
voo=cooo; %cooo中面坐标
doo=zeros(size(cooo,1),1);
xx=xvec(1:5:size(k)-4,mn);
yy=xvec(2:5:size(k)-3,mn);
zz=xvec(3:5:size(k)-2,mn);

for ni=1:jdy1
    for nj=1:jdxp+1
             
     Theta1(ni+(nj-1)*(jdy1+jdy2),1)=0.5*Theta*(ni-1);       

    end
end
for ni=1:jdy2
    for nj=1:jdxp
             
        Theta1(ni+nj*(jdy1)+(nj-1)*(jdy2),1)=Theta*(ni-1);       
     
    end
end
for ni=1:jdy1
    for nj=jdxp+2:jdx+1          
        Theta1(ni+(nj-1)*(jdy1+jdy2),1)=0.5*Theta*(ni-1); 

    end
end

for ni=1:jdy2
    for nj=jdxp+1:jdx
             
        Theta1(ni+nj*(jdy1)+(nj-1)*(jdy2),1)=Theta*(ni-1);       

    end
end
    
% 为划分颜色
for i=1:0.5*size(xx)  
    dis(i,1)=xx(i)*sin(alpha)+(zz(i)*cos(Theta1(i))+yy(i)*sin(Theta1(i)))*cos(alpha);        
end    

for i=0.5*size(xx)+1:size(xx)  
    dis(i,1)=zz(i)*cos(Theta1(i))+yy(i)*sin(Theta1(i));        
end 

% dis(2400,1)=0

CA=0.002;
voo(1:size(voo)-jdy1,1)=voo(1:size(voo)-jdy1,1)+CA*xx;
voo(1:size(voo)-jdy1,2)=voo(1:size(voo)-jdy1,2)+CA*yy;
voo(1:size(voo)-jdy1,3)=voo(1:size(voo)-jdy1,3)+CA*zz;
doo(1:size(doo)-jdy1,1)=doo(1:size(doo)-jdy1,1)+dis;
x1=voo(:,1);
y1=voo(:,2);
z1=voo(:,3);

for i=1:jdx+1
    for j=1:jdy  
     X1(2*i-1,j)=x1((2*j-1)+3*jdy*(i-1));
     Y1(2*i-1,j)=y1((2*j-1)+3*jdy*(i-1));
     Z1(2*i-1,j)=z1((2*j-1)+3*jdy*(i-1));
     DIS(2*i-1,j)=doo((2*j-1)+3*jdy*(i-1));
    end
end
for i=1:jdx
    for j=1:jdy 
X1(2*i,j)=x1(2*jdy+j+3*jdy*(i-1));
Y1(2*i,j)=y1(2*jdy+j+3*jdy*(i-1));
Z1(2*i,j)=z1(2*jdy+j+3*jdy*(i-1));
DIS(2*i,j)=doo(2*jdy+j+3*jdy*(i-1));
    end
end
X1(:,jdy+1)=X1(:,1);
Y1(:,jdy+1)=Y1(:,1);
Z1(:,jdy+1)=Z1(:,1);
DIS(:,jdy+1)=DIS(:,1);

figure(11)
surf(X1,Y1,Z1,abs(DIS));
% surf(X1,Y1,Z1);
daspect([1,1,1])  %控制为笛卡尔坐标系 各坐标轴比例一致

