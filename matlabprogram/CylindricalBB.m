function [KB,MB,T2]=CylindricalBB(nn)
% nn=40;
No_INTpoint_x=4;
No_INTpoint_y=4;
No_INTpoint_z=4;
a=0.01;
A=0.0025;   %圆锥壳半个单元长度
t=0.001;
R=0.07;
Theta=9*pi/180;
% alpha=30*pi/180;
density=2700;
Niu=0.33;
E=70e9;
h=t;

jdx=10;               %number of nodes in x direction
% jdx1=2*jdx+1;
jdy=40;               %number of nodes in y direction
jdy1=2*jdy;     
jdy2=jdy;
    k(1:5*(jdy1*(jdx+1)+jdy*jdx)-5*jdy1,1:5*(jdy1*(jdx+1)+jdy*jdx)-5*jdy1)=0;     %system stiffness matrix
    m(1:5*(jdy1*(jdx+1)+jdy*jdx)-5*jdy1,1:5*(jdy1*(jdx+1)+jdy*jdx)-5*jdy1)=0;    %system mass matrix

%     k(1:5*(jdy1*(jdx+1)+jdy*jdx),1:5*(jdy1*(jdx+1)+jdy*jdx))=0;     %system stiffness matrix
%     m(1:5*(jdy1*(jdx+1)+jdy*jdx),1:5*(jdy1*(jdx+1)+jdy*jdx))=0;    %system mass matrix

    
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
        
        elseif ni==jdy
        en(ni+(nj-1)*(jdy),4)=(2*ni-1)+(nj-1)*jdy1+(nj-1)*jdy2;
        en(ni+(nj-1)*(jdy),7)=2*ni+(nj-1)*jdy1+(nj-1)*jdy2;
        en(ni+(nj-1)*(jdy),3)=1+(nj-1)*jdy1+(nj-1)*jdy2;
        en(ni+(nj-1)*(jdy),1)=(2*ni-1)+(nj)*jdy1+(nj)*jdy2;
        en(ni+(nj-1)*(jdy),5)=2*ni+(nj)*jdy1+(nj)*jdy2;
        en(ni+(nj-1)*(jdy),2)=1+(nj)*jdy1+(nj)*jdy2;
        
      end
   end
end

for ni=1:jdy
    for nj=1:jdx
        if ni~=jdy
        en(ni+(nj-1)*(jdy),8)=ni+(nj)*jdy1+(nj-1)*jdy2;
        en(ni+(nj-1)*(jdy),6)=ni+1+(nj)*jdy1+(nj-1)*jdy2;

        
        elseif ni==jdy
        en(ni+(nj-1)*(jdy),8)=ni+(nj)*jdy1+(nj-1)*jdy2;
        en(ni+(nj-1)*(jdy),6)=1+(nj)*jdy1+(nj-1)*jdy2;
        
        end
    end
end

%圆柱壳上表面节点全局坐标
for ni=1:jdy1
    for nj=1:jdx+1 

        coo(ni+(nj-1)*(jdy1+jdy2),1)=2*a*(nj-1)+2*A*(jdx);       
        coo(ni+(nj-1)*(jdy1+jdy2),2)=(R+0.5*h)*sin(0.5*Theta*(ni-1));   
        coo(ni+(nj-1)*(jdy1+jdy2),3)=(R+0.5*h)*cos(0.5*Theta*(ni-1));

    end
end

for ni=1:jdy2
    for nj=1:jdx 
             
        coo(ni+nj*(jdy1)+(nj-1)*(jdy2),1)=a+2*a*(nj-1)+2*A*(jdx);       
        coo(ni+nj*(jdy1)+(nj-1)*(jdy2),2)=(R+0.5*h)*sin(Theta*(ni-1));   
        coo(ni+nj*(jdy1)+(nj-1)*(jdy2),3)=(R+0.5*h)*cos(Theta*(ni-1));

    end
end
% 圆锥壳的全局坐标下表面
for ni=1:jdy1
    for nj=1:jdx+1 
             
        co(ni+(nj-1)*(jdy1+jdy2),1)=2*a*(nj-1)+2*A*(jdx);     
        co(ni+(nj-1)*(jdy1+jdy2),2)=(R-0.5*h)*sin(0.5*Theta*(ni-1));   
        co(ni+(nj-1)*(jdy1+jdy2),3)=(R-0.5*h)*cos(0.5*Theta*(ni-1));

    end
end

for ni=1:jdy2
    for nj=1:jdx 
             
        co(ni+nj*(jdy1)+(nj-1)*(jdy2),1)=a+2*a*(nj-1)+2*A*(jdx);    
        co(ni+nj*(jdy1)+(nj-1)*(jdy2),2)=(R-0.5*h)*sin(Theta*(ni-1));   
        co(ni+nj*(jdy1)+(nj-1)*(jdy2),3)=(R-0.5*h)*cos(Theta*(ni-1));

    end
end

cooo=(coo+co)/2;
x3=coo(:,1);
y3=coo(:,2);
z3=coo(:,3);
x4=co(:,1);
y4=co(:,2);
z4=co(:,3);
figure(1)
plot3(x3,y3,z3);
hold on
plot3(x4,y4,z4);
hold on

disp(1:(jdy1*(jdx+1)+jdy*jdx),1:5)=1;     % node displacement
% constraints=1:jdy1;  % constraints
% disp(constraints,:)=0;

% disp(1:jdy1,1:3)=0;
disp((jdy1*(jdx+1)+jdy*jdx)-jdy1+1:(jdy1*(jdx+1)+jdy*jdx),1:5)=0;
dof=0;                   %degree of freedom

for ni=1:(jdy1*(jdx+1)+jdy*jdx-jdy1)
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
for loopi=1:jdy*jdx
    dyhm=loopi;
     [ek,theta,xv2i,xv1i,xv3i,zmtemp,v3i,D,jtemp]=shellek(E,Niu,t,dyhm,jdzb,jdzb1,dybh);   
     em=shellem(zmtemp,v3i,density,theta,t,xv2i,xv1i,No_INTpoint_x,No_INTpoint_y,No_INTpoint_z);
 
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
                  k(index(jx),index(jy))=k(index(jx),index(jy))+ek(jx,jy);
                  m(index(jx),index(jy))=m(index(jx),index(jy))+em(jx,jy);
            end
        end
    end
end

% AA=det(k);
% BB=det(m);

% [v,d] = eig(k,m);
% tempd=diag(d);
% [nd,sortindex]=sort(tempd);
% v=v(:,sortindex);
% mode_number=1:15;
% frequency(mode_number)=sqrt(nd(mode_number))/(2*pi);

% 
% up_mtrxA=k;up_massA=m;
% nmodes=25;F=zeros(size(k,1),1);F(1,1)=1;
% [xvec,freq]= Lanczos(up_mtrxA,up_massA,F,nmodes);
% [freq1,sortindex]=sort(freq);
% xvec=xvec(:,sortindex);

% w=884;
mii=m(jdy1*5+1:size(k),jdy1*5+1:size(k));
kii=k(jdy1*5+1:size(k),jdy1*5+1:size(k));
kij=k(jdy1*5+1:size(k),1:jdy1*5);
% mij=m(jdy1*5+1:size(k),1:jdy1*5);
I=eye(5*jdy1);
% [xvec,d] =eigs(k,m,nn,'SM'); 
[xvec,d] =eigs(kii,mii,nn,'SM'); 
tempd=diag(d);
[d,sortindex]=sort(tempd);
% mode_number=1:nn;
% frequency(mode_number)=sqrt(d(mode_number))/(2*pi);
% Frequency=real(frequency);
xvec=xvec(:,sortindex);
v=xvec(:,1:nn);

% for i=1:size(kii)
%     for j=1:nn
%         v(i:j)=v(i,j)/sqrt(mii(i,j));
%     end
% end

T2=zeros(size(k,1),5*jdy1+nn);
% % T2(1:1:5*jdy1,1:5*jdy1)=I;
% % T2(5*jdy1+1:size(k,1),1:5*jdy1)=-inv(kii)*kij;
% % T2(5*jdy1+1:size(k,1),5*jdy1+1:5*jdy1+nn)=v;
% T2(1:size(kii,1),1:nn)=v;
% T2(1:size(kii,1),nn+1:5*jdy1+nn)=-inv(kii)*kij;
% T2(size(kii,1)+1:size(k,1),nn+1:5*jdy1+nn)=I;

T2(1:jdy1*5,1:jdy1*5)=I;
T2(jdy1*5+1:size(k,1),1:jdy1*5)=-inv(kii)*kij;
T2(jdy1*5+1:size(k,1),jdy1*5+1:jdy1*5+nn)=v;

KB=T2'*k*T2;
MB=T2'*m*T2;