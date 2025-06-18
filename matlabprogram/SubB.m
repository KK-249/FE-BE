function [KB,MB,T2]=SubB(cooor,enr,nn,con_z)
%高斯积分的点数
No_INTpoint_x = 2;        
No_INTpoint_y = 2;
No_INTpoint_z = 2;

%密度 泊松比 杨氏模量
density=7900;
Niu=0.3;
E=196e9;
t=0.05;

k(1:size(cooor,1)*5,1:size(cooor,1)*5)=0;        
m(1:size(cooor,1)*5,1:size(cooor,1)*5)=0; 

indices1 = find(cooor(:,3)==con_z);     % 子结构固定界面节点数目
% indices2 = find(cooor(:,3)> 5);

dof = 0;
disp(1:size(cooor,1),1:5) = 1;  

for ni=1:size(cooor,1) 
    for nj=1:5
        if disp(ni,nj)~=0
            dof=dof+1;
            disp(ni,nj)=dof;
        end
    end
end

[jdzb1, jdzb, nodeNormals] = genShellSurfacesFromConn(cooor, enr, t);


dybh=enr;
index(1:40)=0;

for loopi=1:size(enr,1)
    dyhm=loopi;
     [ek,theta,xv2i,xv1i,~,zmtemp,v3i,~,jtemp,~]=shellek(E,Niu,t,dyhm,jdzb,jdzb1,dybh);   
     [em1,~]=shellem(zmtemp,v3i,density,theta,t,xv2i,xv1i,No_INTpoint_x,No_INTpoint_y,No_INTpoint_z,jtemp);
%     [em1]=shellem_corrected(density,new_nodes,dybh,total_area1,total_area2,total_area3);
    if dyhm == 84
        aaa=1;
    end
    for zi=1:8
        index((zi-1)*5+1)=disp(enr(loopi,zi),1);
        index((zi-1)*5+2)=disp(enr(loopi,zi),2);
        index((zi-1)*5+3)=disp(enr(loopi,zi),3);
        index((zi-1)*5+4)=disp(enr(loopi,zi),4);
        index((zi-1)*5+5)=disp(enr(loopi,zi),5);
    end
    
     for jx=1:40
        for jy=1:40
            if(index(jx)*index(jy)~=0)
                  k(index(jx),index(jy))=k(index(jx),index(jy))+ek(jx,jy);
                  m(index(jx),index(jy))=m(index(jx),index(jy))+em1(jx,jy);
            end
        end
     end
end

% iidx = size(indices2,1)*5;
condx = size(indices1,1)*5;
kdx = size(k,1);
mii=m(condx+1:kdx,condx+1:kdx);
kii=k(condx+1:kdx,condx+1:kdx);
kij=k(condx+1:kdx,1:condx);

I=eye(condx);
[xvec,d] =eigs(kii,mii,nn,'SM'); 
tempd=diag(d);
[d,sortindex]=sort(tempd);
xvec=xvec(:,sortindex);
v=xvec(:,1:nn);
mode_number=1:20;
frequency(mode_number)=sqrt(d(mode_number))/(2*pi);
Frequency=real(frequency);
Frequency=Frequency';

T2=zeros(size(k,1),condx+nn);
T2(1:condx,1:condx) = I;
T2(condx+1:kdx,1:condx)=-inv(kii)*kij;
T2(condx+1:kdx,condx+1:condx+nn)=v;


KB=T2'*k*T2;
MB=T2'*m*T2;