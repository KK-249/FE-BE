function [KA,MA,T1]=SubA(coool,enl,nn,con_z)
%高斯积分的点数
No_INTpoint_x = 2;        
No_INTpoint_y = 2;
No_INTpoint_z = 2;

%密度 泊松比 杨氏模量
density=7900;
Niu=0.3;
E=196e9;
t=0.05;

k(1:size(coool,1)*5,1:size(coool,1)*5)=0;        
m(1:size(coool,1)*5,1:size(coool,1)*5)=0; 

indices1 = find(coool(:,3)==con_z);         % 子结构固定界面节点数目
indices2 = find(coool(:,3)< con_z);

dof = 0;
disp(1:size(coool,1),1:5) = 1;  

for ni=1:size(coool,1) 
    for nj=1:5
        if disp(ni,nj)~=0
            dof=dof+1;
            disp(ni,nj)=dof;
        end
    end
end

[jdzb1, jdzb, nodeNormals] = genShellSurfacesFromConn(coool, enl, t);


dybh=enl;
index(1:40)=0;

for loopi=1:size(enl,1)
    dyhm=loopi;
     [ek,theta,xv2i,xv1i,~,zmtemp,v3i,~,jtemp,~]=shellek(E,Niu,t,dyhm,jdzb,jdzb1,dybh);   
     [em1,~]=shellem(zmtemp,v3i,density,theta,t,xv2i,xv1i,No_INTpoint_x,No_INTpoint_y,No_INTpoint_z,jtemp);
%     [em1]=shellem_corrected(density,new_nodes,dybh,total_area1,total_area2,total_area3);

    for zi=1:8
        index((zi-1)*5+1)=disp(enl(loopi,zi),1);
        index((zi-1)*5+2)=disp(enl(loopi,zi),2);
        index((zi-1)*5+3)=disp(enl(loopi,zi),3);
        index((zi-1)*5+4)=disp(enl(loopi,zi),4);
        index((zi-1)*5+5)=disp(enl(loopi,zi),5);
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

% elemIdx = find( any(enl == 1, 2) );

iidx = size(indices2,1)*5;
condx = size(indices1,1)*5;
mii=m(1:iidx,1:iidx);
kii=k(1:iidx,1:iidx);
kij=k(1:iidx,iidx+1:size(k));

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


T1=zeros(size(k,1),condx+nn);
T1(1:size(kii,1),1:nn)=v;
T1(1:size(kii,1),nn+1:condx+nn)=-inv(kii)*kij;
T1(size(kii,1)+1:size(k,1),nn+1:condx+nn)=I;


% data1 = load('matrixA.mat');  
% data2 = load('matrixB.mat');  
% A1    = data1.subk1;               % 重命名为 A1
% A2    = data2.subm1; 
% 
% % 现在就可以比较 A1, A2
% tf = isequal(A1, k);
% fprintf('两矩阵完全相同吗？ %d\n', tf);
% 
% tf = isequal(A2, m);
% fprintf('两矩阵完全相同吗？ %d\n', tf);


KA=T1'*k*T1;
MA=T1'*m*T1;