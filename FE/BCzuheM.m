%% 该程序利用动态子结构法计算组合壳的固有频率和声响应
clc
clear

%% 导入网格信息 no 1 书上经典例子 2 组合壳带球尾模型
no = 2;
[msh] = qcl(no); 

nodes = msh.POS(:,:)./1000;
elements = msh.QUADS(:,1:4);

% cooo 中面坐标 en 8单元节点编号 
[cooo, en] = generateEightNodeElements(nodes, elements);

% 根据纵坐标划分子结构
% indices0 = find(cooo(:,3) > -10 & cooo(:,3) < 6.6 );
indices1 = find(cooo(:,3) > 5 & cooo(:,3) < 6.6);
indices2 = find(cooo(:,3) < 5 );

% 固定界面的坐标索引
indices3 = find(cooo(:,3) == 5);

% 整体结构的约束
indices4 = find(cooo(:,3) >= 6.6);
%% 有限元获得K、M
%高斯积分的点数
No_INTpoint_x = 2;        
No_INTpoint_y = 2;
No_INTpoint_z = 2;

%密度 泊松比 杨氏模量
density=7900;
Niu=0.3;
E=196e9;
t=0.05;

%system stiffness matrix 、 system mass matrix
k(1:(size(cooo,1)-size(indices4,1))*5,1:(size(cooo,1)-size(indices4,1))*5)=0;        
m(1:(size(cooo,1)-size(indices4,1))*5,1:(size(cooo,1)-size(indices4,1))*5)=0;        

dof = 0;
disp(1:size(cooo,1),1:5) = 1;  
disp(indices4,:) = 0;

for ni=1:size(cooo,1) 
    for nj=1:5
        if disp(ni,nj)~=0
            dof=dof+1;
            disp(ni,nj)=dof;
        end
    end
end

[jdzb1, jdzb, nodeNormals] = genShellSurfacesFromConn(cooo, en, t);

dybh=en;
index(1:40)=0;

for loopi=1:size(en,1)
    dyhm=loopi;
     [ek,theta,xv2i,xv1i,xv3i,zmtemp,v3i,D,jtemp,xv3ii]=shellek(E,Niu,t,dyhm,jdzb,jdzb1,dybh);   
     [em1,tt]=shellem(zmtemp,v3i,density,theta,t,xv2i,xv1i,No_INTpoint_x,No_INTpoint_y,No_INTpoint_z,jtemp);
%     [em1]=shellem_corrected(density,new_nodes,dybh,total_area1,total_area2,total_area3);

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
                  m(index(jx),index(jy))=m(index(jx),index(jy))+em1(jx,jy);
            end
        end
    end
end
%% BC 得各结构的k、m、T
nFree = 5;  % 每个节点的自由度
nn = 40;
% 将节点编号转换为自由度索引函数
node2dof =@(indices0) reshape((1:nFree)' + nFree*(indices0(:)' - 1), 1, []);

% 子结构 I
oldDOF1 = node2dof(indices1);
% newDOF1 = [];
newDOF1 = node2dof(indices3);
allDOF1 = [oldDOF1, newDOF1];

subk1 = k(allDOF1, allDOF1);
subm1 = m(allDOF1, allDOF1);
mii1=subm1(1:size(indices1,1)*5,1:size(indices1,1)*5);
kii1=subk1(1:size(indices1,1)*5,1:size(indices1,1)*5);
kij1=subk1(1:size(indices1,1)*5,size(indices1,1)*5+1:(size(indices1,1)+size(indices3,1))*5);
I1=eye(5*size(indices3,1));
[xvec1,d1] =eigs(kii1,mii1,nn,'SM'); 
tempd1=diag(d1);
[d1,sortindex]=sort(tempd1);
xvec1=xvec1(:,sortindex);
v1=xvec1(:,1:nn);
mode_number1=1:nn;
frequency1(mode_number1)=sqrt(d1(mode_number1))/(2*pi);
Frequency1=real(frequency1);
Frequency1=Frequency1';
V1=real(v1);

T1=zeros(size(subk1,1),5*size(indices3,1)+nn);
T1(1:size(kii1,1),1:nn)=v1;
T1(1:size(kii1,1),nn+1:5*size(indices3,1)+nn)=-inv(kii1)*kij1;
T1(size(kii1,1)+1:size(subk1,1),nn+1:5*size(indices3,1)+nn)=I1;
K1=T1'*subk1*T1;
M1=T1'*subm1*T1;

% 子结构II
oldDOF2 = node2dof(indices2);
newDOF2 = node2dof(indices3);
allDOF2 = [newDOF2, oldDOF2];

subk2 = k(allDOF2, allDOF2);
subm2 = m(allDOF2, allDOF2);

mii2=subm2(size(indices3,1)*5+1:size(subm2,1),size(indices3,1)*5+1:size(subm2,1));
kii2=subk2(size(indices3,1)*5+1:size(subm2,1),size(indices3,1)*5+1:size(subm2,1));
kij2=subk2(size(indices3,1)*5+1:size(subm2,1),1:size(indices3,1)*5);
I2=eye(5*size(indices3,1));

eps_k = 1e-5;         % 加扰动强度
eps_m = 1e-5;  
kii2_reg = kii2 + eps_k * speye(size(kii2));
mii2_reg = mii2 + eps_m * speye(size(mii2));
   
[xvec2, d2] = eigs(kii2_reg, mii2_reg, nn,'smallestreal');

tempd2=diag(d2);
[d2,sortindex2]=sort(tempd2);
xvec2=xvec2(:,sortindex2);
v2=xvec2(:,1:nn);
mode_number2=1:nn;
frequency2(mode_number2)=sqrt(d2(mode_number2))/(2*pi);
Frequency2=real(frequency2);
Frequency2=Frequency2';
V2=real(v2);

T2=zeros(size(subk2,1),5*size(indices3,1)+nn);
T2(1:size(indices3,1)*5,1:size(indices3,1)*5)=I2;
T2(size(indices3,1)*5+1:size(subk2,1),1:size(indices3,1)*5)=-inv(kii2_reg)*kij2;
T2(size(indices3,1)*5+1:size(subk2,1),size(indices3,1)*5+1:size(indices3,1)*5+nn)=v2;

K2=T2'*subk2*T2;
M2=T2'*subm2*T2;

%% 第二次坐标变换
KK=zeros((size(indices3,1)*5+nn)*2,(size(indices3,1)*5+nn)*2);
MM=zeros((size(indices3,1)*5+nn)*2,(size(indices3,1)*5+nn)*2);
KK(1:size(indices3,1)*5+nn,1:1:size(indices3,1)*5+nn)=K1;KK(size(indices3,1)*5+1+nn:(size(indices3,1)*5+nn)*2,size(indices3,1)*5+1+nn:(size(indices3,1)*5+nn)*2)=K2;
MM(1:size(indices3,1)*5+nn,1:1:size(indices3,1)*5+nn)=M1;MM(size(indices3,1)*5+1+nn:(size(indices3,1)*5+nn)*2,size(indices3,1)*5+1+nn:(size(indices3,1)*5+nn)*2)=M2;
T3=zeros((size(indices3,1)*5+nn)*2,size(indices3,1)*5+nn*2);
T3(1:nn,1:nn)=eye(nn);
T3(nn+1:size(indices3,1)*5+nn,2*nn+1:size(indices3,1)*5+2*nn)=eye(size(indices3,1)*5);
T3(nn+size(indices3,1)*5+1:size(indices3,1)*5*2+nn,2*nn+1:size(indices3,1)*5+2*nn)=eye(size(indices3,1)*5);
T3(size(indices3,1)*5*2+nn+1:size(indices3,1)*5*2+2*nn,nn+1:2*nn)=eye(nn);
K=T3'*KK*T3;
M=T3'*MM*T3;

nm=20;
[v,d] =eigs(K,M,nm,'SM');  %v特征向量 d特征值
tempd=diag(d);
[d,sortindex]=sort(tempd);
omega=sqrt(d);
v=v(:,sortindex);
mode_number=1:nm;
frequency(mode_number)=sqrt(d(mode_number))/(2*pi);
Frequency=real(frequency);
Frequency=Frequency';
V=real(v);