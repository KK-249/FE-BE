function [K1, M1, T1] = substruct1(cooo, en, indices1, indices4)
clear
clc

%% 获得子结构I的单元编号、节点编号
z_max = 5;
selected_elements = [];

for i = 1:size(en,1)
    z_coords = cooo(en(i,:), 3);
    if all(z_coords >= z_max)
        selected_elements(end+1) = i;
    end
end

en1 = en(selected_elements,:);
cooo1(indices1,:) = cooo(indices1,:);
%% 有限元获得子结构k、m
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
k1(1:size(cooo1,1)*5-size(indices4,1)*5,1:size(cooo1,1)*5-size(indices4,1)*5)=0;        
m1(1:size(cooo1,1)*5-size(indices4,1)*5,1:size(cooo1,1)*5-size(indices4,1)*5)=0;        

dof = 0;
disp(1:size(cooo1,1),1:5) = 1;  
disp(indices4,:) = 0;

for ni=1:size(cooo1,1) 
    for nj=1:5
        if disp(ni,nj)~=0
            dof=dof+1;
            disp(ni,nj)=dof;
        end
    end
end

[jdzb1, jdzb] = genShellSurfacesFromConn(cooo1, en1, t);

dybh=en1;
index(1:40)=0;

for loopi=1:size(en1,1)
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
                  k1(index(jx),index(jy))=k1(index(jx),index(jy))+ek(jx,jy);
                  m1(index(jx),index(jy))=m1(index(jx),index(jy))+em1(jx,jy);
            end
        end
    end
end