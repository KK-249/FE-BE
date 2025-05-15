clear
clc

%% 前处理 导入网格获得八结点单元
no=2;
[msh] = qcl(no);
nodes = msh.POS(:,:);   
elements = msh.QUADS(:,1:4);

% 调用生成八结点单元的函数
[new_nodes,eight_node_elements]=generateEightNodeElements(nodes, elements);
%% FEM部分
No_INTpoint_x=2;        %高斯积分的点数
No_INTpoint_y=2;
No_INTpoint_z=2;

density=7900;           %常数设置
Niu=0.3;
E=196e9;
t=0.05;

en=eight_node_elements;       %element node  
cooo=new_nodes/1000;

new_nodes=cooo;
disp(1:size(cooo,1),1:5)=1;     % node displacement

%找到约束面节点编号
indices = find(cooo(:,3) == 0);
rowsToZero = indices;
disp(rowsToZero, :) = 0;

%激励自由度
indices1 = find(cooo(:,3) == 6.6);

k(1:size(new_nodes,1)*5-size(indices,1)*5,1:size(new_nodes,1)*5-size(indices,1)*5)=0;        %system stiffness matrix
m(1:size(new_nodes,1)*5-size(indices,1)*5,1:size(new_nodes,1)*5-size(indices,1)*5)=0;        %system mass matrix

dof=0;  %degree of freedom

for ni=1:size(cooo,1) %节点总数-
    for nj=1:5
        if disp(ni,nj)~=0
            dof=dof+1;
            disp(ni,nj)=dof;
        end
    end
end

[jdzb1, jdzb, nodeNormals] = genShellSurfacesFromConn(cooo, eight_node_elements, t);

dybh=en;
index(1:40)=0; % vector sontaining system dofs of nodes in each element. 
for loopi=1:size(en,1)
    dyhm=loopi;
     [ek,theta,xv2i,xv1i,xv3i,zmtemp,v3i,D,jtemp,xv3ii]=shellek(E,Niu,t,dyhm,jdzb,jdzb1,dybh);   
     [em1,tt]=shellem(zmtemp,v3i,density,theta,t,xv2i,xv1i,No_INTpoint_x,No_INTpoint_y,No_INTpoint_z,jtemp);
    
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

TT = zeros(3*size(nodeNormals,1),size(nodeNormals,1));   %节点法向量

for i = 1:size(nodeNormals,1)
    TT(3*i-2,i) = nodeNormals(i,1);
    TT(3*i-1,i) = nodeNormals(i,2);
    TT(3*i-0,i) = nodeNormals(i,3);
end
%% 扫频计算
freq0 = 1;
step = 1;
freq = 1000;

for fi = freq0:step:freq
    
    [Ma] = add_mass(fi,new_nodes,eight_node_elements);
    Ma1 = TT*Ma*TT';
    [NewMat] =  expandDiagonalMatrix(Ma1,nodeNormals);

    k0 = k;
    m0 = m;
    for i=1:size(m,1)
        for j=1:size(m,1)
            
        m0(i,j)=m0(i,j)+NewMat(i,j);
        
        end
    end
    
    W0 = 2*pi*fi;
    F = zeros(size(k,1),1);
    F(43,1) = 100; 

    X = inv((k-W0^2*m))*F;
    
    X0 = zeros(size(cooo,1)*5,1);
    X0(:,1) = 1; 
    X0(indices*5-4,1) = 0; 
    X0(indices*5-3,1) = 0; 
    X0(indices*5-2,1) = 0; 
    X0(indices*5-1,1) = 0; 
    X0(indices*5-0,1) = 0; 

    xj = 1;
    for xi = 1:5*size(cooo,1)
       if X0(xi,1) ~= 0   
           
           X0 (xi,1) = X(xj,1);
           xj = xj+1;
           
       end
    end
    V = X0*W0;
    vn = zeros(size(cooo,1),1);
    for vi = 1:size(cooo,1)
        vn(vi,1) = V(5*vi-4,1)*nodeNormals(vi,1)+V(5*vi-3,1)*nodeNormals(vi,2)+V(5*vi-2,1)*nodeNormals(vi,3);
    end
    
    % 关于处理数据文件的部分
    fid = fopen('example4.dat', 'r');
    lines = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    lines = lines{1};

    % 2. 修改第 5 行为 i,0,0
    lines{5} = sprintf('%d,0,0', fi);

    % 3. 写回整个文件（覆盖写）
    fid = fopen('example4.dat', 'w');
    for j = 1:length(lines)
        fprintf(fid, '%s\n', lines{j});
    end
    fclose(fid);
    
    % 参数设置
    start_row = 3566;
    end_row = 6231;

    % ① 读取原始文件为逐行 cell
    fid = fopen('example4.dat', 'r');
    lines = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    lines = lines{1};

    % ② 构造要替换的新行内容
    new_lines = cell(end_row - start_row + 1, 1);
    for i = 1:length(new_lines)
        row_num =  i ;
        new_lines{i} = sprintf('%d,0,0,1,0,%.10f,0', row_num, vn(i,1));
    end

    % ③ 替换目标行
    lines(start_row:end_row) = new_lines;

    % ④ 写回文件
    fid = fopen('example66.dat', 'w');
    for i = 1:length(lines)
        fprintf(fid, '%s\n', lines{i});
    end 
    fclose(fid);
    
    status = system('FinalBEM.exe');    
end


