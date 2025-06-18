%% 该程序利用动态子结构法计算组合壳的固有频率和声响应
clc
clear

%密度 泊松比 杨氏模量 壳体厚度
density = 7900;
Niu = 0.3;
E = 196e9;
t = 0.05;
%% 导入网格信息 no 1 书上经典例子 2 组合壳带球尾模型
no = 2;
[msh] = qcl(no); 

nodes = msh.POS(:,:)./1000;
elements = msh.QUADS(:,1:4);
[cooo, en] = generateEightNodeElements(nodes, elements);
%% 重新排序网格
% [coords_sorted, idx] = sortrows(nodes, [3, 2, 1]);
% % [~, idx] = sort(nodes(:,3), 'ascend');
% % coords_sorted = nodes(idx, :);
% n = size(nodes, 1);
% invIdx = zeros(n,1);
% for newID = 1:n
%     oldID = idx(newID);
%     invIdx(oldID) = newID;
% end
% 
% elems_new = invIdx(elements);
% 
% 
% %子结构划分处
% con_z = 5;
% 
% z_vals = coords_sorted(:,3);  
% low_idx  = find(z_vals <= con_z);
% high_idx = find(z_vals >= con_z);
% 
% coords_low  = coords_sorted(low_idx, :);   % y<5 的节点
% coords_high = coords_sorted(high_idx,:);   % y>=5 的节点
% 
% inv_low = containers.Map(low_idx, 1:numel(low_idx));
% inv_high = containers.Map(high_idx, 1:numel(high_idx));
% 
% elems_low_list  = {};  
% elems_high_list = {};
% 
% for i = 1:size(elems_new,1)
%     nodes = elems_new(i,:);   % 第 i 个单元的 4 个节点新编号（范围 1..n）
% 
%     % 检查这 4 个节点是否都属于 low 组
%     if all( ismember(nodes, low_idx) )
%         % 把 “原来的节点编号” 映射到“coords_low 里对应的行索引”
%         newNodesLow = zeros(1,8);
%         for k = 1:8
%             newNodesLow(k) = inv_low(nodes(k));  
%         end
%         elems_low_list{end+1,1} = newNodesLow;  %#ok<SAGROW>
%     end
% 
%     % 检查这 4 个节点是否都属于 high 组
%     if all( ismember(nodes, high_idx) )
%         newNodesHigh = zeros(1,8);
%         for k = 1:8
%             newNodesHigh(k) = inv_high(nodes(k));
%         end
%         elems_high_list{end+1,1} = newNodesHigh;  %#ok<SAGROW>
%     end
% 
%     % 注意：如果一个单元的 4 个节点有的在 low，有的在 high，那么它不算进任何一组。
% end
% 
% % 把 cell 列表转换成常规数组（行数是符合条件的单元个数，列是 4）
% if isempty(elems_low_list)
%     elems_low = zeros(0,8);
% else
%     elems_low = vertcat(elems_low_list{:});
% end
% 
% if isempty(elems_high_list)
%     elems_high = zeros(0,8);
% else
%     elems_high = vertcat(elems_high_list{:});
% end
% 
% coool = coords_low;
% enl = elems_low;
% cooor = coords_high;
% enr = elems_high;
% 
% [jdzb1l, jdzbl, nodeNormalsl] = genShellSurfacesFromConn(coool, enl, t);
% [jdzb1r, jdzbr, nodeNormalsr] = genShellSurfacesFromConn(cooor, enr, t);
% 
% cooo = coool;
% cooo(1502:2666,:) = cooor(49:1213,:);
% en = elems_new;
%% 根据纵坐标划分子结构
% indices0 = find(cooo(:,3) > -10 & cooo(:,3) < 6.6 );
indices1 = find(cooo(:,3) < 5 );
indices2 = find(cooo(:,3) > 5 );

% 固定界面的坐标索引
indices3 = find(cooo(:,3) == 5);

% 整体结构的约束
indices4 = find(cooo(:,3) == 6.6);

ssss = size(en,1);
selected_elems = false(ssss,1);   % 用逻辑向量标记哪些单元被选中

for si = 1:ssss
    % elems(i,:) 是第 i 个单元的节点编号数组（长度 = p）
    % ismember(elems(i,:), target_nodes) 会返回一个 1×p 的逻辑数组，
    % 如果有任意一个元素为 true，就说明该单元含有目标节点。
    if any( ismember(en(si,:), indices3) )
        selected_elems(si) = true;
    end
end
elem_ids = find(selected_elems);
%% 有限元获得K、M
%高斯积分的点数
No_INTpoint_x = 2;        
No_INTpoint_y = 2;
No_INTpoint_z = 2;

k(1:size(cooo,1)*5,1:size(cooo,1)*5)=0;        
m(1:size(cooo,1)*5,1:size(cooo,1)*5)=0; 
% k = zeros((size(cooo,1)-size(indices4,1))*5);        
% m = zeros((size(cooo,1)-size(indices4,1))*5); 
subkl(1:size(cooo,1)*5,1:size(cooo,1)*5)=0;  
subkr(1:size(cooo,1)*5,1:size(cooo,1)*5)=0;  
subml(1:size(cooo,1)*5,1:size(cooo,1)*5)=0;  
submr(1:size(cooo,1)*5,1:size(cooo,1)*5)=0;  

dof = 0;
disp(1:size(cooo,1),1:5) = 1;  
% disp(indices4,:) = 0;

for ni=1:size(cooo,1) 
    for nj=1:5
        if disp(ni,nj)~=0
            dof=dof+1;
            disp(ni,nj)=dof;
        end
    end
end

subklidx = [];     % 用来收集满足条件的单元行号
subkridx = [];
for i = 1:size(elem_ids,1)  % 遍历每一个单元
    node_ids = en(elem_ids(i,1),:);            % 1×8，这个单元的 8 个节点编号
    z_vals   = cooo(node_ids,3);     % 取出这 8 个节点的 y 坐标，1×8

    if any(z_vals < 5)
        % 只要有一个 y<5，就把单元 i 记录下来
        subklidx(end+1) = elem_ids(i,1);    
    else
        subkridx(end+1) = elem_ids(i,1);
    end
end


% jdzb(1:1501,:) = jdzbl;
% jdzb(1502:2666,:) = jdzbr(49:1213,:);
% jdzb1(1:1501,:) = jdzb1l;
% jdzb1(1502:2666,:) = jdzb1r(49:1213,:);
% nodeNormals(1:1501,:) = nodeNormalsl;
% nodeNormals(1502:2666,:) = nodeNormalsr(49:1213,:);

dybh=en;
index(1:40)=0;
ekl = zeros(15,15);
eml = zeros(15,15);
ekr = zeros(15,15);
emr = zeros(15,15);
jdindex = zeros(3,1);
jdil = zeros(15,1);
jdir = zeros(15,1);

[jdzb1, jdzb, nodeNormals] = genShellSurfacesFromConn(cooo, en, t);

for loopi=1:size(en,1)
    dyhm=loopi;
     [ek,theta,xv2i,xv1i,xv3i,zmtemp,v3i,D,jtemp,xv3ii]=shellek(E,Niu,t,dyhm,jdzb,jdzb1,dybh);   
     [em1,tt]=shellem(zmtemp,v3i,density,theta,t,xv2i,xv1i,No_INTpoint_x,No_INTpoint_y,No_INTpoint_z,jtemp);
  
     if ismember(dyhm, [1, 336, 338, 390])
         aaa = 1;
     end
     
     
    for zi=1:8
        index((zi-1)*5+1)=disp(en(loopi,zi),1);
        index((zi-1)*5+2)=disp(en(loopi,zi),2);
        index((zi-1)*5+3)=disp(en(loopi,zi),3);
        index((zi-1)*5+4)=disp(en(loopi,zi),4);
        index((zi-1)*5+5)=disp(en(loopi,zi),5);
    end
 
    if ismember(dyhm,subklidx)
       ii0 = 0;
        for jdi =1:8
            if cooo(en(dyhm,jdi),3) == 5
                ii0 = ii0+1;
                jdindex(ii0,1) = jdi; 
                jdil(ii0*5-4,1) = jdi*5-4;
                jdil(ii0*5-3,1) = jdi*5-3;
                jdil(ii0*5-2,1) = jdi*5-2;
                jdil(ii0*5-1,1) = jdi*5-1;
                jdil(ii0*5-0,1) = jdi*5-0;
            end
        end
                ekl = ek(jdil,jdil);
                eml = em1(jdil,jdil);
                indexl = index(1,jdil);
                subkl(indexl,indexl) = subkl(indexl,indexl)+ekl;
                subml(indexl,indexl) = subml(indexl,indexl)+eml;
      
    elseif ismember(dyhm,subkridx)
        ii0 = 0;
        for jdi =1:8
            if cooo(en(dyhm,jdi),3) == 5
                ii0 = ii0+1;
                jdindex(ii0,1) = jdi;
                jdir(ii0*5-4,1) = jdi*5-4;
                jdir(ii0*5-3,1) = jdi*5-3;
                jdir(ii0*5-2,1) = jdi*5-2;
                jdir(ii0*5-1,1) = jdi*5-1;
                jdir(ii0*5-0,1) = jdi*5-0;
            end
        end
                ekr = ek(jdir,jdir);
                emr = em1(jdir,jdir);
                indexr = index(1,jdir);
                subkr(indexr,indexr) = subkr(indexr,indexr)+ekr;
                submr(indexr,indexr) = submr(indexr,indexr)+emr;       
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
rows_with_node = find(any(elements == 1, 2));
%% 附加质量计算
%     TT = zeros(3*size(nodeNormals,1),size(nodeNormals,1)); 
    m0(1:size(cooo,1)*5,1:size(cooo,1)*5)=0; 
%     for i = 1:size(nodeNormals,1)
%         TT(3*i-2,i) = nodeNormals(i,1);
%         TT(3*i-1,i) = nodeNormals(i,2);
%         TT(3*i-0,i) = nodeNormals(i,3);
%     end
% 
%     fi = 5.4;
%     [Ma] = add_mass(fi,cooo,en);
%     Ma1 = TT*Ma*TT';
%     [NewMat] =  expandDiagonalMatrix(Ma1,nodeNormals);
%     for i=1:size(m,1)
%         for j=1:size(m,1)
%             
%         m0(i,j)=m0(i,j)+NewMat(i,j);
%         
%         end
%     end
%% BC 得各结构的k、m、T
nFree = 5;  % 每个节点的自由度
nn = 40;
% 将节点编号转换为自由度索引函数
node2dof =@(indices0) reshape((1:nFree)' + nFree*(indices0(:)' - 1), 1, []);   %给定一组节点编号（indices0），返回它们对应的自由度（DOF）全局索引列表。

% conDOF = node2dof(indices4);
% k(conDOF,:) = [];
% k(:,conDOF) = [];
% m(conDOF,:) = [];
% m(:,conDOF) = [];

conDOF = node2dof(indices4);

% 子结构 I
oldDOF1 = node2dof(indices1);
newDOF1 = node2dof(indices3);
allDOF1 = [oldDOF1, newDOF1];
% m0 = m + m0;
kl = k;
kl(newDOF1,newDOF1) = 0;
kl(newDOF1,newDOF1) = subkl(newDOF1,newDOF1);
ml = m;
ml(newDOF1,newDOF1) = 0;
ml(newDOF1,newDOF1) = subml(newDOF1,newDOF1);
subk1 = kl(allDOF1, allDOF1);
subm1 = ml(allDOF1, allDOF1);

% subk1 = k(allDOF1, allDOF1);
% subm1 = m(allDOF1, allDOF1);

mii1 = subm1(1:size(indices1,1)*5,1:size(indices1,1)*5);
kii1 = subk1(1:size(indices1,1)*5,1:size(indices1,1)*5);
kij1 = subk1(1:size(indices1,1)*5,size(indices1,1)*5+1:(size(indices1,1)+size(indices3,1))*5);
I1=eye(5*size(indices3,1));
[xvec1,d1] = eigs(kii1,mii1,nn,'SM'); 
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
oldDOF2 = setdiff(oldDOF2,conDOF);
newDOF2 = node2dof(indices3);
allDOF2 = [newDOF2, oldDOF2];

kr = k;
kr(newDOF2,newDOF2) = 0;
kr(newDOF2,newDOF2) = subkr(newDOF2,newDOF2);
mr = m;
mr(newDOF2,newDOF2) = 0;
mr(newDOF2,newDOF2) = submr(newDOF2,newDOF2);
subk2 = kr(allDOF2, allDOF2);
subm2 = mr(allDOF2, allDOF2);
% 
% subk2 = k(allDOF2, allDOF2);
% subm2 = m(allDOF2, allDOF2);

mii2=subm2(size(indices3,1)*5+1:size(subm2,1),size(indices3,1)*5+1:size(subm2,1));
kii2=subk2(size(indices3,1)*5+1:size(subm2,1),size(indices3,1)*5+1:size(subm2,1));
kij2=subk2(size(indices3,1)*5+1:size(subm2,1),1:size(indices3,1)*5);
I2=eye(5*size(indices3,1));

[xvec2, d2] = eigs(kii2,mii2,nn,'SM'); 
tempd2=diag(d2);
[d2,sortindex2]=sort(tempd2);
xvec2=xvec2(:,sortindex2);
v2=xvec2(:,1:nn);
mode_number2=1:nn;
frequency2(mode_number2)=sqrt(d2(mode_number2))/(2*pi);
Frequency2=real(frequency2);
Frequency2=Frequency2';
V2=real(v2);


% T2=zeros(size(subk2,1),5*size(indices3,1)+nn);
% T2(1:size(indices3,1)*5,1:size(indices3,1)*5)=I2;
% T2(size(indices3,1)*5+1:size(subk2,1),1:size(indices3,1)*5)=-inv(kii2_reg)*kij2;
% T2(size(indices3,1)*5+1:size(subk2,1),size(indices3,1)*5+1:size(indices3,1)*5+nn)=v2;

T2=zeros(size(subk2,1),5*size(indices3,1)+nn);
T2(1:size(indices3,1)*5,1:size(indices3,1)*5)=I2;
T2(size(indices3,1)*5+1:size(subk2,1),1:size(indices3,1)*5)=-inv(kii2)*kij2;
T2(size(indices3,1)*5+1:size(subk2,1),size(indices3,1)*5+1:size(indices3,1)*5+nn)=v2;


K2=T2'*subk2*T2;
M2=T2'*subm2*T2;

%% 第二次坐标变换
condofs = size(indices3,1)*5;
KK=zeros((condofs+nn)*2,(condofs+nn)*2);
MM=zeros((condofs+nn)*2,(condofs+nn)*2);
KK(1:size(indices3,1)*5+nn,1:1:size(indices3,1)*5+nn)=K1;KK(size(indices3,1)*5+1+nn:(size(indices3,1)*5+nn)*2,size(indices3,1)*5+1+nn:(size(indices3,1)*5+nn)*2)=K2;
MM(1:size(indices3,1)*5+nn,1:1:size(indices3,1)*5+nn)=M1;MM(size(indices3,1)*5+1+nn:(size(indices3,1)*5+nn)*2,size(indices3,1)*5+1+nn:(size(indices3,1)*5+nn)*2)=M2;

% T3=zeros((size(indices3,1)*5+nn)*2,size(indices3,1)*5+nn*2);
% T3(1:nn,1:nn)=eye(nn);
% T3(nn+1:size(indices3,1)*5+nn,2*nn+1:size(indices3,1)*5+2*nn)=eye(size(indices3,1)*5);
% T3(nn+size(indices3,1)*5+1:size(indices3,1)*5*2+nn,2*nn+1:size(indices3,1)*5+2*nn)=eye(size(indices3,1)*5);
% T3(size(indices3,1)*5*2+nn+1:size(indices3,1)*5*2+2*nn,nn+1:2*nn)=eye(nn);

T3=zeros((condofs+nn)*2,condofs+nn*2);
T3(1:nn,1:nn)=eye(nn);
T3(nn+1:condofs+nn,2*nn+1:condofs+2*nn)=eye(condofs);
T3(nn+condofs+1:condofs*2+nn,2*nn+1:condofs+2*nn)=eye(condofs);
T3(condofs*2+nn+1:condofs*2+2*nn,nn+1:2*nn)=eye(nn);

K=T3'*KK*T3;
M=T3'*MM*T3;


nm=100;
% sigma = 1e-5;         % 指定一个小的非零移位值
% opts = struct('disp', 0);
% [v,d] = eigs(k, m, nm, sigma, opts);
% 
% sigma = 1e-3;         % 指定一个小的非零移位值
% opts = struct('disp', 0);
% [v,d] = eigs(K, M, nm, sigma, opts);
% [v,d] = eigs(k, m, nm, sigma, opts);
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
