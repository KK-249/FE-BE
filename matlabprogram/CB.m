%% 该程序利用动态子结构法计算组合壳的固有频率和声响应 分别组装
clc
clear
%% 导入网格信息 no 1 书上经典例子 2 组合壳带球尾模型
no = 2;
nn = 40;
[msh] = qcl(no); 

nodes = msh.POS(:,:)./1000;
elements = msh.QUADS(:,1:4);
[nodes, elements] = generateEightNodeElements(nodes, elements);

[coords_sorted, idx] = sortrows(nodes, [3, 2, 1]);
% [~, idx] = sort(nodes(:,3), 'ascend');
% coords_sorted = nodes(idx, :);
n = size(nodes, 1);
invIdx = zeros(n,1);
for newID = 1:n
    oldID = idx(newID);
    invIdx(oldID) = newID;
end

elems_new = invIdx(elements);

%子结构划分处
con_z = 5;

z_vals = coords_sorted(:,3);  
low_idx  = find(z_vals <= con_z);
high_idx = find(z_vals >= con_z);

coords_low  = coords_sorted(low_idx, :);   % y<5 的节点
coords_high = coords_sorted(high_idx,:);   % y>=5 的节点

inv_low = containers.Map(low_idx, 1:numel(low_idx));
inv_high = containers.Map(high_idx, 1:numel(high_idx));

elems_low_list  = {};  
elems_high_list = {};

for i = 1:size(elems_new,1)
    nodes = elems_new(i,:);   % 第 i 个单元的 4 个节点新编号（范围 1..n）

    % 检查这 4 个节点是否都属于 low 组
    if all( ismember(nodes, low_idx) )
        % 把 “原来的节点编号” 映射到“coords_low 里对应的行索引”
        newNodesLow = zeros(1,8);
        for k = 1:8
            newNodesLow(k) = inv_low(nodes(k));  
        end
        elems_low_list{end+1,1} = newNodesLow;  %#ok<SAGROW>
    end

    % 检查这 4 个节点是否都属于 high 组
    if all( ismember(nodes, high_idx) )
        newNodesHigh = zeros(1,8);
        for k = 1:8
            newNodesHigh(k) = inv_high(nodes(k));
        end
        elems_high_list{end+1,1} = newNodesHigh;  %#ok<SAGROW>
    end

    % 注意：如果一个单元的 4 个节点有的在 low，有的在 high，那么它不算进任何一组。
end

% 把 cell 列表转换成常规数组（行数是符合条件的单元个数，列是 4）
if isempty(elems_low_list)
    elems_low = zeros(0,8);
else
    elems_low = vertcat(elems_low_list{:});
end

if isempty(elems_high_list)
    elems_high = zeros(0,8);
else
    elems_high = vertcat(elems_high_list{:});
end

coool = coords_low;
enl = elems_low;
cooor = coords_high;
enr = elems_high;
%% 子结构部分 函数分别计算
[KA,MA,T1] = SubA(coool,enl,nn,con_z);
[KB,MB,T2] = SubB(cooor,enr,nn,con_z);

indices3 = find(coords_sorted(:,3) == con_z);
condofs = size(indices3,1)*5;
KK=zeros((condofs+nn)*2,(condofs+nn)*2);
MM=zeros((condofs+nn)*2,(condofs+nn)*2);
KK(1:size(indices3,1)*5+nn,1:1:size(indices3,1)*5+nn)=KA;KK(size(indices3,1)*5+1+nn:(size(indices3,1)*5+nn)*2,size(indices3,1)*5+1+nn:(size(indices3,1)*5+nn)*2)=KB;
MM(1:size(indices3,1)*5+nn,1:1:size(indices3,1)*5+nn)=MA;MM(size(indices3,1)*5+1+nn:(size(indices3,1)*5+nn)*2,size(indices3,1)*5+1+nn:(size(indices3,1)*5+nn)*2)=MB;

T3=zeros((condofs+nn)*2,condofs+nn*2);
T3(1:nn,1:nn)=eye(nn);
T3(nn+1:condofs+nn,2*nn+1:condofs+2*nn)=eye(condofs);
T3(nn+condofs+1:condofs*2+nn,2*nn+1:condofs+2*nn)=eye(condofs);
T3(condofs*2+nn+1:condofs*2+2*nn,nn+1:2*nn)=eye(nn);

K=T3'*KK*T3;
M=T3'*MM*T3;

nm=30;
sigma = 1e-5;         % 指定一个小的非零移位值
opts = struct('disp', 0);
[v,d] = eigs(K, M, nm, sigma, opts);
% [v,d] =eigs(K,M,nm,'SM');  %v特征向量 d特征值
tempd=diag(d);
[d,sortindex]=sort(tempd);
omega=sqrt(d);
v=v(:,sortindex);
mode_number=1:nm;
frequency(mode_number)=sqrt(d(mode_number))/(2*pi);
Frequency=real(frequency);
Frequency=Frequency';
V=real(v);


% data1 = load('matrixA.mat');  
% data2 = load('matrixB.mat');  
% A1    = data1.K1;               % 重命名为 A1
% A2    = data2.M1; 
% 
% % 现在就可以比较 A1, A2
% tf = isequal(A1, KA);
% fprintf('两矩阵完全相同吗？ %d\n', tf);
% 
% tf = isequal(A2, MA);
% fprintf('两矩阵完全相同吗？ %d\n', tf);
