clear
clc
%% 书上经典例子 1 书上经典例子 2 组合壳带球尾模型 qcl
no=2;
% [msh]=ban(); 
[msh] = qcl(no);
%% 计算中点坐标并重组单元
% 原始节点坐标（四边形单元）
nodes = msh.POS(:,:);   
elements = msh.QUADS(:,1:4);
nodes(abs(nodes) < 0.0001) = 0;
% elements = elements(:, [1, 4, 3, 2]);  % 转换为顺时针排列
% 步骤1: 去除重复节点并获取唯一节点及其索引
% [unique_nodes, ~, idx] = unique(nodes, 'rows', 'stable');

% 步骤2: 更新单元矩阵的节点编号
% updated_elements = idx(elements);

% 调用生成八结点单元的函数
[new_nodes,eight_node_elements]=generateEightNodeElements(nodes, elements);

% 设定阈值
threshold = 6600;

% 找出节点纵坐标 Y > 6.6 的节点编号
nodes_above_threshold = find(nodes(:, 3) >= threshold);

% 初始化存储满足条件的单元号
matching_elements = [];

% 遍历单元
for i = 1:size(elements, 1)
    % 当前单元的节点编号
    current_nodes = elements(i, :);
    
    % 检查是否包含满足条件的节点
    if any(ismember(current_nodes, nodes_above_threshold))
        matching_elements = [matching_elements; i]; % 保存单元号
    end
end
%% 输出文件
% % % new_nodes(:,:)=new_nodes(:,:)/1000;
% new_nodes = new_nodes./1000;
% % 生成行号并拼接到矩阵左侧
% [rowCount, ~] = size(new_nodes);       % 获取矩阵行数
% rowNumbers = (1:rowCount)';    % 生成行号列向量
% A_with_row = [rowNumbers, new_nodes];  % 拼接行号为第一列
% 
% % 打开文件
% fid = fopen('C:\Users\Administrator\Desktop\matlab shuju\shuju1. txt', 'w');
% 
% % 循环逐行写入文件
% for i = 1:rowCount
%     fprintf(fid, '%d,', A_with_row(i, 1));        % 输出行号
%     fprintf(fid, '%.4f,', A_with_row(i, 2:end-1)); % 输出中间的元素
%     fprintf(fid, '%.4f\n', A_with_row(i, end));    % 输出最后一个元素并换行
% end
% 
% % 关闭文件
% fclose(fid);
% 
% [rowCount, ~] = size(eight_node_elements);       % 获取矩阵行数
% rowNumbers = (1:rowCount)';    % 生成行号列向量
% A_with_row = [rowNumbers, eight_node_elements];  % 拼接行号为第一列
% 
% % 打开文件
% fid = fopen('C:\Users\Administrator\Desktop\matlab shuju\shuju2. txt', 'w');
% 
% % 循环逐行写入文件
% for i = 1:rowCount
%     fprintf(fid, '%d,', A_with_row(i, 1));        % 输出行号
%     fprintf(fid, '%d,', A_with_row(i, 2:end-1)); % 输出中间的元素
%     fprintf(fid, '%d\n', A_with_row(i, end));    % 输出最后一个元素并换行
% end
% 
% % 关闭文件
% fclose(fid);
%% 有限元部分
No_INTpoint_x=2;        %高斯积分的点数
No_INTpoint_y=2;
No_INTpoint_z=2;

%密度 泊松比 杨氏模量
density=7900;
Niu=0.3;
E=196e9;

%厚度
t=0.05;

en=eight_node_elements;       %element node  
cooo=new_nodes/1000;
% threshold = 1e-4;
% cooo(abs(cooo) < threshold) = 0;

new_nodes=cooo;
disp(1:size(cooo,1),1:5)=1;     % node displacement

%找到约束面节点编号
indices = find(cooo(:,3) == 0);
rowsToZero = indices;
disp(rowsToZero, :) = 0;
% disp(3:5,:)=0;

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
% jdzb=cooo;
% jdzb1=cooo;
% 
% jdzb(:,3)=cooo(:,3)+0.025;
% jdzb1(:,3)=cooo(:,3)-0.025;

%  [jdzb,jdzb1] = coo_cal(cooo);
[jdzb1, jdzb, nodeNormals] = genShellSurfacesFromConn(cooo, eight_node_elements, t);

% TT=zeros(3*size(nodes,1),size(nodes,1));
dybh=en;
index(1:40)=0; % vector sontaining system dofs of nodes in each element. 
%% 
for loopi=1:size(en,1)
    dyhm=loopi;
     [ek,theta,xv2i,xv1i,xv3i,zmtemp,v3i,D,jtemp,xv3ii]=shellek(E,Niu,t,dyhm,jdzb,jdzb1,dybh);   
     [em1,tt]=shellem(zmtemp,v3i,density,theta,t,xv2i,xv1i,No_INTpoint_x,No_INTpoint_y,No_INTpoint_z,jtemp);
%    [em1]=shellem_corrected(density,new_nodes,dybh,total_area1,total_area2,total_area3);
        
     
    
    
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
%% 附加质量计算
B_A = load('mat_outputa.txt');
B_G = load('mat_outputb.txt');
B_H = load('mat_outputh.txt');
B_CP = load('mat_outputCP.txt');
TT = zeros(3*size(nodeNormals,1),size(nodeNormals,1));

for i = 1:size(nodeNormals,1)
    TT(3*i-2,i) = nodeNormals(i,1);
    TT(3*i-1,i) = nodeNormals(i,2);
    TT(3*i-0,i) = nodeNormals(i,3);
end

Ma = B_A*real(inv(B_H+B_CP))*B_G;
Ma1 = TT*Ma*TT';
[NewMat] =  expandDiagonalMatrix(Ma1,nodeNormals);

for i=1:size(m,1)
    for j=1:size(m,1)
    m(i,j)=m(i,j)+NewMat(i,j);
    end
end
%% 频域响应
freq = 200;
W0 = 2*pi*freq;
% F0 = 1.3083;
F = zeros(size(k,1),1);
% indices11 = find (indices1>890);
% F(indices1*5-2,1) = -1.090278*2;                  % 此处0.8723 是加在r=0.1m处的中点载荷 角点*2
% F(indices1(indices11,1)*5-2,1) = -1.090278;       % 此处1.090278 是加在r=0.1m处的中点载荷 角点*2
% F(indices1*5-2,1) = F0;      % 加激励注意自由度
F(43,1) = 100; 
% [eta,y,omega1,sdof2] = HarmonicRespt(k,m,omega0,a,b,indices1);  %时域响应

X = inv((k-W0^2*m))*F;
% X = inv(k)*F;
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
%%
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

% % 输出边界法向振速 文件shuju3.txt储存
% [rowCount, ~] = size(vn);       % 获取矩阵行数
% rowNumbers = (1:rowCount)';    % 生成行号列向量
% A_with_row = [rowNumbers, vn];  % 拼接行号为第一列
% fid = fopen('C:\Users\Administrator\Desktop\matlab shuju\shuju3. txt', 'w');
% AAA = [0 0 1 0];
% % 循环逐行写入文件
% for i = 1:rowCount
%     fprintf(fid, '%d,', A_with_row(i, 1));        % 输出行号
%     fprintf(fid, '%d,', AAA); % 输出中间的元素
%     fprintf(fid, '%d,', A_with_row(i, 2)); % 输出中间的元素
%     fprintf(fid, '%d\n', 0);    % 输出最后一个元素并换行
% end
% 
% % 关闭文件
% fclose(fid);
%% 自然频率计算
nm=50;

sigma = 1e-3;         % 指定一个小的非零移位值
opts = struct('disp', 0);
[v,d] = eigs(k, m, nm, sigma, opts);
% [v,d] =eigs(k,m,nm,'SM');  %v特征向量 d特征值
tempd=diag(d);
[d,sortindex]=sort(tempd);
omega=sqrt(d);
v=v(:,sortindex);
mode_number=1:nm;
frequency(mode_number)=sqrt(d(mode_number))/(2*pi);
Frequency=real(frequency);
Frequency=Frequency';
%% 绘制模态图
close all
xvec = v;
for i =1:12
    
    mn = 6+2*i;
    xx=xvec(1:5:size(k,1)-4,mn);
    yy=xvec(2:5:size(k,1)-3,mn);
    zz=xvec(3:5:size(k,1)-2,mn);

    ScaleFactor = 55;  % 放大因子（调节可视化效果）
    condition1 = disp(:,1)~=0;
    bx=cooo;
    bx(condition1,:) = cooo(condition1,:) + ScaleFactor * [xx, yy, zz];  % 计算变形后的坐标
    final_bx(:,:) = bx(1:890,:);
    y_coord = final_bx(:, 3); 
    
    figure;
    ax1 = axes('Position', [0.05 0.15 0.6 0.75]);  % [左 下 宽 高]
%     subplot(1, 2, 1);
    patch('Faces', elements, ...
          'Vertices', final_bx, ...
          'FaceVertexCData', y_coord, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'black');
    axis equal; view(3)
    title('主要视图'); 
%     lighting gouraud         % 平滑光照模型
%     camlight('headlight');   % 光源随相机旋转
%     material dull            % 材质属性（可选 smooth / shiny / dull）
    
    ax2 = axes('Position', [0.65 0.4 0.3 0.3]);  % 比例较小
%     subplot(1, 2, 2);  % 1行2列的第2个图
    patch('Faces', elements, ...
          'Vertices', final_bx, ...
          'FaceVertexCData', y_coord, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'black');
    axis equal; view(0, 90);  % 俯视图的视角（方位角0，仰角90）
    title('俯视图'); 
%     lighting gouraud         % 平滑光照模型
%     camlight('headlight');   % 光源随相机旋转
%     material dull            % 材质属性（可选 smooth / shiny / dull）
end
% % P=23;                                %模态阶数
% for i=1:22
% %     mn=50;  
%     mn=6+2*i;                           
%     xx=xvec(1:5:size(k,1)-4,mn);
%     yy=xvec(2:5:size(k,1)-3,mn);
%     zz=xvec(3:5:size(k,1)-2,mn);
% 
%     ScaleFactor = 100;  % 放大因子（调节可视化效果）
%     condition1 = disp(:,1)~=0;
%     bx=cooo;
%     bx(condition1,:) = cooo(condition1,:) + ScaleFactor * [xx, yy, zz];  % 计算变形后的坐标
% 
%     % 定义 3D 网格的 8 个节点（四边形八节点单元）
%     X = cooo(:,1); % X 坐标
%     Y = cooo(:,2); % Y 坐标
%     Z = cooo(:,3); % Z 坐标（初始为平面）
% 
%     % 变形后的坐标（假设一些小的变形）
%     Xd = bx(:,1);
%     Yd = bx(:,2);
%     Zd = bx(:,3);
% 
%     % 原始四边形网格
%     Xq = reshape(X, [2, 1189]);
%     Yq = reshape(Y, [2, 1189]);
%     Zq = reshape(Z, [2, 1189]);
% 
%     % 变形后的网格
%     Xqd = reshape(Xd, [2, 1189]);
%     Yqd = reshape(Yd, [2, 1189]);
%     Zqd = reshape(Zd, [2, 1189]);
% 
%     figure;
%     subplot(1,2,1);
%     surf(Xq, Yq, Zq, 'FaceColor', 'cyan', 'EdgeColor', 'k'); % 原始结构
%     hold on; scatter3(X, Y, Z, 10, 'r', 'filled');           % 标记节点
%     title('原始结构');
%     axis equal; grid on;
% 
%     subplot(1,2,2);
%     surf(Xqd, Yqd, Zqd, 'FaceColor', 'magenta', 'EdgeColor', 'k'); % 变形结构
%     hold on; scatter3(Xd, Yd, Zd, 10, 'g', 'filled');              % 标记节点
%     title('变形后结构');
%     axis equal; grid on;
% end
