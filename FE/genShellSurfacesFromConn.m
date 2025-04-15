function [jdzb1, jdzb, nodeNormals] = genShellSurfacesFromConn(cooo1, en1, t)
% genShellSurfacesFromConn 根据中面节点和单元连接信息生成上下表面的节点坐标，
% 并构造完整带厚度的壳单元连接矩阵
%
% 输入：
%   midSurface - n×3 矩阵，每行表示一个中面节点的 [x, y, z] 坐标
%   elem       - m×8 矩阵，每行存储一个单元8个节点的编号，其中：
%                   节点1-4为角点（构成四边形），节点5-8为中点
%   t          - 壳体厚度
%
% 输出：
%   upperNodes - n×3 矩阵，上表面各节点坐标
%   lowerNodes - n×3 矩阵，下表面各节点坐标
%   fullElem   - m×16 矩阵，完整壳单元的连接矩阵，
%                前8个节点对应上表面，后8个节点对应下表面（下表面节点编号 = 上表面节点编号 + n）

    nNodes = size(cooo1, 1);
    nElem  = size(en1, 1);
    
    % 初始化每个节点的累积法向量及计数器
    nodeNormals = zeros(nNodes, 3);
    count       = zeros(nNodes, 1);
    
    % 对每个单元计算单元法向量，采用四边形角点（节点1-4）构成的四边形
    % 将四边形分为两个三角形计算法向量：三角形1：(1,2,3)；三角形2：(1,3,4)
    for iElem = 1:nElem
        % 取角点编号：节点1-4
        corners = en1(iElem, [1, 2, 3, 4]);
        p1 = cooo1(corners(1), :);
        p2 = cooo1(corners(2), :);
        p3 = cooo1(corners(3), :);
        p4 = cooo1(corners(4), :);
        
        % 计算两个三角形的法向量
        normal1 = cross(p2 - p1, p3 - p1);
        normal2 = cross(p3 - p1, p4 - p1);
        elemNormal = normal1 + normal2;
        normElem = norm(elemNormal);
        if normElem ~= 0
            elemNormal = elemNormal / normElem;
        end
        
        % 将该单元的法向量累加到单元内所有8个节点上
        for j = 1:8
            node = en1(iElem, j);
            nodeNormals(node, :) = nodeNormals(node, :) + elemNormal;
            count(node) = count(node) + 1;
        end
    end
    
    % 对每个节点求平均后归一化，得到节点单位法向量
    for iNode = 1:nNodes
        if count(iNode) > 0
            nodeNormals(iNode, :) = nodeNormals(iNode, :) / count(iNode);
            normN = norm(nodeNormals(iNode, :));
            if normN ~= 0
                nodeNormals(iNode, :) = nodeNormals(iNode, :) / normN;
            end
        else
            warning('节点 %d 没有相邻单元，法向量设为零', iNode);
        end
    end
    
    % 生成上下表面坐标：沿单位法向方向正负偏移 t/2
    jdzb1 = cooo1 + (t/2) * nodeNormals;
    jdzb = cooo1 - (t/2) * nodeNormals;
    
    % 构造完整壳单元的连接矩阵
    % 假定上表面节点编号为1~nNodes，下表面编号为 nNodes+1 ~ 2*nNodes
    % 则每个壳单元的16个节点编号为：
    % [上表面: 原单元的8个节点, 下表面: 原单元的8个节点 + nNodes]
    fullElem = [en1, en1 + nNodes];
end
