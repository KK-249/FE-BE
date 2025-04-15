function [new_nodes, eight_node_elements] = generateEightNodeElements(nodes, elements)
    % generateEightNodeElements 生成八结点四边形单元
    %
    % 输入:
    %   nodes    - N x 2 的节点坐标矩阵
    %   elements - M x 4 的单元连接矩阵（四结点单元）
    %
    % 输出:
    %   new_nodes          - 新的节点坐标矩阵
    %   eight_node_elements - 新的八结点单元连接矩阵

    % 初始化
    new_nodes = nodes;  % 新节点列表，初始为原节点
    edge_midpoint_map = containers.Map('KeyType', 'char', 'ValueType', 'int32');  % 边到中点节点编号的映射
    eight_node_elements = [];  % 新单元列表

    % 遍历每个单元
    for elem_idx = 1:size(elements, 1)
        element = elements(elem_idx, :);  % 当前单元的四个角点节点编号
        midpoints = zeros(1, 4);  % 存储当前单元四条边的中点节点编号

        % 遍历单元的四条边
        for i = 1:4
            n1 = element(i);
            n2 = element(mod(i,4) +1);  % 下一个节点，确保闭合
            edge = sort([n1, n2]);  % 对边的节点编号排序，确保唯一性
            edge_key = sprintf('%d-%d', edge(1), edge(2));  % 边的唯一键

            if isKey(edge_midpoint_map, edge_key)
                % 如果边的中点已存在，直接获取中点节点编号
                midpoints(i) = edge_midpoint_map(edge_key);
            else
                % 计算中点坐标
                midpoint = (nodes(n1, :) + nodes(n2, :)) / 2;
                new_nodes = [new_nodes; midpoint];  % 添加新节点
                mid_idx = size(new_nodes, 1);  % 新节点的编号
                edge_midpoint_map(edge_key) = mid_idx;  % 更新映射
                midpoints(i) = mid_idx;  % 记录中点节点编号
            end
        end

        % 获取原单元的四个角点和四个边中点
        n1 = element(1);
        n2 = element(2);
        n3 = element(3);
        n4 = element(4);
        m1 = midpoints(1);  % 边 n1-n2 的中点
        m2 = midpoints(2);  % 边 n2-n3 的中点
        m3 = midpoints(3);  % 边 n3-n4 的中点
        m4 = midpoints(4);  % 边 n4-n1 的中点

        % 定义八结点单元的节点顺序（角点和边中点）
        eight_node = [n1, n2, n3, n4, m1, m2, m3, m4];
        eight_node_elements = [eight_node_elements; eight_node];
    end
end
