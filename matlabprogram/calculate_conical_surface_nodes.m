function [jdzb, jdzb1] = calculate_conical_surface_nodes(cooo, total_height)
    % mid_height_nodes: 二分之一高度处的节点坐标 (n x 3)
    % total_height: 圆台的总高度 (这里假设 total_height = 0.05)
    
    % 圆台的总高度
    h = total_height;
    
    % 计算给定节点在中截面的半径
    % 假设 mid_height_nodes 中的第二列表示半径（即 x 或 y 坐标的模长）
    r_mid = sqrt(cooo(:,1).^2 + cooo(:,2).^2);
    
    % 定义上下表面半径的变化比例
    % 假设圆台在 z = 0 处为底面，z = h 处为顶面
    r_bottom = r_mid * 1.6;  % 底面半径，假设为中截面半径的一半
    r_top = r_mid * 0.4;  % 顶面半径，假设与中截面相同
    
    % 计算上下表面的节点坐标
    % 1. 计算底面坐标
    jdzb1 = cooo;
    jdzb1(:,3) = 0;  % 底面的 z 坐标为 0（假设底面在 z = 0 处）
    
    % 2. 计算顶面坐标
    jdzb = cooo;
    jdzb(:,3) = h;  % 顶面的 z 坐标为总高度 h
    
    % 根据比例调整坐标
    for i = 1:size(cooo, 1)
        % 计算当前节点的径向比例
        r_mid_node = r_mid(i);
        
        % 计算顶面和底面对应的节点坐标
        bottom_radius = r_bottom * (r_mid_node / max(r_mid));
        top_radius = r_top * (r_mid_node / max(r_mid));
        
        % 底面和顶面坐标 (在圆台内的插值模型)
        jdzb1(i,1:2) = bottom_radius * cooo(i,1:2) / r_mid_node;
        jdzb(i,1:2) = top_radius * cooo(i,1:2) / r_mid_node;
    end
end