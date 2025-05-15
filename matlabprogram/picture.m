clc; clear; close all;

%% 1. 定义 3D 结构的节点坐标 (n × 3) 
NodeCoord = [0, 0, 0;   % 1号节点
             1, 0, 0;   % 2号节点
             1, 1, 0;   % 3号节点
             0, 1, 0;   % 4号节点
             0, 0, 1;   % 5号节点
             1, 0, 1;   % 6号节点
             1, 1, 1;   % 7号节点
             0, 1, 1];  % 8号节点

%% 2. 定义单元连接 (m × 2)，假设是 12 条梁单元（立方框架）
ElementConn = [1 2; 2 3; 3 4; 4 1;   % 底部四边
               5 6; 6 7; 7 8; 8 5;   % 顶部四边
               1 5; 2 6; 3 7; 4 8];  % 竖向四边

%% 3. 定义模态振型位移 (n × 3)
ModeShape = [ 0.1,  0.0,  0.0;  % 1号节点
              0.2,  0.1,  0.0;  % 2号节点
              0.1,  0.2,  0.0;  % 3号节点
              0.0,  0.1,  0.0;  % 4号节点
              0.1,  0.0,  0.1;  % 5号节点
              0.2,  0.1,  0.1;  % 6号节点
              0.1,  0.2,  0.1;  % 7号节点
              0.0,  0.1,  0.1]; % 8号节点

ScaleFactor = 5;  % 放大因子（根据情况调整）
DeformedCoord = NodeCoord + ScaleFactor * ModeShape;  % 计算变形后坐标

%% 4. 绘制 3D 模态振型
figure; hold on; axis equal; grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D 结构模态振型');

% 绘制原始结构 (黑色实线)
for i = 1:size(ElementConn,1)
    nodes = ElementConn(i,:);
    plot3(NodeCoord(nodes,1), NodeCoord(nodes,2), NodeCoord(nodes,3), 'k-', 'LineWidth', 1.5);
end

% 绘制变形后的结构 (红色虚线)
for i = 1:size(ElementConn,1)
    nodes = ElementConn(i,:);
    plot3(DeformedCoord(nodes,1), DeformedCoord(nodes,2), DeformedCoord(nodes,3), 'r--', 'LineWidth', 1.5);
end

% 显示节点位移方向 (蓝色箭头)
quiver3(NodeCoord(:,1), NodeCoord(:,2), NodeCoord(:,3), ...
        ModeShape(:,1), ModeShape(:,2), ModeShape(:,3), 0.5, 'b');

legend('原始结构', '变形结构', '位移方向', 'Location', 'Best');
view(3);  % 设置 3D 视角
