% 清空环境
clc; clear; close all;

% 定义 3D 网格的 8 个节点（四边形八节点单元）
X = [0 1 1 0  0.5 1 0.5 0]; % X 坐标
Y = [0 0 1 1  0 0.5 1 0.5]; % Y 坐标
Z = [0 0 0 0  0 0 0 0];     % Z 坐标（初始为平面）

% 变形后的坐标（假设一些小的变形）
Xd = X + [0 0.1 0.2 -0.1  0  0.1 0  -0.1];
Yd = Y + [0 0 0.1 0  0 -0.1 0.2 0.1];
Zd = Z + [0 0.1 -0.2 0.2  0.1 -0.1 0.2 0.1];

% 原始四边形网格
Xq = reshape(X, [2, 4]);
Yq = reshape(Y, [2, 4]);
Zq = reshape(Z, [2, 4]);

% 变形后的网格
Xqd = reshape(Xd, [2, 4]);
Yqd = reshape(Yd, [2, 4]);
Zqd = reshape(Zd, [2, 4]);

% 画图
figure;
subplot(1,2,1);
surf(Xq, Yq, Zq, 'FaceColor', 'cyan', 'EdgeColor', 'k'); % 原始结构
hold on; scatter3(X, Y, Z, 50, 'r', 'filled'); % 标记节点
title('原始结构');
axis equal; grid on;

subplot(1,2,2);
surf(Xqd, Yqd, Zqd, 'FaceColor', 'magenta', 'EdgeColor', 'k'); % 变形结构
hold on; scatter3(Xd, Yd, Zd, 50, 'b', 'filled'); % 标记节点
title('变形后结构');
axis equal; grid on;
