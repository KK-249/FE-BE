function [jdzb,jdzb1] = coo_cal(cooo)
% 根据中面坐标计算上下表面坐标
% jdzb 下表面 jdzb1 上表面
% 由于各个壳体部分的形状 分三种情况计算

delta_r = 0.025;        % 半径减少值
r1 = 5;               % 原半径计算
r1_new = r1 - delta_r;
r1_new1= r1 + delta_r;
jdzb = cooo;            % 给定初始值
jdzb1 = cooo;

%% 圆柱面坐标映射 
% %（第二套程序） 计算example_book 
% % 条件 1：第三列值小于 10
% condition1 = (cooo(:, 3) <= 60)&(cooo(:,3)>=0) & (sqrt(cooo(:,1).^2 + cooo(:,2).^2) >= 4.975);
% % 映射到新圆柱面的坐标
%     jdzb1(condition1, 1) = cooo(condition1, 1) .* (r1_new ./ r1);    % x 坐标更新
%     jdzb1(condition1, 2) = cooo(condition1, 2) .* (r1_new ./ r1);    % y 坐标更新
%     jdzb(condition1, 1) = cooo(condition1, 1) .* (r1_new1 ./ r1);  % x 坐标更新
%     jdzb(condition1, 2) = cooo(condition1, 2) .* (r1_new1 ./ r1);  % y 坐标更新
%     
% condition2 = sqrt(cooo(:,1).^2 + cooo(:,2).^2) < 4.975;
% % 映射到新圆柱面的坐标
%     jdzb1(condition2, 1) = cooo(condition2, 1);     % x 坐标更新
%     jdzb1(condition2, 2) = cooo(condition2, 2) ;  % y 坐标更新
%     jdzb(condition2, 1) = cooo(condition2, 1)   ;% x 坐标更新
%     jdzb(condition2, 2) = cooo(condition2, 2) ; % y 坐标更新
%     jdzb1(condition2, 3) = cooo(condition2, 3) + delta_r;     % x 坐标更新
%     jdzb(condition2, 3) = cooo(condition2, 3) + delta_r;  % y 坐标更新

%% 第一套程序  用于计算组合壳结构

%% 圆柱面坐标映射

% 条件 1：第三列值小于 10
condition1 = (cooo(:, 3) <= 5)&(cooo(:,3)>=0);
% 映射到新圆柱面的坐标
    jdzb1(condition1, 1) = cooo(condition1, 1) .* (r1_new ./ r1);    % x 坐标更新
    jdzb1(condition1, 2) = cooo(condition1, 2) .* (r1_new ./ r1);    % y 坐标更新
    jdzb(condition1, 1) = cooo(condition1, 1) .* (r1_new1 ./ r1);  % x 坐标更新
    jdzb(condition1, 2) = cooo(condition1, 2) .* (r1_new1 ./ r1);  % y 坐标更新

%% 圆锥壳坐标映射
% 条件 2：第三列值介于 10 和 20 之间
condition2 = (cooo(:, 3) > 5)& (cooo(:,3)<=6.6);
% z0 = 0;         % 圆锥顶点 z0

% 原半径计算
r = sqrt(cooo(:, 1).^2 + cooo(:, 2).^2);
r2_new = r - delta_r;
r2_new1 = r + delta_r;

    R2 = (r2_new ./ r); % 结果是一个列向量
    R21= (r2_new1 ./ r);

    % 映射到新圆锥面的坐标
    jdzb1(condition2, 1) = cooo(condition2, 1) .*  R2(condition2); % x 坐标更新
    jdzb1(condition2, 2) = cooo(condition2, 2) .*  R2(condition2); % y 坐标更新
    jdzb(condition2, 1) = cooo(condition2, 1) .*  R21(condition2); % x 坐标更新
    jdzb(condition2, 2) = cooo(condition2, 2) .*  R21(condition2); % y 坐标更新
% end
%% 球面坐标映射

for i=1:size(cooo,1)
    if cooo(i,3)<0
    % 当前节点坐标
    x0 = cooo(i, 1);
    y0 = cooo(i, 2);
    z0 = cooo(i, 3);
%   for i=1:size(condition3,1)  
    % 计算法向量
    norm_factor = sqrt(x0^2 + y0^2 + z0^2);
    n = [x0, y0, z0] / norm_factor;  % 法向量（单位向量）
   
    % 计算上表面和下表面坐标
    jdzb(i, :) = cooo(i, :) + 0.025 * n;      % 上表面
    jdzb1(i, :) = cooo(i, :) - 0.025 * n;     % 下表面
    elseif cooo(i,3)>6.6
    x0 = cooo(i, 1);
    y0 = cooo(i, 2);
    z0 = cooo(i, 3);
%   for i=1:size(condition3,1)  
    % 计算法向量
    norm_factor = sqrt(x0^2 + y0^2 + z0^2);
    n = [x0, y0, z0] / norm_factor;  % 法向量（单位向量）
   
    % 计算上表面和下表面坐标
    jdzb(i, :) = cooo(i, :) + 0.025 * n;      % 上表面
    jdzb1(i, :) = cooo(i, :) - 0.025 * n;     % 下表面
    end
end
