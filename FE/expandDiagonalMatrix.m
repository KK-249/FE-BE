function [NewMat]= expandDiagonalMatrix(Ma,B_A)
% 假设 n 已定义，且 A 是 3*n x 3*n 的矩阵
n = size(B_A,1); % 示例，n可以任意定义
% Ma = rand(3*n, 3*n);  % 示例数据

% 初始化一个全零的 5*n x 5*n 矩阵
NewMat = zeros(5*n, 5*n);

% 将 A 放入 B 的对应位置（这里选择放在左上角）
for i = 1:n
    for j = 1:n
    
        NewMat(5*i-4,5*j-4) = Ma(3*i-2,3*j-2);
        NewMat(5*i-3,5*j-3) = Ma(3*i-1,3*j-1);
        NewMat(5*i-2,5*j-2) = Ma(3*i,3*j);
        NewMat(5*i-4,5*j-3) = Ma(3*i-2,3*j-1);
        NewMat(5*i-4,5*j-2) = Ma(3*i-2,3*j);
        NewMat(5*i-3,5*j-4) = Ma(3*i-1,3*j-2);
        NewMat(5*i-3,5*j-2) = Ma(3*i-1,3*j);
        NewMat(5*i-2,5*j-4) = Ma(3*i,3*j-2);
        NewMat(5*i-2,5*j-3) = Ma(3*i,3*j-1);
        
    end
end

end