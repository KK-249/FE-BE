function [NewMat]= expandDiagonalMatrix1(Ma,B_A)

    n = size(B_A, 1);              % 获取 Ma 的行数（也就是列数）
    NewMat = zeros(5*n, 5*n);       % 初始化 5n x 5n 的零矩阵

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