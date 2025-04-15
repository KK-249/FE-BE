function [em1]=shellem_corrected(density,new_nodes,dybh,total_area1,total_area2,total_area3)
em1=zeros(40,40);

V1=pi*1.6*(0.425^2-0.375^2)/3;    %圆锥壳体积
V2=5*pi*(0.425^2-0.375^2);        %圆柱壳体积
V3=2*pi*(0.425^3-0.375^3)/3;      %半球壳体积

for i = 1:size(dybh, 1)
    % 获取当前单元的节点坐标
    idx = dybh(i, :);  % 当前单元的节点编号
    P1 = new_nodes(idx(1), :);  % 节点P1坐标
    P2 = new_nodes(idx(2), :);  % 节点P2坐标
    P3 = new_nodes(idx(3), :);  % 节点P3坐标
    P4 = new_nodes(idx(4), :);  % 节点P4坐标
    P=[P1;
        P2;
        P3;
        P4];
    if any(P(:, 3) > 5)
            x1 = P(1,1); y1 = P(1,2); z1=P(1,3);
            x2 = P(2,1); y2 = P(2,2); z2=P(2,3);
            x3 = P(3,1); y3 = P(3,2); z3=P(3,3);
            x4 = P(4,1); y4 = P(4,2); z4=P(4,3);   
            r1=[x1,y1,z1];
            r2=[x2,y2,z2];
            r3=[x3,y3,z3];
            r4=[x4,y4,z4];
            
            % 计算三角形 (1, 2, 3) 的面积
            v1 = r2 - r1; % 向量 (r2 - r1)
            v2 = r3 - r1; % 向量 (r3 - r1)
            A1 = abs(0.5 * norm(cross(v1, v2))); % 计算叉积并求模，得到三角形面积

            % 计算三角形 (1, 3, 4) 的面积
            v3 = r3 - r1; % 向量 (r3 - r1)
            v4 = r4 - r1; % 向量 (r4 - r1)
            A2 = abs(0.5 * norm(cross(v3, v4))); % 计算叉积并求模，得到三角形面积
       
            % 计算四边形单元的面积（使用行列式法）
            area = A1+A2;

            % 计算每个单元的面积占比
            area_ratios = area / total_area1;
            W1=area_ratios*V1*density;
            m=W1/8*[1 0 0 0 0;
                    0 1 0 0 0;
                    0 0 1 0 0;
                    0 0 0 0 0;
                    0 0 0 0 0;];
            em1 = blkdiag(m, m, m, m, m, m, m, m);
            
    elseif any(P(:, 3) < 0)
            x1 = P(1,1); y1 = P(1,2); z1=P(1,3);
            x2 = P(2,1); y2 = P(2,2); z2=P(2,3);
            x3 = P(3,1); y3 = P(3,2); z3=P(3,3);
            x4 = P(4,1); y4 = P(4,2); z4=P(4,3);
            r1=[x1,y1,z1];
            r2=[x2,y2,z2];
            r3=[x3,y3,z3];
            r4=[x4,y4,z4];
            
            % 计算三角形 (1, 2, 3) 的面积
            v1 = r2 - r1; % 向量 (r2 - r1)
            v2 = r3 - r1; % 向量 (r3 - r1)
            A1 = abs(0.5 * norm(cross(v1, v2))); % 计算叉积并求模，得到三角形面积
    
            % 计算三角形 (1, 3, 4) 的面积
            v3 = r3 - r1; % 向量 (r3 - r1)
            v4 = r4 - r1; % 向量 (r4 - r1)
            A2 = abs(0.5 * norm(cross(v3, v4))); % 计算叉积并求模，得到三角形面积
            % 计算四边形单元的面积（使用行列式法）
            area = A1+A2;
            area_ratios = area / total_area2;
            W3=area_ratios*V3*density;
            m=W3/8*[1 0 0 0 0;
                0 1 0 0 0;
                0 0 1 0 0;
                0 0 0 0 0;
                0 0 0 0 0;];
            em1 = blkdiag(m, m, m, m, m, m, m, m);
    else
            x1 = P(1,1); y1 = P(1,2); z1=P(1,3);
            x2 = P(2,1); y2 = P(2,2); z2=P(2,3);
            x3 = P(3,1); y3 = P(3,2); z3=P(3,3);
            x4 = P(4,1); y4 = P(4,2); z4=P(4,3);
            r1=[x1,y1,z1];
            r2=[x2,y2,z2];
            r3=[x3,y3,z3];
            r4=[x4,y4,z4];
            
            % 计算三角形 (1, 2, 3) 的面积
            v1 = r2 - r1; % 向量 (r2 - r1)
            v2 = r3 - r1; % 向量 (r3 - r1)
            A1 = abs(0.5 * norm(cross(v1, v2))); % 计算叉积并求模，得到三角形面积
            % 计算三角形 (1, 3, 4) 的面积
            v3 = r3 - r1; % 向量 (r3 - r1)
            v4 = r4 - r1; % 向量 (r4 - r1)
            A2 = abs(0.5 * norm(cross(v3, v4))); % 计算叉积并求模，得到三角形面积
            
             % 计算四边形单元的面积（使用行列式法）
            area = A1+A2;
        
            area_ratios = area / total_area3;
            W2=area_ratios*V2*density;
            m=W2/8*[1 0 0 0 0;
                0 1 0 0 0;
                0 0 1 0 0;
                0 0 0 0 0;
                0 0 0 0 0;];
            em1 = blkdiag(m, m, m, m, m, m, m, m);
    end  
end
