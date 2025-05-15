function [Ma]=add_mass(FREQ,nodes,elements)

XQ = zeros(8,1);
YQ = zeros(8,1);
ZQ = zeros(8,1);
GG = zeros(8,1);
A_sys = zeros(size(nodes,1));
H_sys = zeros(size(nodes,1));
G_sys = zeros(size(nodes,1));
CP_sys = zeros(size(nodes,1));
IRHOW = 1i*1000*6.283185*FREQ;

for i =1:size(nodes,1)
   XP = nodes(i,1);
   YP = nodes(i,2);
   ZP = nodes(i,3);
   CP = 1;
   for j = 1:size(elements,1)
        for kk = 1:8
            XQ(kk,1) = nodes(elements(j,kk),1);
            YQ(kk,1) = nodes(elements(j,kk),2);
            ZQ(kk,1) = nodes(elements(j,kk),3);
        end   
        
            [H,G,CP,AA,~] = elem8(XP,YP,ZP,XQ,YQ,ZQ,FREQ,CP);
            CP_sys(i,i) = CP_sys(i,i)+CP;
            
            GG(:,1) = -IRHOW*G(:,1)/(1i*6.283185*FREQ);
 
            
            for mm = 1:8
               A_sys(elements(j,mm),elements(j,mm)) = A_sys(elements(j,mm),elements(j,mm))+AA(mm,1);  
               H_sys(i,elements(j,mm)) = H_sys(i,elements(j,mm))+H(mm,1);
               G_sys(i,elements(j,mm)) = G_sys(i,elements(j,mm))+GG(mm,1);          
            end
            
   end
end
Ma = A_sys*real(inv(H_sys+CP_sys))*real(G_sys);
