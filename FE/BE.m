clear
clc

%%
no=1;
[msh]=qcl(no);  
nodes = msh.POS(:,:)./1000;   
elements = msh.QUADS(:,1:4);
XQ = zeros(4,1);
YQ = zeros(4,1);
ZQ = zeros(4,1);
GG = zeros(4,1);
A_sys = zeros(size(nodes,1));
H_sys = zeros(size(nodes,1));
G_sys = zeros(size(nodes,1));
CP_sys = zeros(size(nodes,1));
FREQ = 5.45901;
IRHOW = 1i*1000*6.283185*FREQ;

for i =1:size(nodes,1)
   XP = nodes(i,1);
   YP = nodes(i,2);
   ZP = nodes(i,3);
   
   for j = 1:size(elements,1)
        for kk = 1:4
            XQ(kk,1) = nodes(elements(j,kk),1);
            YQ(kk,1) = nodes(elements(j,kk),2);
            ZQ(kk,1) = nodes(elements(j,kk),3);
        end   
        
            [H,G,CP,AA,QN] = elem4(XP,YP,ZP,XQ,YQ,ZQ);
            CP_sys(i,i) = CP_sys(i,i)+CP;
            
            for l = 1:4
                GG(l,1) = -IRHOW*G(l,1)/(1i*6.283185*FREQ);
            end
            
            for mm = 1:4
               %if i==1
                    A_sys(elements(j,mm),elements(j,mm)) = A_sys(elements(j,mm),elements(j,mm))+AA(mm,1);
               %end   
               H_sys(i,elements(j,mm)) = H_sys(i,elements(j,mm))+H(mm,1);
               G_sys(i,elements(j,mm)) = G_sys(i,elements(j,mm))+GG(mm,1);          
            end
            
   end
end

Ma = A_sys*real(inv(H_sys+CP_sys))*real(G_sys);
