clc
clear
% 第二个子结构自由度少，因为底部加了一圈固定约束

nn=40;

% [KB,MB,T2]=CylindricalB(nn);
[KA,MA,T1]=ConeA(nn);
[KB,MB,T2]=CylindricalBB(nn);

% 400 = 5*80 自由度*固定界面节点数

KK=zeros((400+nn)*2,(400+nn)*2);
MM=zeros((400+nn)*2,(400+nn)*2);
KK(1:400+nn,1:1:400+nn)=KA;KK(401+nn:(400+nn)*2,401+nn:(400+nn)*2)=KB;
MM(1:400+nn,1:1:400+nn)=MA;MM(401+nn:(400+nn)*2,401+nn:(400+nn)*2)=MB;
T3=zeros((400+nn)*2,400+nn*2);
T3(1:nn,1:nn)=eye(nn);
T3(nn+1:400+nn,2*nn+1:400+2*nn)=eye(400);
T3(nn+401:800+nn,2*nn+1:400+2*nn)=eye(400);
T3(800+nn+1:800+2*nn,nn+1:2*nn)=eye(nn);
K=T3'*KK*T3;
M=T3'*MM*T3;

nm=20;
[v,d] =eigs(K,M,nm,'SM');  %v特征向量 d特征值
tempd=diag(d);
[d,sortindex]=sort(tempd);
omega=sqrt(d);
v=v(:,sortindex);
mode_number=1:nm;
frequency(mode_number)=sqrt(d(mode_number))/(2*pi);
Frequency=real(frequency);
Frequency=Frequency';
V=real(v);

% T=zeros(800+2*nn,2*nn+400);
% T(1:nn,1:nn)=eye;
% T(nn+1:nn+400,2*nn+1:2*nn+400)=eye;
% T(nn+400+1:2*nn+400,2*nn+1:2*nn+400)=eye;
% T(800+nn+1:2*nn+800,nn+1:2*nn)=eye;
% 
% k=[KA,zeros(420,nn);zeros(420,nn),KB];
% m=[MA,zeros(420,nn);zeros(420,nn),MB];
% 
% K=T'*k;
% M=T'*m;

% mija=MA(1:nn,nn+1:nn+400);
% mjja=MA(nn+1:nn+400,nn+1:nn+400);
% kiia=KA(1:nn,1:nn);
% kjja=KA(nn+1:nn+400,nn+1:nn+400);
% 
% mijb=MB(401:nn+400,1:400);
% mjjb=MB(1:400,1:400);
% kiib=KB(401:nn+400,401:nn+400);
% kjjb=KB(1:400,1:400);
% 
% M=zeros(2*nn+400);
% M(1:nn,1:nn)=eye;
% M(nn+1:2*nn,nn+1:2*nn)=eye;
% M(1:nn,2*nn+1:2*nn+400)=mija;
% M(nn+1:2*nn,2*nn+1:2*nn+400)=mijb;
% M(2*nn+1:2*nn+400,1:nn)=mija';
% M(2*nn+1:2*nn+400,nn+1:2*nn)=mijb';
% M(2*nn+1:2*nn+400,2*nn+1:2*nn+400)=mjja+mjjb;
% 
% K=zeros(2*nn+400);
% K(1:nn,1:nn)=kiia;
% K(nn+1:2*nn,nn+1:2*nn)=kiib;
% K(2*nn+1:2*nn+400,2*nn+1:2*nn+400)=kjja+kjjb;








