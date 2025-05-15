function [eta,y,omega1,sdof2] = HarmonicRespt(kk,mm,omega0,a,b,indices1)
%% Solve the eigenvalue problem
T = linspace(0, 5, 501);
F = 1e6;
fd = zeros(size(kk,1),1);
fd(indices1,1) = F;
T = T';
nm=50;
q0 = 0;
dq0 = 0;

sigma = 1e-3;         % 指定一个小的非零移位值
opts = struct('disp', 0);
[sdof,n1] = size(kk);
[nstep,n2] = size(T);
[V,D] = eigs(kk, mm, nm, sigma, opts);
[lambda,ki] = sort(diag(D));
omega = sqrt(lambda);
omega1 = sqrt(lambda)/(2*pi);

V=V(:,ki);
%% Check the eigenvalues
jk=0;

for ii = 1:sdof
    
    check=omega(ii);
    if check>1.0e12
        jk=jk+1;
        omi(jk)=ii;
    end
end

sdof2 = sdof-jk;
V1 = V(:,1:sdof2);
%% Normalize the eigenvectors and compute parameters
Factor = diag(V1'*mm*V1);
Vnorm = V1*inv(sqrt(diag(Factor)));
omega2 = diag(sqrt(Vnorm'*kk*Vnorm));

Fnorm = Vnorm'*fd;
%% compuete madal damping matrix from the proportional damping matrix
Modamp=Vnorm'*(a*mm+b*kk)*Vnorm;
zeta=diag((1/2)*Modamp*inv(diag(omega2)));

if( max(zeta) >= 1)
    disp('Warning -Your maximum damping ratio is grater than or equal to 1')
    disp('You have to reselect a and b')
    pause
    disp('If you want to continue, type return key')
end
%% Find out harmonic response of each modal coordinate analytically
eta0 = Vnorm'*mm*q0;
deta0 = Vnorm'*mm*dq0;
    %initial conditions for modal coordinates both displacement and velocity
eta=zeros(nstep,sdof2);
phase0=omega0*T;
    
for i=1:sdof2 	%responses are obtained for n modes
    gama = omega0/omega(i);
    omegad = omega(i)*sqrt(1-zeta(i)^2);
    phase = omegad*T;
    Exx = exp(-zeta(i)*omega(i)*T);
    C1 = eta0(i);
    C2 = (deta0(i)+eta0(i)*zeta(i)*omega(i))/omegad;

    X0 = sqrt((1-gama^2)^2+(2*zeta(i)*gama)^2);
    XX = Fnorm(i)/(omega(i)^2*X0);
    XP = atan((2*zeta(i)*gama)/(1-gama^2));
    D1 = (zeta(i)*omega(i)*cos(XP)+omega0*sin(XP))/omegad;
    D2 = cos(XP);

    eta(:,i) = C1*Exx.*cos(phase)+C2*Exx.*sin(phase)...
                -XX*Exx.*(D1*sin(phase)+D2*cos(phase))...
                +XX*cos(phase0-XP);
end
%% Convert modal coordinate responses to physical coordinate response
eta = eta';
y = Vnorm*eta;

if (a+b) == 0
    disp('The response results of the undamping system')
else
    disp('The response results of damping system')
end

disp('The excitation is harmonic force')