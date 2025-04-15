%--------------------------------------------------------------------------
function [point1,weight1]=GaussPoint1(No_points)
%--------------------------------------------------------------------------
%  Purpose:
%     determine the integration points and weighting coefficients of Gauss quadrature
%  Variable Description:
%     No_points - number of integration points
%     point1  - vector containing integration points
%     weight1 - vector containing weighting coefficients
%-------------------------------------------------------------------
%------------------------------------------
%  (1) initialization of the vectors
%------------------------------------------
   points=zeros(No_points,1);
   weights=zeros(No_points,1);
%----------------------------------------------------------
%  (2) find corresponding integration points and weights
%----------------------------------------------------------
if No_points==1                 % 1-point quadrature rule
    point1(1)=0.0;
    weight1(1)=2.0;
elseif No_points==2             % 2-point quadrature rule
    point1(1)=-0.577350269189626;
    point1(2)=-point1(1);
    
    weight1(1)=1.0;
    weight1(2)=weight1(1);
elseif No_points==3             % 3-point quadrature rule
    point1(1)=-0.774596669241483;
    point1(2)=0.0;
    point1(3)=-point1(1);
    
    weight1(1)=0.555555555555556;
    weight1(2)=0.888888888888889;
    weight1(3)=weight1(1);
elseif No_points==4             % 4-point quadrature rule
    point1(1)=-0.861136311594053;
    point1(2)=-0.339981043584856;
    point1(3)=-point1(2);
    point1(4)=-point1(1);
    
    weight1(1)=0.347854845137454;
    weight1(2)=0.652145154862546;
    weight1(3)=weight1(2);
    weight1(4)=weight1(1);
elseif No_points==5             % 5-point quadrature rule
    point1(1)=-0.906179845938664;
    point1(2)=-0.538469310105683;
    point1(3)=0.0;
    point1(4)=-point1(2);
    point1(5)=-point1(1);
    
    weight1(1)=0.236926885056189;
    weight1(2)=0.478628670499366;
    weight1(3)=0.568888888888889;
    weight1(4)=weight1(2);
    weight1(5)=weight1(1);
elseif  No_points==6            % 6-point quadrature rule
    point1(1)=-0.932469514203152;
    point1(2)=-0.661209386466265;
    point1(3)=-0.238619186083197;
    point1(4)=-point1(3);
    point1(5)=-point1(2);
    point1(6)=-point1(1);
    
    weight1(1)=0.171324492379170;
    weight1(2)=0.360761573048139;
    weight1(3)=0.467913934572691;
    weight1(4)=weight1(3);
    weight1(5)=weight1(2);
    weight1(6)=weight1(1);
else
   disp('number of integration points should be')
   disp('> 0 or < 7')
end
%-------------------------------------------------------------------
%    The end
%-------------------------------------------------------------------