% OED problem of a model of microbial growth, Muñoz-Tamayo 2022

function [Fv] = moedse(t,X,param)


k    = 2; 
kI   = 50;
D = 0.1;

%u = param(1); 
tt = [0 6 12 24]; % time of piece wise linear function
u = interp1(tt,param(1:4),t);

Np= 2;  
Nx= 2; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x1 = X(1);
x2 = X(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	 Fv(1,:)=(x1*x2)/(k + x2 + x2^2/kI) - D*x1;
	 Fv(2,:)=D*(u - x2) - (x1*x2)/(k + x2 + x2^2/kI);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative of F with respect to the parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	dFdpv(1,1)=-(x1*x2)/(k + x2 + x2^2/kI)^2;
	dFdpv(1,2)=(x1*x2^3)/(kI^2*(k + x2 + x2^2/kI)^2);
	dFdpv(2,1)=(x1*x2)/(k + x2 + x2^2/kI)^2;
	dFdpv(2,2)=-(x1*x2^3)/(kI^2*(k + x2 + x2^2/kI)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative of F with respect to the state 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	dFdxv(1,1)=x2/(k + x2 + x2^2/kI) - D;
	dFdxv(1,2)=x1/(k + x2 + x2^2/kI) - (x1*x2*((2*x2)/kI + 1))/(k + x2 + x2^2/kI)^2;
	dFdxv(2,1)=-x2/(k + x2 + x2^2/kI);
	dFdxv(2,2)=(x1*x2*((2*x2)/kI + 1))/(k + x2 + x2^2/kI)^2 - x1/(k + x2 + x2^2/kI) - D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	 j =Nx+1; 
	   for k=1:Np; 
		 Fv(j:j+Nx-1,:)=dFdxv*X(j:j+Nx-1)+dFdpv(:,k); 
		       j=j+Nx; 
	 end 
