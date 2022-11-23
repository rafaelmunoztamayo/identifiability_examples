% OED problem of a model of microbial growth, Muñoz-Tamayo 2022

function [Fv] = moed(t,X,param)


k    = 2; 
kI   = 50;
D = 0.1; 

tt = [0 6 12 24]; % time of piece wise linear function
u = interp1(tt,param(1:4),t);

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
