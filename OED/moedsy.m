% OED problem of a model of microbial growth, Muñoz-Tamayo 2022

function [sy1, sy2] = moedsy(X,sx1,sx2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x1 = X(1);
x2 = X(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative of H with respect to the parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	dHdpv(1,1)=0;
	dHdpv(1,2)=0;
	dHdpv(2,1)=0;
	dHdpv(2,2)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivative of H with respect to the state  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	dHdxv(1,1)=1;
	dHdxv(1,2)=0;
	dHdxv(2,1)=0;
	dHdxv(2,2)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity equations of the output 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		 sy1 = 		 dHdxv*sx1 + dHdpv(:,1); 
		 sy2 = 		 dHdxv*sx2 + dHdpv(:,2); 
