% OED problem of a model of microbial growth, Muñoz-Tamayo 2022

function [J] = moedcost(paramesti)

st = 0:1.0:24; 

% Initial conditions of the state variables 

 Cinit = [  10  0  ];
 
 param = round(10*paramesti)/10;


Np= 2;  
Nx= 2; 

ceros=zeros(1,Np*Nx); 

Cin = [Cinit  ceros]; 
[ti,Ci] = ode45(@(ti,Ci) moedse(ti,Ci,param), st, Cin);

Xm=Ci(2:end,1:Nx); 
Ym = Xm'; 

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of sensitivities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivities of the state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	 s_1=(Ci(2:end,3:4))';
	 s_2=(Ci(2:end,5:6))';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivities of the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 clear k
 
for k=1:length(ti)-1
[dymds_1(:,k), dymds_2(:,k)] = moedsy(Xm(k,:),s_1(:,k),s_2(:,k));
end

[rowy,coly]=size(Ym);
nt=length(st);
     Sigma = 0.1; 
	 IS=inv(Sigma);
	 F = zeros(Np,Np); % FIM
	    for  i=1:coly 
	    dymdp = [	 dymds_1(:,i) 	 dymds_2(:,i) 	];
	    F  = dymdp'*IS*dymdp +F;
	   
	 end 
J = -det(F);  % determinant of the FIM 
