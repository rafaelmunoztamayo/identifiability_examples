% OED problem of a model of microbial growth, Muñoz-Tamayo 2022


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Optimization 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all


	 Np= 2;  
	 Ny= 2; 

alt_optim= 1; 

p0 = 10; 

load moedparam.mat  
load moedparam0.mat

%paramMi = p0; 
pcoptim = 465.5; 
pwl = pcoptim*ones(1,4); % PWL 
%paramMi = pcoptim; 
%paramMi = pwl; 
paramMi = moedparam; 

	 % paramMi = log(paramMi); 
tic;
disp('                                       ');
	 % disp('   Type 1 i you are using normal scale, or 2 if it is in log scale : '); 
	 % olog = input('    '); 
	 olog = 1; 

clc ; 
	 if olog==1
	 LBp = paramMi*0.9;
	 UBp = paramMi*1.1;
	 else 
	 LBp = log(exp(paramMi*0.8));
	 UBp = log(exp(paramMi*1.2));
	 end 

	 N = 1; % Number of runs of the optimization algorithm 
	 options1 = optimset('LargeScale','off', 'Display','iter', 'Diagnostics','on','MaxIter',500);
	 options2 = optimset('LargeScale','on', 'GradObj','on','Hessian', 'on', 'Display','iter', 'Diagnostics','on');

	 if alt_optim ==1 
	 options = options1; 
	 else 
	 options = options2; 
	 end 
	 clear i; 
disp('                                       ');
disp('           *********************************** ');
disp('                  Finding the estimates');
disp('           *********************************** ');
	 for i = 1:N 
		  aleat = rand(1,Np);
 		  if i==1
 		  paramo(1,:) = paramMi; 
 		  else 
 		   for j=1:Np  
 		    paramo(i,j) = LBp(j)+aleat(j)*(UBp(j)-LBp(j));  
 		  end 
 		  end 
 	  [x(i,:),fval(i), exitflag(i), output] = fminsearch(@moedcost,paramo(i,:),options);
disp('                                       ');
	 fvalues(i) = fval(i);
	 estimates(i,:) = x(i,:);
		  end 
 
	 F = sort(fvalues);
	 minfval = min(fvalues); 
	 indexfval = find(fvalues==minfval); 
	 teta = estimates(indexfval,:);
	 Jc = minfval; 
	 teta_estimates = teta;
	 % teta_estimates = exp(teta); % For the case of change in the parameterization 
	 zero_teta = length(find(teta_estimates>=0)); 
	 J_cost = Jc; 
     teta_estimates = round(10*teta_estimates)/10;

     moedparam = teta_estimates; 
     moedparam0 = paramMi; 

    save moedparam.mat moedparam
    save moedparam0.mat moedparam0
