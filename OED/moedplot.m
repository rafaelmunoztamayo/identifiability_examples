% OED problem of a model of microbial growth, Muñoz-Tamayo 2022


clear all
st = 0:1.0:24; 

% Initial conditions of the state variables 

 moedparamPWL = [ 286.0000  600.9000  663.4000  836.5000]; % optimal solution with dynamic input as a piecewise linear function 
 tt = [0 6 12 24]; % time of piece wise linear function
 u1 = interp1(tt,moedparamPWL,st);
 moedparamC = 465.5*ones(1,4); % optimal solution with input as constant 
 u2 = interp1(tt,moedparamC,st);
 X0 = [  10  0  ]; % iniitial conditions 
 [t1,X1] = ode45(@(t1,X1) moed(t1,X1,moedparamPWL), st, X0);
 
 [F1] = moedFIM(moedparamPWL); % FIM
 dF1 = det(F1);                % determinant of the FIM
 P1 = inv(F1);                 % covariance matrix of the parameter estimates
 sP1 = P1.^0.5;                
 sD1 = diag(sP1);              % standand deviation of the parameter estimates
  
 [t2,X2] = ode45(@(t2,X2) moed(t2,X2,moedparamC), st, X0);
 [F2] = moedFIM(moedparamC);   % FIM 
 dF2 = det(F2);                % determinant of the FIM
 P2 = inv(F2);                 % covariance matrix of the parameter estimates
 sP2 = P2.^0.5;
 sD2 = diag(sP2);              % standand deviation of the parameter estimates

 
 Results = [2 50; sD1'; sD2'];
 

figure; 
subplot(131),plot(st,X1(:,1),'linewidth',1.5); 
hold on 
plot(st,X2(:,1),'--r','linewidth',1.5); 


set(gca,'Fontsize',18);
	 xlabel(' Time (h) ','fontsize',16);
	 ylabel(' x_1 ','fontsize',16);
subplot(132), plot(st,X1(:,2),'linewidth',1.5); 
hold on 
plot(st,X2(:,2),'--r','linewidth',1.5);  
set(gca,'Fontsize',18);
	 xlabel(' Time (h) ','fontsize',16);
	 ylabel(' x_2 ','fontsize',16);
 subplot(133), plot(st,u1,'linewidth',1.5); 
 set(gca,'Fontsize',18);
 hold on 
 plot(st,u2,'r--','linewidth',1.5); 
	 xlabel(' Time (h) ','fontsize',16);
	 ylabel(' u ','fontsize',16);
     %legend('PWL', 'Constant')
