% Muñoz-Tamayo, R., M. Popova, M. Tillier, D. P. Morgavi, J. P. Morel, G. Fonty, and N. Morel-Desrosiers. 2019. 
% Hydrogenotrophic methanogens of the mammalian gut: Functionally similar, thermodynamically different—A modelling approach. 
% PLoS One. 14:e0226243. doi:10.1371/journal.pone.0226243.

function model = methanogenesis()
    % Mathematical model of methanogenesis 
    
    % Reference article PLOS ONE 
    %    

    % Symbolic variables
     syms x_h2 ng_h2 ng_co2 ng_ch4 s_co2
  
%x_h2 = X(1);   % Biomass concentration, mol/L
%ng_h2 = X(2);  % moles of H2 in gas phase
%ng_co2 = X(3); % moles of CO2 in gas phase
%ng_ch4 = X(4); % moles of CH4 in gas phase
%s_co2 = X(5);  % CO2 concentration in liquid phase, mol/L

    % Known parameters
  
     syms kLa KH_co2  R T V_l  V_g 

  
    % Unknown Parameters 
    
      syms muMax Y Ks kd Ych4 Yco2 
      model.sym.p = [muMax;Y;Ks;kd;Ych4;Yco2];
     
     %syms muMax Y Ks kd 
     %model.sym.p = [muMax;Y;Ks;kd];
  
    
    % Dependent parameters
    
    km_h2 = muMax/Y;           % maximum specific substrate rate constant (mol H2/(mol biomass*h))
    
    % using knowledge on the stoichiometry
    
    %fh2 = 1-10*Y;              % fraction of H2 utilized in catabolism
    %Ych4 = fh2*(1/4);          % Yield factor of methane production (mol CH4/molH2)
    %Yco2 = ((1/4)*fh2 + 5*Y);  % Yield factor of CO2 consumption (mol CO2/molH2). 

    % State variables
    model.sym.x = [x_h2;ng_h2;ng_co2;ng_ch4;s_co2];
    
    % Functions 
     roh   = km_h2*exp(-Ks*V_g/(ng_h2))*x_h2*V_l;    % Kinetic rate of H2 utilization                
     rho_xh2 = kd*x_h2;                              % Microbial decay rate
     rhoT_co2 = kLa*(s_co2 - KH_co2*R*T*ng_co2/V_g); % Liquid-gas transfer rate

    % Control vectors (g)
    model.sym.g = [];
    
 
    % ODE model
    
       model.sym.xdot=[ Y*roh/V_l - rho_xh2
                        -roh 
                        V_l*rhoT_co2 
                        Ych4*roh 
                        -Yco2*roh/V_l - rhoT_co2];

    % Initial conditions
    
    model.sym.x0 = [ 0.0012;    0.0014;    0.0008;    0;  0.0268]; 
     
    % Observables
     model.sym.y = [x_h2;ng_h2];
     % model.sym.y = [x_h2;ng_h2;ng_co2;ng_ch4];
end
