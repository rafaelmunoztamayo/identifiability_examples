% Muñoz-Tamayo, R., S. Giger-Reverdin, and D. Sauvant. 2016. 
% Mechanistic modelling of in vitro fermentation and methane production by rumen microbiota. 
% Anim. Feed Sci. Technol. 220:1–21. doi:10.1016/j.anifeedsci.2016.07.005.

clear;

% states
syms z_ndf z_nsc z_pro s_su s_aa s_ac s_bu s_pr s_IN s_IC s_h2 s_ch4 x_su x_aa x_h2 ng_co2 ng_h2 ng_ch4 

x = [z_ndf; z_nsc; z_pro; s_su; s_aa; s_ac; s_bu; s_pr; s_IN; s_IC; s_h2; s_ch4; x_su; x_aa; x_h2; ng_co2; ng_h2; ng_ch4];

% outputs
h = [z_ndf; z_nsc; z_pro; s_su; s_aa; s_ac; s_bu; s_pr;s_IN; ng_co2; ng_h2; ng_ch4 ]; 


% no input:
u    = [];


% Unknown parameters 


syms  khyd_ndf khyd_nsc khyd_pro km_su Ks_su Ysu km_aa Ks_aa Yaa km_h2 Ks_h2 Yh2  lambda_1 lambda_2 
     
p = [khyd_ndf; khyd_nsc; khyd_pro; km_su; Ks_su; Ysu; km_aa; Ks_aa; Yaa; km_h2; Ks_h2; Yh2;  lambda_1; lambda_2];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Physicochemical parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Trumen = 312.15;    % Temperature, K
V_l    = 0.060;     % Volume in the liquid phase, L 
V_g    = 0.265;     % Volume in the gas phase, L 
Ptot   = 1.01325;   % Pressure, bar. 
R = 8.314*1e-2;     % Ideal gas constant, bar*L/(mol*K)
pgas_h2o = 0.08274; % Partial pressure water vapour, bar
kLa      = 8.3333;      % Liquid–gas transfer constant, 1/h

 deltaH0_KH_co2 = -19410; % J
 deltaH0_KH_ch4 = -14240; % J
 deltaH0_KH_h2 = -4180;   % J


% Henrys constant  M/bar, 

% at T = 25C (298.15K)

KH_co2_s = 0.035;
KH_ch4_s = 0.0014;
KH_h2_s =  7.8e-4; 

KH_co2 = KH_co2_s*exp(-(19410/(R*100))*(1/298.15-1/Trumen));
KH_ch4 = KH_ch4_s*exp(-(14240/(R*100))*(1/298.15-1/Trumen));
KH_h2  = KH_h2_s*exp(-(4180/(R*100))*(1/298.15-1/Trumen));
 

% Equilibrium constants 
deltaH0_Ka_w = 55900;   
deltaH0_Ka_co2 = 7646;
deltaH0_Ka_nh4 = 51965;

% Acid-base constants mol/L

K_w = exp(deltaH0_Ka_w/(R*100)*(1/298.15-1/Trumen))*1e-14;
K_a_ac = 10^(-4.76);
K_a_bu = 10^(-4.82);
K_a_pr = 10^(-4.88);
K_a_vfa = K_a_ac;
K_a_co2 = 10^(-6.35)*exp(deltaH0_Ka_co2/(R*100)*(1/298.15-1/Trumen));
K_a_nh4 =  10^(-9.25)*exp(deltaH0_Ka_nh4/(R*100)*(1/298.15-1/Trumen));

% Molecular weights (g/mol)

w_h2 = 2;
w_nh3 = 17;   
w_vfa = 60.05; 
w_ch4 = 16.0425;
w_co2 = 44.0095; 
w_mb=   113; 

N_mb = 1;
w_su = 180.16;   
w_aa = 134;    

N_aa = 2;

w_ac = 60.05; 
w_bu = 88.1051; 
w_pr = 74.08;


     

    pH = 6.6;         % Average pH; 
    ionH = 10^(-pH); 
    s_co2 = s_IC - K_a_co2*s_IC/(K_a_co2+ionH);  % Concentration of carbon dioxide in the liquid phase, mol/L
    
    fac_aa   = 0.67;        % Stoichiometric coefficient of acetate production from amino acids fermentation 
    fpr_aa   = 0.062;       % Stoichiometric coefficient of propionate production from amino acids fermentation 
    fbu_aa   = 0.24;        % Stoichiometric coefficient of butyrate production from amino acids fermentation  
    fh2_aa   = 0.82;        % Stoichiometric coefficient of hydrogen production from amino acids fermentation 
    fIC_aa   = 0.88;        % Stoichiometric coefficient of inorganic carbon production from amino acids fermentation 
    fch_x    = 0.20;        % Fraction of carbohydrates of the biomass, g/g
    fpro_x   = 0.55;        % Fraction of proteins of the biomass, g/g 

    K_S_IN = 2.0e-4;   % Nitrogen limitation constant, mol/L 
    k_d = 8.33e-04;    % Microbial decay rate constant, 1/h
  
    
% Liquid-gas transfer phenomena

  pgas_co2 =  R*Trumen*ng_co2/V_g;  % Partial pressure of CO2, bars
  pgas_ch4 =  R*Trumen*ng_ch4/V_g;  % Partial pressure of CH4, bars
  pgas_h2  =  R*Trumen*ng_h2/V_g;   % Partial pressure of H2, bars
  
  rhoT_co2 = kLa*(s_co2 - KH_co2*R*Trumen*ng_co2/V_g);   % Liquid–gas transfer rate of CO2, mol/(L h)
  rhoT_ch4 = kLa*(s_ch4 - KH_ch4*R*Trumen*ng_ch4/V_g);   % Liquid–gas transfer rate of CH4, mol/(L h)
  rhoT_h2  = kLa*(s_h2  - KH_h2*R*Trumen*ng_h2/V_g);     % Liquid–gas transfer rate of H2, mol/(L h)

% Microbial fermentation 

     % First order kinetics of hydrolysis
     
     rho_ndf  = khyd_ndf*z_ndf; % Hydrolysis rate of NDF crabohydrates, g/(L h)
     rho_nsc  = khyd_nsc*z_nsc; % Hydrolysis rate of NSC crabohydrates, g/(L h)
     rho_pro  = khyd_pro*z_pro; % Hydrolysis rate of proteins, g/(L h)
         
     
     % Stoichiometry of the rumen fermentation 

% Glucose utilization  

%  C6H12O6 + 2H2O  -> 2CH3COOH +2CO2 + 4H2                   (R1)
%  3C6H12O6        -> 2CH3COOH  + 4CH3CH2COOH + 2CO2 + 2H2O  (R2)
%  C6H12O6         -> CH3CH2CH2COOH + 2CO2 + 2H2             (R3)
%  5C6H12O6 + 6NH3 -> 6C5H7O2N + 18H2O                       (R4)

fsu = 1-(5/6)*Ysu;                                         % Fraction of glucose utilized in catabolism%
lambda_3  = 1- (lambda_1 + lambda_2);                      % Molar fraction of the sugars utilized via reaction 3

Yac_su  = fsu*(2*lambda_1 + (2/3)*lambda_2);               % Yield factor of acetate during sugars utilization, mol/mol                        
Ypr_su  = fsu*((4/3)*lambda_2);                            % Yield factor of propionate during sugars utilization, mol/mol
Ybu_su  = fsu*(1*lambda_3);                                % Yield factor of butyrate during sugars utilization, mol/mol
Yh2_su  = fsu*(4*lambda_1 + 2*lambda_3);                   % Yield factor of hydrogen during sugars utilization, mol/mol
YIC_su  = fsu*(2*lambda_1 + (2/3)*lambda_2 + 2*lambda_3);  % Yield factor of inorganic carbon during sugars utilization, mol/mol
YIN_su = -Ysu;                                             % Yield factor of inorganic nitrogen during sugars utilization, mol/mol

% Amino acids utilization

Yac_aa  = (1-Yaa)*fac_aa; % Yield factor of acetate during amino acids utilization, mol/mol
Ypr_aa  = (1-Yaa)*fpr_aa; % Yield factor of propionate during amino acids utilization, mol/mol
Ybu_aa  = (1-Yaa)*fbu_aa; % Yield factor of butyrate during amino acids utilization, mol/mol
Yh2_aa  = (1-Yaa)*fh2_aa; % Yield factor of hydrogen during amino acids utilization, mol/mol
YIC_aa = (1-Yaa)*fIC_aa;  % Yield factor of inorganic carbon during amino acids utilization, mol/mol
YIN_aa =  N_aa -Yaa*N_mb; % Yield factor of inorganic nitrogen during amino acids utilization, mol/mol


% H2 utilization 

% 4H2 + CO2         -> CH4 + 2H2O      (R5)
% 10H2 + 5CO2 + NH3 -> C5H7O2N + 8H2O  (R6)

fh2 = 1-10*Yh2; % Fraction of H2 utilized in catabolism

Ych4_h2 = fh2*(1/4);                  % Yield factor of methane during hydrogen utilization, mol/mol
YIC_h2  = -((1/4)*fh2 + 5*Yh2);       % Yield factor of inorganic carbon during hydrogen utilization, mol/mol
YIN_h2  = -Yh2;                       % Yield factor of inorganic nitrogen during hydrogen utilization, mol/mol

% Microbial kinetic rates
     
     % Nitrogen limitation 
     I_IN_lim = 1/(1 + K_S_IN/s_IN);  % Nitrogen limitation factor 
     
     rho_su   = km_su*s_su*x_su*I_IN_lim/(Ks_su + s_su);   % Utilization rate of sugars with H2 regulation, mol/(L h)
     rho_aa   = km_aa*s_aa*x_aa/(Ks_aa + s_aa);                 % Utilization rate of amino acids, mol/(L h)
     rho_h2   = km_h2*s_h2*x_h2*I_IN_lim/(Ks_h2 + s_h2);   % Utilization rate of hydrogen with bromoform inhibition, mol/(L h) 
     rho_xsu  = k_d*x_su;                                       % Cell death rate of sugars utilizers, mol/(L h)
     rho_xaa  = k_d*x_aa;                                       % Cell death rate of amino acids utilizers, mol/(L h)
     rho_xh2  = k_d*x_h2;                                       % Cell death rate of hydrogen utilizers, mol/(L h)
    



% dynamic equations
f = [ -rho_ndf;  
      -rho_nsc + (fch_x*w_mb)*(rho_xsu + rho_xaa + rho_xh2);
      -rho_pro + (fpro_x*w_mb)*(rho_xsu + rho_xaa + rho_xh2); 
      (rho_ndf/w_su + rho_nsc/w_su  - rho_su); 
      (rho_pro/w_aa - rho_aa); 
      (Yac_su*rho_su + Yac_aa*rho_aa); 
      (Ybu_su*rho_su + Ybu_aa*rho_aa); 
      (Ypr_su*rho_su + Ypr_aa*rho_aa); 
      (YIN_su*rho_su + YIN_aa*rho_aa + YIN_h2*rho_h2); 
      (YIC_su*rho_su + YIC_aa*rho_aa + YIC_h2*rho_h2 - rhoT_co2); 
      (Yh2_su*rho_su + Yh2_aa*rho_aa - rho_h2 - rhoT_h2); 
      (Ych4_h2*rho_h2 - rhoT_ch4); 
      (Ysu*rho_su - rho_xsu);  
      (Yaa*rho_aa - rho_xaa);  
      (Yh2*rho_h2 - rho_xh2); 
      (V_l*rhoT_co2); 
      (V_l*rhoT_h2); 
      (V_l*rhoT_ch4)];


% initial conditions:    
ics  = [    2.4920    3.4939    1.2933    0.0007         0    0.0263    0.0051    0.0074    0.0075    0.1404    0.0000    0.0007    0.0910    0.0048    0.0010 0 0 0];

% which initial conditions are known:
known_ics = ones(1,18); %[1 1 1 1 1];

save('rumen','x','h','u','p','f','ics','known_ics');


