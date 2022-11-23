 % Muñoz-Tamayo, R., de Groot, J., Bakx, E., Wierenga, P.A., Gruppen, H., Zwietering, M.H., Sijtsma, L., 2011. 
 % Hydrolysis of B-casein by the cell-envelope-located PI-type protease of Lactococcus lactis: A modelling approach.
 % Int. Dairy J. 21. https://doi.org/10.1016/j.idairyj.2011.03.012
    
clear;

% states
syms x1 
x = [x1];

% 1 output
h = x1;

% no input:
u    = [];

% parameters 
syms b1 b2 
p = [b1; b2];

% functions 
 x0 = 10; % Initial concentrationg of B-casein
 
 roh   = b1*x1/(b2 - x1); % Competitive inhibition kinetics

% dynamic equations
f = [-roh];


% initial conditions:    
ics  = [10];  

% which initial conditions are known:
known_ics = [1];

save('Bcaseinp','x','h','u','p','f','ics','known_ics');


