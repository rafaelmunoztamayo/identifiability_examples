 % Mu?oz-Tamayo, R., de Groot, J., Bakx, E., Wierenga, P.A., Gruppen, H., Zwietering, M.H., Sijtsma, L., 2011. 
 % Hydrolysis of B-casein by the cell-envelope-located PI-type protease of Lactococcus lactis: A modelling approach.
 % Int. Dairy J. 21. https://doi.org/10.1016/j.idairyj.2011.03.012
    

% Step 1: create symbols for the dynamic states, unknown parameters,
% measured input(s) and unmeasured input(s) of the system
syms x1  p1 p2


% Step 2: write the system functions 
F=[-p1*x1/(p2 - x1)];

% Step 3: write the output functions
h=[x1];

% Step 4: write the vector of dynamic states
X=[x1];

% Step 5: write the vector of unknown parameters
Theta=[ p1 p2];

% Step 6: write the vector of measured inputs U=[u1,...,ur]; U=[], if no
% measured inputs are applied
U=[];

% Step 7: write the vector of unmeasured inputs W=[w1,...,wnw]; W=[], if no
% unmeasured inputs are applied
W=[];

% Step 8: choose the maximum order of the time derivative of unmeasured
% inputs considered; normally, it is recommended that kmax=2^i-1>n+l+nw 
% (i=1,2,3,..., n is the number of dynamic states, l is the number of
% parameters, nw is the number of unmeasured inputs); if no unmeasured
% inputs are applied, the default kmax is n+l-1.

kmax=4; %i is chosen as 4

% Step 9: input to the function

NonlinearObservabilityTest(F,h,X,Theta,U,W,kmax)