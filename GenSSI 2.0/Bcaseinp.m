    % Muñoz-Tamayo, R., de Groot, J., Bakx, E., Wierenga, P.A., Gruppen, H., Zwietering, M.H., Sijtsma, L., 2011. 
    % Hydrolysis of B-casein by the cell-envelope-located PI-type protease of Lactococcus lactis: A modelling approach.
    % Int. Dairy J. 21. https://doi.org/10.1016/j.idairyj.2011.03.012
    

function model = Baseinp()


    % Symbolic variables
     syms x1 
     syms b1 b2
        
     x0 = 10; % Initial concentrationg of B-casein
      
    
    % Parameters
     model.sym.p = [b1; b2];
 
     % State variables
    model.sym.x = [x1];
    
    % Functions 
     
     roh   = b1*x1/(b2-x1); % parameterization
               
     
    % Control vectors (g)
    model.sym.g = [];
    
    % Autonomous dynamics (f)
                 
                         
       model.sym.xdot=[ -roh];

    % Initial conditions
     model.sym.x0 = [x0]; 
    
     % Observables
     model.sym.y = [x1];
end
