    # Muñoz-Tamayo, R., de Groot, J., Bakx, E., Wierenga, P.A., Gruppen, H., Zwietering, M.H., Sijtsma, L., 2011. 
    # Hydrolysis of β-casein by the cell-envelope-located PI-type protease of Lactococcus lactis: A modelling approach.
    # Int. Dairy J. 21. https://doi.org/10.1016/j.idairyj.2011.03.012
    
using SIAN

ode = @ODEmodel(
   x1'(t)  = -k*x1(t)/(Km*(1+(10-x1(t))/kI)+x1(t)),
   y1(t)   = x1(t)
   );

output = identifiability_ode(ode, get_parameters(ode));