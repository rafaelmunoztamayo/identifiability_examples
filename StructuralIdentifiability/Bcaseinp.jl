    # Muñoz-Tamayo, R., de Groot, J., Bakx, E., Wierenga, P.A., Gruppen, H., Zwietering, M.H., Sijtsma, L., 2011. 
    # Hydrolysis of β-casein by the cell-envelope-located PI-type protease of Lactococcus lactis: A modelling approach.
    # Int. Dairy J. 21. https://doi.org/10.1016/j.idairyj.2011.03.012

using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)



ode = @ODEmodel(
x1'(t)  = -b1*x1(t)/(b2-x1(t)),
   y1(t)   = x1(t)
   )
   
assess_identifiability(ode)