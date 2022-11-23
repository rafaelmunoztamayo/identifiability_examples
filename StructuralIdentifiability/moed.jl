    # Model microbial growth with inhibition 

using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)



ode = @ODEmodel(
   x1'(t)  = x1(t)*(x2(t)/(x2(t) + k + x2(t)^2/kI)) -D*x1(t),
   x2'(t)  = D*(u(t)-x2(t))-x1(t)*(x2(t)/(x2(t) + k + x2(t)^2/kI)),
   y1(t)   = x1(t),
   y2(t)   = x2(t)
   )
   
assess_identifiability(ode,[k,kI])
