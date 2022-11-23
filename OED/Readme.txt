This folder contains the files to perform an optimal experiment design (OED) for a model of microbial growth.

moed.m: file with the ordinary differential equations (ODE) of the model 
moedse.m: file with ODE of the model and sensitivity equations 
moedsy.m: auxiliary file to compute the sensitivity equations 
moedFIM.m: file to calculate the Fisher Information Matrix (FIM)
moedcost.m: file with the cost function to optimize for the OED problem (maximization of the determinant of the FIM)
moedoptim: file to perform the optimization 
moedplot: file to plot the responses perform the optimization 
