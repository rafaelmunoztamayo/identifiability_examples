# Muñoz-Tamayo, R., S. Giger-Reverdin, and D. Sauvant. 2016. 
# Mechanistic modelling of in vitro fermentation and methane production by rumen microbiota. 
# Anim. Feed Sci. Technol. 220:1–21. doi:10.1016/j.anifeedsci.2016.07.005.


# parameters for identifiability testing: khyd_ndf khyd_nsc khyd_pro km_su Ks_su Ysu km_aa Ks_aa Yaa km_h2 Ks_h2 Yh2  lambda_1 lambda_2 

using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)


ode = @ODEmodel(
z_ndf'(t) =   -khyd_ndf*z_ndf(t),
z_nsc'(t) =   fch_x*w_mb*(k_d*x_aa(t) + k_d*x_h2(t) + k_d*x_su(t)) - khyd_nsc*z_nsc(t),
z_pro'(t) =   fpro_x*w_mb*(k_d*x_aa(t) + k_d*x_h2(t) + k_d*x_su(t)) - khyd_pro*z_pro(t),
s_su'(t)  =   (khyd_ndf*z_ndf(t))/w_su + (khyd_nsc*z_nsc(t))/w_su - (km_su*s_su(t)*x_su(t))/((Ks_su + s_su(t))*(K_s_IN/s_IN(t) + 1)),     
s_aa'(t)  =   (khyd_pro*z_pro(t))/w_aa - (km_aa*s_aa(t)*x_aa(t))/(Ks_aa + s_aa(t)),
s_ac'(t) =  - (fac_aa*km_aa*s_aa(t)*x_aa(t)*(Yaa - 1))/(Ks_aa + s_aa(t)) - (km_su*s_su(t)*x_su(t)*((5*Ysu)/6 - 1)*(2*lambda_1 + (2*lambda_2)/3))/((Ks_su + s_su(t))*(K_s_IN/s_IN(t) + 1)),
s_bu'(t) =  (km_su*s_su(t)*x_su(t)*((5*Ysu)/6 - 1)*(lambda_1 + lambda_2 - 1))/((Ks_su + s_su(t))*(K_s_IN/s_IN(t) + 1)) - (fbu_aa*km_aa*s_aa(t)*x_aa(t)*(Yaa - 1))/(Ks_aa + s_aa(t)),  
s_pr'(t) = 	- (fpr_aa*km_aa*s_aa(t)*x_aa(t)*(Yaa - 1))/(Ks_aa + s_aa(t)) - (4*km_su*lambda_2*s_su(t)*x_su(t)*((5*Ysu)/6 - 1))/(3*(Ks_su + s_su(t))*(K_s_IN/s_IN(t) + 1)), 
s_IN'(t) =  (km_aa*s_aa(t)*x_aa(t)*(N_aa - N_mb*Yaa))/(Ks_aa + s_aa(t)) - (Yh2*km_h2*s_h2(t)*x_h2(t))/((Ks_h2 + s_h2(t))*(K_s_IN/s_IN(t) + 1)) - (Ysu*km_su*s_su(t)*x_su(t))/((Ks_su + s_su(t))*(K_s_IN/s_IN(t) + 1)),  
s_IC'(t) = kLa*((K_a_co2*s_IC(t))/(K_a_co2 + ionH) - s_IC(t) + (KH_co2*R*Trumen*ng_co2(t))/V_g) - (fIC_aa*km_aa*s_aa(t)*x_aa(t)*(Yaa - 1))/(Ks_aa + s_aa(t)) - (km_h2*s_h2(t)*x_h2(t)*((5*Yh2)/2 + 1/4))/((Ks_h2 + s_h2(t))*(K_s_IN/s_IN(t) + 1)) + (km_su*s_su(t)*x_su(t)*((5*Ysu)/6 - 1)*((4*lambda_2)/3 - 2))/((Ks_su + s_su(t))*(K_s_IN/s_IN(t) + 1)),   
s_h2'(t) = - kLa*(s_h2(t) - (KH_h2*R*Trumen*ng_h2(t))/V_g) - (km_h2*s_h2(t)*x_h2(t))/((Ks_h2 + s_h2(t))*(K_s_IN/s_IN(t) + 1)) - (fh2_aa*km_aa*s_aa(t)*x_aa(t)*(Yaa - 1))/(Ks_aa + s_aa(t)) - (km_su*s_su(t)*x_su(t)*((5*Ysu)/6 - 1)*(2*lambda_1 - 2*lambda_2 + 2))/((Ks_su + s_su(t))*(K_s_IN/s_IN(t) + 1)),
s_ch4'(t) = - kLa*(s_ch4(t) - (KH_ch4*R*Trumen*ng_ch4(t))/V_g) - (km_h2*s_h2(t)*x_h2(t)*((5*Yh2)/2 - 1/4))/((Ks_h2 + s_h2(t))*(K_s_IN/s_IN(t) + 1)),
x_su'(t) = (Ysu*km_su*s_su(t)*x_su(t))/((Ks_su + s_su(t))*(K_s_IN/s_IN(t) + 1)) - k_d*x_su(t),
x_aa'(t) = (Yaa*km_aa*s_aa(t)*x_aa(t))/(Ks_aa + s_aa(t)) - k_d*x_aa(t),
x_h2'(t) = (Yh2*km_h2*s_h2(t)*x_h2(t))/((Ks_h2 + s_h2(t))*(K_s_IN/s_IN(t) + 1)) - k_d*x_h2(t),
ng_co2'(t) = -V_l*kLa*((K_a_co2*s_IC(t))/(K_a_co2 + ionH) - s_IC(t) + (KH_co2*R*Trumen*ng_co2(t))/V_g),
ng_h2'(t) = V_l*kLa*(s_h2(t) - (KH_h2*R*Trumen*ng_h2(t))/V_g),
ng_ch4'(t) = V_l*kLa*(s_ch4(t) - (KH_ch4*R*Trumen*ng_ch4(t))/V_g),
   y1(t)   = z_ndf(t),
   y2(t)   = z_nsc(t),
   y3(t)   = z_pro(t),
   y4(t)   = s_su(t),
   y5(t)   = s_aa(t),
   y6(t)   = s_ac(t),
   y7(t)   = s_bu(t),
   y8(t)   = s_pr(t),
   y9(t)   = s_IN(t),
   y10(t)  = ng_co2(t),
   y11(t)  = ng_h2(t),
   y12(t)  = ng_ch4(t)   
   );
   
QQ = StructuralIdentifiability.Nemo.QQ

ode = set_parameter_values(ode, Dict(
    ionH => QQ(1186204829293589,4722366482869645213696),
	Trumen => QQ(6243,20),
	R => QQ(4157,50000),
	V_l  => QQ(3,5),
	V_g  => QQ(53,200),
	kLa => QQ(83333,10000),
	KH_co2 => QQ(1775116163948307,72057594037927936),
	KH_ch4 => QQ(311868769628171,288230376151711744),
	KH_h2 => QQ(6670194274691929,9223372036854775808),
	K_a_co2  => QQ(2422372177661713,4722366482869645213696),
 	fac_aa => QQ(67,100),
	 fpr_aa=> QQ(62,1000),
	fbu_aa=> QQ(24,100),
	 fh2_aa=> QQ(82,100),
	fIC_aa => QQ(88,100),
	 fch_x => QQ(20,100),
	fpro_x  => QQ(55,100),
	K_s_IN => QQ(1,5000), 
    k_d => QQ(3843056309736095,4611686018427387904),
	 w_mb => QQ(113),
	 w_su=> QQ(4504,25),
	 w_aa => QQ(134),
	 N_aa => QQ(2),
	 N_mb  => QQ(1) 
	 ))


assess_local_identifiability(ode,[khyd_ndf,khyd_nsc,khyd_pro,km_su,Ks_su,Ysu,km_aa,Ks_aa,Yaa,km_h2,Ks_h2,Yh2,lambda_1,lambda_2 ])
#assess_global_identifiability(ode,[khyd_ndf,khyd_nsc,khyd_pro,km_su,Ks_su,Ysu,km_aa,Ks_aa,Yaa,km_h2,Ks_h2,Yh2,lambda_1,lambda_2])
#assess_local_identifiability(ode)


 