
# PARAMETERS

#time horizon 
param j default 0;
param n default 1;
param T:=300-j; 
param m:=100;
set f;
let f:={0..T+m};
set i;

#let i:={1..6};
#set o;
#let o:={1..6};
        

# time invariant parameters
param gamma; #:=0.4;  # capital elasticity in production function				#OK
param rho; #:=0.015; #time preference rate
param eta; #:=1.45; #intertemporal elasticity of marginal utility of consumption
param deltaK;
param Tbar;

#FaIR Parameters
param q1:= 0.33;
param q2:= 0.41;
param d1:= 239;
param d2:= 4.1;
param r0:= 35;
param rc:= 0.019;
param rT:= 4.165;
param F2x:=3.71;
param tau1:=10^6;
param tau2:=394.4;
param tau3:=36.54;
param tau4:=4.304;
param a1:=0.2173;
param a2:=0.2240;
param a3:=0.2824;
param a4:=0.2763;

# Geoffroy Parameters
param lambdag:= 1.18;
param upperc:= 8.2; 
param lowerc:= 109.0;
param betag:=0.67;
param effg:= 1.28; 
param dtg:= 1;
param cdeepp:=lowerc*effg;
param gammap:=betag*effg;
param g1:= (lambdag+gammap)/upperc;
param g2:= gammap/cdeepp;
param g = g1 + g2;
param gstar = g1 - g2;
param delsqrt = (g*g - 4*g2*lambdag/upperc)^0.5;
param afast = (g + delsqrt)/2;
param aslow = (g - delsqrt)/2;
param cc = 0.5/(upperc*delsqrt);
param amixf = cc*(gstar + delsqrt);
param amixs = -cc*(gstar - delsqrt);
param adeepf = -gammap/(upperc*cdeepp*delsqrt);
param adeeps = -adeepf;
param adf = 1/(afast*1);
param ads = 1/(aslow*1);
param expf = exp(-1/adf);
param exps = exp(-1/ads);

#Damage parameters
param Psi; #:=0; #these parameters can be left unspecified and specified in the .run file.
param Psi2;
param Psi3;
param Theta:=2.6; # Exponent of control cost function



#params defining initial values
param Qgross0>=0; # Initial world gross output (trill 2010 USD)
param K0>=0; 		#initial capital value (trillions 2005 USD)
param A0>=0; 		#initial level of total factor productivity
param I0>=0;		# initial level of Investment
param R10>=0; 		#:=57.59868223;
param R20>=0; 		#:=42.34107373;
param R30>=0; 		#:=17.11210069;
param R40>=0; 		#:=2.785685633;
param CCA0>=0;		#:=307.9811494;
param alpha0>=0; 	#:=0.491714646;
param dalpha0>=0; 	#:=19.6448;
param TMF0; 		#:=0.858580961;
param TMS0;		#:=0.105895355;
param TDF0;		#:=-0.021543931;
param TDS0;		#:=0.248033649;



#param Ecum0:=90; #Initial cumulative emissions (GtCO2)
    

# time variant parameters
#param R {t in 0..T}:=1/(1+rho)^(t+j); # use this specification for the restart solution
param R {t in 0..T}:=1/(1+rho)^(t);


# discount factor
#let R[0]:=1/(;
#let {t in 1..T} R[t] := R[t-1]/((1+rho));

# time dependent parameters
param year{f}>=0;
param muo{f};

# time and scenario dependent parameters

param L{f,i}>=0;
param A{f,i}>=0;
param sigma{f,i}>=0;
param pback{f,i}>=0;
param phead{t in 0..T}=pback[t+j,n]*sigma[t+j,n]/Theta/1000;
#param phead {f,i}=pback[f,n]*sigma[f,n]/Theta/1000;
param Fex{f,i};
param Eland{f,i};
#param deltaK{f,i}; #:=0.0344030175444404; #0.07852405; #depreciation rate on capital (per year)				#OK

# VARIABLES

# capital
var K {t in 0..T}>=10;

# aggregate consumption
var C {t in 0..T}>=0;

# Investment
var I {t in 0..T}>=0;

# per capita consumption
var c {t in 0..T}>=0.01 default 5;

# Gross output
var Qgross {t in 0..T}>=0;
# = A[t+j]*K[t]^gamma*(L[t])^(1-gamma); # gross output (trillion PPP dollars) 


# carbon stocks
var R1 {t in 0..T};
var R2 {t in 0..T};
var R3 {t in 0..T};
var R4 {t in 0..T};
var Catm{t in 0..T}=278+R1[t]+R2[t]+R3[t]+R4[t];
var CCA{t in 0..T};
var Catmstock{t in 0..T}=Catm[t]*2.12;

# total radiative forcing
#var F {t in 0..T}=kappa*((log(MAT[t]/MATEQ))/log(2))+Fex[t]; 
var F {t in 0..T}=(F2x/log(2))*log(Catm[t]/278)+Fex[t+j,n];

#276.60278
#281.9352715993539

var intf{t in 1..T} = (F[t-1]*adf + F[t]*(1 - adf) - expf*(F[t-1]*(1 + adf) - F[t]*adf))/afast;
var ints{t in 1..T} = (F[t-1]*ads + F[t]*(1 - ads) - exps*(F[t-1]*(1 + ads) - F[t]*ads))/aslow;

#feedback
var alpha {t in 0..T} >=0.01, <=9;
var dalpha{t in 0..T} default 25.035;


#  temperature
#var T1 {t in 0..T}>=0, <=10;
#var T2 {t in 0..T}>=0, <=10;
#var Temp {t in 0..T}=T1[t]+T2[t] default 0.61;

var TMF {t in 0..T}; #>=0, <=10;
var TMS {t in 0..T}; #>=0, <=10;
var TDF {t in 0..T}; #>=0, <=10;
var TDS {t in 0..T}; #>=0, <=10;
var Temp {t in 0..T}=TMF[t]+TMS[t]; #default 0.61;
var TempIm {t in 0..T}=max(0,Temp[t]-0.0485899); #note that the damage function is calibrated for an temp increase against average 1850-1899

# damage fraction
var Omega {t in 0..T}=Psi*(TempIm[t])^2+Psi2*(TempIm[t])^2+Psi3*(1-(1/((1+(TempIm[t]/20.5847)^2)+(TempIm[t]/6.081)^(6.754))));
#var Omega {t in 0..T}=1-(1/((1+(TempIm[t]/20.5847)^2)+(TempIm[t]/6.081)^(6.754)));


# damages 
var damage {t in 0..T}=Omega[t]*Qgross[t];

# output net of damages and abatement and CDR
var Q {t in 0..T};

# emission control
var mu {t in 0..T}>=0;


# abatement costs (fraction of output)
var Lambda {t in 0..T}=Qgross[t]*phead[t]*(mu[t]^Theta);

# industrial emissions
var EInd {t in 0..T}=sigma[t+j,n]*Qgross[t]*(1-mu[t]);


# total emissions
var E {t in 0..T}=EInd[t]+Eland[t+j,n]*3.664;
#var E {t in 0..T}=EInd[t];
var Ec {t in 0..T}=E[t]/3.664;
var Ecp{t in 0..T}=Ec[t]/2.1293970467500096;


# maximum cumulative extraction fossil fuels (GtC)
var Ecum {t in 0..T}<=6000;  

# Marginal cost of abatement (carbon price)
var cprice {t in 0..T}=pback[t+j,n]*mu[t]^(Theta-1);


# utility
var U {t in 0..T}=L[t+j,n]*(((c[t]^(1-eta))-1)/(1-eta)-1);


# total period utility
#var U_period {t in 0..T}=U[t]*R[t];


# welfare/objective function
var W=sum{t in 0..T} U[t]*R[t];

maximize objective_function: W;


subject to constr_captital_dynamics {t in 0..T-1}: K[t+1]=(1-deltaK)*K[t]+I[t];					
#subject to constr_captital_dynamics {t in 0..T-1}: K[t+1]=(1-deltaK[t+j,n])*K[t]+I[t];					
subject to constr_output_gross {t in 0..T}: Qgross[t]=A[t+j,n]*((L[t+j,n]/1000)^(1-gamma))*(K[t]^gamma); 		
subject to constr_output_net {t in 0..T}: Q[t]=(Qgross[t]*(1-Omega[t]))-Lambda[t];				


subject to constr_accounting {t in 0..T}: C[t]=Q[t]-I[t];								
subject to constr_consumtionpercapita {t in 0..T}: c[t]= 1000*C[t]/L[t+j,n];					

#subject to constr_cumulativeemissions {t in 1..T}: Ecum[t]=Ecum[t-1]+(Ec[t-1]*1/3.666); 			


subject to R1step {t in 0..T-1}: R1[t+1]=R1[t]+a1*Ecp[t]-R1[t]/(alpha[t]*tau1);
subject to R2step {t in 0..T-1}: R2[t+1]=R2[t]+a2*Ecp[t]-R2[t]/(alpha[t]*tau2);
subject to R3step {t in 0..T-1}: R3[t+1]=R3[t]+a3*Ecp[t]-R3[t]/(alpha[t]*tau3);
subject to R4step {t in 0..T-1}: R4[t+1]=R4[t]+a4*Ecp[t]-R4[t]/(alpha[t]*tau4);
subject to CCAstep {t in 0..T-1}: CCA[t+1]=CCA[t]+Ec[t]-(Catm[t+1]-Catm[t])*2.12;

#subject to T1step {t in 1..T}: T1[t]=T1[t-1]+(q1*F[t]-T1[t-1])/d1;
#subject to T2step {t in 1..T}: T2[t]=T2[t-1]+(q2*F[t]-T2[t-1])/d2;

subject to TMFstep {t in 1..T}: TMF[t]=expf*TMF[t-1]+amixf*intf[t];
subject to TMSstep {t in 1..T}: TMS[t]=exps*TMS[t-1]+amixs*ints[t];
subject to TDFstep {t in 1..T}: TDF[t]=expf*TDF[t-1]+adeepf*intf[t];
subject to TDSstep {t in 1..T}: TDS[t]=exps*TDS[t-1]+adeeps*ints[t];

##subject to dalphastep {t in 1..T}: dalpha[t]=(a1*alpha[t]*tau1-a1*exp (-100/(alpha[t]*tau1))*(100 + alpha[t]*tau1)+a2*alpha[t]*tau2-a2*exp (-100/(alpha[t]#*tau2))*###(100 ##+ alpha[t]*tau2) + a3*alpha[t]*tau3-a3*exp (-100/(alpha[t]*tau3))*(100 + alpha[t]*tau3)+a4*alpha[t]*tau4-a4*exp (-100/(alpha[t]*tau4))*(100 + #alpha[t]#*tau4))/alpha
##[t];

subject to dalphastep {t in 1..T}: dalpha[t]=-a1*100*exp(-100/(alpha[t-1]*tau1))/alpha[t-1]-a2*100*exp(-100/(alpha[t-1]*tau2))/alpha[t-1]-a3*100*exp(-100/(alpha[t-1]*tau3))/alpha[t-1]-a4*100*exp(-100/(alpha[t-1]*tau4))/alpha[t-1]+a1*(1-exp(-100/(alpha[t-1]*tau1)))*tau1+a2*(1-exp(-100/(alpha[t-1]*tau2)))*tau2+a3*(1-exp(-100/(alpha[t-1]*tau3)))*tau3+a4*(1-exp(-100/(alpha[t-1]*tau4)))*tau4;

subject to alphastep {t in 1..T}: alpha[t]=alpha[t-1]+(rT*(Temp[t]-Temp[t-1])+rc*(CCA[t]-CCA[t-1]))/dalpha[t];

#subject to dalphastep {t in 0..T-1}: dalpha[t+1]=-a1*100*exp(-100/(alpha[t]*tau1))/alpha[t]-a2*100*exp(-100/(alpha[t]*tau2))/alpha[t]-a3*100*exp(-100/(alpha[t]#*tau3))/alpha[t]-a4*100*exp(-100/(alpha[t]*tau4))/alpha[t]+a1*(1-exp(-100/(alpha[t]*tau1)))*tau1+a2*(1-exp(-100/(alpha[t]*tau2)))*tau2+a3*(1-exp(-100/(alpha[t]#*tau3)))*tau3+a4*(1-exp(-100/(alpha[t]*tau4)))*tau4;

#subject to alphastep {t in 1..T}: alpha[t]=alpha[t-1]+(rT*(Temp[t]-Temp[t-1])+rc*(CCA[t]-CCA[t-1]))/dalpha[t];



# Initial conditions
subject to initial_capital: K[0] = K0;
subject to initial_investment: I[0]=I0;
#subject to initial_output: Qgross[0] = Qgross0; 
#subject to initial_Ecum: Ecum[0]=Ecum0;
subject to initial_R1: R1[0]=R10;
subject to initial_R2: R2[0]=R20;
subject to initial_R3: R3[0]=R30;
subject to initial_R4: R4[0]=R40;
subject to initial_CCA: CCA[0]=CCA0;
subject to initial_TMF: TMF[0]=TMF0;
subject to initial_TMS: TMS[0]=TMS0;
subject to initial_TDF: TDF[0]=TDF0;
subject to initial_TDS: TDS[0]=TDS0;
subject to initial_Talpha: alpha[0]=alpha0;
subject to initial_dalpha: dalpha[0]=dalpha0;
subject to controlrestrict1 {t in 0..99}: mu[t]<=1;
subject to controlrestrict2 {t in 100..T}: mu[t]<=1.2;
subject to ceiling  {t in 0..T}: Temp[t]<=Tbar;



