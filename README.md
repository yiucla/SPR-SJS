# SPR-SJS

Liu, Y., Xu, J., & Li, G. Sure joint feature screening in nonparametric transformation model with right censored data.

Description:

This program implements the SPR-SJS joint screening method of the paper.

Main function: SPR-SJS(X,Y,delta,beta_ini,k)

X: Design matrix

Y: Time to event or censoring

delta: 0=censoring; 1 otherwise

beta_ini: Initial value for IHT algorithm

k: Screened model size

S: Active index set  

beta_S: Parameter estimates of active variables

Output: Active index set and parameter estimates of active variables.
