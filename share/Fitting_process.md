# Fitting process in `fitDR()`

a `trace` level is available in `control` argument to display fitting process's traces

## distributions mmbefd(a,b) 

- initial value set via log/exp link function
- prefit is performed
- if non NA, initial value is updated for the two domains
- two optimizations are carried out on the two domains via `mledist()` for MLE and `constrOptim.nl()` directly
- select the optimal value
- try to compute optimale Hessian value

## distribution MMBEFD(g,b)

- initial value set via log/exp link function for b only
- try to update b via a second-order equation
- prefit is performed
- if non NA, initial value is updated for the two domains
- two optimizations are carried out on the two domains via `mledist()` for MLE and `constrOptim.nl()` directly
- select the optimal value
- try to compute optimale Hessian value

## distribution oiunif

- initial value set via `etl()`
- optimisation is performed via `fitdist()` directly

## distributions oistpareto, oibeta, oigbeta


- initial value set : a=1, MME, Theil coefficient respectively
- prefit is performed for MLE only
- optimisation is performed via `fitdist()` for MLE and `optim()` for TLMME

