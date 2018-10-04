# IVQR-GMM-R-codes

R codes for exact computation of the IVQR GMM estimator. Detailed description of the estimator and computation can be found in the paper:

Chen, Le-Yu and Lee, Sokbae (September 2017), "Exact computation of GMM estimators for instrumental variable quantile regression models".

The paper has been published at Journal of Applied Econometrics. See https://onlinelibrary.wiley.com/doi/full/10.1002/jae.2619.

Le-Yu Chen also provided a MATLAB version of this model here: https://github.com/LeyuChen/IVQR-GMM-computation-codes. The Python and R codes are inspired by the MATLAB codes.

Both the R function named IVQR_GMM can be used to calculate the GMM estimator for instrumental variable quantile regression models. The implementation involves Gurobi solver and some additional function defined in the same file. The Gurobi solver is freely available for academic purposes.

There is also a Python version of the codes: https://github.com/jordanxzz/IVQR-GMM-Python-codes.

## Download (require the devtools package)
```
install_github("jordanxzz/ivqrgmm-r")
```
## Main Functions
- IVQR_GMM:
  Used to do the instrumental variable quantile regression using the GMM estimator.

## Explanation
IVQR_GMM(y, w, z, tau, intercept=False, T=0, abgap=0, bnd = None)

function input :

    y        : vector of outcomes
    
    w        : (n by k) matrix of the covariate dataset, include exogeneous variables
    
    z        : (n by p ) matrix of the instrument variable dataset, include exogeneous variables
    
    tau      : quantile index
    
    intercept: False ==> The function will NOT add intercept term automatically, include it in w and z if needed
      
               True  ==> The function will ADD intercept term to w and z automatically
               
    T        : scalar. T=0 by default. If T>0, then T is the time limit specified for early termination
    
               of the MIO solver. Otherwise, the MIO solver keeps running until convergence.
               
    abgap    : the absolute gap specified for early termination of the MIO solver. abgap=0 by default.
    
    bnd      : (k by 2) matrix where the first and second columns respectively store the 
               
               lower and upper bounds of the unknown coefficients. No boundary by default.
                
function output :

    theta_hat   : the vector of the coefficient estimates
    
    s_hat       : the estimated asymptotic standard errors
    
    obj_v       : the value of the GMM objective function
    
    gap         : the MIO optimization gap value in case of early termination
    
    rtime       : the time used by the MIO solver in the estimation procedure
    
    ncount      : the number of nodes already explored by the MIO solver 

## Requirements
Requires R 3.5.1 or above and Gurobi solver R API (available free for academic purposes).

## Acknowledgments
* Le-Yu Chen (lychen@econ.sinica.edu.tw)
* Sokbae Lee (sl3841@columbia.edu)
The R codes in this package are developed based on the Matlab implementation by Le-Yu Chen (https://github.com/LeyuChen/IVQR-GMM-computation-codes).
