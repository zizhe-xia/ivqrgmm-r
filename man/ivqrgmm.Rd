\name{ivqrgmm}
\alias{ivqrgmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Instrument Variable Quantile Regression with GMM Estimators
%%  ~~function to do ... ~~
}
\description{
This is a package for Instrumental Variable Quantile Regression (IVQR) using Generalized Methods of Moments (GMM) estimation. It is based on the paper: Chen, Le-Yu and Lee, Sokbae (September 2017), "Exact computation of GMM estimators for instrumental variable quantile regression models".
}
\usage{
IVQR_GMM(y, w, z, tau, intercept, tlimit, abgap, bnd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{vector of outcomes.}
  \item{w}{(n by k) matrix of the covariate dataset.}
  \item{z}{(n by p ) matrix of the instrument variable dataset.}
  \item{tau}{quantile index}
  \item{intercept}{FALSE by default, and the function will NOT add intercept term automatically, include it in w and z if needed.}
  \item{tlimit}{scalar, 0 by default. If tlimit>0, then tlimit is the time limit specified for early termination of the MIO solver. Otherwise, the MIO solver keeps running until convergence.}
  \item{abgap}{scalar, 0 by default. The absolute gap specified for early termination of the MIO solver.}
  \item{bnd}{(k by 2) matrix, NA by default. The first and second columns of  the matrix respectively store the lower and upper bounds of the  unknown coefficients.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{theta_hat}{the vector of the coefficient estimates}
  \item{s_hat}{the estimated asymptotic standard errors}
  \item{obj_v}{the value of the GMM objective function}
  \item{gap}{the MIO optimization gap value in case of early termination gap = 0 ==> optimal solution is found within the time limit}
  \item{rtime}{the time used by the MIO solver in the estimation procedure}
  \item{ncount}{the number of nodes already explored by the MIO solver}
}
\references{
\emph{Exact computation of GMM estimators for instrumentalvariable quantile regression models} \url{https://onlinelibrary.wiley.com/doi/full/10.1002/jae.2619} by Chen, Le-Yu and Lee, Sokbae (September 2017). The R codes are based on the MATLAB codes available at: \url{https://github.com/LeyuChen/IVQR-GMM-computation-codes}
}
\author{
Zizhe Xia
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Original MATLAB implementation of the model \url{https://github.com/LeyuChen/IVQR-GMM-computation-codes}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
\keyword{ models }% __ONLY ONE__ keyword per line
