\name{dynNetModelSamp}
\alias{dynNetModelSamp}
\title{
dynNetModelSamp
}
\description{
Usage:\cr
posteriorSamples <- dynNetModelSamp(y, X, l, k, tau, a, nu, rho, B, timeIdxs, nBurn, nSamp, nThin)\cr
\cr
Where:\cr
y is a vector of length n data-points for the response variable.\cr
X is a matrix of dimension n*p, with p columns of predictor variables, each with n data-points.\cr
l is the `sparsity within time' hyperparameter.\cr
k is the `sparsity across time' hyperparameter.\cr
tau is the initial value of the model precision parameter. Suggested value is 1.\cr
a is the initial value of the model intercept parameter. Suggested value is 0.\cr
nu is a vector of length p of initial values for the model `precision scaling' parameter nu. Suggested value is a vector of ones.\cr
rho is a vector of length p of the model `correlation across time' parameter rho. Suggested value is a vector of ones.\cr
B is a matrix of dimension T*p of regression coefficients (T is the number of time-points). Suggested value is a matrix of zeros.\cr
timeIdxs is a matrix of dimension T*2, with each row giving the start and finish index of the data for each time-point. These indices correspond to the elements of y, and to the rows of X.\cr
nBurn is the number of burn-in samples taken, before samples start being saved.\cr
nSamp is the number of samples to be taken, after nBurn is completed.\cr
nThin is the number of samples to thin by (i.e, save only every nThin-th sample).\cr
\cr
posteriorSamples is a list containing elements:\cr
aSamp - a vector of length nSamp/nThin, containing the posterior samples of a.\cr
tauSamp - a vector of length nSamp/nThin, containing the posterior samples of tau.\cr
nuSamp - a matrix of dimension (nSamp/nThin) * p, with each column containing the posterior samples for the nu parameter for the corresponding column of X.\cr
rhoSamp - a matrix of dimension (nSamp/nThin) * p, with each column containing the posterior samples for the rho parameter for the corresponding column of X.\cr
Bsamp - an array of dimension T * p * (nSamp/nThin), with columns corresponding to the columns of X, rows corresponding to the T time-points, and slices corresponding to posterior samples.
}
\details{
}
\author{
Thomas E. Bartlett, Ioannis Kosmidis, and Ricardo Silva.
}
\references{
Thomas E. Bartlett, Ioannis Kosmidis, and Ricardo Silva (2020). Two-way sparsity for time-varying networks, with applications in genomics. Preprint available: http://arxiv.org/abs/1802.08114
}
\keyword{}
\seealso{
}
\examples{
}
