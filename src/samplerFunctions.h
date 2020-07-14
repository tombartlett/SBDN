#ifndef GUARD_samplerFunctions
#define GUARD_samplerFunctions

#include <RcppArmadillo.h>

NumericMatrix precMat(double rho,
                      int nPath) {
    
    double adj, invDetTemp;
    NumericMatrix precTemp(nPath,nPath);
    
    adj=1e-15;
    
    if(rho>(1-adj)||rho<adj){
        rho=adj+(1-2*adj)*rho;
    }
    
    invDetTemp=1/(1-pow(rho,2));
    precTemp[0]=invDetTemp;
    precTemp[nPath*(nPath-1)+nPath-1]=invDetTemp;
    for (int iPath = 2; iPath != nPath; ++iPath) {
        precTemp[(iPath-1)*nPath+iPath-1]=(1+pow(rho,2))*invDetTemp;
    }
    for (int iPath = 1; iPath != nPath; ++iPath) {
        precTemp[iPath*nPath+iPath-1]=(-rho*invDetTemp);
        precTemp[(iPath-1)*nPath+iPath]=(-rho*invDetTemp);
    }
    return precTemp;
    
}

double SSE(NumericVector y,
           arma::mat X,
           double a,
           arma::mat B,
           NumericMatrix dataTimeLocs){
    int n = X.n_rows, nPath = dataTimeLocs.nrow(), nParents = B.n_cols,tempStart, tempFinish;
    Rcpp::NumericVector tempSum;
    double sseTemp=0;
    for(int iPath = 0; iPath != nPath; ++iPath) {
        tempStart=dataTimeLocs[iPath];
        tempFinish=dataTimeLocs[nPath+iPath];
        for(int iLoc = tempStart; iLoc != (tempFinish+1); ++iLoc) {
            tempSum=0;
            for(int iParent = 0; iParent != nParents; ++iParent) {
                tempSum=tempSum+X[iParent*n+iLoc-1]*B[iParent*nPath+iPath];
            }
            sseTemp=sseTemp+sum(pow(y[iLoc-1]-tempSum-a,2));
        }
    }
    return(sseTemp);
}

double gRhoLog(double rhoTemp,
               double nuTemp,
               arma::vec bTemp,
               double kTemp,
               int nPath) {

    NumericMatrix precTemp;
    arma::mat precArma;
    double logF;
    
    precTemp=precMat(rhoTemp,nPath);
    precArma=arma::mat(precTemp.begin(), precTemp.nrow(), precTemp.ncol(), false);
    logF=sum(abs(arma::log_det(precArma))/2-bTemp.t()*precArma*bTemp*nuTemp/2+kTemp*rhoTemp);
    
    return(logF);

}

double rinvgauss(const double mu, const double lambda)
{

    double v,y,x,z;
    
    v=R::rnorm(0,1);
    y=pow(v,2);
    x=mu+pow(mu,2)*y/(2.0*lambda)-(mu/(2.0*lambda))*sqrt(4.0*mu*lambda*y+pow(mu,2)*pow(y,2));
    
    z=R::runif(0,1);
    if(z<=mu/(mu+x)) return x;
    else return pow(mu,2)/x;
    
}

#endif
