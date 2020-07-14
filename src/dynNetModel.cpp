# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

#include "slice_samp.h"
#include "samplerFunctions.h"
#define minVal (DOUBLE_EPS*1.1)

// [[Rcpp::export()]]
Rcpp::List dynNetModelSamp(NumericVector y,
                           arma::mat X,
                           double lambda,
                           double k,
                           double tau,
                           double a,
                           arma::vec nu,
                           arma::vec rho,
                           arma::mat B,
                           NumericMatrix dataTimeLocs,
                           int nBurn,
                           int nSamp,
                           int nThin) {
    
    RNGScope scope; /* ensure RNG gets set / reset */

    
    /* define parameters and variables */
    int n = X.n_rows, nPath=dataTimeLocs.nrow(), nParents = B.n_cols, tempStart, tempFinish;
    double tempSum, rhoTemp, nuTemp, psiTemp, nRep=nBurn+nSamp, sampSize=nSamp/nThin;
    NumericMatrix precTemp(nPath,nPath);
    arma::mat precArma(nPath,nPath), nuSamp(sampSize,nParents), rhoSamp(sampSize,nParents);
    arma::vec bTemp(nPath), aSamp(sampSize), tauSamp(sampSize);
    arma::cube Bsamp(nPath,nParents,sampSize);

    
    /* summary stats for calculating posteriors */
    NumericMatrix sumX2j(nPath,nParents);
    arma::mat sumX2jArma;
    for (int iPath = 0; iPath != nPath; ++iPath) {
        tempStart=dataTimeLocs[iPath];
        tempFinish=dataTimeLocs[nPath+iPath];
        for (int iParent = 0; iParent != nParents; ++iParent) {
            tempSum=0;
            for(int iLoc = tempStart; iLoc != (tempFinish+1); ++iLoc) {
                tempSum=tempSum+pow(X[iLoc-1+n*iParent],2);
            }
            sumX2j[iPath+nPath*iParent]=tempSum;
        }
    }
    sumX2jArma=arma::mat(sumX2j.begin(), sumX2j.nrow(), sumX2j.ncol(), false);
 
    
    /* upper bound rho for slice-sampling */
    for(int j = 0; j != nParents; ++j) {
        
        if(rho[j]>1-1e-6){
            rho[j]=1-1e-6;
        }
        
    }
    
    
    /* main sampling loop */
    for (int iRep = 0; iRep != nRep; ++iRep) {
       
  
        /* sample tau */
        tau=R::rgamma(n/2+1,1/(1+SSE(y,X,a,B,dataTimeLocs)/2));
        if(tau<minVal){
            tau=minVal;
        }
        
        
        /* sample a */
        double tempAmean=0;
        for(int iPath = 0; iPath != nPath; ++iPath) {
            tempStart=dataTimeLocs[iPath];
            tempFinish=dataTimeLocs[nPath+iPath];
            for(int iLoc = tempStart; iLoc != (tempFinish+1); ++iLoc) {
                tempSum=0;
                for(int iParent = 0; iParent != nParents; ++iParent) {
                    tempSum=tempSum+X[iParent*n+iLoc-1]*B[iParent*nPath+iPath];
                }
                tempAmean=tempAmean+(y[iLoc-1]-tempSum);
            }
        }
        a=R::rnorm(tempAmean*tau/(1+n*tau),1/sqrt(1+n*tau));

        
        /* sample rho */
        for(int j = 0; j != nParents; ++j) {
            rhoTemp=rho[j];
            nuTemp=nu[j];
            bTemp=B.col(j);
            rho[j]=slice_samp(gRhoLog,nuTemp,bTemp,k,nPath,rhoTemp,INFINITY,1,0,1-1e-6);
        }
        
        
        /* sample nu */
        for(int j = 0; j != nParents; ++j) {
            rhoTemp=rho[j];
            bTemp=B.col(j);
            precTemp=precMat(rhoTemp,nPath);
            precArma=arma::mat(precTemp.begin(), precTemp.nrow(), precTemp.ncol(), false);
            psiTemp=sum(bTemp.t()*precArma*bTemp);
            if(psiTemp<minVal){
                nu[j]=1/R::rgamma(0.5,2/pow(lambda,2));
            }else{
                nu[j]=rinvgauss(lambda/sqrt(psiTemp),pow(lambda,2));
            }
            if(nu[j]<minVal){
                nu[j]=minVal;
            }
        }

        
        /* sample B */
        arma::mat VinvArma(nPath,nPath);
        arma::mat precTildArma(nPath,nPath);
        arma::vec mTemp(nPath);
        arma::vec mTildTemp(nPath);
        arma::mat tempMat(nPath,nPath);
        arma::vec tempVec(nPath);
        arma::mat Bsub;
        arma::mat Xsub;

        for(int j = 0; j != nParents; ++j) {

            precTemp=nu[j]*precMat(rho[j],nPath);
            precArma=arma::mat(precTemp.begin(), precTemp.nrow(), precTemp.ncol(), false);
            VinvArma=arma::diagmat(sumX2jArma.col(j))*tau;
            precTildArma=precArma+VinvArma;
            Xsub=X;
            Xsub.shed_col(j);
            Bsub=B;
            Bsub.shed_col(j);
            for (int iPath = 0; iPath != nPath; ++iPath) {
                tempStart=dataTimeLocs[iPath];
                tempFinish=dataTimeLocs[nPath+iPath];
                tempSum=0;
                for(int iLoc = tempStart; iLoc != (tempFinish+1); ++iLoc) {
                    tempSum=tempSum+sum(X[n*j+iLoc-1]*(y[iLoc-1]-Xsub.row(iLoc-1)*Bsub.row(iPath).t()-a));
                }
                mTemp[iPath]=tau*tempSum;
            }
            mTildTemp=solve(precTildArma,mTemp);
            
            for(int iPath=0; iPath!=nPath; ++iPath){
                tempVec[iPath]=R::rnorm(0,1);
            }
            tempMat=arma::chol(precTildArma);
            tempVec=mTildTemp+arma::solve(tempMat,tempVec);
            B.col(j)=tempVec;
        
        }

        
        /* store samples */
        if(iRep>=nBurn&&(iRep-nBurn)%nThin==0){
            aSamp[(iRep-nBurn)/nThin]=a;
            tauSamp[(iRep-nBurn)/nThin]=tau;
            nuSamp.row((iRep-nBurn)/nThin)=nu.t();
            rhoSamp.row((iRep-nBurn)/nThin)=rho.t();
            Bsamp.slice((iRep-nBurn)/nThin)=B;
        }

        
    }

    
    
    return Rcpp::List::create(Rcpp::Named("aSamp")=aSamp,Rcpp::Named("tauSamp")=tauSamp,Rcpp::Named("nuSamp")=nuSamp,Rcpp::Named("rhoSamp")=rhoSamp,Rcpp::Named("Bsamp")=Bsamp);

    
}
