#ifndef GUARD_slice_samp
#define GUARD_slice_samp


/* Code adapted from slice-R-prog.R as downloaded from         */
/* www.cs.toronto.edu/~radford/slice.software.html on 29/6/17  */


#include <Rcpp.h>
using namespace Rcpp;


double slice_samp(double (*g)(double, double, arma::vec, double, int),
                  double sTemp,
                  arma::vec bTemp,
                  double kTemp,
                  int nPath,
                  double x,
                  int m = INFINITY,
                  double w = 1,
                  double lower = -INFINITY,
                  double upper = INFINITY) {
    
    double x0=x, L=x, R=x, gx0, logy, u, x1, gx1;
    int J, K;
    
    gx0 = g(x0,sTemp,bTemp,kTemp,nPath);
    
    /* Determine the slice level, in log terms */
    
    logy=gx0-R::rexp(1);
    
    /* Find the initial interval to sample from */
    
    u=R::runif(0,w);
    L=x0-u;
    R=x0+(w-u);
    
    /* Expand the interval until its ends are outside the slice, or until the limit on steps is reached */
    
    if(isinf(m)){

        while(L>lower&&g(L,sTemp,bTemp,kTemp,nPath)>logy){
            L=L-w;
        }
        
        while(R<upper&&g(R,sTemp,bTemp,kTemp,nPath)>logy){
            R=R+w;
        }
        
    }else{
        
        if(m>1){
            
            J=static_cast <int> (std::floor(R::runif(0,m)));
            K=(m-1)-J;
        
            while(J>0&&L>lower&&g(L,sTemp,bTemp,kTemp,nPath)>logy){
                L=L-w;
                J=J-1;
            }
            
            while(K>0&&R<upper&&g(R,sTemp,bTemp,kTemp,nPath)>logy){
                R=R+w;
                K=K-1;
                
           }
            
        }
        
    }
    
    /* Shrink interval to lower and upper bounds */
    
    if(L<lower){
        L=lower;
    }
    if(R>upper){
        R=upper;
    }
    
    /* Sample from the interval, shrinking it on each rejection */
    

    x1=R::runif(L,R);
    gx1=g(x1,sTemp,bTemp,kTemp,nPath);
        
    while(gx1<logy){
        
        if(x1>x0){
            R=x1;
        }else{
            L=x1;
        }
        
        x1=R::runif(L,R);
        gx1=g(x1,sTemp,bTemp,kTemp,nPath);
        
    }
    
    return(x1);
    
}

#endif
