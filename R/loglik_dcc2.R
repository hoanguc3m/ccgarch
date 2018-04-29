# the objective function in the 2nd step DCC estimation. This is to be minimised!
loglik.dcc2 <- function(param, dvar){        # dvar is the standardised residuals
   nobs <- dim(dvar)[1]
   ndim <- dim(dvar)[2]
   DCC <- dcc.est(dvar, param)$DCC

#   lf <- numeric(ndim)  
   lf <- numeric(nobs)  # bug fixed on 2013.08.18
   for( i in 1:nobs){                        
      R <- matrix(DCC[i,], ndim, ndim)
      invR <- solve(R)
      lf[i] <- 0.5*(log(det(R)) +sum(dvar[i,]*crossprod(invR,dvar[i,])))  
   }
   sum(lf)
}


loglik.dcc2.std <- function(param, dvar){        # dvar is the standardised residuals
    nobs <- dim(dvar)[1]
    ndim <- dim(dvar)[2]
    DCC <- dcc.est(dvar, param)$DCC
    nu = param[3]
    
    #   lf <- numeric(ndim)  
    lf <- numeric(nobs)  # bug fixed on 2013.08.18
    for( i in 1:nobs){                        
        R <- matrix(DCC[i,], ndim, ndim)
        invR <- solve(R)
        lf[i] <- 0.5*(log(det(R)) ) + 0.5*(nu + ndim) * log( 1 + sum(dvar[i,]*crossprod(invR,dvar[i,]))/nu ) +
                    0.5 * ndim * log( pi * nu) + lgamma(nu*0.5) - lgamma( 0.5 * (nu + ndim) )
        #  dmvt(dvar[i,], delta = rep(0, ndim), sigma = R, df = nu, log = T)
    }
    sum(lf)
}

loglik.dcc2.ghst <- function(param, dvar){        # dvar is the standardised residuals
    nobs <- dim(dvar)[1]
    ndim <- dim(dvar)[2]
    DCC <- dcc.est(dvar, param)$DCC
    nu = param[3]
    gamma = rep(param[4],ndim)
    
    #   lf <- numeric(ndim)  
    lf <- numeric(nobs)  # bug fixed on 2013.08.18
    for( i in 1:nobs){                        
        R <- matrix(DCC[i,], ndim, ndim)
        invR <- solve(R)
        Q_x <- sum(dvar[i,]*crossprod(invR,dvar[i,]))
        Q_gamma <- sum(gamma*crossprod(invR,gamma)) 
        
        const <- log(2) * (1 - (nu + ndim)*0.5) - lgamma(nu*0.5) - 0.5*ndim* log(pi * nu) - 0.5* log( det(R) )
        const = const + sum(dvar[i,]*crossprod(invR,gamma)) - 0.5*(nu + ndim) * log( 1 + Q_x/nu)
        
        
        bessel_out <- besselK(x = sqrt( (nu + Q_x) * Q_gamma), nu = (nu+ndim)*0.5)
        
        
        
        lf[i] <- - const - log(bessel_out)  - 0.25 * (nu+ndim) * log( (nu + Q_x) * Q_gamma )
    }
    sum(lf)
}


# the log-likelihood function for the 2nd step DCC estimation
#loglik.dcc2 <- function(param, dvar){        # dvar is the standardised residuals
#   nobs <- dim(dvar)[1]
#   ndim <- dim(dvar)[2]
#   
#   if(sum(param)>1|sum(param)<0){
#      param <- c(0, 0)
#   }
#   
#   DCC <- dcc.est(dvar, param)$DCC
#
#   lf <- numeric(ndim)
#   for( i in 1:nobs){                        
#      R <- matrix(DCC[i,], ndim, ndim)
#      invR <- solve(R)
#      lf[i] <- -0.5*(log(det(R)) +sum(dvar[i,]*crossprod(invR,dvar[i,])))
#   }
#   
#   sum(-lf)
#}





