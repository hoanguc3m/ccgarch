#*****************************************************************************************************************
# The 2nd step DCC estimation.
dcc.estimation2 <- function(dvar, para, gradient=0, dist = "gauss"){ # dvar must be standardised residuals
   if (dist == "gauss"){
       
       resta <- rbind(c(-1, -1), diag(2))
       restb <- c(-1, 0, 0)
       
       if(gradient!=0){
           step2 <- constrOptim(theta=para, f=loglik.dcc2, grad=grad.dcc2, ui=resta, ci=restb, mu=1e-5, dvar=dvar)
       } else {
           step2 <- constrOptim(theta=para, f=loglik.dcc2, grad=NULL, ui=resta, ci=restb, mu=1e-5, dvar=dvar)
       }
   }
   
   if (dist == "std"){
       resta <- rbind(c(-1, -1,0), diag(3), c(0, 0, -1))
       restb <- c(-1, 0, 0, 2, -30)
           step2 <- constrOptim(theta=para, f=loglik.dcc2.std, grad=NULL, ui=resta, ci=restb, mu=1e-5, dvar=dvar)
   
   }
   
   step2
}
#para = c(0.01, 0.98, 5)
