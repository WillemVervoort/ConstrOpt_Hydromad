## simple constrained optimisation
## to be implemented later using fitByOptim in hydromad

## what we want:
# minimise RMSE s.t. all ||e_i|| < beta
# the question is how do we write ||e_i|| < beta as a simple condition
# My first thought was minimise beta for sum(||e_i|| - beta)
# but that seems very similar to RMSE
# or to write a scaled penalty function if any(||e_i|| > beta)
# here the scaling would be ||e_i|| - beta and then sum?
# but isn't that same as using LCC?


# what about a likelihood function where we control the sd of the distribution
# seems like the first possible solution
require(hydromad)

foo <- function(parsToFit,Hmodel) {
  # parsToFit are the parameters to fit in the overall scheme
  # Hmodel is the defined hydromad model
  # put parameters in Hmodel
 # browser()
  n <- length(coef(Hmodel))
  Nmodel <- update(Hmodel,newpars=parsToFit[1:n])
  # simulate values
  pred <- Nmodel$fitted.values
  #browser()
  # calculate objective function
  R <- dnorm(Nmodel$data$Q - pred, parsToFit[n+1],parsToFit[n+2])
  objfun <-  -  sum(log(R),na.rm=T)
  return(objfun)
}




data(Cotter)

data.cal <- window(Cotter, start = "1990-01-01",
                   end = "1998-12-31")

# Define the model, important to define return_state=T
Cotter_mod <- hydromad(DATA=data.cal,
                       sma = "gr4j", routing = "gr4jrouting", 
                       x1 = c(100,1500), x2 = c(-30,20), x3 = c(5,500), 
                       x4 = c(0.5,10), etmult=c(0.01,0.5), 
                       return_state=TRUE)

# Fit traditional fit
# Using Optim algorithm for fitting
Cotter_fit_o<- fitByOptim(Cotter_mod,  
                        objective=~hmadstat("r.squared")(Q,X))
# Extract the coefficients and the summary
summary(Cotter_fit_o)
xyplot(Cotter_fit_o)


inpars <- suppressWarnings(sapply(coef(Cotter_mod), mean))
Cotter_in <- update(Cotter_mod, newpars=inpars)

inpars$mu <- 0
inpars$sd <- 0.5
lb <- c(sapply(coef(Cotter_mod),min),-1,0.01)
ub <- c(sapply(coef(Cotter_mod),max),1,2)

fit <- optim(do.call(c,inpars),foo,Hmodel=Cotter_in, method = "L-BFGS-B",
             lower = lb, upper = ub)
str(fit)
xyplot(update(Cotter_in,newpars=fit$par[1:5]))
Cotter_fit <-update(Cotter_in,newpars=fit$par[1:5])
hist(residuals(Cotter_fit))

nseStat(data.cal$Q,Cotter_fit$fitted.values)
