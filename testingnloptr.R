# using nloptr to do a constrained optimisation of a hydromad moddel
# following the same lines as the first attempt
require(hydromad)
require(nloptr)
require(epiR)

# data and model
# use the satellite hydromad example data
load("C:/Users/rver4657/ownCloud/working/SatelliteHydromad/Examples/Cotter.rdata")
load("C:/Users/rver4657/ownCloud/working/SatelliteHydromad/Examples/CotterMODISET.rdata")
#data(Cotter)


data.cal <- window(Cotter, start = "2000-01-01",
                   end = "2005-12-31")

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
plot(residuals(Cotter_fit_o))
hist(residuals(Cotter_fit_o),breaks=50)
linsccc <- epi.ccc(data.cal$Q[101:length(data.cal$Q)],
                   fitted(Cotter_fit_o))
linsccc$rho.c$est

# Now define a function to fit with nloptr
# following: https://cran.r-project.org/web/packages/nloptr/vignettes/nloptr.pdf
# but our function will not use a gradient
# for this we need to have the objectivefunction

objfun <- function(par, beta = 0.5, in_data,mod) {
  # par is vector of initial values of parameters to be optimised
  # this includes beta
  # in_data are the input data
  # mod is a hydromad model
  mod_run <- update(mod, x1 = par[1], x2 = par[2], x3 = par[3], 
                    x4 = par[4], etmult=par[5])
  #browser()
  error <- fitted(mod_run) - in_data$Q
  
  SSE <- sqrt(sum((error)^2,na.rm=T))/length(error)
  
  
  ValToMinimise <- SSE
  return(ValToMinimise)
}

x01 <- objfun(600, 0, 100, 2, 0.15,
               mod=Cotter_mod,in_data=data.cal)

# Constraint function
Constr_Mod <- function(par, beta, in_data, mod) {
  # par is vector of initial values of parameters to be optimised
  # this includes beta
  # in_data are the input data
  # mod is a hydromad model
  mod_run <- update(mod, x1 = par[1], x2 = par[2], x3 = par[3], 
                    x4 = par[4], etmult=par[5])
  #browser()
  error <- fitted(mod_run) - in_data$Q
  minvalue <- sqrt(error^2) - beta
  return(coredata(minvalue))
}

x02 <- Constr_Mod(par=c(600, 0, 100, 2, 0.15), beta=0.5,
              mod=Cotter_mod,in_data=data.cal)
plot(x02)
# Now invoke nloptr:
res1 <- nloptr( x0=c(600, 0, 100, 2, 0.15),
                  eval_f=objfun,
                  lb = c(100, -30, 5, 0.6 , 0.01),
                  ub = c(1500, 20, 500, 10, 0.5),
                  eval_g_ineq = Constr_Mod,
                  # algorithm simply what is suggested in the vignette
                  opts = list("algorithm"="NLOPT_LN_COBYLA",
                              "xtol_rel"=1.0e-7, print_level = 1,
                              maxeval=1500),
                  beta = 5,
                in_data = data.cal,
                mod = Cotter_mod
                )

#str(res1)
# plot the solution:
model_fit <- update(Cotter_mod, x1 = res1$solution[1], 
                    x2 = res1$solution[2], x3 = res1$solution[3], 
                    x4 = res1$solution[4], etmult=res1$solution[5])
# make a plot
plot(fitted(model_fit))
lines(data.cal$Q,col="red")
hist(data.cal$Q[101:length(data.cal$Q)]-fitted(model_fit),breaks=50)


# Changing beta makes a big difference
# Now we want to run a loop over beta values
# plot the resulting r.sq.sqrt and rel.bias on a figure together
# with the result of the unconstrained optimisation

beta1 <- seq(2,4,length=10)
n <- length(beta1)

result <- data.frame(beta = beta1, mean_error = rep(0,n),
                     var_error = rep(0,n),
                     r.sq.srqt = rep(0,n), 
                     Linccc = rep(0,n))

unconstrResults <- rep(0,4)

unc_error <- data.cal$Q[101:length(data.cal$Q)]-fitted(model_fit)
unconstrResults[1:2] <- c(mean(unc_error),var(unc_error))

for (i in 1:length(beta1)) {
  res1 <- nloptr( x0=c(600, 0.1, 100, 2, 0.15),
                  eval_f=objfun,
                  lb = c(100, -30, 5, 0.6 , 0.01),
                  ub = c(1500, 20, 500, 10, 0.5),
                  eval_g_ineq = Constr_Mod,
                  opts = list("algorithm"="NLOPT_LN_COBYLA",
                              "xtol_rel"=1.0e-8, print_level = 1,
                              maxeval=1500),
                  beta = beta1[i],
                  in_data = data.cal,
                  mod = Cotter_mod
  )
  
  model_fit <- update(Cotter_mod, x1 = res1$solution[1], 
                        x2 = res1$solution[2], x3 = res1$solution[3], 
                        x4 = res1$solution[4], etmult=res1$solution[5])
  # make a plot
  par(mfrow=c(1,2))
  plot(fitted(model_fit), xlab = "Date", ylab="Flow in mm",
       main = paste("Fitted and observed"))
  lines(data.cal$Q,col="red")
  hist(data.cal$Q-fitted(model_fit),breaks=50, 
       main = paste("beta =", beta1[i]) )
  par(mfrow=c(1,1))
  
  
  # calculate statistics
  result[i,2] <- mean(data.cal$Q[101:length(data.cal$Q)]-
                                      fitted(model_fit))
  result[i,3] <- var(data.cal$Q[101:length(data.cal$Q)]-
                        fitted(model_fit))
  result[i,4] <- hmadstat("r.sq.sqrt")(data.cal$Q[101:length(data.cal$Q)],
                                       fitted(model_fit))
  result[i,5] <- epi.ccc(data.cal$Q[101:length(data.cal$Q)],
                         fitted(model_fit))$rho.c$est
}

unconstrResults[3:4] <- summary(model_fit)[c("rel.bias","r.sq.sqrt")]

# now make plots of the resulting r.sq.sqrt relative to the unconstrained fit
par(mfrow=c(2,2))
plot(result[,1], result[,2], xlab="beta", ylab="mean error",
     ylim=c(-0.5,1), col="blue", lwd=5)
lines(beta1,rep(unconstrResults[1],n), col="red",lty=2)
plot(result[,1], result[,3], xlab="beta", ylab="var error",
     ylim=c(-0.2,1), col="blue", lwd=5)
lines(beta1,rep(unconstrResults[2],n), col="red",lty=2)
legend("topright",c("constrained","standard"),lty=c(NA,2),pch=c(1,NA),
       lwd=c(5,1),col=c("blue","red"))
plot(result[,1], result[,4], xlab="beta", ylab="r.sq.sqrt flow",
     ylim=c(-0.5,1), col="blue", lwd=5)
lines(beta1,rep(unconstrResults[4],n), col="red",lty=2)
plot(result[,1], result[,5], xlab="beta", ylab="Lin's ccc flow",
     ylim=c(0,1), col="blue", lwd=5)
lines(beta1,rep(linsccc$rho.c$est,n), col="red",lty=2)
par(mfrow=c(1,1))


####################################
# Including ET in the nloptr
# again start of with the earlier work on hydromad
setwd("C:/Users/rver4657/ownCloud/working/SatelliteHydromad/Examples")
source("setup.R")

# using ETa.merge()
Flow.Modis.zoo <- ETa.merge(Flowdata=Cotter,ETdata=Cot_MODISET)

# remake the calibration data
data.modis.cal <- window(Flow.Modis.zoo, start = "2000-01-01",end = "2005-12-31")

# Because we have rebuilt data.cal, redefine the model
Cotter_mod_M <- hydromad(DATA=data.modis.cal,
                         sma = "gr4j", routing = "gr4jrouting", 
                         x1 = c(100,1500), x2 = c(-30,20), 
                         x3 = c(5,500), x4 = c(0.5,10), 
                         etmult=c(0.01,0.5), 
                         return_state=TRUE)
hydromad.options(trace=TRUE)
options(warn=1)

Cotter_fit_q_ET <- FitMODbySCE(mod=Cotter_mod_M, FIT_Q=T, 
                               FIT_Q_ET = T, FIT_ET = T)
summary(Cotter_fit_q_ET, items = c("rel.bias", "r.squared","r.sq.sqrt", "r.sq.log"))

linsccc <- epi.ccc(data.modis.cal$Q[101:length(data.modis.cal$Q)],
                   fitted(Cotter_fit_q_ET[["Fit ET and Q"]]))
unconstrResults[4] <- linsccc$rho.c$est

# essentially the same solution as earlier, ET calibration has no effect
# plot
xyplot(Cotter_fit_q_ET)
aET_obs <- aggregate(data.modis.cal$aET,
                     list(date=data.modis.cal$et.period),sum)
ET_pred <- aggregate(Cotter_fit_q_ET[["Fit ET and Q"]]$U$ET,
                     list(date=data.modis.cal$et.period),
                     sum)

nseStat(coredata(aET_obs),coredata(ET_pred))
plot.ET(data.modis.cal,Cotter_fit_q_ET[["Fit ET and Q"]])


# Do the same using constrained optimisation using nloptr
# Rewrite the optimisation function to include ET
# This requires using some of the functions from the satellite hydromad
# define a function to fit with Optim
# follow the earlier nloptr example

objfun_ET <- function(par, beta = 0.5, in_data, mod) {
  # par is vector of initial values of parameters to be optimised
  # this includes beta
  # in_data are the input data, which includes MODIS ET
  # mod is a hydromad model
  mod_run <- update(mod, x1 = par[1], x2 = par[2], x3 = par[3], 
                    x4 = par[4], etmult=par[5])
  #browser()
  # Q objective functions
  error_Q <- fitted(mod_run) - in_data$Q
  
  SSE_Q <- sqrt(sum((error_Q)^2,na.rm=T))/length(error_Q)
  # ET objective function
  #browser()
  aET_fin <- aggregate(in_data$aET,list(date=in_data$et.period),sum)
  ET_fin <- aggregate(mod_run$U$ET,
                      list(date=in_data$et.period),
                      sum)
  #browser()
  error_ET <- coredata(ET_fin) - coredata(aET_fin)
  
  SSE_ET <- sqrt(sum((error_ET)^2,na.rm=T))/length(error_ET)
  
  w <- 0.5
  
  ValToMinimise <- w*(SSE_Q) +
    (1-w)*(SSE_ET)
  return(ValToMinimise)
}

test_QET <- objfun_ET(par=c(600, 0, 100, 2, 0.15, 5),
                   mod=Cotter_mod_M,in_data=data.modis.cal)
test_QET
x01

# Constraint function
Constr_Mod_ET <- function(par, beta, in_data, mod) {
  # par is vector of initial values of parameters to be optimised
  # this includes beta
  # in_data are the input data
  # mod is a hydromad model
  mod_run <- update(mod, x1 = par[1], x2 = par[2], x3 = par[3], 
                    x4 = par[4], etmult=par[5])
  #browser()
  aET_fin <- aggregate(in_data$aET,list(date=in_data$et.period),sum)
  ET_fin <- aggregate(mod_run$U$ET,
                      list(date=in_data$et.period),
                      sum)
    #browser()

  error_Q <- coredata(fitted(mod_run)) - coredata(in_data$Q[101:length(data.modis.cal$Q)])
  error_ET <- coredata(ET_fin) - coredata(aET_fin)
  
  adj_error_ET <- zoo(sqrt(error_ET^2) - beta*8,
                      order.by = time(time(aET_fin)))
  adj_error_Q <- zoo(sqrt(error_Q^2) - beta,
                     order.by = time(in_data))
  adj_error_fin <- merge(adj_error_Q,adj_error_ET[13:length(adj_error_ET)], all=T)
  out <- apply(adj_error_fin, 1, sum, na.rm=T)
  
  return(na.omit(out))
}

x02 <- Constr_Mod_ET(par=c(600, 0, 100, 2, 0.15), beta=0.5,
                  mod=Cotter_mod,in_data=data.modis.cal)
plot(x02)
# Now invoke nloptr:
res1 <- nloptr( x0=c(600, 0.1, 100, 2, 0.15),
                eval_f=objfun_ET,
                lb = c(100, -30, 5, 0.6 , 0.01),
                ub = c(1500, 20, 500, 10, 0.5),
                eval_g_ineq = Constr_Mod_ET,
                opts = list("algorithm"="NLOPT_LN_COBYLA",
                            "xtol_rel"=1.0e-8, print_level = 1,
                            maxeval=1500),
                beta = 5,
                in_data = data.modis.cal,
                mod = Cotter_mod_M
)

#str(res1)
# plot the solution:
model_fit_M <- update(Cotter_mod_M, x1 = res1$solution[1], 
                    x2 = res1$solution[2], x3 = res1$solution[3], 
                    x4 = res1$solution[4], etmult=res1$solution[5])
# make a plot
plot(fitted(model_fit_M))
lines(data.cal$Q,col="red")
hist(data.cal$Q-fitted(model_fit_M),breaks=50)



# Show the ET calibration
aET_obs <- aggregate(data.modis.cal$aET,
                     list(date=data.modis.cal$et.period),sum)
ET_pred <- aggregate(model_fit_M$U$ET,
                     list(date=data.modis.cal$et.period),
                     sum)

epi.ccc(coredata(aET_obs),coredata(ET_pred))

plot.ET(caldata=data.modis.cal,model_fit_M)

# calculate statistics
nseStat(data.modis.cal$Q,fitted(model_fit_M))
hmadstat("rel.bias")(data.modis.cal$Q[101:length(data.modis.cal$Q)],
                     fitted(model_fit_M))
hmadstat("r.squared")(data.modis.cal$Q[101:length(data.modis.cal$Q)],
                      fitted(model_fit_M))
hmadstat("r.sq.sqrt")(data.modis.cal$Q[101:length(data.modis.cal$Q)],
                      fitted(model_fit_M))
hmadstat("r.sq.log")(data.modis.cal$Q[101:length(data.modis.cal$Q)],
                     fitted(model_fit_M))

# Now we want to run a loop over beta values
# plot the resulting r.sq.sqrt and rel.bias on a figure together
# with the result of the unconstrained optimisation

beta1 <- seq(1,10,by=1)
n <- length(beta1)

result <- data.frame(beta = beta1, mean_error = rep(0,n),
                     var_error = rep(0,n),
                     r.sq.srqt = rep(0,n), 
                     Linccc = rep(0,n),
                     mean_error_ET = rep(0,n),
                     var_error_ET = rep(0,n),
                     r.sq.srqt_ET = rep(0,n), 
                     Linccc_ET = rep(0,n))

unconstrResults <- rep(0,8)

unc_error <- data.modis.cal$Q[101:length(data.modis.cal$Q)] - 
  fitted(Cotter_fit_q_ET[["Fit ET and Q"]])
unconstrResults[1:2] <- c(mean(unc_error),var(unc_error))
unc_errorET <- coredata(aET_obs) - coredata(ET_pred)
unconstrResults[5:6] <- c(mean(unc_errorET),var(unc_errorET))



for (i in 1:length(beta1)) {
  res1 <- nloptr( x0=c(600, 0.1, 100, 2, 0.15),
                  eval_f=objfun_ET,
                  lb = c(100, -30, 5, 0.6 , 0.01),
                  ub = c(1500, 20, 500, 10, 0.5),
                  eval_g_ineq = Constr_Mod_ET,
                  opts = list("algorithm"="NLOPT_LN_COBYLA",
                              "xtol_rel"=1.0e-8, print_level = 1,
                              maxeval=5000),
                  beta = beta1[i],
                  in_data = data.modis.cal,
                  mod = Cotter_mod_M
  )
  
  model_fit_M <- update(Cotter_mod_M, x1 = res1$solution[1], 
                        x2 = res1$solution[2], x3 = res1$solution[3], 
                        x4 = res1$solution[4], etmult=res1$solution[5])
  # make a plot
  par(mfrow=c(1,2))
  plot(fitted(model_fit_M), xlab = "Date", ylab="Flow in mm",
       main = paste("Fitted and observed"))
  lines(data.cal$Q,col="red")
  hist(data.cal$Q-fitted(model_fit_M),breaks=50, 
       main = paste("beta =", beta1[i]) )
  par(mfrow=c(1,1))
  
  
  # ET calibration
  aET_obs <- aggregate(data.modis.cal$aET,
                       list(date=data.modis.cal$et.period),sum)
  ET_pred <- aggregate(model_fit_M$U$ET,
                       list(date=data.modis.cal$et.period),
                       sum)
  
  result[i,2] <- mean((data.modis.cal$Q[101:length(data.modis.cal$Q)]-
                       fitted(model_fit_M)))
  result[i,3] <- var((data.modis.cal$Q[101:length(data.modis.cal$Q)]-
                      fitted(model_fit_M)))
  result[i,4] <- hmadstat("r.sq.sqrt")(data.modis.cal$Q[101:length(data.modis.cal$Q)],
                                       fitted(model_fit_M))
  result[i,5] <- epi.ccc(data.modis.cal$Q[101:length(data.modis.cal$Q)],
                          fitted(model_fit_M))$rho.c$est
  # calculate statistics
  result[i,6] <- mean(coredata(aET_obs) - coredata(ET_pred))
  result[i,7] <- var(coredata(aET_obs) - coredata(ET_pred))
  result[i,8] <- hmadstat("r.sq.sqrt")(coredata(aET_obs),
                                       coredata(ET_pred))
  result[i,9] <- epi.ccc(coredata(aET_obs),
                         coredata(ET_pred))$rho.c$est
  
}

unconstrResults[3] <- summary(Cotter_fit_q_ET[["Fit ET and Q"]], 
                           items = "r.sq.sqrt")[9]
ET_pred_unconstr <- aggregate(Cotter_fit_q_ET[["Fit ET and Q"]]$U$ET,
                              list(date=data.modis.cal$et.period),
                              sum)
unconstrResults[7] <- hmadstat("r.sq.sqrt")(coredata(aET_obs),
                                            coredata(ET_pred_unconstr))
linsET <- epi.ccc(coredata(aET_obs),coredata(ET_pred))
unconstrResults[8] <- LinsET$rho.c$est


# now make plots of the resulting r.sq.sqrt relative to the unconstrained fit
par(mfrow=c(2,2))
plot(result[,1], result[,2], xlab="beta", ylab="mean error",
     ylim=c(-0.5,1), col="blue", lwd=5)
lines(beta1,rep(unconstrResults[1],n), col="red",lty=2)
plot(result[,1], result[,3], xlab="beta", ylab="var error",
     ylim=c(-0.2,1), col="blue", lwd=5)
lines(beta1,rep(unconstrResults[2],n), col="red",lty=2)
legend("topright",c("constrained","standard"),lty=c(NA,2),pch=c(1,NA),
       lwd=c(5,1),col=c("blue","red"))
plot(result[,1], result[,4], xlab="beta", ylab="r.sq.sqrt flow",
     ylim=c(-0.5,1), col="blue", lwd=5)
lines(beta1,rep(unconstrResults[3],n), col="red",lty=2)
plot(result[,1], result[,5], xlab="beta", ylab="Lin's ccc flow",
     ylim=c(0,1), col="blue", lwd=5)
lines(beta1,rep(linsccc$rho.c$est,n), col="red",lty=2)
par(mfrow=c(1,1))

# Do again for ET
par(mfrow=c(2,2))
plot(result[,1], result[,6], xlab="beta", ylab="mean error ET",
     ylim=c(-3,1), col="blue", lwd=5)
lines(beta1,rep(unconstrResults[5],n), col="red",lty=2)
plot(result[,1], result[,7], xlab="beta", ylab="var error ET",
     ylim=c(5,20), col="blue", lwd=5)
lines(beta1,rep(unconstrResults[6],n), col="red",lty=2)
legend("topright",c("constrained","standard"),lty=c(NA,2),pch=c(1,NA),
       lwd=c(5,1),col=c("blue","red"))
plot(result[,1], result[,8], xlab="beta", ylab="r.sq.sqrt ET",
     ylim=c(-2,1), col="blue", lwd=5)
lines(beta1,rep(unconstrResults[7],n), col="red",lty=2)
plot(result[,1], result[,9], xlab="beta", ylab="Lin's ccc ET",
     ylim=c(0,1), col="blue", lwd=5)
lines(beta1,rep(linsET$rho.c$est,n), col="red",lty=2)
par(mfrow=c(1,1))
