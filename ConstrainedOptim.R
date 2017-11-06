# testing constrained optimisation
# using norm(error - beta)
# automated version


require(hydromad)

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

# define a function to fit with Optim
# for this we need to have the objectivefunction being 
# a function of beta and the error

objfun <- function(par, in_data, mod, ...) {
  # par is vector of initial values of parameters to be optimised
  # this includes beta
  # in_data are the input data
  # mod is a hydromad model
  mod_run <- update(mod, x1 = par[1], x2 = par[2], x3 = par[3], 
                         x4 = par[4], etmult=par[5])
  #browser()
  error <- fitted(mod_run) - in_data$Q
  
  adj_error <- sqrt((sqrt(error^2) - par[6])^2)
  SSE <- sqrt(sum((error)^2,na.rm=T))/length(error)
  
  
  ValToMinimise <- SSE + sum((adj_error),na.rm=T)
  return(ValToMinimise)
}

test <- objfun(par=c(600, 0, 100, 2, 0.15, 5),
               mod=Cotter_mod,in_data=data.cal)

# Now write an optimisation routine that includes this
# Optimisation routine

# lower boundaries
# x1, x2, x3, x4, etmult, beta
lb <- c(100, -30, 5, 0.6 , 0.01, 0.01)
# upper boundaries
up <- c(1500, 20, 500, 10, 0.5, 0.05)

# initial values
par_in <- c(600, 0.1, 100, 2, 0.15, 0.02)

# use optim and L-BFGS-B
sol <- optim(par_in, objfun,method="L-BFGS-B", 
             in_data=data.cal, mod=Cotter_mod,
             lower = lb, upper = up)

# apply solution to model and plot
model_fit <- update(Cotter_mod, x1 = sol$par[1], x2 = sol$par[2], x3 = sol$par[3], 
                    x4 = sol$par[4], etmult=sol$par[5])
# make a plot
plot(fitted(model_fit))
lines(data.cal$Q,col="red")

sol$par
# beta is 0.1127...

# calculate statistics
nseStat(data.cal$Q,fitted(model_fit))
hmadstat("rel.bias")(data.cal$Q[101:length(data.cal$Q)],fitted(model_fit))
hmadstat("r.squared")(data.cal$Q[101:length(data.cal$Q)],fitted(model_fit))
hmadstat("r.sq.sqrt")(data.cal$Q[101:length(data.cal$Q)],fitted(model_fit))
hmadstat("r.sq.log")(data.cal$Q[101:length(data.cal$Q)],fitted(model_fit))

# this is slightly better in low flow fit and similar in r sq sqrt and has a lower bias
# That would be correct given the way we fit beta

# Same, but now also fit MODIS ET
# load functions from Satellite calibration
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


# Do the same using constrained optimisations
# for this we need to rewrite the whole optimisation function to include ET
# This requires using some of the functions from the satellite hydromad
# define a function to fit with Optim
# for this we need to have the objectivefunction being 
# a function of beta and the error

objfun <- function(par, in_data, mod, ...) {
  # par is vector of initial values of parameters to be optimised
  # this includes beta
  # in_data are the input data, which includes MODIS ET
  # mod is a hydromad model
  mod_run <- update(mod, x1 = par[1], x2 = par[2], x3 = par[3], 
                    x4 = par[4], etmult=par[5])
  #browser()
  # Q objective functions
  error_Q <- fitted(mod_run) - in_data$Q
  
  adj_error_Q <- abs(abs(error_Q) - par[6])
  SSE_Q <- sqrt(sum((error_Q)^2,na.rm=T))/length(error_Q)
  # ET objective function
  #browser()
  aET_fin <- aggregate(in_data$aET,list(date=in_data$et.period),sum)
  ET_fin <- aggregate(mod_run$U$ET,
                      list(date=in_data$et.period),
                      sum)
  #browser()
  error_ET <- coredata(ET_fin) - coredata(aET_fin)
  
  adj_error_ET <- sqrt((sqrt(error_ET^2) - par[6]))
  SSE_ET <- sqrt(sum((error_ET)^2,na.rm=T))/length(error_ET)
  
  w <- 0.5
  
  ValToMinimise <- w*(SSE_Q + sum((adj_error_Q),na.rm=T)) +
                   (1-w)*(SSE_ET + sum((adj_error_ET),na.rm=T))
  return(ValToMinimise)
}

test_QET <- objfun(par=c(600, 0, 100, 2, 0.15, 5),
               mod=Cotter_mod_M,in_data=data.modis.cal)

# Now write an optimisation routine that includes this
# Optimisation routine

# use optim and L-BFGS-B
sol <- optim(par_in, objfun,method="L-BFGS-B", 
             in_data=data.modis.cal, mod=Cotter_mod_M,
             lower = lb, upper = up)

# apply solution to model and plot
model_fit_M <- update(Cotter_mod_M, x1 = sol$par[1], x2 = sol$par[2], x3 = sol$par[3], 
                    x4 = sol$par[4], etmult=sol$par[5])
# make a plot
plot(fitted(model_fit_M))
lines(data.modis.cal$Q,col="red")

sol$par
# beta is 0.13254...

# calculate statistics
nseStat(data.modis.cal$Q,fitted(model_fit_M))
hmadstat("rel.bias")(data.modis.cal$Q[101:length(data.modis.cal$Q)],fitted(model_fit_M))
hmadstat("r.squared")(data.modis.cal$Q[101:length(data.modis.cal$Q)],fitted(model_fit_M))
hmadstat("r.sq.sqrt")(data.modis.cal$Q[101:length(data.modis.cal$Q)],fitted(model_fit_M))
hmadstat("r.sq.log")(data.modis.cal$Q[101:length(data.modis.cal$Q)],fitted(model_fit_M))

# Show the ET calibration
aET_obs <- aggregate(data.modis.cal$aET,
                     list(date=data.modis.cal$et.period),sum)
ET_pred <- aggregate(model_fit_M$U$ET,
                     list(date=data.modis.cal$et.period),
                     sum)

nseStat(coredata(aET_obs),coredata(ET_pred))

plot.ET(caldata=data.modis.cal,model_fit_M)

# improved the low flow and general flow calibration, reduced the peaks
# reduced the ET calibration
# this model is not spactial, so not really taking advantage of the ET data