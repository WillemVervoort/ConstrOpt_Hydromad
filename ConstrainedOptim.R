

require(hydromad)

# data and model
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
  
  adj_error <- abs(abs(error) - par[6])
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
lb <- c(100, -30, 5, 0.6 , 0.01, 0.1)
# upper boundaries
up <- c(1500, 20, 500, 10, 0.5, 10)

# initial values
par_in <- c(600, 0.1, 100, 2, 0.15, 5)

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

# calculate statistics
nseStat(data.cal$Q,fitted(model_fit))
hmadstat("rel.bias")(data.cal$Q[101:length(data.cal$Q)],fitted(model_fit))
hmadstat("r.squared")(data.cal$Q[101:length(data.cal$Q)],fitted(model_fit))
hmadstat("r.sq.sqrt")(data.cal$Q[101:length(data.cal$Q)],fitted(model_fit))
hmadstat("r.sq.log")(data.cal$Q[101:length(data.cal$Q)],fitted(model_fit))
