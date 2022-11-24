library("pomp")
library("deSolve")



### ODE integrator
closed.sir.model <- function (t, x, params) {
  S = x[1]
  I = x[2]
  R = x[3]
  ## now extract the parameters
  beta <- params["beta"]
  gamma <- params["gamma"]
  N <- S+I+R
  ## now code the model equations
  dSdt <- -beta*S*I/N
  dIdt <- beta*S*I/N-gamma*I
  dRdt <- gamma*I
  ## combine results into a single vector
  dxdt <- c(dSdt,dIdt,dRdt)
  ## return result as a list!
  list(dxdt)
}


### Transition equation

tran <- discrete_time(
  function(beta, S, I, R, sdlog, ...){
    
    # Run integration
    ode_result <- ode(
      func = closed.sir.model,
      y = c(S, I, R),
      times = 1,
      parms = beta
    )
    
    
    
    c(beta = beta_next, S = S_next, I = S_next, R = R_next)
  },
  
  delta.t = 1
)



### Observation equation

meas_pois_interval <- function(t, beta, S, I, R, T, Y, y, g_shape, g_scale, ..., log){
  if(t > 1){
    lambda = beta[t]*S[t]*I[t]
    y*log(lambda)-lambda
  }else{
    dpois(round(y), lambda = 1)
  }
}


### Run particle filter

run_pf <- function(d, sdlog){
  init_vals = c(sdlog = sdlog,
                g_shape = sid_covid_shape, g_scale = sid_covid_scale)
  
  pf <- pfilter(
    Np=1000,
    times = "idx",
    t0 = 1,
    data = d,
    params = init_vals,
    rinit = function(x0, ...){c(x=1.5, s = 0.1)},
    rprocess = tran,
    dmeasure = meas_pois_interval,
    T=t,
    Y = d$y,
    statenames = c("x"),
    filter.mean = TRUE,
    filter.traj = TRUE,
  )
  
  return(result)
}


