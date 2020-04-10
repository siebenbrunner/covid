######################################################################
######################################################################
# Support functions for SIR and SIRD Model
######################################################################
######################################################################


sir_1 <- function(beta_0, beta_L, gamma, delta, S0, I0, R0, D0, times, lockdown) {
  require(deSolve) # for the "ode" function
  
  # the differential equations:
  sir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      
      #beta <- ((1-lockdown[time]) * beta_0 + lockdown[time] * beta_L)
      #R >= 0
      #D >= 0
      #S >= 0
      dS <- -((1-lockdown[time]) * beta_0 + lockdown[time] * beta_L) * I * S
      dI <-  ((1-lockdown[time]) * beta_0 + lockdown[time] * beta_L)  * I * S - gamma * I - delta * I
      dR <-  gamma * I
      dD <-  delta * I
      return(list(c(dS, dI, dR, dD)))
    })
  }
  
  # the parameters values:
  parameters_values <- c(beta_0  = beta_0, beta_L = beta_L, gamma = gamma, delta = delta)
  
  # the initial values of variables:
  initial_values <- c(S = S0, I = I0, R = R0, D = D0)
  
  # solving
  out <- ode(y = initial_values, 
             times = times, 
             func = sir_equations, 
             parms = parameters_values,
             method = c("lsode"))
  
  # returning the output:
  as.data.frame(out)
}

ss_sir <- function(beta_0, 
                   beta_L, 
                   gamma, 
                   delta, 
                   data,
                   I0,
                   times,
                   lockdown_var,
                   N = 1){
 
  data = data
  
  I0 <- data[1, I0]
  
  times <- data[,times]
  
  lockdown <- data[, lockdown_var]
  
  N = N
  
  predictions <- sir_1(beta_0 = beta_0, beta_L = beta_L, gamma = gamma, delta = delta,   # parameters
                       S0 = N - I0, I0 = I0, R0 = 0,
                       D0 = 0,# variables' intial values
                       times = times, lockdown = lockdown)                # time points
  return(sum((predictions$I - data$Confirmed)^2) + 
         sum((predictions$R - data$Recovered)^2) +
         sum((predictions$D - data$Deaths)^2))
}


ss_SIR <- function(x,
                   data,
                   I0,
                   times,
                   lockdown_var,
                   N) {
  ss_sir(beta_0 = x[1], beta_L = x[2], gamma = x[3], delta = x[4], data = data, I0 =I0, times = times, lockdown_var = lockdown_var, N = N)
}

