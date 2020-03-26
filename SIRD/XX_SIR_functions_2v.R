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
      
      dS <- -((1-lockdown[time]) * beta_0 + lockdown[time] * beta_L) * I * S
      dI <-  ((1-lockdown[time]) * beta_0 + lockdown[time] * beta_L) * I * S - gamma * I - delta * I
      dR <-  gamma * I
      dD <-  delta * I
      return(list(c(dS, dI, dR, dD)))
    })
  }
  
  # the parameters values:
  parameters_values <- c(beta_0  = beta, beta_L = beta_L, gamma = gamma, delta = delta)
  
  # the initial values of variables:
  initial_values <- c(S = S0, I = I0, R = R0, D = D0)
  
  # solving
  out <- ode(initial_values, times, sir_equations, parameters_values)
  
  # returning the output:
  as.data.frame(out)
}


ss_sir <- function(beta_0, beta_L, gamma, delta, data = corona, N = corona$Population[1]) {
  I0 <- data$Confirmed[1]
  times <- data$Day
  lockdown <- data$Lockdown
  predictions <- sir_1(beta_0 = beta_0, beta_L, beta_L, gamma = gamma, delta = delta,   # parameters
                       S0 = N - I0, I0 = I0, R0 = 0,
                       D0 = 0,# variables' intial values
                       times = times, lockdown = lockdown)                # time points
 
  sum((predictions$I - data$Confirmed)^2) + 
    sum((predictions$R - data$Recovered)^2) +
    sum((predictions$D - data$Deaths)^2)
}


ss_SIR <- function(x) {
  ss_sir(beta = x[1], gamma = x[2], delta = x[3])
}

# Optimization with nls: not completed
ss_sir_nls <- function(beta, gamma, data = corona, N = corona$Population[1]) {
  I0 <- data$Confirmed[1]
  times <- data$Day
  predictions <- sir_1(beta = beta, gamma = gamma,   # parameters
                       S0 = N - I0, I0 = I0, R0 = 0, # variables' intial values
                       times = times)                # time points
  
  sum((predictions$I - data$Confirmed)^2) + sum((predictions$R - data$Recovered)^2)
}

ss_SIR_nls <- function(x) {
  ss_sir(beta = x[1], gamma = x[2], delta = x[3])
}


mLL <- function(beta, gamma, delta, sigma_I, sigma_R, sigma_D, Day, Confirmed, Recovered, Deaths, Population) {
  beta <- exp(beta) # to make sure that the parameters are positive
  gamma <- exp(gamma)
  delta <- exp(delta)
  sigma_I <- exp(sigma_I)
  sigma_R <- exp(sigma_R)
  simga_D <- exp(sigma_D)
  I0 <- Confirmed[1] # initial number of infectious
  observations <- Confirmed[-1] # the fit is done on the other data points
  recovered <- Recovered[-1]
  death <- Deaths[-1]
  N <- Population[1]
  predictions <- sir_1(beta = beta, gamma = gamma, delta = delta,
                       S0 = N - I0, I0 = I0, R0 = 0, D0 = 0, times = Day)
  pred_infected <- predictions$I[-1] # removing the first point too
  pred_recovered <- predictions$R[-1]
  pred_death <- predictions$D[-1]
  # returning minus log-likelihood:
  # this double sum might be incorrect?
  - (sum(dnorm(x = observations, mean = pred_infected, sd = sigma_I, log = TRUE)) +
     sum(dnorm(x = recovered, mean = pred_recovered, sd = sigma_R, log = TRUE)) +
     sum(dnorm(x = death, mean = pred_death, sd = sigma_D, log = TRUE)) )   
}



