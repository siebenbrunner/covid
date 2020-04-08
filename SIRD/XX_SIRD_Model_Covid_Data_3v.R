######################################################################
######################################################################
# SIRD Model 
######################################################################
######################################################################

# based on: 
# https://rpubs.com/choisy/sir

pacman::p_load(deSolve)
pacman::p_load(optimx)
pacman::p_load(bbmle)

# Save all three files in the same folder:
Corona.path.to.R.Files  <-  paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/")


#####################################################################
# Step 1: Load Covid Data from Github:
#####################################################################

source(paste(Corona.path.to.R.Files, "XX_Git_Hub_Dowload_Data.R", sep = ""))

#####################################################################
# Step 2: Load Support Functions:
#####################################################################

source(paste(Corona.path.to.R.Files, "XX_SIR_functions_2v.R", sep = ""))

#####################################################################
# Step 3: Load data for Austria
#####################################################################

corona <- subset(covid, covid$Country.Region == c("Austria"))



#####################################################################
# Step 4: Writing differential equations:
#####################################################################

predictions <- sir_1(beta_0 = 0.00000005, 
                     beta_L = 0.00000005, 
                     gamma = 0.0003,
                     delta = 0.0001,
                     S0 = corona$Population[1], 
                     I0 = 1,
                     D0 = 0,
                     R0 = 0, 
                     times = corona$Day,
                     lockdown = corona$Lockdown)

#####################################################################
# Use function: Optim:
#####################################################################

starting_param_val <- c(0.00000005, 0.00000005, 0.00000005, 0.00000005)

ss_optim_sir <- optim(par = starting_param_val, 
                      fn = ss_SIR,
                      method = c("Nelder-Mead"),
                      control = list(maxit = 1000000, pgtol = 1e-10)
                      )
ss_optim_sir

predictions <- sir_1(beta_0 = ss_optim_sir$par[1], 
                     beta_L = ss_optim_sir$par[2], 
                     gamma = ss_optim_sir$par[3],
                     delta = ss_optim_sir$par[4],
                     S0 = corona$Population[1], 
                     I0 = corona$Confirmed[1], 
                     R0 = 0,
                     D0 = 0,
                     times = corona$Day,
                     lockdown = corona$Lockdown)

predictions


#######################################################################
# Use optimx with lower and upper bound.
######################################################################
starting_param_val <- c(0.0005,0.0005,0,0)

starting_param_val <- ss_optim_sir$par


ss_optim_sir_lower_bound <- optimx(par = starting_param_val, 
                                  fn = ss_SIR,
                                  method = c("L-BFGS-B"),
                                  #all.methods = TRUE,
                                  lower = 0,
                                  upper = 0.005,
                                  control = list(maxit = 1000000, pgtol = 1e-09)
)

ss_optim_sir_lower_bound


predictions_lower_bound <- sir_1(beta_0 = ss_optim_sir_lower_bound$p1, 
                     beta_L = ss_optim_sir_lower_bound$p2, 
                     gamma = ss_optim_sir_lower_bound$p3,
                     delta = ss_optim_sir_lower_bound$p4,
                     S0 = corona$Population[1], 
                     I0 = corona$Confirmed[1], 
                     R0 = 0,
                     D0 = 0,
                     times = corona$Day,
                     lockdown = corona$Lockdown)

predictions_lower_bound


#####################################################################
# Use function: Optimx:
#####################################################################

ss_optim_sir2 <- optimx(par = starting_param_val, 
                        fn = ss_SIR)

predictions <- sir_1(beta = ss_optim_sir2$p1[1], 
                     gamma = ss_optim_sir2$p2[1],
                     delta = ss_optim_sir2$p3[1],
                     S0 = corona$Population[1], 
                     I0 = 1, 
                     R0 = 0,
                     D0 = 0,
                     times = corona$Day)

#####################################################################
# Maximum likelihood estimation with the bbmle package
#####################################################################

# does not work very well.

# starting_param_val <- list(beta = 0.04, 
#                            gamma = 0.01,
#                            delta = 0.04,
#                            sigma_I = 1, 
#                            sigma_R = 1,
#                            sigma_D = 1)
# 
# mLL(beta = 0.004, 
#     gamma =0.5,
#     delta = 0.004,
#     sigma_I = 1,
#     sigma_R = 1,
#     sigma_D = 1,
#     Day = corona$Day,
#     Confirmed = corona$Confirmed,
#     Recovered = corona$Recovered,
#     Deaths = corona$Deaths,
#     Population = corona$Population[1])
# 
# estimates <- mle2(minuslogl = mLL, 
#                   start = lapply(starting_param_val, log),
#                   method = "Nelder-Mead", 
#                   data = corona)
# 
# # The point estimates (we need to back transform):
# maxlik_estimation <- exp(coef(estimates))
# 
# predictions <- sir_1(beta = maxlik_estimation[1], 
#                      gamma = maxlik_estimation[2], 
#                      S0 = corona$Population[1], 
#                      I0 = 1, 
#                      R0 = 0, 
#                      times = corona$Day)
