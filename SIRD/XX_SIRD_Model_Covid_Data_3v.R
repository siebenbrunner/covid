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

predictions <- sir_1(beta_0 = 0.05, 
                     beta_L = 0.05, 
                     gamma = 0.003,
                     delta = 0.01,
                     S0 = corona$Population[1], 
                     I0 = 1,
                     D0 = 0,
                     R0 = 0, 
                     times = corona$Day)

#####################################################################
# Use function: Optim:
#####################################################################

starting_param_val <- c(0.0, 0.0, 0.0)

ss_optim_sir <- optim(starting_param_val, 
                      ss_SIR)
ss_optim_sir


predictions <- sir_1(beta = ss_optim_sir$par[1], 
                     gamma = ss_optim_sir$par[2],
                     delta = ss_optim_sir$par[3],
                     S0 = corona$Population[1], 
                     I0 = 1, 
                     R0 = 0,
                     D0 = 0,
                     times = corona$Day)

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

starting_param_val <- list(beta = 0.04, 
                           gamma = 0.01,
                           delta = 0.04,
                           sigma_I = 1, 
                           sigma_R = 1,
                           sigma_D = 1)

mLL(beta = 0.004, 
    gamma =0.5,
    delta = 0.004,
    sigma_I = 1,
    sigma_R = 1,
    sigma_D = 1,
    Day = corona$Day,
    Confirmed = corona$Confirmed,
    Recovered = corona$Recovered,
    Deaths = corona$Deaths,
    Population = corona$Population[1])

estimates <- mle2(minuslogl = mLL, 
                  start = lapply(starting_param_val, log),
                  method = "Nelder-Mead", 
                  data = corona)

# The point estimates (we need to back transform):
maxlik_estimation <- exp(coef(estimates))

predictions <- sir_1(beta = maxlik_estimation[1], 
                     gamma = maxlik_estimation[2], 
                     S0 = corona$Population[1], 
                     I0 = 1, 
                     R0 = 0, 
                     times = corona$Day)





