######################################################################
######################################################################
# SIRD Model 
######################################################################
######################################################################

# based on: 
# https://rpubs.com/choisy/sir

# Remove all objects:
rm(list=ls(all=TRUE))

# Server is true/false
sysinfo <- Sys.info()

R.Server <- sysinfo[1] == "Linux"

pacman::p_load(deSolve)
pacman::p_load(bbmle)


# Save all three files in the same folder:
Corona.path.to.R.Files  <-  paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/")


#####################################################################
# Step 1: Load Covid Data from Github:
#####################################################################

source(paste(Corona.path.to.R.Files, "XX_Git_Hub_Dowload_Data.R", sep = ""))

covid$Country.Region <- as.character(covid$Country.Region)

#####################################################################
# Step 2: Load Support Functions:
#####################################################################

source(paste(Corona.path.to.R.Files, "XX_SIR_functions_3v.R", sep = ""))

#####################################################################
# Definie countries for the analysis.
#####################################################################

# not all countries have useful data:
countries_to_fit <- c("Albania",
                      "Algeria",
                      "Andorra",
                      "Argentina",
                      "Armenia",
                      "Australia",
                      "Austria",
                      "Bahrain",
                      "Belgium",
                      "Bosnia and Herzegovina",
                      "Brazil",
                      "Bulgaria",
                      "Canada",
                      "Chile",
                      "Colombia",
                      "Costa Rica",
                      "Croatia",
                      "Cyprus",
                      "Denmark",
                      "Ecuador",
                      "Estonia",
                      "Finland",
                      "France",
                      "Germany",
                      "Greece",
                      "Hungary",
                      "Iceland",
                      "India",
                      "Indonesia",
                      "Iran",
                      "Ireland",
                      "Israel",
                      "Italy",
                      # "Japan",
                      # "Latvia",
                      # "Lithuania",
                      # "Luxembourg",
                      # "Netherlands",
                      # "New Zealand",
                      # "Norway",
                      # "Poland",
                      # "Portugal",
                      # "Romania",
                      # "Saudi Arabia",
                      # "Serbia",
                      # "Singapore",
                      # "Slovenia",
                      # "Spain",
                      # "Sweden",
                      # "Switzerland",
                      # "Turkey",
                      # "Ukraine",
                      "United Kingdom",
                      "US")

if (R.Server == TRUE){
  
  pacman::p_load("doParallel")
  doParallel::registerDoParallel(length(countries_to_fit))
  
}
# Make parallelization possible for windows
if (R.Server == FALSE){
  
  pacman::p_load(doSNOW)
  cl<-makeCluster(2)
  registerDoSNOW(cl)
}

corona_models_SIRD <- list()


#####################################################################
# Parallelization:
#####################################################################

corona_models_SIRD <- foreach (r = 1:length(countries_to_fit),
                                  .errorhandling = 'pass') %dopar% {
                                    
  corona <- subset(covid, covid$Country.Region == countries_to_fit[r])
  
  corona[,c("Confirmed",
            "Deaths",
            "Recovered")] <-  corona[,c("Confirmed",
                                        "Deaths",
                                        "Recovered")] / corona$Population 
  
  pacman::p_load(optimx)
  #####################################################################
  # Use function: Optim:
  #####################################################################
  
  starting_param_val <- c(0.25, 0.1, 0.07, 0.0007)
  
  ss_optim_sir <- optim(par = starting_param_val, 
                        fn = ss_SIR,
                        method = c("Nelder-Mead"),
                        control = list(maxit = 1000000, pgtol = 1e-10)
  )
  
  predictions <- sir_1(beta_0 = ss_optim_sir$par[1], 
                       beta_L = ss_optim_sir$par[2], 
                       gamma = ss_optim_sir$par[3],
                       delta = ss_optim_sir$par[4],
                       S0 = 1, 
                       I0 = corona$Confirmed[1], 
                       R0 = 0,
                       D0 = 0,
                       times = corona$Day,
                       lockdown = corona$Lockdown_Lag)
  
  predictions[,c("I", "R", "D")] <- predictions[,c("I", "R", "D")] * corona$Population
                                    
  predictions
  
  #######################################################################
  # Use optimx with lower and upper bound.
  ######################################################################
  
  ss_optim_sir_lower_bound <- optimx(par = starting_param_val, 
                                     fn = ss_SIR,
                                     method = c("L-BFGS-B"),
                                     #all.methods = TRUE,
                                     lower = 0,
                                     upper = 1,
                                     control = list(maxit = 1000000, pgtol = 1e-09)
  )
  
  
  
  
  predictions_lower_bound <- sir_1(beta_0 = ss_optim_sir_lower_bound$p1, 
                                   beta_L = ss_optim_sir_lower_bound$p2, 
                                   gamma = ss_optim_sir_lower_bound$p3,
                                   delta = ss_optim_sir_lower_bound$p4,
                                   S0 = 1, 
                                   I0 = corona$Confirmed[1], 
                                   R0 = 0,
                                   D0 = 0,
                                   times = corona$Day,
                                   lockdown = corona$Lockdown)
  
  predictions_lower_bound[,c("I", "R", "D")] <- predictions_lower_bound[,c("I", "R", "D")] * corona$Population[1]
  
  
  corona_models_SIRD[[r]] <- list(ss_optim_sir,
                                  predictions,
                                  ss_optim_sir_lower_bound,
                                  predictions_lower_bound)
  
  
                                    
}

names(corona_models_SIRD) <- as.character(countries_to_fit)


