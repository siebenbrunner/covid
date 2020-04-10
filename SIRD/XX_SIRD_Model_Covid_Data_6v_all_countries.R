######################################################################
######################################################################
# SIRD Model 
######################################################################
######################################################################

# based on: 
# https://rpubs.com/choisy/sir

# Remove all objects:
# rm(list=ls(all=TRUE))

# Server is true/false
sysinfo <- Sys.info()

Linux <- sysinfo[1] == "Linux"

pacman::p_load(deSolve)
pacman::p_load(bbmle)
pacman::p_load(optimx)


# Save all three files in the same folder:
Corona.path.to.R.Files  <-  paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/")


#####################################################################
# Step 1: Load Covid Data from Github:
#####################################################################

# source(paste(Corona.path.to.R.Files, "XX_Git_Hub_Dowload_Data.R", sep = ""))

covid$Country.Region <- as.character(covid$Country.Region)

#####################################################################
# Step 2: Load Support Functions:
#####################################################################

source(paste(Corona.path.to.R.Files, "XX_SIR_functions_4v.R", sep = ""))

#####################################################################
# Definie countries for the analysis.
#####################################################################

# not all countries have useful data:
countries_to_fit <- c(#"Albania",
                      # "Algeria",
                      # "Andorra",
                      # "Argentina",
                      # "Armenia",
                      # "Australia",
                      "Austria",
                      # "Bahrain",
                      # "Belgium",
                      # "Bosnia and Herzegovina",
                      # "Brazil",
                      # "Bulgaria",
                      # "Canada",
                      # "Chile",
                      # "Colombia",
                      # "Costa Rica",
                      # "Croatia",
                      # "Cyprus",
                      # "Denmark",
                      # "Ecuador",
                      # "Estonia",
                      # "Finland",
                      # "France",
                      # "Germany",
                      # "Greece",
                      # "Hungary",
                      # "Iceland",
                      # "India",
                      # "Indonesia",
                      # "Iran",
                      # "Ireland",
                      # "Israel",
                      # "Italy",
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
                      # "United Kingdom",
                      "US")

#####################################################################
# Define search grid
#####################################################################
vec_beta0 = seq(0.1,0.3,0.02)
vec_betaL = seq(0.03,0.9,0.02)
vec_gamma = seq(0.04,0.1,0.01)
vec_delta = seq(0.01,0.04,0.01)

vec_beta0 = c(0.3,0.25)
vec_betaL = 0.07
vec_gamma = 0.03132798
vec_delta = 0.03132343


# Make parallelization possible for windows

if (Linux == TRUE){
  pacman::p_load("doParallel")
  doParallel::registerDoParallel(length(countries_to_fit))
  
} else {
  pacman::p_load(doSNOW)
  cl<-makeCluster(2)
  registerDoSNOW(cl)
}

corona_models_SIRD <- list()


#####################################################################
# Parallelization:
#####################################################################

# corona_models_SIRD <- foreach (r = 1:length(countries_to_fit),
#                                   .errorhandling = 'pass') %dopar% 
  
  for (r in 1:length(countries_to_fit)) {
                                    
  corona <- subset(covid, covid$Country.Region == countries_to_fit[r])
  
  corona[,c("Confirmed",
            "Deaths",
            "Recovered")] <-  corona[,c("Confirmed",
                                        "Deaths",
                                        "Recovered")] / corona$Population 
  
  bln_model_found <- FALSE
  
  #####################################################################
  # Grid search:
  #####################################################################
  
  for (i_beta0 in vec_beta0) {
    for (i_betaL in vec_betaL) {
      for (i_gamma in vec_gamma) {
        for (i_delta in vec_delta) {
          
          #####################################################################
          # Use function: Optim:
          #####################################################################
          
          # starting_param_val <- c(0.3, 0.07, 0.03132798, 0.03132343)
          starting_param_val <- c(i_beta0, i_betaL, i_gamma, i_delta)
          
          ss_optim_sir <- try({optim(par = starting_param_val, 
                                fn = ss_SIR,
                                gr = NULL,
                                method = c("Nelder-Mead"),
                                control = list(maxit = 1000000,  pgtol = 1e-10),
                                data = corona,
                                lockdown_var = c("Lockdown_Lag_Rescale"),
                                I0 = c("Confirmed"),
                                times = c("Day"),
                                N = 1)},silent=TRUE)
          
          #####################################################################
          # Model selection: rules:
          # Store at least one model that converged
          # Always prefer a model that satisfies parameter constraints and has a gradient
          # Among those models, choose the one with the best fit
          #####################################################################
          if (!bln_model_found) {
            corona_models_SIRD[[r]] <- ss_optim_sir
            if (length(ss_optim_sir)>1) {
              bln_model_found <- TRUE
              } 
            }
          if (bln_model_found) {
              if ((any(corona_models_SIRD[[r]]$par < 0) && all(ss_optim_sir$par >= 0)) || 
                  (is.na(corona_models_SIRD[[r]]$counts["Gradient"]) && !is.na(ss_optim_sir$counts["Gradient"]))) {
                corona_models_SIRD[[r]] <- ss_optim_sir
              } else {
                if (ss_optim_sir$counts["function"] < corona_models_SIRD[[r]]$counts["function"]) {
                  corona_models_SIRD[[r]] <- ss_optim_sir
                }
              }
              
            }
        }
      }
    }
  }
  
  
  
                                    
}

names(corona_models_SIRD) <- as.character(countries_to_fit)






