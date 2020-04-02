#############################################################
# read data
#############################################################
require(reshape2)
require(dplyr)
require(wbstats)
require(ggplot2)
require(ggforce)
require(gridExtra)

set.url.data <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
covid = read.csv(url(set.url.data))
set.url.data.us <- c("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv")
covid.us = read.csv(url(set.url.data.us))
covid.us <- covid.us[,7:ncol(covid.us)]
covid.us <- covid.us %>% dplyr::rename(Province.State = Province_State, Country.Region = Country_Region) %>%
  dplyr::select(-Lat,-Long_,-Combined_Key)
covid <- dplyr::select(covid,-Lat,-Long)
covid <- reshape2::melt(covid, id.vars=c("Province.State","Country.Region"),variable.name="Day",value.name="Confirmed")
covid.us <- reshape2::melt(covid.us, id.vars=c("Province.State","Country.Region"),variable.name="Day",value.name="Confirmed")
covid <- rbind(covid,covid.us)
covid$Day <- as.numeric(covid$Day)
covid$Country.Region <- as.factor(covid$Country.Region)

# aggregate provinces
covid <- covid %>% dplyr::group_by(Country.Region,Day) %>% summarize(Confirmed = sum(Confirmed, na.rm = TRUE))

# unbalance panel: start at first case for each country
covid <- covid[covid$Confirmed > 0,]
covid <- droplevels(covid)

# renumber days
for (c in levels(covid$Country.Region)) {
  covid[covid$Country.Region==c,"Day"] <- 1:sum(covid$Country.Region==c)
}

# remove countries with less than x observations
covid <- dplyr::filter(covid,n()>20)

# remove countries with too little variation
to_keep <- covid %>% group_by(Country.Region) %>% 
  summarize(min_confirmed = min(Confirmed), max_confirmed = max(Confirmed)) %>% 
  filter(min_confirmed + 20 < max_confirmed)

covid <- covid[covid$Country.Region %in% to_keep$Country.Region,]

covid <- droplevels(covid)

# download population data (for logistic models)
population <- wb(indicator = "SP.POP.TOTL",country = "countries_only", startdate = 2015, enddate = 2020)
population <- filter(population,date==max(date)) %>% select("Country.Region"="country","Population"="value")
population[population$Country.Region=="United States",1] <- "US"
population$Country.Region[population$Country.Region == "Iran, Islamic Rep."] <- c("Iran")

#############################################################
# Fit exponential models
#############################################################

exponential_models <- data.frame(row.names = levels(covid$Country.Region))
exponential_models$Intercept = 0
exponential_models$Rate = 0
exponential_models$R2 = 0

for (c in levels(covid$Country.Region)) {
  exp_model <- lm(log(Confirmed) ~ Day,filter(covid,Country.Region==c))
  exponential_models[c,c("Intercept","Rate")] <- exp_model$coefficients
  exponential_models[c,"R2"] <- summary(exp_model)$r.squared
}

#############################################################
# Fit logistic models
#############################################################

# merge population data
covid <- inner_join(covid, population)
covid$Country.Region <- as.factor(covid$Country.Region)
covid <- droplevels(covid)

logistic_models <- data.frame(row.names = levels(covid$Country.Region))
logistic_models$Intercept = 0
logistic_models$Rate = 0
logistic_models$N = 0
logistic_models$current_cases = 0
logistic_models$Max_Infection_Rate = 0
logistic_models$R2 = 0

# use binary search to find the best N (warning: assuming single peak in the range 2*current_cases:population)
get_fit <- function(N,covid,c) {summary(lm(log(Confirmed/N/(1-Confirmed/N)) ~ Day,filter(covid,Country.Region==c)))$r.squared}

for (c in levels(covid$Country.Region)) {
  current_cases = max(covid[covid$Country.Region==c,"Confirmed"])
  population_size = max(covid[covid$Country.Region==c,"Population"])
  
  smaller <- current_cases*1.1
  larger <- population_size
  middle <- (current_cases+population_size)/2
  
  smaller_fit <- get_fit(smaller*2,covid,c)
  larger_fit <- get_fit(larger,covid,c)
  
  for (i in 1:20) {
    middle_fit <- get_fit(middle,covid,c)
    if (smaller_fit < larger_fit) {
      smaller <- middle
      smaller_fit <- middle_fit
      middle <- (middle+larger)/2
    } else {
      larger <- middle
      larger_fit <- middle_fit
      middle <- (smaller+middle)/2
    }
  }
  
  if (smaller_fit > middle_fit) {N = smaller} else {
    if (middle_fit > larger_fit) {N = middle} else {
      N = larger
    }
  }
  
  logistic_model <- lm(log(Confirmed/N/(1-Confirmed/N)) ~ Day,filter(covid,Country.Region==c))
  logistic_models[c,c("Intercept","Rate")] <- logistic_model$coefficients
  logistic_models[c,"R2"] <- summary(logistic_model)$r.squared
  logistic_models[c,"N"] <- N
  logistic_models[c,"current_cases"] <- current_cases
  logistic_models[c,"Max_Infection_Rate"] <- N / population_size
}

#############################################################
# Make plots
#############################################################

data_plot <- ggplot(data = covid) + 
  ggtitle("Covid cases") +
  facet_wrap(~ Country.Region, scales='free', ncol = 4) + 
  geom_point(aes(x=Day, y=Confirmed))

ggsave("Raw Data.png")

data_plot <- ggplot(data = covid) + 
  ggtitle("Covid cases (log scale)") +
  facet_wrap(~ Country.Region, scales='free', ncol = 4) + 
  geom_point(aes(x=Day, y=Confirmed)) + scale_y_log10()

ggsave("Log Data.png")

covid$Forecast <- 0

for (c in levels(covid$Country.Region)) {
  days <- max(covid[covid$Country.Region == c,"Day"])
  if (logistic_models[c,"Max_Infection_Rate"] > 0.9) {
    function_to_plot <- function(x) exp(exponential_models[c,"Intercept"] + exponential_models[c,"Rate"]*x)
    maxT <- round(days * 1.25)
    
  } else {
    function_to_plot <- function(x) logistic_models[c,"N"]/(1+exp(-(logistic_models[c,"Intercept"]+logistic_models[c,"Rate"]*x)))
    maxT <- days * 2
  }
  
  # model_plot <- ggplot(data = filter(covid,Country.Region==c)) + 
  #   ggtitle(c) + 
  #   geom_point( aes(x = Day, y = Confirmed)) +
  #   stat_function(fun = function_to_plot) + xlim(0,maxT)
  # 
  # ggsave(paste0("./plots/",c,".png"))
  
  covid[covid$Country.Region==c,"Forecast"] <- function_to_plot(1:days)
  covid = bind_rows(covid,data.frame(Country.Region = toString(c), Day = (days+1):maxT, Forecast = function_to_plot((days+1):maxT)))
}


pl <- lapply(1:ceiling(nrow(logistic_models)/6), function(x) {
  ggplot(data = covid) + 
    facet_wrap_paginate(~ Country.Region, scales='free', ncol = 2, nrow = 3,page=x) + 
    geom_point(aes(x=Day, y=Confirmed)) +
    geom_line(aes(x=Day,y=Forecast))
})

ml <- marrangeGrob(pl, nrow=1, ncol=1)
ggsave("Model Report.pdf",ml)

