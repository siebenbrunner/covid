#############################################################
# read data
#############################################################
require(reshape2)
require(dplyr)
require(wbstats)

set.url.data <- c("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv")
covid = read.csv(url(set.url.data))
colnames(covid) <- c(colnames(covid)[1:4],1:(ncol(covid)-4))
covid <- reshape2::melt(covid, id.vars=c("Province.State","Country.Region","Lat","Long"),variable.name="Day",value.name="Confirmed")
covid$Day <- as.numeric(covid$Day)
covid$Country.Region <- as.factor(covid$Country.Region)

# aggregate provinces
covid <- covid %>% dplyr::group_by(Country.Region,Day) %>% summarize(Confirmed = sum(Confirmed, na.rm = TRUE))

# unbalance panel: start at first case for each country
covid <- covid[covid$Confirmed > 0,]

# renumber days
for (c in levels(covid$Country.Region)) {
  covid[covid$Country.Region==c,"Day"] <- 1:sum(covid$Country.Region==c)
}

# remove countries with less than x observations
covid <- dplyr::filter(covid,n()>10)

# remove countries with no variation
to_keep <- covid %>% group_by(Country.Region) %>% 
  summarize(min_confirmed = min(Confirmed), max_confirmed = max(Confirmed)) %>% 
  filter(min_confirmed < max_confirmed)

covid <- covid[covid$Country.Region %in% to_keep$Country.Region,]

covid <- droplevels(covid)

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

# download population data
population <- wb(indicator = "SP.POP.TOTL",country = "countries_only", startdate = 2015, enddate = 2020)
population <- filter(population,date==max(date)) %>% select("Country.Region"="country","Population"="value")
population[population$Country.Region=="United States",1] <- "US"

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
  
  smaller <- current_cases*2
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

