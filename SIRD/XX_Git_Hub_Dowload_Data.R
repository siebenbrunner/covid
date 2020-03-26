#############################################################
# read data
#############################################################
pacman::p_load(reshape2)
pacman::p_load(dplyr)
pacman::p_load(wbstats)
pacman::p_load(ggplot2)
pacman::p_load(ggforce)
pacman::p_load(gridExtra)
pacman::p_load(tvReg)

#############################################################
# Covid 19 Confirmed:
#############################################################

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

#############################################################
# Covid 19 Deaths:
#############################################################

set.url.data <- c("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Deaths.csv")
covid_deaths = read.csv(url(set.url.data))
colnames(covid_deaths) <- c(colnames(covid_deaths)[1:4],1:(ncol(covid_deaths)-4))
covid_deaths <- reshape2::melt(covid_deaths, 
                               id.vars=c("Province.State","Country.Region","Lat","Long"),
                               variable.name="Day",
                               value.name="Deaths")

covid_deaths$Day <- as.numeric(covid_deaths$Day)
covid_deaths$Country.Region <- as.factor(covid_deaths$Country.Region)

# aggregate provinces
covid_deaths <- covid_deaths %>% dplyr::group_by(Country.Region,Day) %>% summarize(Deaths = sum(Deaths, na.rm = TRUE))

#############################################################
# Covid 19 Recovered:
#############################################################

set.url.data <- c("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Recovered.csv")
covid_recovered = read.csv(url(set.url.data))
colnames(covid_recovered) <- c(colnames(covid_recovered)[1:4],1:(ncol(covid_recovered)-4))
covid_recovered <- reshape2::melt(covid_recovered, 
                               id.vars=c("Province.State","Country.Region","Lat","Long"),
                               variable.name="Day",
                               value.name="Recovered")

covid_recovered$Day <- as.numeric(covid_recovered$Day)
covid_recovered$Country.Region <- as.factor(covid_recovered$Country.Region)

# aggregate provinces
covid_recovered <- covid_recovered %>% dplyr::group_by(Country.Region,Day) %>% summarize(Recovered = sum(Recovered, na.rm = TRUE))

covid_2 <- left_join(covid, covid_deaths)

covid_3 <- left_join(covid_2, covid_recovered)

covid <- covid_3


# renumber days
for (c in levels(covid$Country.Region)) {
  covid[covid$Country.Region==c,"Day"] <- 1:sum(covid$Country.Region==c)
}

# remove countries with less than x observations
covid <- dplyr::filter(covid,n()>20)

# remove countries with no variation
to_keep <- covid %>% group_by(Country.Region) %>% 
  summarize(min_confirmed = min(Confirmed), max_confirmed = max(Confirmed)) %>% 
  filter(min_confirmed < max_confirmed)

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

# merge population data
covid <- inner_join(covid, population)
covid$Country.Region <- as.factor(covid$Country.Region)
covid <- droplevels(covid)

