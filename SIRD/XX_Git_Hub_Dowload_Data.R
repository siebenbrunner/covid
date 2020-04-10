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

set.url.data <- c("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
covid = read.csv(url(set.url.data))
covid <- reshape2::melt(covid, id.vars=c("Province.State","Country.Region","Lat","Long"),variable.name="Day",value.name="Confirmed")
covid$Country.Region <- as.factor(covid$Country.Region)

# aggregate provinces
covid <- covid %>% dplyr::group_by(Country.Region,Day) %>% summarize(Confirmed = sum(Confirmed, na.rm = TRUE))

#############################################################
# Covid 19 Deaths:
#############################################################

set.url.data <- c("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv")
covid_deaths = read.csv(url(set.url.data))
covid_deaths <- reshape2::melt(covid_deaths, 
                               id.vars=c("Province.State","Country.Region","Lat","Long"),
                               variable.name="Day",
                               value.name="Deaths")

covid_deaths$Country.Region <- as.factor(covid_deaths$Country.Region)

# aggregate provinces
covid_deaths <- covid_deaths %>% dplyr::group_by(Country.Region,Day) %>% summarize(Deaths = sum(Deaths, na.rm = TRUE))

#############################################################
# Covid 19 Recovered:
#############################################################

set.url.data <- c("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv")
covid_recovered = read.csv(url(set.url.data))
colnames(covid_recovered)[1] <- "Province"
covid_recovered <- reshape2::melt(covid_recovered, 
                                  id.vars=c("Province","Country.Region","Lat","Long"),
                                  variable.name="Day",
                                  value.name="Recovered")

covid_recovered$Country.Region <- as.factor(covid_recovered$Country.Region)

# aggregate provinces
covid_recovered <- covid_recovered %>% dplyr::group_by(Country.Region,Day) %>% summarize(Recovered = sum(Recovered, na.rm = TRUE))

#############################################################
# Merge Covid Data Sets
#############################################################

covid_2 <- left_join(covid, covid_deaths)

covid_3 <- left_join(covid_2, covid_recovered)

covid <- covid_3

#############################################################
# merge population and lockdown data
#############################################################

# download population data
population <- wb(indicator = "SP.POP.TOTL",country = "countries_only", startdate = 2015, enddate = 2020)
population <- filter(population,date==max(date)) %>% 
  select("CountryCode" = "iso3c", "Country.Region"="country","Population"="value")

population[population$Country.Region=="United States",1] <- "US"
population$Country.Region[population$Country.Region == "Iran, Islamic Rep."] <- c("Iran")
population[population$Country.Region=="Korea, Rep.",1] <- "Korea, South"

# download lockdown data
lockdowns <- read.csv(url("https://ocgptweb.azurewebsites.net/CSVDownload")) %>% select("CountryCode", "Date", "Lockdown" = "StringencyIndexForDisplay")

# merge lockdown and population data
lockdowns <- inner_join(population,lockdowns) %>%
  select("Country.Region", "Day" = "Date", "Population", "Lockdown")

# rename dates
for (i in 1:nrow(lockdowns)) {
  new_date <- paste0("X",gsub("(^|[^0-9])0+", "\\1", substr(lockdowns[i,"Day"],5,6), perl = TRUE),".") # Month
  new_date <- paste0(new_date,gsub("(^|[^0-9])0+", "\\1", substr(lockdowns[i,"Day"],7,8), perl = TRUE),".") # Day
  new_date <- paste0(new_date,substr(lockdowns[i,"Day"],3,4)) # Year
  lockdowns[i,"Day"] <- new_date
}

# merge with covid data
covid <- left_join(covid, lockdowns) %>% filter(!is.na(Population))
covid$Country.Region <- as.factor(covid$Country.Region)
covid <- droplevels(covid)

# remove NAs and add lag
covid$Lockdown_Lag <- 0
lag_length <- 10
covid <- as.data.frame(covid)
for (c in levels(covid$Country.Region)) {
  pos <- covid$Country.Region==c
  for (i in 1:sum(pos)) {
    if (is.na(covid[pos,"Lockdown"][i])) {
      if (i==1) {
        covid[pos,"Lockdown"][i] <- 0
      } else {
        covid[pos,"Lockdown"][i] <- covid[pos,"Lockdown"][i-1]
      }
    }
  }
  covid[pos,"Lockdown_Lag"] <- c(rep(0,lag_length),covid[pos,"Lockdown"][1:(sum(pos)-lag_length)])
}

# rescale lockdown per country:
max_measures_country <- covid %>% group_by(Country.Region) %>%
                        summarise_at(c("Lockdown"),
                                    list(~max(., na.rm = TRUE)))

for (c in levels(covid$Country.Region)){

  covid[covid$Country.Region == c, c("Lockdown_Rescale")] <-
      covid[covid$Country.Region == c, c("Lockdown")] /
        as.double(max_measures_country[max_measures_country$Country.Region == c, c("Lockdown")])

  covid[covid$Country.Region == c, c("Lockdown_Lag_Rescale")] <-
    covid[covid$Country.Region == c, c("Lockdown_Lag")] /
    as.double(max_measures_country[max_measures_country$Country.Region == c, c("Lockdown")])
}



#############################################################
# data filtering and formatting
#############################################################

# unbalance panel: start at first case for each country
covid <- covid[covid$Confirmed > 0,]

# renumber days
covid$Country.Region <- as.factor(covid$Country.Region)
for (c in levels(covid$Country.Region)) {
  covid[covid$Country.Region==c,"Day"] <- 1:sum(covid$Country.Region==c)
}

# remove countries with less than x observations
covid <- dplyr::filter(covid,n()>20)

# remove countries with not enough variation
to_keep <- covid %>% group_by(Country.Region) %>% 
  summarize(min_confirmed = min(Confirmed), max_confirmed = max(Confirmed)) %>% 
  filter(min_confirmed + 20 < max_confirmed)

covid <- covid[covid$Country.Region %in% to_keep$Country.Region,]

covid <- droplevels(covid)

covid$Day <- as.numeric(covid$Day)

