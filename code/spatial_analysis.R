#### READ IN STATE DATA ####
library(tidyverse)
path = "/Users/paigepark/Desktop/us-mortality/data"
states = c("AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DE", "FL", "GA", "HI", "ID", "IL", "IN", "IA", "KS", 
           "KY", "LA", "ME", "MD", "MA", "MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY", 
           "NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VT", "VA", "WA", 
           "WV", "WI", "WY")
state_data = list()
for (i in 1:length(states)){
  file = paste(states[i], "bltper_1x1.csv", sep = "_")
  state_data[[i]] = read.csv(paste(path, file, sep = "/"))
}

# create dataframe for all states life expectancy at birth
s_e0 <- lapply(state_data, function(x) x$ex[x$Age == 0])
full_state_names <- tolower(state.name[match(states, state.abb)]) # full state names to match map
names(s_e0) = full_state_names
state_e0 <- do.call(rbind, s_e0)
colnames(state_e0) <- c(seq(1959,2020,1))
state_e0 <- as.data.frame(state_e0)

# create long version of dataframe
e0_df_long <- state_e0 %>%
  gather(key = "Year", value = "Life_Expectancy") 

e0_df_long <- cbind(states, e0_df_long)
e0_df_long$Year <- as.numeric(e0_df_long$Year)


#### CREATING NEIGHBOR OBJECTS ####
library(sf)
library(tigris)
library(spdep) # spatial stats package
library(rgdal) # supporter package for spdep

options(tigris_use_cache = TRUE)

# get US states shapefile from tigris
states_sf = states(class = "sf")

# merge spatial and life expectancy data
data_sf = left_join(states_sf, e0_df_long, by = c("STUSPS" = "states"), multiple = "all")

# remove alaska and hawaii for the purposes of clearer spatial analysis
data_sf = data_sf %>%
  filter(!(STUSPS %in% c("AK", "HI")))

# subset to 1 year to get neighbors 
states_unique_sf = data_sf %>% filter(Year == 1959)

# define neighbors for each state
nb = poly2nb(states_unique_sf)

# assign spatial weights to each relationship
weights = nb2listw(nb, style = "B")

# run tests for spatial autocorrelation in life expectancy for each year
moran_results = list()

for(y in unique(data_sf$Year)) {
  
  # get specific year
  year_data = data_sf %>% filter(Year == y)
  
  # compute moran's I
  moran_res = moran.test(year_data$Life_Expectancy, weights)
  
  # store results in list
  moran_results[[as.character(y)]] <- moran_res
}









