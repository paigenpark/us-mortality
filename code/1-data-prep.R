###### READ IN DATA - SET UP E0 DATASETS - FIT PRELIMINARY MODELS #####
library(here)
path = here("data")

#### ENTIRE US ####
# read in life tables for US
us_data = read.csv(paste(path, "USA_bltper_1x1.csv", sep = "/"))
us_e0 = us_data$ex[us_data$Age == 0]


#### CENSUS DIVSIONS ####
# read in life tables for Census divisions
divisions = c(1:9)
div_data = list()
for (i in 1:length(divisions)){
  file = paste("Div", divisions[i], "_bltper_1x1.csv", sep = "")
  div_data[[i]] = read.csv(paste(path, file, sep = "/"))
}

# create dataframe for census divisions life expectancy at birth
d_e0 <- lapply(div_data, function(x) x$ex[x$Age == 0])
full_div_names <- c("new england", "middle atlantic", "east north central", "west north central", 
                    "south atlantic","east south central", "west south central", "mountain", 
                    "pacific")
names(d_e0) = full_div_names
div_e0 <- do.call(rbind, d_e0)
colnames(div_e0) <- c(seq(1959,2020,1))

#### STATES ####
# read in life tables for all 50 states 
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
names(s_e0) = states
state_e0 <- do.call(rbind, s_e0)
colnames(state_e0) <- c(seq(1959,2020,1))

#### SAVE COMBINED STATE AND US DATASET - WHAT IS USED IN THE PAPER ####
combined_e0 <- rbind(USA = us_e0, state_e0)
write.csv(combined_e0, paste(path, "combined_e0.csv", sep = "/"))