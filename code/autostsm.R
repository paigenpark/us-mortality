library(autostsm)
#library(data.table)

# read in life tables for all 50 states 
path = "/Users/paigepark/Desktop/Summer_2023"
states = c("AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DE", "FL", "GA", "HI", "ID", "IL", "IN", "IA", "KS", 
           "KY", "LA", "ME", "MD", "MA", "MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY", 
           "NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VT", "VA", "WA", 
           "WV", "WI", "WY")
data = list()
for (i in 1:length(states)){
  file = paste(states[i], "bltper_1x1.csv", sep = "_")
  data[[i]] = read.csv(paste(path, file, sep = "/"))
}

# create dataframe for all states life expectancy at birth
e0 <- lapply(data, function(x) x$ex[x$Age == 0])
full_state_names <- tolower(state.name[match(states, state.abb)]) # full state names to match map
names(e0) = full_state_names
e0_df <- do.call(rbind, e0)
colnames(e0_df) <- c(seq(1959,2020,1))
head(e0_df)

# creating state-specific structural time series models
CA <- as.data.frame(cbind(seq(1959,2020,1), e0_df[5,]))
colnames(CA) = c("year", "e0")
CA$year <- as.Date(paste(CA$year, 1, 1, sep = "-"), format = "%Y-%m-%d")
stsm = stsm_estimate(CA, verbose = TRUE, det_drift = TRUE)
stsm_fc = stsm_forecast(stsm, y = CA, n.ahead = floor(stsm$freq)*3, plot = TRUE)
stsm_fc = merge(stsm_fc, 
                stsm_detect_anomalies(stsm, y = CA, plot = TRUE), 
                by = "date", all = TRUE)
stsm_fc = merge(stsm_fc, 
                stsm_detect_breaks(stsm, y = CA, plot = TRUE, show_progress = TRUE), 
                by = "date", all = TRUE)
