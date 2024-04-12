### SET-UP ###
# some packages
library(tidyverse) # version 
library(whereami) # version
library(ggExtra)
library(gridExtra)

# path
script_dir <- whereami()
path <- file.path(script_dir, "../../data") # assumes this R script is in a code folder, 
                                            # this goes up two levels from current script location  
                                            # into the main directory, then back down into data folder
path = normalizePath(path)




### ADDING WRAPPER FUNCTIONS ###

### STRUCTTS ###
# random walk with drift (RWD) specification
rwd <- function(y){
  sts = StructTS(x = y, type = "trend", fixed = c(NA, 0, 0))
  d.hat = sts$model$a[2]
  var_innov.hat = sts$coef["level"]
  return(list(d = d.hat,
              var_innov = var_innov.hat))
}

# structural time series (STS) specification = RWD + observation error
rwd_with_obs_error <- function(y, ..., maxit = 5000){
  sts = StructTS(x = y, type = "trend", fixed = c(NA, 0, NA), ...,
                 optim.control = list(maxit = maxit))
  d.hat = sts$model$a[2]
  var_innov.hat = sts$coef["level"]
  var_obs.hat = sts$coef["epsilon"]
  names(d.hat) <- names(var_innov.hat) <- names(var_obs.hat) <- ""
  
  return(c(d = d.hat,
           var_innov = var_innov.hat,
           var_obs = var_obs.hat))
}

# function that estimates drift, variance, and AIC and BIC for RWD
rwd_aic_bic <- function(y){
  sts = StructTS(x = y, type = "trend", fixed = c(NA, 0, 0))
  d.hat = sts$model$a[2]
  var_innov.hat = sts$coef["level"]
  
  # get AIC and BIC (have to do this manually)
  aic_val = -2*sts$loglik + 2
  bic_val = -2*sts$loglik + log(length(y))
  
  return(c(d = d.hat,
           var_innov = var_innov.hat,
           aic = aic_val,
           bic = bic_val))
}

# function that estimates drift, variances, AIC and BIC for STS 
rwd_with_obs_aic_bic <- function(y, ..., maxit = 4000){
  sts = StructTS(x = y, type = "trend", fixed = c(NA, 0, NA), ...,
                 optim.control = list(maxit = maxit))
  d.hat = sts$model$a[2]
  var_innov.hat = sts$coef["level"]
  var_obs.hat = sts$coef["epsilon"]
  names(d.hat) <- names(var_innov.hat) <- names(var_obs.hat) <- ""
  
  # get AIC and BIC (have to do this manually)
  aic_val = -2*sts$loglik + 2*2
  bic_val = -2*sts$loglik + 2*log(length(y))
  
  return(c(d = d.hat,
           var_innov = var_innov.hat,
           var_obs = var_obs.hat,
           aic = aic_val,
           bic = bic_val))
}

# from my research, it seems that StructTS is the least flexible method for structural 
# time series estimation. If there is some kind of issue with the optimization leading to 
# the 0s, it is hard to adjust anything substantial in the StructTS specifications.
# also, based on my software simulations, it produces more estimates of 0 than it should
# therefore, we also set up function that use the MARSS() package




### MARSS ###
## BFGS ##
rwd_with_obs_error_bfgs <- function(y, ...)
{
  library(MARSS)
  mod.list.2 <- list(B = matrix(1),
                     U = matrix("d"),
                     Q = matrix("q"),
                     Z = matrix(1),
                     A = matrix(0),
                     R = matrix("r"),
                     x0 = matrix("mu"),
                     tinitx = 0)
  out2.marss = MARSS(y, model = mod.list.2, control=list(maxit=1000), method= "BFGS", ...)
  
  if (out2.marss$convergence == 0) {
    d = coef(out2.marss)$U[1,1]
    var_innov = coef(out2.marss)$Q[1,1]
    var_obs = coef(out2.marss)$R[1,1]
  } else {
    d = NA
    var_innov = NA
    var_obs = NA
  }
  
  names(d) <- names(var_innov) <- names(var_obs) <- ""
  return(c(d= d,
           var_innov= var_innov,
           var_obs = var_obs))
}

# for use with later AIC/AICc function - RWD specification
rwd_bfgs <- function(y, ...)
{
  library(MARSS)
  mod.list.2 <- list(B = matrix(1),
                     U = matrix("d"),
                     Q = matrix("q"),
                     Z = matrix(1),
                     A = matrix(0),
                     R = matrix(0), # set the observation var to 0 for rwd
                     x0 = matrix("mu"),
                     tinitx = 0)
  out2.marss = MARSS(y, model = mod.list.2, control=list(maxit=3000), method= "BFGS", ...)
  
  return(out2.marss)
}


# for use with later AIC/AICc function - STS specification
rwd_obs_error_bfgs_out <- function(y, ...)
{
  library(MARSS)
  mod.list.2 <- list(B = matrix(1),
                     U = matrix("d"),
                     Q = matrix("q"),
                     Z = matrix(1),
                     A = matrix(0),
                     R = matrix("r"),
                     x0 = matrix("mu"),
                     tinitx = 0)
  out2.marss = MARSS(y, model = mod.list.2, control=list(maxit=3000), method= "BFGS", ...)

  return(out2.marss)
}

# AIC/AICc function
get_aic_aicc = function(y, ...) {
  rwd_result = rwd_bfgs(y)
  obserr_result = rwd_obs_error_bfgs_out(y)
  
  return(c(rwd_result$AIC, obserr_result$AIC, rwd_result$AICc, obserr_result$AICc))
}

# this function needs work - I can get AIC and AICc from the MARSS output, but not the bootstrap AICs
# that are supposedly better for time series. This function can help me do that but its not working 
# properly right now
get_aicbp_aicbb <- function(y, ...) {
  rwd_result = rwd_bfgs(y)
  obserr_result = rwd_obs_error_bfgs_out(y)
  rwd_aic = MARSSaic(rwd_result, output = c("AICbp", "AICbb"))
  obserr_aic = MARSSaic(obserr_result, output = c("AICbp", "AICbb"))
  
  return(c(rwd_aic, obserr_aic))
}



## KEM ##
rwd_with_obs_error_kem <- function(y, ...)
{
  library(MARSS)
  mod.list.2 <- list(B = matrix(1),
                     U = matrix("d"),
                     Q = matrix("q"),
                     Z = matrix(1),
                     A = matrix(0),
                     R = matrix("r"),
                     x0 = matrix("mu"),
                     tinitx = 0)
  out2.marss = MARSS(y, model = mod.list.2, control=list(maxit=1000), method= "kem", ...)
  
  if (out2.marss$convergence == 0) {
    d = coef(out2.marss)$U[1,1]
    var_innov = coef(out2.marss)$Q[1,1]
    var_obs = coef(out2.marss)$R[1,1]
  } else {
    d = NA
    var_innov = NA
    var_obs = NA
  }
  
  names(d) <- names(var_innov) <- names(var_obs) <- ""
  return(c(d= d,
           var_innov= var_innov,
           var_obs = var_obs))
}

# No model selection criteria functions created yet for KEM method, could do this later if needed 






###### READ IN DATA - SET UP E0 DATASETS - FIT PRELIMINARY MODELS #####

#### ENTIRE US ####
# read in life tables for US
us_data = read.csv(paste(path, "USA_bltper_1x1.csv", sep = "/"))
us_e0 = us_data$ex[us_data$Age == 0]
no2020 = us_e0[1:61]
plot(us_e0)

# fitting the structural time series models to aggregate US data
# with StructTS
us_sts = rwd_with_obs_error(us_e0)
# with MARSS (BFGS and kem)
us_bfgs = rwd_with_obs_error_bfgs(us_e0) 
us_kem = rwd_with_obs_error_kem(us_e0)
# this is not converging
    # user guide to MARSS suggests that if one of the elements on the diagonal 
    # of Q or R are going to 0 (are degenerate) then it will take the EM algorithm forever
    # to get to 0 
rbind(us_sts, us_bfgs, us_kem)




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

# MARSS BFGS
div_bfgs = matrix(NA, nrow=9, ncol=3)
div_bfgs = as.data.frame(apply(div_e0, 1, rwd_with_obs_error_bfgs))

# MARSS KEM
div_kem = matrix(NA, nrow=9, ncol=3)
div_kem = as.data.frame(apply(div_e0, 1, rwd_with_obs_error_kem))

# StructTS
div_sts = matrix(NA, nrow=9, ncol=3)
div_sts = as.data.frame(apply(div_e0, 1, rwd_with_obs_error))


# get names of divisions with var_obs = 0
div_no_obs = which(div_sts[3,]==0)
div_no_obs = colnames(div_sts)[div_no_obs]
div_no_obs



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
full_state_names <- tolower(state.name[match(states, state.abb)]) # full state names to match map
names(s_e0) = full_state_names
state_e0 <- do.call(rbind, s_e0)
colnames(state_e0) <- c(seq(1959,2020,1))

# create data with 2020 removed for testing impact of 2020 mortality on models 
s_no2020 = state_e0[, -62]

# MARSS kem
state_kem = matrix(NA, nrow=50, ncol=3)
state_kem = as.data.frame(apply(state_e0, 1, rwd_with_obs_error_kem))
  # results in 6 states that don't converge when used with "kem" method
  # results in 8 states that don't converge when run on data excluding 2020

# MARSS BFGS
state_bfgs = matrix(NA, nrow=50, ncol=3)
state_bfgs = as.data.frame(apply(state_e0, 1, rwd_with_obs_error_bfgs))
  # results in 6 states with near 0 observation variance 
  # results in 6 states with 0 observation variance when 2020 is excluded

# StructTS
state_sts = matrix(NA, nrow=50, ncol=3)
state_sts = as.data.frame(apply(state_e0, 1, rwd_with_obs_error))
    # results in 13 states with 0 observation variance
    # results in 10 states with 0 observation variance when using data excluding 2020
    # help page for StructTS says that it's not uncommon for 0's to occur
   

# get names of states with var_obs = 0 
state_no_obs_sts = which(state_sts[3,]==0)
state_no_obs_sts = colnames(state_sts)[state_no_obs_sts]
state_no_obs_sts

state_no_obs_kem = which(is.na(state_kem[3,]))
state_no_obs_kem = colnames(state_kem)[state_no_obs_kem]
state_no_obs_kem

state_no_obs_bfgs = which(state_bfgs[3,]<=0.0000000001)
state_no_obs_bfgs = colnames(state_bfgs)[state_no_obs_bfgs]
state_no_obs_bfgs

## BAR PLOT OF 0 OBS_VAR BY METHOD ##


barchart_data <- data.frame(
  state = full_state_names,
  MARSS_BFGS = full_state_names %in% state_no_obs_bfgs,
  StructTS = full_state_names %in% state_no_obs_sts
)

# Reshape data to long format
long_data <- pivot_longer(barchart_data, cols = c(MARSS_BFGS, StructTS), names_to = "method", values_to = "has_no_obs")

# Filter to keep only TRUE values
long_data <- long_data[long_data$has_no_obs == TRUE, ]

# Create a stacked bar chart
ggplot(long_data, aes(x = method, fill = state)) +
  geom_bar() +
  labs(y = "Count of States with Estimate of Zero for Observation Variance", fill = "State", x=NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust the x-axis labels if necessary

##### GENERATING TABLES OF DRIFT, INNOV SD, OBS SD, SAMPLING SD, & SHOCK SD #####
library(knitr)
library(kableExtra)

# get estimate for N to use in sampling SD approximation: I use the population estimates from the mid-point year of the time series (1989) 
# if population growth is roughly linear during this time span that would be a reasonable thing to do

# load in population data from the Census (data includes pop by state from 1985 to 1989)
pop_89 = read.csv(paste(path, "89_pop.csv", sep = "/"), header = TRUE, row.names = 1)
pop_89 = pop_89["Jul.89"] # filter to just year 1989

pop_samp_sd = 150/sqrt(pop_89) # create sampling standard deviation estimates using N from 89 (approximation given in Hanley 2022)
samp_sd = subset(samp_sd, !(rownames(samp_sd) %in% c("US", "DC"))) # get rid of observations not in USMDB 

# construct StructTS table 
state_sts_table = t(state_sts)
state_sts_table = as.data.frame(cbind(state_sts_table, samp_sd))
colnames(state_sts_table)[4] = "samp_sd"
state_sts_table = state_sts_table %>%
  mutate(shock_sd = var_obs - samp_sd^2) %>% # get shock sd through additive property of var
  mutate(var_innov = sqrt(var_innov), var_obs = sqrt(var_obs), shock_sd = sqrt(shock_sd)) # changing all of the variances to SDs
  
colnames(state_sts_table) = c("drift", "process_SD", "observation_SD", "sampling_SD", "shock_SD")
state_sts_table$shock_SD <- ifelse(state_sts_table$shock_SD < 0, NA, state_sts_table$shock_SD)

# round 
library(scales)
state_sts_table_round <- data.frame(lapply(state_sts_table, function(x) if(is.numeric(x)) number(x, accuracy = 0.00000001) else x))

# convert cells with 0 to red 
state_sts_table_red = state_sts_table_round
format_cell = function(cell_value) {
  if (is.na(cell_value)) {
    return(cell_spec("NA", "latex", background = "gray"))
  } else {
    cell_spec(cell_value, "latex", background = ifelse(cell_value == "0.00000000", "#FF9999", "white"))
  }
}
  
state_sts_table_red[] = as.data.frame(sapply(state_sts_table_round, function(col) {
  sapply(col, format_cell) 
  }))

rownames(state_sts_table_red) <- rownames(state_sts_table)

# get latex code for StructTS
kable(state_sts_table_red, caption = "StructTS Results", format="latex", escape = FALSE) %>%
  kable_styling()


# construct MARSS BFGS table
state_bfgs_table = t(state_bfgs)
state_bfgs_table = as.data.frame(cbind(state_bfgs_table, samp_sd))
state_bfgs_table = state_bfgs_table %>%
  mutate(shock_sd = var_obs - Jul.89^2) %>% # get shock sd through additive property of var
  mutate(var_innov = sqrt(var_innov), var_obs = sqrt(var_obs), shock_sd = sqrt(shock_sd)) # changing all of the variances to SDs

colnames(state_bfgs_table) = c("drift", "process_SD", "observation_SD", "sampling_SD", "shock_SD")
state_bfgs_table$shock_SD <- ifelse(state_bfgs_table$shock_SD < 0, NA, state_bfgs_table$shock_SD)

# round 
library(scales)
state_bfgs_table_round <- data.frame(lapply(state_bfgs_table, function(x) if(is.numeric(x)) number(x, accuracy = 0.00000001) else x))

# convert cells with 0 to red 
state_bfgs_table_red = state_bfgs_table_round
format_cell = function(cell_value) {
  if (is.na(cell_value)) {
    return(cell_spec("NA", "latex", background = "gray"))
  } else {
    # cell_spec(cell_value, "latex", background = ifelse(cell_value <= "0.00010000", "#FF9999", "white")) 
    cell_spec(cell_value, "latex", background = "white")# changing this code to capture near 0 results in red 
  }
}

state_bfgs_table_red[] = as.data.frame(sapply(state_bfgs_table_round, function(col) {
  sapply(col, format_cell) 
}))

rownames(state_bfgs_table_red) <- rownames(state_bfgs_table)

# get latex code for MARSS BFGS
kable(state_bfgs_table_red, caption = "MARSS BFGS Results", format="latex", escape = FALSE) %>%
  kable_styling()


# construct MARSS kem table
state_kem_table = t(state_kem)
state_kem_table = as.data.frame(cbind(state_kem_table, samp_sd))
state_kem_table = state_kem_table %>%
  mutate(shock_sd = var_obs - Jul.89^2) %>% # get shock sd through additive property of var
  mutate(var_innov = sqrt(var_innov), var_obs = sqrt(var_obs), shock_sd = sqrt(shock_sd)) # changing all of the variances to SDs

colnames(state_kem_table) = c("drift", "process_SD", "observation_SD", "sampling_SD", "shock_SD")
state_kem_table$shock_SD <- ifelse(state_kem_table$shock_SD < 0, NA, state_kem_table$shock_SD)

# round 
library(scales)
state_kem_table_round <- data.frame(lapply(state_kem_table, function(x) if(is.numeric(x)) number(x, accuracy = 0.00000001) else x))

# convert cells with 0 to red 
state_kem_table_red = state_kem_table_round
format_cell = function(cell_value) {
  if (is.na(cell_value)) {
    return(cell_spec("NA", "latex", background = "gray"))
  } else {
    cell_spec(cell_value, "latex", background = ifelse(cell_value == "0.00000000", "#FF9999", "white"))
  }
}

state_kem_table_red[] = as.data.frame(sapply(state_kem_table_round, function(col) {
  sapply(col, format_cell) 
}))

rownames(state_kem_table_red) <- rownames(state_kem_table)

# get latex code for MARSS kem
kable(state_kem_table_red, caption = "MARSS KEM Results", format="latex", escape = FALSE) %>%
  kable_styling()



# VARIANCE TABLES 

# construct StructTS table with variances 
state_sts_table = t(state_sts)
state_sts_table = as.data.frame(cbind(state_sts_table, samp_sd))
state_sts_table = state_sts_table %>%
  mutate(var_samp = Jul.89^2) %>%
  mutate(var_shock = var_obs - var_samp)  # get shock variance through additive property of var
state_sts_table = select(state_sts_table, -Jul.89)  

colnames(state_sts_table) = c("drift", "process_var", "observation_var", "sampling_var", "shock_var")
state_sts_table$shock_var <- ifelse(state_sts_table$shock_var < 0, NA, state_sts_table$shock_var)

# round 
library(scales)
state_sts_table_round <- data.frame(lapply(state_sts_table, function(x) if(is.numeric(x)) number(x, accuracy = 0.001) else x))

# convert cells with 0 to red 
state_sts_table_red = state_sts_table_round
format_cell = function(cell_value) {
  if (is.na(cell_value)) {
    return(cell_spec("NA", "latex", background = "gray"))
  } else {
    cell_spec(cell_value, "latex", background = ifelse(cell_value == "0.000", "#FF9999", "white"))
  }
}

state_sts_table_red[] = as.data.frame(sapply(state_sts_table_round, function(col) {
  sapply(col, format_cell) 
}))

rownames(state_sts_table_red) <- rownames(state_sts_table)

# get latex code for StructTS
kable(state_sts_table_red, caption = "StructTS Results", format="latex", escape = FALSE) %>%
  kable_styling()


# construct MARSS BFGS table with variances instead of SD
state_bfgs_table = t(state_bfgs)
state_bfgs_table = as.data.frame(cbind(state_bfgs_table, samp_sd))
state_bfgs_table = state_bfgs_table %>%
  mutate(var_samp = Jul.89^2) %>%
  mutate(var_shock = var_obs - var_samp)  # get shock variance through additive property of var
state_bfgs_table = select(state_bfgs_table, -Jul.89)  

colnames(state_bfgs_table) = c("drift", "process_var", "observation_var", "sampling_var", "shock_var")
state_bfgs_table$shock_var <- ifelse(state_bfgs_table$shock_var < 0, NA, state_bfgs_table$shock_var)

# round 
library(scales)
state_bfgs_table_round <- data.frame(lapply(state_bfgs_table, function(x) if(is.numeric(x)) number(x, accuracy = 0.001) else x))

# convert cells with 0 to red 
state_bfgs_table_red = state_bfgs_table_round
format_cell = function(cell_value) {
  if (is.na(cell_value)) {
    return(cell_spec("NA", "latex", background = "gray"))
  } else {
    cell_spec(cell_value, "latex", background = ifelse(cell_value <= "0.00010000", "#FF9999", "white")) # changing this code to capture near 0 results in red 
  }
}

state_bfgs_table_red[] = as.data.frame(sapply(state_bfgs_table_round, function(col) {
  sapply(col, format_cell) 
}))

rownames(state_bfgs_table_red) <- rownames(state_bfgs_table)

# get latex code for MARSS BFGS
kable(state_bfgs_table_red, caption = "MARSS BFGS Results", format="latex", escape = FALSE) %>%
  kable_styling()


# construct MARSS kem table with variances
state_kem_table = t(state_kem)
state_kem_table = as.data.frame(cbind(state_kem_table, samp_sd))
state_kem_table = state_kem_table %>%
  mutate(var_samp = Jul.89^2) %>%
  mutate(var_shock = var_obs - var_samp)  # get shock variance through additive property of var
state_kem_table = select(state_kem_table, -Jul.89)  

colnames(state_kem_table) = c("drift", "process_var", "observation_var", "sampling_var", "shock_var")
state_kem_table$shock_var <- ifelse(state_kem_table$shock_var < 0, NA, state_kem_table$shock_var)

# round 
library(scales)
state_kem_table_round <- data.frame(lapply(state_kem_table, function(x) if(is.numeric(x)) number(x, accuracy = 0.001) else x))

# convert cells with 0 to red 
state_kem_table_red = state_kem_table_round
format_cell = function(cell_value) {
  if (is.na(cell_value)) {
    return(cell_spec("NA", "latex", background = "gray"))
  } else {
    cell_spec(cell_value, "latex", background = ifelse(cell_value == "0.00000000", "#FF9999", "white"))
  }
}

state_kem_table_red[] = as.data.frame(sapply(state_kem_table_round, function(col) {
  sapply(col, format_cell) 
}))

rownames(state_kem_table_red) <- rownames(state_kem_table)

# get latex code for MARSS kem
kable(state_kem_table_red, caption = "MARSS KEM Results", format="latex", escape = FALSE) %>%
  kable_styling()






#### TABLES WITH BOOTSTRAP SAMPLING SD AND CONFIDENCE INTERVALS ####
# load in sampling standard deviation estimates (and CIs) created using bootstrap approach
# code for creating of these estimates can be found at in code file: 
# bootstrap_variance_of_variance_of_e0.RMD
boot_results = read.csv(paste(path, "boot_results.csv", sep = "/"), header = TRUE, row.names = 1)

# the e0.var column in boot_results should replace the samp_var column in the prior table
# the sd column in boot_results can be used to calculate confidence intervals for the samp_var column

# to get sd/ci for first three columns of the table we're constructing, we need to get these from the 
# MARSS output 
set.seed(10)
marss_out <- apply(state_e0, 1, rwd_obs_error_bfgs_out)
bfgs_with_cis <- lapply(marss_out, MARSSparamCIs)
obs_var_new <- unlist(lapply(marss_out, function(x) coef(x)$R[1,1]))

# to get the sd/ci for the shock_sd column, we need to subtract cis from obs_var and samp_var (I think
# this will work anyway)
# to produce a conservative interval for shock variance, we need to subtract the lower bound obs var interval 
# from the upper bound samp var interval and vice versa 
# this gives the widest and most conservative interval of the shock variance 
# R param = observation variance 
obs_lb <- unlist(lapply(bfgs_with_cis, function(x) x[[26]][['R']]))
obs_ub <- unlist(lapply(bfgs_with_cis, function(x) x[[25]][['R']]))

samp_lb <- boot_results$e0.var - 2 * boot_results$se_of_var_e0
samp_ub <- boot_results$e0.var + 2 * boot_results$se_of_var_e0

# shock_lb <- obs_lb - samp_ub
# shock_ub <- obs_ub - samp_lb
shock_mean <- obs_var_new - boot_results$e0.var

### trying suggestion from Josh 
# assuming observation variance as a given (conditional on obs variance)
# shock variance should have same se as samping variance 
shock_lb <- shock_mean - 2 * boot_results$se_of_var_e0
shock_ub <- shock_mean + 2 * boot_results$se_of_var_e0

library(ggplot2)

# prepare data
data <- data.frame(
  State = factor(names(shock_lb)),
  Mean = shock_mean,
  Lower = shock_lb,
  Upper = shock_ub
)

data$State <- reorder(data$State, data$Mean)

# Dot chart with confidence intervals
ggplot(data, aes(x = State, y = Mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  coord_flip() +
  xlab("") +
  ylab("Shock Variance") +
  theme_minimal()


# construct MARSS BFGS table with bootstrap variance
state_bfgs_table = t(state_bfgs)
samp_var = boot_results$e0.var
state_bfgs_table = as.data.frame(cbind(state_bfgs_table, samp_var))
state_bfgs_table = state_bfgs_table %>%
  mutate(shock_sd = var_obs - samp_var) %>% # get shock sd through additive property of var
  mutate(var_innov = sqrt(var_innov), var_obs = sqrt(var_obs), samp_var = sqrt(samp_var), shock_sd = sqrt(shock_sd)) # changing all of the variances to SDs

colnames(state_bfgs_table) = c("drift", "process SD", "observation SD", "sampling SD", "shock SD")
# state_bfgs_table$shock_SD <- ifelse(state_bfgs_table$shock_SD < 0, NA, state_bfgs_table$shock_SD)

# round 
library(scales)
library(kableExtra)
state_bfgs_table_round <- data.frame(lapply(state_bfgs_table, function(x) if(is.numeric(x)) number(x, accuracy = 0.001) else x))

# convert cells with 0 to red 
state_bfgs_table_red = state_bfgs_table_round
format_cell = function(cell_value) {
  if (is.na(cell_value)) {
    return(cell_spec("sim 0", "latex"))
  } else {
    cell_spec(cell_value, "latex", background = ifelse(cell_value <= "0.00010000", 
                                                       "#FF9999", "white")) # changing this code to capture near 0 results in red 
  }
}

state_bfgs_table_red[] = as.data.frame(sapply(state_bfgs_table_round, function(col) {
  sapply(col, format_cell) 
}))

rownames(state_bfgs_table_red) <- rownames(state_bfgs_table)

# get latex code for MARSS BFGS
kable(state_bfgs_table_red, 
      caption = "Structural Time Series Model Results for Life Expectancy at Birth", 
      format="latex", escape = FALSE) %>%
  kable_styling()

### FIGURE VERSION OF TABLE ###
pop_89_for_figure = subset(pop_89, !(rownames(pop_89) %in% c("US", "DC")))
pop_89_for_figure = pop_89_for_figure / 1e6 # get population in millions
state_table_for_figure = cbind(state_bfgs_table_round, pop_89_for_figure)
state_table_for_figure = as.data.frame(state_table_for_figure)
state_table_for_figure[] = lapply(state_table_for_figure, function(x) as.numeric(x))

# remove states with 0 observation variance 
state_table_for_figure = subset(state_table_for_figure, state_table_for_figure$observation.SD >= 0.000000001)

# Create the scatterplot for process SD
process <- ggplot(state_table_for_figure, aes(x = Jul.89, y = as.numeric(process.SD))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  theme_minimal() +
  ylim(0, 0.5) +
  labs(x = "Population Size (In Millions)", y = "Standard Deviation in Years", 
       title = "(a) Process Standard Deviation")

# Add the marginal histogram to the y-axis
process <- ggExtra::ggMarginal(process, type = "histogram", margins = "y")

# Display the plot
print(process)


# Create the scatterplot for observation variance
obs <- ggplot(state_table_for_figure, aes(x = Jul.89, y = as.numeric(observation.SD))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  theme_minimal() +
  ylim(0, 0.5) +
  labs(x = "Population Size (In Millions)", y = "Standard Deviation in Years", 
       title = "(b) Obs. Standard Deviation")

# Add the marginal histogram to the y-axis
obs <- ggExtra::ggMarginal(obs, type = "histogram", margins = "y")

# Display the plot
print(obs)

# Create the scatterplot for sampling SD
sampling <- ggplot(state_table_for_figure, aes(x = Jul.89, y = as.numeric(sampling.SD))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  theme_minimal() +
  ylim(0, 0.5) +
  labs(x = "Population Size (In Millions)", y = "Standard Deviation in Years", 
       title = "(c) Sampling Standard Deviation")

# Add the marginal histogram to the y-axis
sampling <- ggExtra::ggMarginal(sampling, type = "histogram", margins = "y")

# Display the plot
print(sampling)

# Create the scatterplot for shock
shock <- ggplot(state_table_for_figure, aes(x = Jul.89, y = as.numeric(shock.SD))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  theme_minimal() +
  ylim(0, 0.5) +
  labs(x = "Population Size (In Millions)", y = "Standard Deviation in Years", 
       title = "(d) Shock Standard Deviation")

# Add the marginal histogram to the y-axis
shock <- ggExtra::ggMarginal(shock, type = "histogram", margins = "y")

# Display the plot
print(shock)

combined_plot <- grid.arrange(process, obs, nrow=1)





### COMPARING MODEL SELECTION CRITERIA BETWEEN RWD AND STS ###
# STS with AIC and BIC 
state_aic_bic_sts = matrix(NA, nrow=5, ncol=50)
state_aic_bic_sts = as.data.frame(apply(state_e0, 1, rwd_with_obs_aic_bic))

# RWD with AIC and BIC
state_aic_bic_rwd = matrix(NA, nrow=4, ncol=50)
state_aic_bic_rwd = as.data.frame(apply(state_e0, 1, rwd_aic_bic))

# list of states with lower AIC for rwd model (is this a similar list to 0 obs var list?)
aic_smaller_indices = which(state_aic_bic_rwd['aic', ] < state_aic_bic_sts['aic', ])
aic_smaller_states = colnames(state_aic_bic_rwd)[aic_smaller_indices]

# list of states with lower BIC for rwd model (similar to 0 obs var list?)
bic_smaller_indices = which(state_aic_bic_rwd['bic', ] < state_aic_bic_sts['bic', ])
bic_smaller_states = colnames(state_aic_bic_rwd)[bic_smaller_indices]

# Print the lists
aic_smaller_states
bic_smaller_states

# we can see that more than half of the states favor the simpler random walk with drift models 
# over the STS trend models considering AIC and BIC 
# does this suggest something disappointing about the goodness-of-fit of the STS models? 
# or do we think AIC and BIC aren't the best criterion to use here? 

### COMPARING MODEL SELECTION CRITERIA USING MARSS BFGS ###
aic_aicc = as.data.frame(t((apply(state_e0, 1, get_aic_aicc))))
colnames(aic_aicc) = c("rwd_aic", "sts_aic", "rwd_aicc", "sts_aicc")
aic_aicc$rwd_aic <- round(aic_aicc$rwd_aic, 1)
aic_aicc$sts_aic <- round(aic_aicc$sts_aic, 1)
aic_aicc$rwd_aicc <- round(aic_aicc$rwd_aicc, 1)
aic_aicc$sts_aicc <- round(aic_aicc$sts_aicc, 1)

# prep barchart
aic = aic_aicc[,1:2]
aic = rownames_to_column(aic, "state")
aic$diff = aic$sts_aic - aic$rwd_aic

# pure AIC value plot 
aic_long = pivot_longer(aic, cols=c("rwd_aic", "sts_aic"), names_to="model", values_to="AIC")

ggplot(aic_long, aes(x = state, y = AIC, fill = model)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Adjust text angle and justification for readability
  labs(title = "Comparison of AIC Values by State",
       x = "State",
       y = "AIC Value",
       fill = "Model")

# difference plot 
ggplot(aic, aes(x = state, y = diff, fill = diff > 0)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = c("TRUE" = "lightblue", "FALSE" = "darkblue"), 
                    name = "Model Better", 
                    labels = c("TRUE" = "RWD Better", "FALSE" = "STS Better")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "State",
       y = "AIC Difference")

# prep table
aic_smaller_for_rwd = which(aic_aicc[ ,'rwd_aic'] < aic_aicc[ ,'sts_aic'])
aic_smaller_for_rwd = rownames(aic_aicc)[aic_smaller_for_rwd]

aicc_smaller_for_rwd = which(aic_aicc[ ,'rwd_aicc'] < aic_aicc[ ,'sts_aicc'])
aicc_smaller_for_rwd = rownames(aic_aicc)[aicc_smaller_for_rwd]

# Print the lists
aic_smaller_for_rwd
aicc_smaller_for_rwd

### CREATING TABLE 4: TABLES OF AIC AND AICc VALUES FOR EACH STATE AND EACH MODEL ### 

aic_aicc$rwd_aic = cell_spec(aic_aicc$rwd_aic, "latex", 
                             background = ifelse(rownames(aic_aicc) %in% aic_smaller_for_rwd, "#FF9999", "white"))

aic_aicc$rwd_aicc <- cell_spec(aic_aicc$rwd_aicc, "latex", 
                               background = ifelse(rownames(aic_aicc) %in% aicc_smaller_for_rwd, "#FF9999", "white"))

kable(aic_aicc, "latex", escape = FALSE) %>%
  kable_styling()




### PLOTTING DATA AND PREDICTIONS FROM RWD AND STS ###

# test with whole US data
rwd_us = rwd_bfgs(us_e0)
sts_us = rwd_obs_error_bfgs_out(us_e0)

fit_rwd = fitted(rwd_us)
fit_sts = fitted(sts_us)

plot(us_e0)
lines(fit_rwd[,4], col = "red")
lines(fit_sts[,4], col = "blue")


# The rwd and sts predictions look identical here, which makes sense given that 
# the sts model estimates an obs variance of 0
rwd_ca = rwd_bfgs(state_e0["california",])
sts_ca = rwd_obs_error_bfgs_out(state_e0["california",])

fit_rwd = fitted(rwd_ca)
fit_sts = fitted(sts_ca)

plot(state_e0["california",])
lines(fit_rwd[,4], col = "red")
lines(fit_sts[,4], col = "blue")

# testing on a state with a smaller population but still 0 obs variance 
rwd_ct = rwd_bfgs(state_e0["connecticut",])
sts_ct = rwd_obs_error_bfgs_out(state_e0["connecticut",])

fit_rwd = fitted(rwd_ct)
fit_sts = fitted(sts_ct)

plot(state_e0["connecticut",])
lines(fit_rwd[,4], col = "red")
lines(fit_sts[,4], col = "blue")

rwd_nc = rwd_bfgs(state_e0["north carolina",])
sts_nc = rwd_obs_error_bfgs_out(state_e0["north carolina",])

fit_rwd = fitted(rwd_nc)
fit_sts = fitted(sts_nc)

plot(state_e0["north carolina",])
lines(fit_rwd[,4], col = "red")
lines(fit_sts[,4], col = "blue")


# an example of states with a high observation variance value and a stronger STS 
# indicated by AIC values 
rwd_ak = rwd_bfgs(state_e0["alaska",])
sts_ak = rwd_obs_error_bfgs_out(state_e0["alaska",])

fit_rwd = fitted(rwd_ak)
fit_sts = fitted(sts_ak)

plot(state_e0["alaska",])
lines(fit_rwd[,4], col = "red")
lines(fit_sts[,4], col = "blue")

rwd_ak = rwd_bfgs(state_e0["massachusetts",])
sts_ak = rwd_obs_error_bfgs_out(state_e0["massachusetts",])

fit_rwd = fitted(rwd_ak)
fit_sts = fitted(sts_ak)

plot(state_e0["massachusetts",])
lines(fit_rwd[,4], col = "red")
lines(fit_sts[,4], col = "blue")







### ADDITIONAL ANALYSES ###


### RUNNING MODELS WITH BAYESIAN DATA ###
bayes_data = readRDS(paste(path, "USA_b_county_lt.rds", sep="/"))
e0_bayes = bayes_data %>% 
  filter(age == 0) %>%
  select(fips, year, ex)
# reshape data with "spread"
e0_bayes = e0_bayes %>%
  spread(key = year, value = ex)
e0_bayes = e0_bayes[,-1]

states_dc = c("AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DE", "DC", "FL", "GA", "HI", "ID", "IL", "IN", "IA", "KS", 
              "KY", "LA", "ME", "MD", "MA", "MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY", 
              "NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VT", "VA", "WA", 
              "WV", "WI", "WY")
rownames(e0_bayes) = states_dc

# MARSS kem
bayes_kem = matrix(NA, nrow=51, ncol=3)
bayes_kem = as.data.frame(apply(e0_bayes, 1, rwd_with_obs_error_kem))
# results in 8 states that don't converge when used with "kem" method

# MARSS BFGS
bayes_bfgs = matrix(NA, nrow=51, ncol=3)
bayes_bfgs = as.data.frame(apply(e0_bayes, 1, rwd_with_obs_error_bfgs))
# results in 10 states with 0 observation variance 

# StructTS
bayes_sts = matrix(NA, nrow=51, ncol=3)
bayes_sts = as.data.frame(apply(e0_bayes, 1, rwd_with_obs_error))
# results in 13 states with 0 observation variance

# get names of states with var_obs = 0 
bayes_no_obs_kem = which(is.na(state_kem[3,]))
bayes_no_obs_kem = colnames(state_kem)[state_no_obs_kem]
bayes_no_obs_kem

bayes_no_obs_bfgs = which(bayes_bfgs[3,]<=0.0000000001)
bayes_no_obs_bfgs = colnames(bayes_bfgs)[bayes_no_obs_bfgs]
bayes_no_obs_bfgs

bayes_no_obs_sts = which(bayes_sts[3,]==0)
bayes_no_obs_sts = colnames(bayes_sts)[bayes_no_obs_sts]
bayes_no_obs_sts


## EXAMINING DISTURBANCE DISTRIBUTIONS ##
# create resid function
marss_kem_resid <- function(y, ...)
{
  library(MARSS)
  mod.list.2 <- list(B = matrix(1),
                     U = matrix("d"),
                     Q = matrix("q"),
                     Z = matrix(1),
                     A = matrix(0),
                     R = matrix("r"),
                     x0 = matrix("mu"),
                     tinitx = 0)
  out2.marss = MARSS(y, model = mod.list.2, control=list(maxit=1000), method= "kem", ...)
  
  resid = residuals(out2.marss, type = "tT")
  model <- subset(resid, name=="model")
  state <- subset(resid, name=="state")
  return(list(model, state))
}

# examining distribution of disturbances with no2020 state data 
kem_resid = marss_kem_resid(no2020)

# model residuals
model_resid = kem_resid[[1]]$.resids
# histogram
hist(model_resid, breaks = 50)
# qqplots
qqnorm(model_resid)
qqline(model_resid)
# Shapiro-Wilk normality test
shapiro.test(model_resid)

# state residuals
state_resid = kem_resid[[2]]$.resids
# histogram
hist(state_resid, breaks = 50)
# qqplots
qqnorm(state_resid)
qqline(state_resid)
# Shapiro-Wilk normality test
shapiro.test(state_resid)


### EXAMINING SPECIFIC STATE DISTURBANCE DISTRIBUTIONS ###
# make residual plots using kem method
s_resid = as.data.frame(apply(state_e0, 1, marss_kem_resid))

# looking at state residuals for states with and without an estimated 0 
# observation variance 

# states with 0 obs variance examined below 
  # Florida 
  # New York

# states with non-0 obs variance examined below
  # Utah
  # Georgia 

# Florida model residuals 
# histogram
hist(s_resid$florida..resids, breaks = 50)
# qqplots
qqnorm(s_resid$florida..resids)
qqline(s_resid$florida..resids)
# Shapiro-Wilk normality test
shapiro.test(s_resid$florida..resids)

# Florida state residuals 
# histogram
hist(s_resid$florida..resids.1, breaks = 50)
# qqplots
qqnorm(s_resid$florida..resids.1)
qqline(s_resid$florida..resids.1)
# Shapiro-Wilk normality test
shapiro.test(s_resid$florida..resids.1)

# NY model residuals 
# histogram
hist(s_resid$new.york..resids, breaks = 50)
# qqplots
qqnorm(s_resid$new.york..resids)
qqline(s_resid$new.york..resids)
# Shapiro-Wilk normality test
shapiro.test(s_resid$new.york..resids)

# NY state residuals 
# histogram
hist(s_resid$new.york..resids.1, breaks = 50)
# qqplots
qqnorm(s_resid$new.york..resids.1)
qqline(s_resid$new.york..resids.1)
# Shapiro-Wilk normality test
shapiro.test(s_resid$new.york..resids.1)

# Utah model residuals 
# histogram
hist(s_resid$utah..resids, breaks = 50)
# qqplots
qqnorm(s_resid$utah..resids)
qqline(s_resid$utah..resids)
# Shapiro-Wilk normality test
shapiro.test(s_resid$utah..resids)

# Utah state residuals 
# histogram
hist(s_resid$utah..resids.1, breaks = 50)
# qqplots
qqnorm(s_resid$utah..resids.1)
qqline(s_resid$utah..resids.1)
# Shapiro-Wilk normality test
shapiro.test(s_resid$utah..resids.1)

# Georgia model residuals 
# histogram
hist(s_resid$georgia..resids, breaks = 50)
# qqplots
qqnorm(s_resid$georgia..resids)
qqline(s_resid$georgia..resids)
# Shapiro-Wilk normality test
shapiro.test(s_resid$georgia..resids)

# Georgia state residuals 
# histogram
hist(s_resid$georgia..resids.1, breaks = 50)
# qqplots
qqnorm(s_resid$georgia..resids.1)
qqline(s_resid$georgia..resids.1)
# Shapiro-Wilk normality test
shapiro.test(s_resid$georgia..resids.1)


