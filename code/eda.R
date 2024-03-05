library(dtwclust)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(maps)
library(sf)
library(ggthemes)
library(maptools)
library(mapproj)
library(GGally) # for correlation plot
library(zoo) # for moving average
library(corrr)
library(ggcorrplot)
library(FactoMineR)
library(PerformanceAnalytics)
library(ggrepel)
library(ggplot2)
library(ggfortify)
library(geosphere) # for calculating distances between states 
library(reshape2) # for melting data frames into long format 
library(gganimate)

### STATE DATA ###
# read in life tables for all 50 states 
path = "/Users/paigepark/repos/us-mortality/data"
states = c("AK", "AL", "AZ", "AR", "CA", "CO", "CT", "DE", "FL", "GA", "HI", "ID", "IL", "IN", "IA", "KS", 
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
e0_df <- as.data.frame(e0_df)

# create long version of dataframe
e0_df_long <- e0_df %>%
  gather(key = "Year", value = "Life_Expectancy") 

e0_df_long <- cbind(states, e0_df_long)
e0_df_long$Year <- as.numeric(e0_df_long$Year)

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
div_e0 = as.data.frame(div_e0)

# create long version of dataframe
e0_df_long_div <- div_e0 %>%
  gather(key = "Year", value = "Life_Expectancy") 

e0_df_long_div <- cbind(full_div_names, e0_df_long_div)
e0_df_long_div$Year <- as.numeric(e0_df_long_div$Year)

### PLOTTING E0 TO SEE OVERALL PATTERNS AND REVERSALS ###
### STATES - APPENDIX? ###
# estimate of sd based on mid point pop size 
pop_89 = read.csv(paste(path, "89_pop.csv", sep = "/"), header = TRUE, row.names = 1)
pop_89 = pop_89["Jul.89"]

# calculating and doubling the standard deviation to allow for a form of significance testing 
# we want to be able to identify the reversals of magnitude greater than 2*standard deviation
# to insure the reversals aren't simply due to random error
samp_sd = 150/sqrt(pop_89)
samp_sd = subset(samp_sd, !(rownames(samp_sd) %in% c("US", "DC")))
double_samp_sd = 2*samp_sd
double_samp_sd = double_samp_sd %>%
  rownames_to_column(var = "states")

# merge sd into long e0 dataframe for states
e0_df_long = e0_df_long %>%
  left_join(double_samp_sd, by = "states")

# set up dataframe for reversal figure
e0_df_long = e0_df_long %>%
  arrange(states, Year) %>%
  group_by(states) %>%
  mutate(
    
    # designates a 1 year decline in life expectancy (1 is decline, 0 is no decline)
    decrease1 = ifelse(lag(Life_Expectancy) > Life_Expectancy, 1, 0),
    
    # designates a 2 year decline in life expectancy 
    decrease2 = ifelse(lag(Life_Expectancy,2) > lag(Life_Expectancy) & lag(Life_Expectancy) > Life_Expectancy, 1, 0),
    
    #designates a 3 year decline in life expectancy
    decrease3 = ifelse(lag(Life_Expectancy,3) > lag(Life_Expectancy,2) & lag(Life_Expectancy,2) > lag(Life_Expectancy) & lag(Life_Expectancy) > Life_Expectancy, 1, 0),
    
    # calculate the size of the decline
    decline_size = ifelse(lag(Life_Expectancy) > Life_Expectancy, lag(Life_Expectancy) - Life_Expectancy, NA),
    
    # designate whether or not the decline was significant (whether decline size is > 2*sd)
    sig = ifelse(decline_size > Jul.89, 1, 0)
  )
 
# replica of Josh and Ron's figure 
# black and grey = 1 year reversals
# blue and light blue = 2 year reversals
# red and light red = 3 year reversals
# darker are statistically significant declines (larger than 2 standard deviations based on pop size)

ggplot(e0_df_long, aes(x = Year, y = Life_Expectancy)) +
    geom_line(aes(group = states)) +
    geom_vline(data = subset(e0_df_long, decrease1 == 1 & sig == 1), aes(xintercept = Year), color = "black") +
    geom_vline(data = subset(e0_df_long, decrease2 == 1 & sig == 1), aes(xintercept = Year), color = "blue") +
    geom_vline(data = subset(e0_df_long, decrease3 == 1 & sig == 1), aes(xintercept = Year), color = "red") +
    facet_wrap(~ states, ncol = 10) +
    theme_minimal() +
    labs(x = "Year", y = "Life Expectancy at Birth") +
    scale_x_continuous(breaks = c(min(e0_df_long$Year), max(e0_df_long$Year)))

### CENSUS DIVSION REVERSALS FIGURE ###
# Creating a lookup dataframe with State and their corresponding Census Division
lookup_table <- data.frame(
  state = c("CT", "ME", "MA", "NH", "RI", "VT", 
            "NJ", "NY", "PA", "IL", "IN", "MI", "OH", "WI", 
            "IA", "KS", "MN", "MO", "NE", "ND", "SD", 
            "DE", "FL", "GA", "MD", "NC", "SC", "VA", "WV", 
            "AL", "KY", "MS", "TN", "AR", "LA", "OK", "TX", 
            "AZ", "CO", "ID", "MT", "NV", "NM", "UT", "WY", 
            "AK", "CA", "HI", "OR", "WA"),
  full_div_names = c(rep("new england", 6), rep("middle atlantic", 3), 
               rep("east north central", 5), rep("west north central", 7), 
               rep("south atlantic", 8), rep("east south central", 4), 
               rep("west south central", 4), rep("mountain", 8), 
               rep("pacific", 5))
)

pop_89_w_states = rownames_to_column(pop_89, var = "state")
pop_89_w_div = merge(lookup_table, pop_89_w_states, by = "state")

pop_89_div = pop_89_w_div %>%
  group_by(full_div_names) %>% 
  summarize(div_pop = sum(Jul.89, na.rm = TRUE))

pop_89_w_sd = pop_89_div %>% 
  mutate(sd = 150/sqrt(div_pop))
pop_89_w_2sd = pop_89_w_sd %>% 
  mutate(double_sd = sd*2)
#double_samp_sd_div = double_samp_sd_div %>%
 # rownames_to_column(var = "full_div_names")


# merge sd into other dataframe
e0_df_long_div_sd = e0_df_long_div %>%
  left_join(pop_89_w_2sd, by = "full_div_names")

# set up dataframe for census division reversal figure
e0_reversals_div = e0_df_long_div_sd %>%
  arrange(full_div_names, Year) %>%
  group_by(full_div_names) %>%
  mutate(
    decrease1 = ifelse(lag(Life_Expectancy) > Life_Expectancy, 1, 0),
    decrease2 = ifelse(lag(Life_Expectancy,2) > lag(Life_Expectancy) & lag(Life_Expectancy) > Life_Expectancy, 1, 0),
    decrease3 = ifelse(lag(Life_Expectancy,3) > lag(Life_Expectancy,2) & lag(Life_Expectancy,2) > lag(Life_Expectancy) & lag(Life_Expectancy) > Life_Expectancy, 1, 0),
    decline_size = ifelse(lag(Life_Expectancy) > Life_Expectancy, lag(Life_Expectancy) - Life_Expectancy, NA),
    sig = ifelse(decline_size > double_sd, 1, 0)
  )

# plot
 
ggplot(e0_reversals_div, aes(x = Year, y = Life_Expectancy)) +
  geom_line(aes(group = full_div_names)) +
  geom_vline(data = subset(e0_reversals_div, decrease1 == 1), aes(xintercept = Year), color = "black") +
  geom_vline(data = subset(e0_reversals_div, decrease2 == 1), aes(xintercept = Year), color = "blue") +
  geom_vline(data = subset(e0_reversals_div, decrease3 == 1), aes(xintercept = Year), color = "red") +
  facet_wrap(~ full_div_names, ncol = 3) +
  theme_minimal() +
  labs(x = "Year", y = "Life Expectancy at Birth") +
  scale_x_continuous(breaks = c(min(e0_reversals_div$Year), 1970, 1980, 1990, 2000, 2010, max(e0_reversals_div$Year)))






### CLUSTER TIME SERIES PLOTS - CLUSTERING OF TIME SERIES, NOT CLUSTERING OF SHOCKS ###

# clustering and validating 
e0_clust <- tsclust(e0_df, type="partitional", k=2L:17L, distance="dtw", centroid="pam")
cvi_results <- lapply(e0_clust, cvi)
names(cvi_results) = paste("k", 2:17, sep="_")
cvi_results <- do.call(cbind, cvi_results)
cvi_results
# could perform more robust checks to determine ideal K but looking at rough heuristics:
# K=2 or k=4 seem like good choices, as it has relatively high Sil, CH, and D values and low DB, DBstar, and COP values.

# final clusters 
set.seed(123)
e0_clust <- tsclust(e0_df, type="partitional", k=2L, distance="dtw", centroid="pam", seed=123)
plot(e0_clust, type = "sc", xaxt = "n")
years <- 1959:2020
# axis(side=1, at=x, labels=years)

clusters = as.data.frame(cbind(region = full_state_names, cluster = e0_clust@cluster))

### CLUSTER MAP - CLUSTERING OF TS NOT SHOCKS ### 
us <- map_data("state")

merged_data <- merge(us, clusters, by = "region", all = F)


# Gen map
gg <- ggplot()
gg <- gg + geom_map(data=us, map=us,
                    aes(x=long, y=lat, map_id=region),
                    fill="#ffffff", color="#ffffff", size=0.15)
gg <- gg + geom_map(data=clusters, map=us,
                    aes(fill=as.factor(cluster), map_id=region),
                    color="#ffffff", size=0.15)
gg <- gg + scale_fill_discrete(guide='legend') # Change this line
gg <- gg + labs(x=NULL, y=NULL)
gg <- gg + coord_map("albers", lat0 = 39, lat1 = 45) 
gg <- gg + theme(panel.border = element_blank())
gg <- gg + theme(panel.background = element_blank())
gg <- gg + theme(axis.ticks = element_blank())
gg <- gg + theme(axis.text = element_blank())
gg

### CORRELATION PLOTS OF SHOCKS ###
### define shocks ### 
e0_df = as.matrix(e0_df)
moving_average = matrix(NA, nrow = 50, ncol = 62)
moving_average_left = matrix(NA, nrow = 50, ncol = 62)
moving_average_right = matrix(NA, nrow = 50, ncol = 62)
for (i in 1:nrow(e0_df)) {
  zoo_object <- zoo(e0_df[i,])

  moving_average[i,] <- rollmean(zoo_object, k = 7, fill = NA, align = "center")
  moving_average_left[i,] <- rollmean(zoo_object, k = 7, fill = NA, align = "left")
  moving_average_right[i,] <- rollmean(zoo_object, k = 7, fill = NA, align = "right")
}


# adding in moving average estimates for first and last 3 years, 
# which cannot be calculated with the "center" alignment
moving_average[,1:3] <- moving_average_left[,1:3]
moving_average[,60:62] <- moving_average_right[,60:62]

shock_ts <- moving_average - e0_df

### MAP OF SHOCK GEOGRAPHY IN BAD FLU YEAR (2009 SWINE FLU) ###
shocks_2009 <- as.data.frame(shock_ts[,"2009"])
shocks_2009 <- rownames_to_column(shocks_2009, var = "region")
colnames(shocks_2009) <- c("region", "shocks")
us <- map_data("state")
#merged_2009 <- merge(us, shocks_2009, by = "region", all = F)

ggplot() + 
  geom_map(data=us, map=us,
          aes(x=long, y=lat, map_id=region)) +
  geom_map(data=shocks_2009, map=us,
            aes(fill=shocks, map_id=region)) +
  scale_fill_gradient(low = "darkred", high = "white") +
  theme_minimal() +
  labs(fill = "Shock Magnitude")

### MAP OF SHOCK GEOGRAPHY IN BAD HEAT WAVE YEAR (1988 HEAT WAVE) ###
shocks_1988 <- as.data.frame(shock_ts[,"1988"])
shocks_1988 <- rownames_to_column(shocks_1988, var = "region")
colnames(shocks_1988) <- c("region", "shocks")
us <- map_data("state")
#merged_2009 <- merge(us, shocks_2009, by = "region", all = F)

ggplot() + 
  geom_map(data=us, map=us,
           aes(x=long, y=lat, map_id=region)) +
  geom_map(data=shocks_1988, map=us,
           aes(fill=shocks, map_id=region)) +
  scale_fill_gradient(low = "darkred", high = "white") +
  theme_minimal() +
  labs(fill = "Shock Magnitude")

### MAP ANIMATION ###
# adjust data first 
shock_ts_long <- rownames_to_column(as.data.frame(t(shock_ts)), var = "year")
shock_ts_long <- pivot_longer(shock_ts_long,
                              cols = -year,
                              names_to = "region",
                              values_to = "shocks")
shock_ts_long$year <- as.integer(shock_ts_long$year)
shock_ts_no2020 <- shock_ts_long[shock_ts_long$year != 2020, ]

base_map <- ggplot() + 
  geom_map(data=us, map=us,
           aes(x=long, y=lat, map_id=region)) +
  geom_map(data=shock_ts_no2020, map=us,
           aes(fill=shocks, map_id=region, group=year)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(fill = "Shock Magnitude")

map_with_animation <- base_map +
  transition_time(year) +
  ggtitle('Year: {frame_time}')

## save as gif
animation = animate(map_with_animation,
                    end_pause = 50,
                    renderer = gifski_renderer(),
                    width = 800,
                    height = 450,
                    fps = 1
                    )
animation

anim_save(
  filename = paste(path, "animation.gif", sep="/")
)
  
### PCA PLOT - NOT THAT USEFUL ###
# changing labels on shock_ts for better looking plot : not good code, fix later 
state_abb <- c(
  "alabama" = "AL",
  "alaska" = "AK",
  "arizona" = "AZ",
  "arkansas" = "AR",
  "california" = "CA",
  "colorado" = "CO",
  "connecticut" = "CT",
  "delaware" = "DE",
  "florida" = "FL",
  "georgia" = "GA",
  "hawaii" = "HI",
  "idaho" = "ID",
  "illinois" = "IL",
  "indiana" = "IN",
  "iowa" = "IA",
  "kansas" = "KS",
  "kentucky" = "KY",
  "louisiana" = "LA",
  "maine" = "ME",
  "maryland" = "MD",
  "massachusetts" = "MA",
  "michigan" = "MI",
  "minnesota" = "MN",
  "mississippi" = "MS",
  "missouri" = "MO",
  "montana" = "MT",
  "nebraska" = "NE",
  "nevada" = "NV",
  "new hampshire" = "NH",
  "new jersey" = "NJ",
  "new mexico" = "NM",
  "new york" = "NY",
  "north carolina" = "NC",
  "north dakota" = "ND",
  "ohio" = "OH",
  "oklahoma" = "OK",
  "oregon" = "OR",
  "pennsylvania" = "PA",
  "rhode island" = "RI",
  "south carolina" = "SC",
  "south dakota" = "SD",
  "tennessee" = "TN",
  "texas" = "TX",
  "utah" = "UT",
  "vermont" = "VT",
  "virginia" = "VA",
  "washington" = "WA",
  "west virginia" = "WV",
  "wisconsin" = "WI",
  "wyoming" = "WY"
)
names(shock_ts) <- state_abb[names(shock_ts)]

shock_ts <- t(shock_ts)
corr_matrix <- cor(shock_ts)
data.pca <- princomp(corr_matrix)
summary(data.pca)
autoplot(data.pca, label = TRUE, shape = FALSE)
#biplot(data.pca)

### LOADINGS PLOT - NOT THAT USEFUL ###
loadings <- sort(data.pca$loadings[, 1], decreasing = TRUE)
barplot(loadings, main="Loadings of the first principal component", xlab="Variables", ylab="Loadings")
l_plot = autoplot(data.pca, loadings = TRUE)
loadings_df <- as.data.frame(data.pca$loadings[, 1:2])
l_plot + 
  geom_text_repel(data = loadings_df, max.overlaps = 20,
                  aes(x = Comp.1, y = Comp.2, label = rownames(loadings_df)))

### CORRELATION PLOT WITH ALL 50 STATES - TOO CROWDED ###
# reorder data according to factor loadings on PC1
# The order you want
new_order <- names(loadings)

# Rearrange the dataframe
shock_ts <- shock_ts[, new_order]

# rename columns for easier reading in plot 
#colnames(shock_ts) <- c("IN", "MS", "IL", "MO", "LA", "MI", "KS", 
                     #   "TN", "AL", "NJ", "TX", "NY", "CT", "OH", 
                      #  "KY", "WV", "NC", "IA", "AZ", "CA", "AR", 
                       # "NM", "CO", "GA", "WI", "MN", "SC", "VA", 
                        #"OK", "MD", "NE", "ND", "PA", "SD", "FL", 
                        #"DE", "NV", "RI", "WY", "UT", "ID", "MT", 
                        #"WA", "MA", "OR", "ME", "AK", "NH", "VT", "HI")
  # change ^^ to not be hard coded

#corr_matrix <- cor(shock_ts)
ggcorrplot(corr_matrix, insig = "blank") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

### CORRELATION PLOT WITH SRS OF STATES - DOESN'T ACTUALLY SAY MUCH ###
plots = list()
for (i in 1:4) {
  sampled_indices <- sort(sample(ncol(shock_ts), 15))
  sampled_columns <- shock_ts[, sampled_indices]
  sub_corr = cor(sampled_columns)
  plots[[i]] = ggcorrplot(sub_corr, insig = "blank") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}
grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], nrow=2, ncol=2)

### CORRELATION PLOT WITH STRATIFIED SAMPLE ###
# correlation plots with random sample stratified by census division (2 per division)
# list of divisions given here: https://www2.census.gov/geo/pdfs/maps-data/maps/reference/us_regdiv.pdf
# used R sample command to select 2 states from each division 
# full list of states 
  # northeast
    # div 1: MA & VT
    # div 2: PA & NJ
  # midwest
    # div 3: WI & MI
    # div 4: MN & IA
  # south
    # div 5: VA & GA
    # div 6: MS & TN
    # div 7: OK & AR
  # west
    # div 8: NM & AZ
    # div 9: OR & WA
states_to_include <- c("MA", "VT", "PA", "NJ", "WI", "MI", "MN", "IA", 
                       "VA", "GA", "MS", "TN", "OK", "AR", "NM", "AZ", "OR", "WA")
subset_shock <- shock_ts[, names(shock_ts) %in% states_to_include]
#subset_shock <- shock_ts[, states_to_include]
sub_corr = cor(subset_shock)
ggcorrplot(sub_corr, insig = "blank") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# sampling code   
# > sample(1:6, 2, replace=FALSE)
# [1] 3 6
# > sample(1:3, 2, replace=FALSE)
# [1] 3 1
# > sample(1:5, 2, replace=FALSE)
# [1] 5 3
# > sample(1:7, 2, replace=FALSE)
# [1] 3 1
# > sample(1:8, 2, replace=FALSE)
# [1] 7 3
# > sample(1:4, 2, replace=FALSE)
# [1] 4 3
# > sample(1:4, 2, replace=FALSE)
# [1] 3 1
# > sample(1:8, 2, replace=FALSE)
# [1] 4 1
# > sample(1:3, 2, replace=FALSE)
# [1] 3 2



### CORRELATION BY DISTANCE SCATTER PLOT ###
# steps
# 1. calculate centroids of each state
us <- map_data("state")
centroids <- us %>% 
  group_by(region) %>%
  group_modify(~ data.frame(centroid(cbind(.x$long, .x$lat))))

# 2. calculate distance between each state - have a dataframe representing distances between all pairwise combos of states 
# Calculate pairwise distances
distance_matrix <- as.data.frame(distm(centroids[, c("lon", "lat")], fun = distVincentySphere)) / 1000  # in kilometers

# Convert to miles
distance_matrix <- distance_matrix * 0.621371

# Add state names as row and column names
rownames(distance_matrix) <- centroids$region
colnames(distance_matrix) <- centroids$region

# 3. calculate correlation coefficients for all pairwise combinations of states 
  # this has already been calculated above (corr_matrix)

# 4. create a plot with distance on x and corr coef on y (each point is a pairwise combo of states)

# melt the two matrices into long format 
distance_matrix <- as.matrix(distance_matrix)
distance_long <- melt(replace(distance_matrix, lower.tri(distance_matrix, TRUE), NA), na.rm = TRUE, 
                      varnames = c("state1", "state2"), value.name = "distance")

corr_long <- melt(replace(corr_matrix, lower.tri(corr_matrix, TRUE), NA), na.rm = TRUE, 
                  varnames = c("state1", "state2"), value.name = "shock_correlation")

combined_long <- distance_long %>%
  inner_join(corr_long, by = c("state1", "state2"))


ggplot(data = combined_long, x = distance, y = shock_correlation)


reg = lm(shock_correlation ~ distance, data = combined_long)
summary(reg)

ggplot(combined_long, aes(x = distance, y = shock_correlation)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal(base_size = 15) + # white background
  labs(
    x = "Distance between state centers (in miles)",
    y = "Correlation in mortality shocks"
  )

### CORRELATION PLOT WITH CENSUS DIVISIONS ###
### DEFINE SHOCKS ### 
div_e0 = as.matrix(div_e0)
moving_average = matrix(NA, nrow = nrow(div_e0), ncol = ncol(div_e0))
moving_average_left = matrix(NA, nrow = nrow(div_e0), ncol = ncol(div_e0))
moving_average_right = matrix(NA, nrow = nrow(div_e0), ncol = ncol(div_e0))
for (i in 1:nrow(div_e0)) {
  zoo_object <- zoo(div_e0[i,])
  
  moving_average[i,] <- rollmean(zoo_object, k = 7, fill = NA, align = "center")
  moving_average_left[i,] <- rollmean(zoo_object, k = 7, fill = NA, align = "left")
  moving_average_right[i,] <- rollmean(zoo_object, k = 7, fill = NA, align = "right")
}


# adding in moving average estimates for first and last 3 years, 
# which cannot be calculated with the "center" alignment
moving_average[,1:3] <- moving_average_left[,1:3]
moving_average[,60:62] <- moving_average_right[,60:62]

shock_div <- moving_average - div_e0
shock_div <- as.data.frame(t(shock_div))

corr_matrix <- cor(shock_div)
ggcorrplot(corr_matrix, insig = "blank", hc.order = TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


#### EXTRA STUFF THAT PROBABLY WON'T GO IN APPENDIX OR PAPER####
#### AGE PATTERNS ####
# look at the age pattern of shocks and distinguish short and longer-lasting effects by age pattern
# start by creating a similar figure to Ron and Josh's paper (I'm not sure exactly how they created 
# theirs so I'm doing something a little different that will hopefully give a similar type of info)

  # 1. Separate years of mortality improvement and decline (by state)
    # add column to long data indicating whether life expectancy improved or not 
    df_improved <- e0_df_long %>%
      group_by(states) %>%
      arrange(Year) %>%
      mutate(improved = Life_Expectancy > lag(Life_Expectancy))

    # create list of vectors containing the years of improvement for each state
    improvement_years <- df_improved %>%
      filter(improved == TRUE) %>%
      group_by(states) %>%
      summarize(years = list(Year))
    
    # create list of vectors containing the years of decline for each state
    decline_years <- df_improved %>%
      filter(improved == FALSE) %>%
      group_by(states) %>%
      summarize(years = list(Year))
    
    # Convert the lists to named lists for easier access
    improvement_years <- setNames(improvement_years$years, improvement_years$states)
    decline_years <- setNames(decline_years$years, decline_years$states)
    
  # 2. Get average age-specific mortality rate for good and bad years 
    improv_mxs = matrix(NA, nrow=110, ncol=50)
    for (i in 1:length(data)) {
      for (j in 1:110) { # add 110+ later
        improv_mxs[j,i] = mean(subset(data[[i]], Age == j-1 & Year %in% improvement_years[[i]])$mx)
      }
    }
    
    improv_mxs = cbind(seq(0, 109, by = 1), improv_mxs)
    colnames(improv_mxs) = c("Age", names(improvement_years))
    #rownames(improv_mxs) = seq(0, 109, by = 1)
    
    improv_mxs = as.data.frame(improv_mxs)
    
    
    decline_mxs = matrix(NA, nrow=110, ncol=50)
    for (i in 1:length(data)) {
      for (j in 1:110) { # add 110+ later
        decline_mxs[j,i] = mean(subset(data[[i]], Age == j-1 & Year %in% decline_years[[i]])$mx)
      }
    }
    
    decline_mxs = cbind(seq(0, 109, by = 1), decline_mxs)
    colnames(decline_mxs) = c("Age", names(decline_years))
    #rownames(decline_mxs) = seq(0, 109, by = 1)
    decline_mxs = as.data.frame(decline_mxs)
    
  # 3. Plot it (can change out the state names to see different states)
    # this wasn't that interesting
    # can maybe come back to this later, but I'm not putting more time into it now
    ggplot() +
      geom_line(data=improv_mxs, aes(x=Age, y=UT), color = "blue") +
      geom_line(data=decline_mxs, aes(x=Age, y=UT), color = "red")


#### FIGURE ?: COVID E0 DROPS BY STATE ####

# Use the e0_df_long dataset 
# Filter for years 2019 and 2020, then calculate the differences in life expectancy
df_diff <- e0_df_long %>%
  filter(Year %in% c(2019, 2020)) %>%
  group_by(states) %>%
  summarize(difference = -diff(Life_Expectancy)) %>%
  arrange(difference)

# Then, create a barplot of the differences
    ggplot(df_diff, aes(x = reorder(states, difference), y = difference)) +
      geom_bar(stat = "identity", fill = "black") +
      coord_flip() +
      theme_minimal() +
      labs(x = "State",
           y = "Difference in Life Expectancy in Years (2020 - 2019)") +
      theme(legend.position = "none")
    
    

#### FIGURE/TABLE ?: HIGHEST AND LOWEST PERFORMERS 1959 TO 2019 ####
# could be cool to do a animation here of the US map with top performers in green 
# and bottom performers in red and show how this has changed throughout the years 

# Filter for the years 1959 and 2019
df_filtered <- e0_df_long %>%
  filter(Year %in% c(1959, 2019))

    # uncomment if running the map code
#state_name_mapping <- setNames(tolower(state.name), state.abb)
#df_filtered$states <- state_name_mapping[df_filtered$states]
#df_filtered <- rename(df_filtered, region = states)



# For each year, select top 10 and bottom 10 states in terms of life expectancy
df_selected <- df_filtered %>%
  group_by(Year) %>%
  top_n(5, Life_Expectancy) %>%
  bind_rows(df_filtered %>%
              group_by(Year) %>%
              top_n(-5, Life_Expectancy))

# Print the tables
df_selected_1959 <- df_selected %>%
  filter(Year == 1959) %>%
  arrange(desc(Life_Expectancy))

df_59rank_2019 <- df_filtered %>%
  filter(Year == 2019, states %in% df_selected_1959$states)

df_59rank <- bind_rows(df_selected_1959, df_59rank_2019)

df_selected_2019 <- df_selected %>%
  filter(Year == 2019) %>%
  arrange(desc(Life_Expectancy))

df_19rank_1959 <- df_filtered %>%
  filter(Year == 1959, states %in% df_selected_2019$states)

df_19rank <- bind_rows(df_selected_2019, df_19rank_1959)

# we can see a shift from more rural states having the highest e0 in 1959, to the states with larger pops
# in 2019 


df_1959 <- df_filtered %>%
  filter(Year == 1959) %>%
  arrange(desc(Life_Expectancy)) %>%
  rename(region = states)

df_2019 <- df_filtered %>%
  filter(Year == 2019) %>%
  arrange(desc(Life_Expectancy))

## DATASETS CURRENTLY NOT COMPATIBLE SO THIS CODE NOT WORKING
# Join the map and e0 data 
df_1959_merged = merge(us, df_1959, by = "region", all = F)
df_2019_merged = merge(us, df_2019, by = "region", all = F)

# plot 
gg <- ggplot()
gg <- gg + geom_map(data=us, map=us,
                    aes(x=long, y=lat, map_id=region),
                    fill="#ffffff", color="#ffffff", size=0.15)
gg <- gg + geom_map(data=df_1959_merged, map=us,
                    aes(fill=Life_Expectancy, map_id=region),
                    color="#ffffff", size=0.15)
gg <- gg + scale_fill_continuous(low = "lightgrey", high = "black")
gg <- gg + labs(x=NULL, y=NULL)
gg <- gg + coord_map("albers", lat0 = 39, lat1 = 45) 
gg <- gg + theme(panel.border = element_blank())
gg <- gg + theme(panel.background = element_blank())
gg <- gg + theme(axis.ticks = element_blank())
gg <- gg + theme(axis.text = element_blank()) 
gg

gg <- ggplot()
gg <- gg + geom_map(data=us, map=us,
                    aes(x=long, y=lat, map_id=region),
                    fill="#ffffff", color="#ffffff", size=0.15)
gg <- gg + geom_map(data=df_2019_merged, map=us,
                    aes(fill=Life_Expectancy, map_id=region),
                    color="#ffffff", size=0.15)
gg <- gg + scale_fill_continuous(low = "lightgrey", high = "black")
gg <- gg + labs(x=NULL, y=NULL)
gg <- gg + coord_map("albers", lat0 = 39, lat1 = 45) 
gg <- gg + theme(panel.border = element_blank())
gg <- gg + theme(panel.background = element_blank())
gg <- gg + theme(axis.ticks = element_blank())
gg <- gg + theme(axis.text = element_blank())
gg

# create slope graph of this data
library(RColorBrewer)
p1 = ggplot(df_59rank, aes(x = factor(Year), y = Life_Expectancy, group = states)) +
  geom_line(aes(color = states), size = 1) +
  geom_point(aes(color = states), size = 3) +
  scale_color_brewer(palette = "Set3") +
  #scale_color_manual(values = colorRampPalette(c("lightblue", "darkblue"))(10)) + 
  geom_text_repel(aes(label = ifelse(Year == 1959, paste(states, round(Life_Expectancy, 1)), NA)),
                  nudge_x = -0.5, # nudges starting position of label to left or right
                  direction = "y",
                  hjust = 1,
                  box.padding = 0.5,
                  segment.size = 0.2) +
  geom_text_repel(aes(label = ifelse(Year == 2019, paste(states, round(Life_Expectancy, 1)), NA)),
                  nudge_x = 0.5, # nudges starting position of label to left or right
                  direction = "y",
                  hjust = 0,
                  box.padding = 0.5,
                  segment.size = 0.2) +
  scale_x_discrete(expand = c(0.2, 0.5)) +  # Extra space for labels
  theme_minimal() +
  labs(x = "", y = "Life Expectancy", title = "Highest and Lowest Five Performers in 1959 Carried Forward to 2019") +
  theme(legend.position = "none")


p2 = ggplot(df_19rank, aes(x = factor(Year), y = Life_Expectancy, group = states)) +
  geom_line(aes(color = states), size = 1) +
  geom_point(aes(color = states), size = 3) +
  scale_color_brewer(palette = "Set3") +
  #scale_color_manual(values = colorRampPalette(c("lightblue", "darkblue"))(10)) + 
  geom_text_repel(aes(label = ifelse(Year == 1959, paste(states, round(Life_Expectancy, 1)), NA)),
                  nudge_x = -0.5, # nudges starting position of label to left or right
                  direction = "y",
                  hjust = 1,
                  box.padding = 0.5,
                  segment.size = 0.2) +
  geom_text_repel(aes(label = ifelse(Year == 2019, paste(states, round(Life_Expectancy, 1)), NA)),
                  nudge_x = 0.5, # nudges starting position of label to left or right
                  direction = "y",
                  hjust = 0,
                  box.padding = 0.5,
                  segment.size = 0.2) +
  scale_x_discrete(expand = c(0.2, 0.5)) +  # Extra space for labels
  theme_minimal() +
  labs(x = "", y = "Life Expectancy", title = "Highest and Lowest Five Performers in 2019 Linked Back to 1959") +
  theme(legend.position = "none")


grid.arrange(p1, p2, nrow = 1)
