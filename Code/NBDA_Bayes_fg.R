## --- Bayesian Multi-network Diffusion Analysis --- ##

####=== Fixed-gear Foraging ===####

# load packages
#if(!require(devtools)){install.packages('devtools'); library(devtools)} # To load NBDA
#if(!require(remotes)){install.packages('remotes'); library(remotes)} 
#remotes::install_github("stan-dev/cmdstanr") # If STBayes doesn't download
# Install NBDA package
#devtools::install_github("whoppitt/NBDA")
# install devtools if not already
#devtools::install_github("michaelchimento/STbayes")
#cmdstanr::set_cmdstan_path(cmdstanr::install_cmdstan())
if(!require(pacman)){install.packages('pacman'); library(pacman)} # to load all packages
pacman::p_load(stringr, tidyr, abind, STbayes, ggplot2,
               dplyr, posterior, tidyverse, asnipe, sf,
               sp, adehabitatHR, kinship2, doParallel)

# Set relative path working directory
setwd("../Data") 

#### PART 1: Create Networks ####
# Horizontal network -----------------------------------------

# Read in filtered data
filtered_data <- read.csv("filtered_data.csv")

# Add individual data
# Read ILVs
ILV_all <- read.csv("Individuals_Residency.csv", header=TRUE, sep=",")
ILV_all <- ILV_all[, c("Alias", "HI_Indiv", "Mom", "Sex", "BirthYear")]

# Get rid of ids that aren't in filtered_data
filtered_data <- subset(filtered_data, Code %in% unique(ILV_all$Alias))

# Group each individual by date and sighting
group_data <- filtered_data[,c("Date", "Sighting", "Code", "Year", 
                               "MonthIndex", "DiffHI")]
group_data$Group <- cumsum(!duplicated(group_data[1:2])) # Create sequential group # by date
group_data <- group_data[,c(1, 3:7)] # Subset ID and group #

#' Begin creating the acquisition data
# Find the first time globally where Confirmed_HI == HC
HC <- "FG" # name of target behavior
first_month <- min(filtered_data$MonthIndex[filtered_data$DiffHI == HC])

# Get individuals who had DiffHI == HC in that first month
first_month_individuals <- unique(filtered_data$Code[
  filtered_data$DiffHI == HC & filtered_data$MonthIndex == first_month
])

# Assign yes/no based on that
ILV_all$Demons_HI_forage <- ifelse(ILV_all$Alias %in% 
                                     first_month_individuals, "yes", "no")

#' Subset the data to include observations only from before each 
#' individual's acquisition event
# 1. Exclude demonstrators
exclude_codes <- ILV_all$Alias[ILV_all$Demons_HI_forage == 'yes']

# 2. Identify first observation of HC for each Code
first_HI_index <- tapply(seq_len(nrow(group_data)), group_data$Code, function(idx) {
  hi_rows <- idx[group_data$Confirmed_HI[idx] == 1]
  if (length(hi_rows) > 0) hi_rows[1] else Inf
})

# 3. Split rows by Code
split_rows <- split(seq_len(nrow(group_data)), group_data$Code)

# 4. For excluded Codes: keep all rows
# For others: keep up to first HC, or all if no HC (cutoff = Inf)
keep_rows <- mapply(function(idx, cutoff, code) {
  if (code %in% exclude_codes || is.infinite(cutoff)) {
    idx  # keep everything for excluded Codes or codes with no HC
  } else {
    idx[idx <= cutoff]  # keep up to first HC for others
  }
}, split_rows, first_HI_index[names(split_rows)], names(split_rows))

# 5. Unlist the rows and subset the group_data
keep_rows <- unlist(keep_rows)
group_data <- group_data[keep_rows, ]


#' Create a column for event intervals
# Ensure Date is Date type
group_data <- group_data %>%
  mutate(Date = as.Date(Date))

# 1. Get unique months when HC occurs
event_months <- group_data %>%
  filter(DiffHI == HC) %>%
  arrange(Date) %>%
  group_by(Code) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(month_period = format(Date, "%Y-%m")) %>%
  distinct(month_period) %>%
  pull(month_period)

event_months <- as.Date(paste0(event_months, "-01")) # add first to each month
event_months <- event_months[!is.na(event_months)]  # remove bad values
event_months <- sort(unique(event_months))         # enforce order

# 2. Add event_period to group_data
group_data <- group_data %>%
  mutate(
    month_date = as.Date(format(Date, "%Y-%m-01")),
    event_period = findInterval(month_date, event_months)
  )

# 3. Get rid of time before event
group_data <- subset(group_data, event_period > 0)

# Save data
write.csv(group_data, "group_data_sd.csv")
write.csv(ILV_all, "ILV_FG_subset.csv")

# Now create a list for each time period to calculate association for inter-events
group_data_list <- split(group_data, group_data$event_period)

# Calculate Gambit of the group
create_gbi <- function(list_years) {
  gbi <- list()
  for (i in seq_along(list_years)) {
    
    # Gambit of the group index
    gbi[[i]] <- get_group_by_individual(list_years[[i]][,c("Code", "Group")], data_format = "individuals")
  }
  return(gbi)                                      
}

gbi_fg <- create_gbi(group_data_list) # Run the function
saveRDS(gbi_fg, "gbi_fg.RData") # Save data

# Calculate association matrix
source("../Code/functions.R") # to extract SRI.func
create_nxn <- function(gbi) {
  n.cores <- detectCores()
  system.time({
    registerDoParallel(n.cores)
    nxn <- list()
    for (i in seq_along(gbi)) {
      nxn[[i]] <- as.matrix(SRI.func(gbi[[i]]))
    }                                 
    # End parallel processing
    stopImplicitCluster()
  })
  return(nxn)
}

nxn_fg <- create_nxn(gbi_fg) # Run the function
saveRDS(nxn_fg, "nxn_fg.RData") # Save data

#### PART 2: Calculate acquisition and predictor data ####

# Read in all network data
nxn <- readRDS("nxn_fg.RData")
SRI_vert_all <- readRDS("SRI_vert_all.RData")
ecol_all <- as.array(readRDS("ecol_all.RData"))

#' Duplicate ecol data for inter-event periods
# Read in group data
group_data <- read.csv("group_data_fg.csv")

# Get unique event periods in order
names(ecol_all) <- sort(unique(group_data$Year))

# Group inter-event periods based on their year
event_periods <- group_data %>%
  arrange(event_period, Year) %>%
  group_by(event_period) %>%
  slice(1) %>%   # keeps first occurrence
  ungroup()

# Get the extend of events
T_event <- length(unique(event_periods$event_period))

# Expand ecol_all to match each event period
ecol_all <- lapply(seq_len(nrow(event_periods)), function(i) {
  yr <- event_periods$Year[i]
  ecol_all[[as.character(yr)]]
})

# Get rid of IDs without data in nxn from ecol data
for (i in seq_along(ecol_all)) {
  target_names <- rownames(nxn[[i]])
  
  rn <- rownames(ecol_all[[i]])
  cn <- colnames(ecol_all[[i]])
  
  keep_r <- target_names[target_names %in% rn]
  keep_c <- target_names[target_names %in% cn]
  
  ecol_all[[i]] <- ecol_all[[i]][keep_r, keep_c, drop = FALSE]
}

#' Add zeros to rows that don't have all individuals so that each network 
#' has every individual
# Get all unique IDs across all matrices
total_ids <- unique(unlist(lapply(nxn, colnames)))

# Update each matrix to include all IDs, filling missing rows/columns with zeros
nxn_full <- lapply(nxn, function(mat) {
  # Current IDs in this matrix
  current_ids <- rownames(mat)
  
  # Create a full zero matrix with all IDs
  full_mat <- matrix(0, nrow = length(total_ids), ncol = length(total_ids),
                     dimnames = list(total_ids, total_ids))
  
  # Fill in existing values
  full_mat[current_ids, current_ids] <- mat
  
  return(full_mat)
})


#' Do the same for the vert matrix
# Current IDs in this matrix
current_ids <- rownames(SRI_vert_all)

# Keep only IDs that appear in total_ids
common_ids <- intersect(total_ids, current_ids)

# Create full zero matrix (395 x 395)
full_mat <- matrix(
  0,
  nrow = length(total_ids),
  ncol = length(total_ids),
  dimnames = list(total_ids, total_ids)
)

# Insert existing values (aligned by name)
full_mat[common_ids, common_ids] <- SRI_vert_all[common_ids, common_ids]

# Turn vertical matrix into list
SRI_vert_all <- replicate(T_event, full_mat, simplify = FALSE)
SRI_vert_all <- as.array(SRI_vert_all)

#' Also do this for the ecol matrix
ecol_all <- lapply(ecol_all, function(mat) {
  # Current IDs in this matrix
  current_ids <- rownames(mat)
  
  # Create a full zero matrix with all IDs
  full_mat <- matrix(0, nrow = length(total_ids), ncol = length(total_ids),
                     dimnames = list(total_ids, total_ids))
  
  # Fill in existing values
  full_mat[current_ids, current_ids] <- mat
  
  return(full_mat)
})


#' Put all networks into data frame called edge_list
edge_list <- do.call(rbind, lapply(seq_along(nxn_full), function(t) {
  mat <- nxn_full[[t]]
  ids <- colnames(mat)  # or rownames(mat), assuming square
  upper_idx <- which(upper.tri(mat, diag = TRUE), arr.ind = TRUE)
  
  data.frame(
    focal = ids[upper_idx[, 1]],
    other = ids[upper_idx[, 2]],
    trial = 1,
    assoc = mat[upper_idx],
    time = t
  )
}))


#' Add order of acquisition data
# Read in demographic data
ILV_all <- read.csv("ILV_FG_subset.csv")

# Subset event data to include only edge_list ids
ILV_all <- subset(ILV_all, Alias %in% unique(edge_list$focal))

# Create acquisition data
# 1. Filter group_data for confirmed HC behavior after first aquisition date
HC <- "FG"
first_month <- min(group_data$event_period[group_data$DiffHI == HC])
hi_data <- group_data[group_data$DiffHI == HC & group_data$event_period != first_month, ]

# 2. Get the first time each ID showed the behavior
first_hi_mon <- aggregate(event_period ~ Code, data = hi_data, FUN = min)

# 3. Create a new column for order of acquisition
first_hi_mon$HI_Order_acquisition <- as.numeric(first_hi_mon$event_period)

# 4. Merge this info back into ILV_all
ILV_all <- merge(ILV_all, first_hi_mon[, c("Code", "HI_Order_acquisition")],
                 by.x = "Alias", by.y = "Code", all.x = TRUE)

# 5. Replace NA with 0 for individuals who had the behavior in the first month
ILV_all$HI_Order_acquisition[ILV_all$Demons_HI_forage == "yes"] <- 0

# Change individuals who didn't require behavior to t_end +1
ILV_all$time <- ifelse(is.na(ILV_all$HI_Order_acquisition), 
                       max(na.omit(ILV_all$HI_Order_acquisition)) + 1, 
                       ILV_all$HI_Order_acquisition)

# Add end time
ILV_all$t_end <- max(na.omit(ILV_all$time)) - 1

# Add trial and id
ILV_all$id <- ILV_all$Alias
ILV_all$trial <- 1

# Get rid of unnecessary columns
event_data <- ILV_all[, c("trial", "id", "time", "t_end")]
write.csv(event_data, "event_data_fg.csv") # Save data

# Remove unused individuals
edge_list <- subset(edge_list, focal %in% event_data$id)
edge_list <- subset(edge_list, other %in% event_data$id)

#' Create static and dynamic individual level variables
# Make sex binary
ILV_all$Sex <- ifelse(ILV_all$Sex == "Female", 1, 0)
ILV_all$Sex <- ifelse(is.na(ILV_all$Sex), 1, 0) # Fix NAs

# Fix age
ILV_all$BirthYear <- as.numeric(ILV_all$BirthYear)
ILV_all$BirthYear <- ifelse(is.na(ILV_all$BirthYear), 1985, ILV_all$BirthYear)

# Create constant ILV dataframe
ILV_c <- data.frame(id = ILV_all$Alias,
                    sex = ILV_all$Sex)
write.csv(ILV_c, "ILV_c_fg.csv") # Save data

# Create time varying ILV dataframe
ILV_tv <- data.frame(
  trial = 1,
  id = rep(ILV_all$Alias, each = T_event),
  time = rep(1:T_event, times = length(ILV_all$Alias)),  
  age = (rep(event_periods$Year, times = length(ILV_all$Alias)) - 
           rep(ILV_all$BirthYear, each = T_event)) - 1
)

# Change age to age groups
ILV_tv$age_group <- ifelse(ILV_tv$age >= 10, "adult", 
                           ifelse(ILV_tv$age > 4, "juvenile",
                                  ifelse(ILV_tv$age > 0, "calf", "unborn")))
ILV_tv <- ILV_tv[, -4] # Get rid of unnecessary row
write.csv(ILV_tv, "ILV_tv_fg.csv") # Save data

#' Create static weights based on the proportion of time each individual 
#' engaged in the target behavior
# Catergorize HC into count data
rawHI_diff <- as.data.frame(table(group_data$Code,
                                  group_data$DiffHI))
colnames(rawHI_diff) <- c("Code", "DiffHI", "Freq")

# Categorize ID to Sightings
IDbehav <- as.data.frame(table(group_data$Code))
colnames(IDbehav) <- c("Code", "Sightings")

# Create a frequency count for each HC behavior
get_IDHI <- function(HI, IDbehav_data, rawHI_diff_data) {
  df <- IDbehav_data
  HI_freq <- rawHI_diff_data$Freq[rawHI_diff_data$DiffHI %in% HI]
  df$Behav <- HI_freq[match(df$Code, rawHI_diff_data$Code)]
  colnames(df) <- c("Code", "Sightings", "Behav")
  return(df)
}
IDbehav_HI <- get_IDHI(HC, IDbehav, rawHI_diff) # Run function

# Proportion of Sightings spent in HC
Prop_HI <- function(df, raw_data, HI) {
  
  year_count <- aggregate(Year ~ Code,
                          data = raw_data[raw_data$DiffHI %in% HI, ],
                          FUN = function(x) length(unique(x))
  )
  
  # Create full list of Codes
  all_codes <- unique(raw_data$Code)
  
  # Merge and fill NAs with 0
  year_count <- merge(data.frame(Code = all_codes), year_count, by = "Code", all.x = TRUE)
  year_count$Year[is.na(year_count$Year)] <- 0
  
  # Ensure order matches df
  year_count <- year_count[match(df$Code, year_count$Code), ]
  
  df$HIprop <- as.numeric(df$Behav) / year_count$Year
  df$HIprop[is.na(df$HIprop)] <- 0
  # Keep only 'Code' and 'HIprop' columns
  prop_df <- df[, c('Code', 'HIprop')]
  return(prop_df)
}
prob_HI <- Prop_HI(IDbehav_HI, group_data, HC) # Run function

# Get ids and HC proportions
ids <- unique(prob_HI$Code)
HIProp <- prob_HI$HIprop

# Create dataframe for HC prop weights
HI_matrix <- data.frame(
  trial = 1,
  id = ids,
  t_weight = HIProp)

# Get rid of unused individuals
HI_matrix <- HI_matrix[HI_matrix$id %in% event_data$id, ]
write.csv(HI_matrix, "HI_matrix_fg.csv") # Save data

# Add vertical network to edge_list
edge_list_vert <- do.call(rbind, lapply(seq_along(SRI_vert_all), function(t) {
  mat <- SRI_vert_all[[t]]
  ids <- colnames(mat)  # or rownames(mat), assuming square
  upper_idx <- which(upper.tri(mat, diag = TRUE), arr.ind = TRUE)
  
  data.frame(
    focal = ids[upper_idx[, 1]],
    other = ids[upper_idx[, 2]],
    trial = 1,
    vert = mat[upper_idx],
    time = t
  )
}))

edge_list$vert <- edge_list_vert$vert

# Add ecol network to edge_list
edge_list_ecol <- do.call(rbind, lapply(seq_along(ecol_all), function(t) {
  mat <- ecol_all[[t]]
  ids <- colnames(mat)  # or rownames(mat), assuming square
  upper_idx <- which(upper.tri(mat, diag = TRUE), arr.ind = TRUE)
  
  data.frame(
    focal = ids[upper_idx[, 1]],
    other = ids[upper_idx[, 2]],
    trial = 1,
    ecol = mat[upper_idx],
    time = t
  )
}))

edge_list$ecol <- edge_list_ecol$ecol

# Rearrange
edge_list <- edge_list[, c("focal", "other", "trial", 
                           "assoc", "vert", "ecol", "time")]

saveRDS(edge_list, "edge_list_fg.RData") # Save data


#### PART 3: Aggregate data for model input ####

# Read in event data
event_data <- read.csv("event_data_fg.csv")

# Read in ILV data
ILV_tv <- read.csv("ILV_tv_fg.csv")
ILV_c <- read.csv("ILV_c_fg.csv")

# Read in weighted data
HI_matrix <- read.csv("HI_matrix_fg.csv")

# Input data
data_list <- import_user_STb(
  event_data = event_data,
  networks = edge_list,
  ILV_c = ILV_c,
  ILV_tv = ILV_tv,
  ILVi = c("age_group", "sex"),
  ILVs = c("age_group", "sex"),
  t_weights = HI_matrix
)

saveRDS(data_list, "data_list_fg.RData") # Save data


#### PART 4: Run the model ####

# Input data_list
data_list <- readRDS("data_list_fg.RData")

# Generate the model
model_full <- generate_STb_model(data_list, 
                                 data_type = "continuous_time",
                                 gq = T, 
                                 est_acqTime = T)

# Diagnoses notes
#' Detected N_veff = 0, likihood provides no information
#' Parameters that must be positive / bounded
#' Weak or default priors
#' Hierarchical or hazard-based structure
#' Therefore we need to change priors and baseline learning rate must be positive

# Do this in the model.STAN file
# writeLines(model_full, "model_strong_priors_fg.stan")

# Precompute E for faster computing
# Extract needed pieces
K <- data_list$K
P <- data_list$P
T_max <- data_list$T_max
N_networks <- data_list$N_networks
A <- data_list$A
Z <- data_list$Z

# Precompute E
E <- array(0, dim = c(N_networks, K, T_max, P))

for (network in 1:N_networks) {
  for (trial in 1:K) {
    for (t in 1:T_max) {
      
      A_mat <- A[network, trial, t, , ]
      Z_vec <- Z[trial, t, ]
      
      E[network, trial, t, ] <- A_mat %*% Z_vec
    }
  }
}

# Add to data list
data_list$E <- E

# REMOVE A (important for memory + speed)
data_list$A <- NULL


# Test run

# fit_test <- fit_STb(
#   data_list,
#   "model_strong_priors_fg.stan",
#   chains = 1,
#   iter_warmup = 1000,
#   iter_sampling = 500,
#   adapt_delta = 0.99,
#   refresh = 50
# )

# Fit the model
fit <- fit_STb(
  data_list,
  "model_strong_priors_fg.stan",
  chains = 4,
  parallel_chains = 4,
  cores = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = 0.99,
  refresh = 100
)

fit$cmdstan_diagnose()

STb_save(fit, output_dir = "cmdstan_saves_fg", name="my_first_fit") # Save model


#### PART 5: Summarize Outputs ####

# Read in model output
fit <- readRDS('cmdstan_saves_fg/my_first_fit.rds')

# View output
results <- STb_summary(fit, digits = 3)

#' The most important output are the intrinsic rate (lambda_0), 
#' and the relative strength of social transmission (s), whose 
#' interpretations are the same as the NBDA package. The relative 
#' strength of social transmission (s = s_prime / lambda_0) is generally 
#' what we’re after. %ST for network n is reported as percent_ST[n]. This 
#' is a single-network model, thus percent_ST[1] is the estimated percentage 
#' of events that occurred through social transmission. The [1] refers to 
#' the “assoc” network, as we’ve only given a single network. If you fit a 
#' multi-network model, all networks will have an estimate. For a number of 
#' reasons, STbayes actually fits lambda_0 and social transmission 
#' rate (s_prime) on the log scale. The linear transformation of s_prime 
#' itself usually isn’t reported and is excluded from the output, but you 
#' could calculate it yourself from the fit.

#' Posterior Predictive Checks
# Read in event data
event_data <- read.csv("event_data_fg.csv")

# Cumulative distribution curve
plot_data_obs <- event_data %>%
  filter(time > 0, time <= t_end) %>% # exclude demonstrators (time == 0) and censored (time > t_end)
  group_by(trial) %>%
  arrange(time, .by_group = TRUE) %>%
  mutate(
    cum_prop = row_number() / n(),
    type = "observed"
  ) %>%
  select(trial, time, cum_prop, type) %>%
  ungroup()

# add in 0,0 starting point
plot_data_obs <- bind_rows(
  plot_data_obs,
  plot_data_obs %>%
    distinct(trial) %>%
    mutate(time = 0, cum_prop = 0, type = "observed")
) %>%
  arrange(trial, time)

draws_df <- as_draws_df(fit$draws(variables = "acquisition_time", inc_warmup = FALSE))

# pivot longer
ppc_long <- draws_df %>%
  select(starts_with("acquisition_time[")) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("trial", "ind"),
    names_pattern = "acquisition_time\\[(\\d+),(\\d+)\\]",
    values_to = "time"
  ) %>%
  mutate(
    trial = as.integer(trial),
    ind = as.integer(ind),
    draw = rep(1:(nrow(draws_df)), 
               each = length(unique(.$trial)) * length(unique(.$ind)))
  )

# thin sample for plotting
sample_idx <- sample(c(1:max(ppc_long$draw)), 100)
ppc_long <- ppc_long %>% filter(draw %in% sample_idx)

# build cumulative curves per draw
plot_data_ppc <- ppc_long %>%
  group_by(draw, trial, time) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(draw, trial) %>%
  arrange(time) %>%
  mutate(cum_prop = cumsum(n) / data_list$Q)

# add in 0,0 starting point
plot_data_ppc <- bind_rows(
  plot_data_ppc,
  plot_data_ppc %>%
    distinct(trial, draw) %>%
    mutate(time = 0, cum_prop = 0, type = "ppc")
) %>%
  arrange(trial, time)

# plot it
ggplot() +
  geom_line(data = plot_data_ppc, 
            aes(x = time, y = cum_prop, 
                group = interaction(draw, trial)), alpha = .1) +
  geom_line(data = plot_data_obs, aes(x = time, y = cum_prop), linewidth = 1) +
  labs(x = "Time", y = "Cumulative proportion informed", color = "Trial") +
  theme_minimal()
# The model reproduces the overall learning dynamics at the population level.

# Estimated versus observed 
acqdata = extract_acqTime(fit, data_list)

ggplot(acqdata, aes(x = observed_time, y = median_time)) +
  geom_segment(
    aes(x = observed_time, xend = observed_time, 
        y = median_time, yend = observed_time),
    color = "red",
    alpha = 0.3) +
  geom_point(size = 2) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  labs(x = "Observed time", y = "Estimated time") +
  theme_minimal()
#' This is exactly what to expect in a hazard‑based model:
#' Early learners are well constrained
#' Later learners accumulate uncertainty across time

# Make a clean summary of results
results_clean <- results %>%
  mutate(
    # Classify parameter type
    type = case_when(
      str_detect(Parameter, "^beta_ILVi") ~ "Individual",
      str_detect(Parameter, "^beta_ILVs") ~ "Social",
      str_detect(Parameter, "^s\\[") ~ "Network strength",
      str_detect(Parameter, "^percent_ST") ~ "Social contribution",
      str_detect(Parameter, "lambda") ~ "Baseline",
      TRUE ~ "Other"
    ),
    
    # Extract variable names
    variable = case_when(
      str_detect(Parameter, "sex") ~ "Sex",
      str_detect(Parameter, "HAB") ~ "HAB",
      str_detect(Parameter, "age_group") ~ "Age",
      str_detect(Parameter, "^s\\[") ~ "Network",
      str_detect(Parameter, "^percent_ST") ~ "Network",
      TRUE ~ Parameter
    ),
    
    # Extract level (for age groups or networks)
    level = str_extract(Parameter, "\\d+") %>% as.numeric(),
    
    # Clean labels for plotting
    label = case_when(
      str_detect(Parameter, "beta_ILVi") ~ paste0("Ind_", variable, "_", level %||% ""),
      str_detect(Parameter, "beta_ILVs") ~ paste0("Soc_", variable, "_", level %||% ""),
      str_detect(Parameter, "^s\\[") ~ paste0("s_net", level),
      str_detect(Parameter, "^percent_ST") ~ paste0("ST_net", level),
      TRUE ~ Parameter
    )
  )

results_clean <- results_clean %>%
  mutate(label = case_when(
    str_detect(Parameter, "beta_ILVi_sex") ~ "Ind: Sex",
    str_detect(Parameter, "beta_ILVs_sex") ~ "Soc: Sex",
    str_detect(Parameter, "beta_ILVi_age_group") ~ paste0("Ind: Age ", level),
    str_detect(Parameter, "beta_ILVs_age_group") ~ paste0("Soc: Age ", level),
    TRUE ~ label
  ))

# Only network variables
Neffects_only <- results_clean %>%
  filter(variable %in% c("Network"))

# Only ILV variables
ILVeffects_only <- results_clean %>%
  filter(variable %in% c("Sex", "Age"))


#' Visual Plot of effect sizes
# Networks
ggplot(Neffects_only, aes(x = label, y = Median)) +
  geom_col(fill = "lightcoral", width = 0.7) +
  
  geom_errorbar(
    aes(ymin = CI_Lower, ymax = CI_Upper),
    width = 0.15,
    size = 0.6
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  theme_classic(base_size = 14) +   # removes gridlines
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    
    # Clean look
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  
  labs(
    x = NULL,
    y = "Effect size",
    title = "Effects of Network Variables"
  )

# ILV
ggplot(ILVeffects_only, aes(x = label, y = Median)) +
  geom_col(fill = "lightcoral", width = 0.7) +
  
  geom_errorbar(
    aes(ymin = CI_Lower, ymax = CI_Upper),
    width = 0.15,
    size = 0.6
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  theme_classic(base_size = 14) +   # removes gridlines
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    
    # Clean look
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  
  labs(
    x = NULL,
    y = "Effect size",
    title = "Effects of Individual and Social Variables"
  )


#' Plot learning hazard over time
# Read in data_list
data_list <- readRDS("data_list_sd.RData")

# 1. Extract posterior draws
draws_df <- fit$draws(format = "df")

lambda_0   <- draws_df$lambda_0
s_prime_1  <- draws_df$`s_prime[1]`
s_prime_2  <- draws_df$`s_prime[2]`
s_prime_3  <- draws_df$`s_prime[3]`

# 2. Supply your trial sequence and informed counts per network 
trials <- 1:data_list$T_max   

# Social input for network k at time t = mean over naive focals of: sum_others(A[k,1,t,focal,other] * Zn[t,other])
n_networks <- data_list$N_networks  # 3 networks
T_max      <- data_list$T_max       # 42 inter-event periods
A          <- data_list$A           # [3, 1, 42, 1021, 1021]
Zn         <- data_list$Zn[1, , ]  # [42, 1021]

# Mean social input per network per trial
social_input <- matrix(NA, nrow = T_max, ncol = n_networks)

for (t in 1:T_max) {
  zn_t <- Zn[t, ]          # informed status of all individuals at time t
  naive_t <- which(zn_t == 0)  # only naive individuals can learn
  
  for (k in 1:n_networks) {
    A_kt <- A[k, 1, t, , ]   # [1021 focal x 1021 other] edge weights at time t
    
    # For each individual, social input = sum of edge weights to informed others
    social_input_per_ind <- A_kt %*% zn_t   # [1021 x 1] vector
    
    # Average over naive individuals only (they're the ones at risk)
    if (length(naive_t) > 0) {
      social_input[t, k] <- mean(social_input_per_ind[naive_t])
    } else {
      social_input[t, k] <- 0
    }
  }
}

social_input_df <- as.data.frame(social_input)
colnames(social_input_df) <- data_list$network_names  # "assoc", "vert", "ecol"
social_input_df$trial <- 1:T_max

# 3. Compute hazard rates per trial, across all posterior draws
n_draws <- nrow(draws_df)

il_hazard_mat <- matrix(rep(draws_df$lambda_0, T_max), nrow = n_draws, ncol = T_max)

sl_hazard_mat <- matrix(NA, n_draws, T_max)
for (t in 1:T_max) {
  sl_hazard_mat[, t] <- draws_df$`s_prime[1]` * social_input_df$assoc[t] +
    draws_df$`s_prime[2]` * social_input_df$vert[t]  +
    draws_df$`s_prime[3]` * social_input_df$ecol[t]
}

# 4. Summarise posterior (median + 90% CI)
summarise_hazard <- function(mat, trials) {
  tibble(
    trial  = trials,
    median = apply(mat, 2, median),
    lo     = apply(mat, 2, quantile, 0.055),
    hi     = apply(mat, 2, quantile, 0.945)
  )
}

il_summary <- summarise_hazard(il_hazard_mat, trials) %>% mutate(type = "Individual Learning")
sl_summary <- summarise_hazard(sl_hazard_mat, trials) %>% mutate(type = "Social Learning")

plot_df <- bind_rows(il_summary, sl_summary)

# Social Learning - Individual Learning difference for crossover shading
diff_df <- tibble(
  trial      = trials,
  diff_med   = apply(sl_hazard_mat - il_hazard_mat, 2, median),
  diff_lo    = apply(sl_hazard_mat - il_hazard_mat, 2, quantile, 0.055),
  diff_hi    = apply(sl_hazard_mat - il_hazard_mat, 2, quantile, 0.945)
)

crossover_trial <- diff_df %>%
  filter(diff_med >= 0) %>%
  slice(1) %>%
  pull(trial)

# 5. Plot
pal <- c("Individual Learning" = "#E07B54", "Social Learning" = "#4A90D9")

ggplot(plot_df, aes(x = trial, colour = type, fill = type)) +
  # Crossover shading
  geom_vline(xintercept = crossover_trial, linetype = "dashed",
             colour = "grey40", linewidth = 0.6) +
  annotate("text", x = crossover_trial + 0.5, y = Inf,
           label = paste0("Crossover\n(trial ", crossover_trial, ")"),
           vjust = 1.5, hjust = 0, size = 3.2, colour = "grey30") +
  # CI ribbons
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
  # Median lines
  geom_line(aes(y = median), linewidth = 1) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(
    title    = "Social vs Individual Learning Hazard Over Trials",
    subtitle = "Shaded bands = 90% posterior credible interval",
    x        = "Trial",
    y        = "Hazard rate",
    colour   = NULL, fill = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position  = c(0.15, 0.88),
    legend.background = element_blank(),
    plot.title       = element_text(face = "bold")
  )

# Input group data to determine year of crossover
group_data <- read.csv("group_data_fg.csv")
unique(group_data$Year[group_data$event_period == 26])
