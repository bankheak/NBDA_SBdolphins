## --- Social Plasticity Fitness --- ##

# Set relative path working directory
setwd("../Data") # set working directory

#### PART 1: Wrangling data ####

# Load package
if(!require(dplyr)){install.packages('dplyr'); library(dplyr)} 
if(!require(tidyr)){install.packages('tidyr'); library(tidyr)} 

# Read in data
filtered_data <- read.csv("filtered_data.csv")
# Find the demographics of the population
ILV_all <- read.csv("Individuals_Residency.csv", header=TRUE, sep=",")

# Fix birthyear
ILV_all$BirthYear <- as.numeric(ILV_all$BirthYear)
ILV_all$BirthYear <- ifelse(is.na(ILV_all$BirthYear), 1985, ILV_all$BirthYear)

# Subset data for those in filtered data
Codes <- unique(filtered_data$Code)
ILV_pat <- subset(ILV_all, Alias %in% Codes)

# Create HAB column
filtered_data$HAB <- ifelse(filtered_data$Year < 2001, "Before",
                            ifelse(filtered_data$Year < 2007, "During",
                                   "After"))

# Create calving success

## Subset data to only include females
female_data <- subset(ILV_pat, ILV_pat$Sex == "Female")
calv_succ <- subset(filtered_data, filtered_data$Code %in% ILV_pat$Alias)
calv_succ <- na.omit(calv_succ)

## Years lived as a reproducer

# 1) Align each Code to the corresponding BirthYear via Alias
idx <- match(calv_succ$Code, ILV_pat$Alias)

# 2) Pull BirthYear for each row in calv_succ
birth_year <- ILV_pat$BirthYear[idx]

# 3) Compute age
calv_succ$age <- calv_succ$Year - birth_year

# 4) Estimate number of years lived as a reproducer
calv_succ$repro_liv <- ifelse(calv_succ$age > 5, 
                              calv_succ$Year - (calv_succ$Year - (calv_succ$age - 6)),
                                  NA)
calv_succ$repro_liv[is.na(calv_succ$repro_liv)] <- 0L

## Number of calves produced

# Count calves per Mom per BirthYear from ILV_pat
calves_by_mom_year <- ILV_pat %>%
  # Ensure types are consistent
  mutate(
    Mom = as.character(Mom),
    BirthYear = as.integer(BirthYear)
  ) %>%
  group_by(Mom, BirthYear) %>%
  summarise(n_calves = n(), .groups = "drop")

# Subset years within study period
calves_by_mom_year <- subset(calves_by_mom_year, BirthYear > 1992)

# Add HAB
calves_by_mom_year$HAB <- ifelse(calves_by_mom_year$BirthYear < 2001, "Before",
                            ifelse(calves_by_mom_year$BirthYear < 2007, "During",
                                   "After"))

# Aggregate by HAB
calv_dat <- aggregate(n_calves ~ HAB + Mom, data = calves_by_mom_year, FUN = sum, na.rm = TRUE)

# Merge onto calv_succ; fill NAs with 0
calv_succ$num_calves <- calv_dat$n_calves[
  match(
    paste(calv_succ$HAB, calv_succ$Code),
    paste(calv_dat$HAB,  calv_dat$Mom)
  )]

calv_succ$num_calves[is.na(calv_succ$num_calves)] <- 0L

## Calving success
repro_dat <- aggregate(repro_liv ~ HAB + Code, data = calv_succ, FUN = max, na.rm = TRUE)

# Compute success per HAB + Code
calv_succ$repro_liv <- repro_dat$repro_liv[
  match(
    paste(calv_succ$HAB, calv_succ$Code),
    paste(repro_dat$HAB, repro_dat$Code)
  )]

# Create a column for Behavioral plasticity

# aggregate with any() then coerce to integer
bg_hab <- aggregate(DiffHI == "BG" ~ Code + HAB, data = calv_succ, FUN = any)
bg_hab$BG_flag <- as.integer(bg_hab$`DiffHI == "BG"`)

fg_hab <- aggregate(DiffHI == "FG" ~ Code + HAB, data = calv_succ, FUN = any)
fg_hab$FG_flag <- as.integer(fg_hab$`DiffHI == "FG"`)

sd_hab <- aggregate(DiffHI == "SD" ~ Code + HAB, data = calv_succ, FUN = any)
sd_hab$SD_flag <- as.integer(sd_hab$`DiffHI == "SD"`)

HC_hab <- merge(bg_hab, fg_hab, all = T)
HC_hab <- merge(HC_hab, sd_hab, all = T)

# Attach to calv_succ and fill missing flags with 0
calv_succ <- merge(calv_succ, HC_hab, by = c("Code", "HAB"), all.x = TRUE)
calv_succ$BG_flag[is.na(calv_succ$BG_flag)] <- 0L
calv_succ$FG_flag[is.na(calv_succ$FG_flag)] <- 0L
calv_succ$SD_flag[is.na(calv_succ$SD_flag)] <- 0L

# Subset calv_succ
write.csv(calv_succ, "calv_succ.csv")

#### PART 2: Run brms ####

## Run packages
if(!require(abind)){install.packages('abind'); library(abind)} # array
if(!require(brms)){install.packages('brms'); library(brms)} # For brm modellibrary(coda)
if(!require(bayesplot)){install.packages('bayesplot'); library(bayesplot)} # plot parameters
if(!require(doParallel)){install.packages('doParallel'); library(doParallel)} # Run parallel processing
if(!require(rstan)){install.packages('rstan'); library(rstan)} # To make STAN run faster
if(!require(tidybayes)){install.packages('tidybayes'); library(tidybayes)} # get_variables
if(!require(dplyr)){install.packages('dplyr'); library(dplyr)} 
if(!require(tidyr)){install.packages('tidyr'); library(tidyr)} 

# Help STAN run faster
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Read in file
calv_succ <- read.csv("calv_succ.csv")

# Look at calving success
hist(calv_succ$num_calves)

# Prepare dataset
fitness_df <- data.frame(Success = calv_succ$num_calves,
                         Total_opp = calv_succ$repro_liv,
                         ID = calv_succ$Code,
                         HAB = calv_succ$HAB,
                         Behav_plas_bg = calv_succ$BG_flag,
                         Behav_plas_fg = calv_succ$FG_flag,
                         Behav_plas_sd = calv_succ$SD_flag)

# Set priors
full_priors <- c(
  # Prior for Prop_BG
  set_prior("normal(0, 1)", class = "b", coef = "Behav_plas_bg"),
  # Prior for Prop_FG
  set_prior("normal(0, 1)", class = "b", coef = "Behav_plas_bg"),
  # Prior for Prop_SD
  set_prior("normal(0, 1)", class = "b", coef = "Behav_plas_bg")
)

# Run model
fit_sc <- brm(Success | trials(Total_opp) ~
                Behav_plas_bg * HAB + 
                Behav_plas_fg * HAB +
                Behav_plas_sd * HAB + 
                (1 | ID),
              chains = 4, iter = 4000, warmup = 3000, 
              family = binomial(), data = fitness_df)

saveRDS(fit_sc, "fit_sc.RData")
fit_sc <- readRDS("fit_sc.RData")
summary(fit_sc)

# Check for model convergence
model <- fit_sc
plot(model)
pp_check(model) # check to make sure they line up

# Check estimates
mcmc_intervals(
  as.array(model), 
  pars = c("b_Prop_BG:Period2MDuring_HAB", "b_Prop_BG:Period3MAfter_HAB",
           "b_Prop_BG"),
  prob = 0.95, # 90% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
) +
  theme_minimal() + # Use a minimal theme
  theme(
    text = element_text(family = "sans"), # Set text family
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_blank(), # Remove panel background
    axis.line = element_line(color = "black") # Add axis lines
  )
