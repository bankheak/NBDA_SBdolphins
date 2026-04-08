## --- Social Plasticity Fitness --- ##

# Set relative path working directory
setwd("../Data") # set working directory

#### PART 1: Wrangling data for reproductive output ####

# Load package
if(!require(dplyr)){install.packages('dplyr'); library(dplyr)} 
if(!require(tidyr)){install.packages('tidyr'); library(tidyr)} 
if(!require(ggplot2)){install.packages('ggplot2'); library(ggplot2)} 

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

# Create a column for Behavioral plasticity

# Counts per HI
hi_counts_by_hab_long <- filtered_data %>%
  filter(!is.na(DiffHI), DiffHI != "None") %>%
  group_by(Code, HAB, DiffHI) %>%
  summarise(
    n_HI_category = n(),
    .groups = "drop"
  )

# First engagement
first_hi_by_code_cat_hab <- filtered_data %>%
  filter(!is.na(DiffHI), DiffHI != "None") %>%
  group_by(Code, DiffHI, HAB) %>%
  summarise(
    first_ts = min(Date, na.rm = TRUE),
    .groups = "drop"
  )


counts_after_first_HI_by_HAB <- filtered_data %>%
  inner_join(
    first_hi_by_code_cat_hab %>%
      rename(HI_start = DiffHI),
    by = c("Code", "HAB")
  ) %>%
  filter(Date >= first_ts) %>%
  group_by(Code, HAB, HI_start) %>%
  summarise(
    n_sightings_after_first = n(),
    .groups = "drop"
  ) %>%
  arrange(Code, HAB, HI_start)

# Now divide the number of HI by the number of sightings
HC_data <- counts_after_first_HI_by_HAB %>%
  left_join(
    hi_counts_by_hab_long,
    by = c(
      "Code"     = "Code",
      "HAB"      = "HAB",
      "HI_start" = "DiffHI"
    )
  )

HC_data$HI_prop <- HC_data$n_HI_category/HC_data$n_sightings_after_first
HC_data$HI_cat <- ifelse(HC_data$HI_prop < 0.2, "Low", 
                         ifelse(HC_data$HI_prop < 0.6, "Med", "High"))

# Years lived as a reproducer

# Estimate number of years lived as a reproducer
## Subset data to only include females
female_data <- subset(ILV_pat, ILV_pat$Sex == "Female")
calv_succ <- subset(filtered_data, filtered_data$Code %in% female_data$Alias)
calv_succ <- na.omit(calv_succ)

## Calculate first born calf to mother
first_calf <- aggregate(BirthYear ~ Mom, data = ILV_pat, FUN = min)
names(first_calf) <- c("Mom", "First_birth_year")

## Merge the files
calv_succ <- calv_succ %>%
  left_join(first_calf, by = c("Code" = "Mom"))

# Calculate reproductive years
calv_succ$repro_years <- ifelse(calv_succ$Year >=
  calv_succ$First_birth_year, calv_succ$Year -
    calv_succ$First_birth_year, 0)


## Number of calves produced

# Count calves per Mom per BirthYear from ILV_pat
calves_by_year <- ILV_pat %>%
  group_by(Mom, BirthYear) %>%
  summarize(
    calves_born = n(),
    .groups = "drop"
  )

calves_cumulative <- calves_by_year %>%
  arrange(Mom, BirthYear) %>%
  group_by(Mom) %>%
  mutate(
    total_calves_to_date = cumsum(calves_born)
  ) %>%
  ungroup()

# Merge the datasets
calv_succ <- calv_succ %>%
  left_join(
    calves_cumulative,
    by = c("Code" = "Mom", "Year" = "BirthYear")
  )
calv_succ$total_calves_to_date[is.na(calv_succ$total_calves_to_date)] <- 0

## Calving success

# Survival until age 4/Total calves
ILV_pat <- ILV_pat %>%
  mutate(
    survived_4 = DeathYear >= (BirthYear + 4) | is.na(DeathYear)
  )

# Moms by surv calves and birthyear
calves_survived4_by_birthyear <- ILV_pat %>%
  filter(survived_4) %>%
  group_by(Mom, BirthYear) %>%
  summarize(
    calves_survived4 = n(),
    .groups = "drop"
  )

calves_survived4_by_year <- calves_survived4_by_birthyear %>%
  mutate(
    Year = BirthYear + 4
  ) %>%
  select(Mom, Year, calves_survived4)

# Cumulate the number of calves survived
calves_survived4_cumulative <- calves_survived4_by_year %>%
  arrange(Mom, Year) %>%
  group_by(Mom) %>%
  mutate(
    n_calves_survived_age4 = cumsum(calves_survived4)
  ) %>%
  ungroup()

# Match this to calv_succ
calv_succ <- calv_succ %>%
  left_join(
    calves_survived4_cumulative,
    by = c("Code" = "Mom", "Year" = "Year")
  )

# Fill for none birthyears
calv_succ$n_calves_survived_age4[
  is.na(calv_succ$n_calves_survived_age4)
] <- 0


# Calculate survival of calves
calv_succ$calv_surv <- calv_succ$n_calves_survived_age4/calv_succ$total_calves_to_date

# Get survival death age

# Align each Code to the corresponding BirthYear via Alias
idx <- match(surv_succ$Code, ILV_pat$Alias)

# Pull DeathYear for each row in surv_succ
death_year <- ILV_pat$DeathYear[idx]

# Assign death year
surv_succ$Death <- death_year

# Now reformat data
calv_data <- calv_succ %>%
  group_by(Code, HAB) %>%
  summarize(
    max_total_calves = max(total_calves_to_date, na.rm = TRUE),
    max_surv4_calves = max(n_calves_survived_age4, na.rm = TRUE),
    max_repro_years = {
      x <- repro_years
      if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
    }
    ,
    .groups = "drop"
  )

# Align each Code to the corresponding IDs for HC Behavior
calv_data <- calv_data %>%
  left_join(
    HC_data %>% select(Code, HI_start),
    by = "Code"
  ) %>%
  rename(HC = HI_start)

# Fix the NAs
calv_data$max_repro_years[is.na(calv_data$max_repro_years)] <- 0
calv_data$HC[is.na(calv_data$HC)] <- "None"

calv_data_no_zero <- subset(calv_data, max_repro_years > 0)

# Visual
ggplot(calv_data_no_zero, aes(x = HC, y = max_surv4_calves, fill = HC)) +
  geom_boxplot(outlier.alpha = 0.5, width = 0.7) +
  geom_jitter(width = 0.12, alpha = 0.5, size = 1.6) +
  facet_wrap(~ HAB, nrow = 1) +
  labs(
    title = "Surviving calves per mom (age ≥ 4) by HC within HAB periods",
    x = "HC category",
    y = "Surviving calves per mom"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank())


# Subset calv_succ
write.csv(calv_succ, "calv_succ.csv")
write.csv(surv_succ, "surv_succ.csv")


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

# Reproductive Success Model

# Read in file
calv_succ <- read.csv("calv_succ.csv") # Reproductive success

# Prepare reproductive dataset
repro_succ_df <- data.frame(Repro_Success = calv_succ$num_calves,
                         Total_Repro_opp = calv_succ$repro_liv,
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
              family = binomial(), data = repro_succ_df)

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


# Survival Model

# Read in data
surv_succ <- read.csv("surv_succ.csv") # Survival

# Prepare reproductive dataset
surv_df <- data.frame(Death = surv_succ$Death,
                      Age  = surv_succ$age,
                      ID = surv_succ$Code,
                      HAB = surv_succ$HAB,
                      Behav_plas_bg = surv_succ$BG_flag,
                      Behav_plas_fg = surv_succ$FG_flag,
                      Behav_plas_sd = surv_succ$SD_flag)

# Survival Opportunities
surv_df$Year <- ifelse(surv_df$HAB == "Before", 2001, 
                       ifelse(surv_df$HAB == "During", 2007, 2014))
surv_df$Surv_success <- ifelse(surv_df$Death > surv_df$Year, 1, 0)

# Get rid of data without survival
surv_df <- na.omit(surv_df)

# Run model
fit_sc <- brm(Surv_success | trials(Age) ~
                Behav_plas_bg * HAB + 
                Behav_plas_fg * HAB +
                Behav_plas_sd * HAB + 
                (1 | ID),
              chains = 4, iter = 4000, warmup = 3000, 
              family = binomial(), data = surv_df)

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
