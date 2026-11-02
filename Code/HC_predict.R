## --- Bayesian Multi-network Diffusion Analysis --- ##

####=== Measure HC predictability and abundance ===####

# load packages
if(!require(pacman)){install.packages('pacman'); library(pacman)} # to load all packages
pacman::p_load(stringr, tidyr, abind, ggmap, ggplot2,
               dplyr, tidyverse, sf, brms,
               sp, adehabitatHR, lubridate)

# Set relative path working directory
setwd("../Data") 

# Create fake data
set.seed(123)

n <- 2031

boat_data <- data.frame(
  X.1 = 1:n,
  
  X = sample(1:10000, n, replace = FALSE),
  
  SightingFID = sample(40000:50000, n, replace = TRUE),
  
  Date = sample(
    seq(as.Date("2012-01-01"),
        as.Date("2023-12-31"),
        by = "day"),
    n,
    replace = TRUE
  ),
  
  SurveyNum = sample(7000:12000, n, replace = TRUE),
  
  Sighting = sample(1:100, n, replace = TRUE),
  
  SurveyType = sample(
    c("Photo-id (Pop Mon)",
      "Behavior Survey",
      "Focal Follow"),
    n,
    replace = TRUE,
    prob = c(0.8, 0.1, 0.1)
  ),
  
  Effort = sample(
    c("Photo-id (Pop Mon)",
      "Focal Follow",
      "Behavior Survey"),
    n,
    replace = TRUE
  ),
  
  StartTime = sprintf(
    "%d:%02d",
    sample(7:18, n, replace = TRUE),
    sample(0:59, n, replace = TRUE)
  ),
  
  EndTime = sprintf(
    "%d:%02d",
    sample(7:19, n, replace = TRUE),
    sample(0:59, n, replace = TRUE)
  ),
  
  StartLat = round(runif(n, 27.0, 27.8), 5),
  
  StartLon = round(runif(n, -82.9, -82.2), 5),
  
  HumanInteraction = sample(
    c("0", "F", "I", "L", "F,I", "F,L", "I,L"),
    n,
    replace = TRUE,
    prob = c(0.70, 0.08, 0.06, 0.06, 0.05, 0.03, 0.02)
  ),
  
  SightingGroup = replicate(
    n,
    paste(
      sample(
        c(
          paste0("F", sample(100:300, sample(2:6,1))),
          paste0("FB", sample(1:99, sample(1:3,1))),
          paste0("C", sample(100:999, sample(0:2,1)))
        ),
        replace = FALSE
      ),
      collapse = ","
    )
  ),
  
  SightingNotes = sample(
    c(
      "Several recreational boats observed.",
      "Fishing activity nearby.",
      "No apparent human interaction.",
      "One vessel transiting through area.",
      "Crab traps present.",
      "Anglers actively fishing.",
      ""
    ),
    n,
    replace = TRUE
  ),
  
  X.Boats = pmin(rpois(n, lambda = 4), 25),
  
  X.Lines = pmin(rpois(n, lambda = 0.5), 10),
  
  X.CrabPots = pmin(rpois(n, lambda = 0.3), 20),
  
  Categories.Defined. = sample(
    c("0", "1"),
    n,
    replace = TRUE,
    prob = c(0.3, 0.7)
  ),
  
  stringsAsFactors = FALSE
)


# Diagnose data before model
hist(boat_data$X.CrabPots)
hist(boat_data$X.Lines)
## Negative binomial

# Look at abundance
pot_data <- boat_data[boat_data$X.CrabPots != 0,]
line_data <- boat_data[boat_data$X.Lines != 0,]

abund_data <- bind_rows(
  data.frame(Abundance = pot_data$X.CrabPots, Gear = "Pots"),
  data.frame(Abundance = line_data$X.Lines, Gear = "Lines")
)

# Get average
mean(abund_data$Abundance[abund_data$Gear == "Pots"])
sample_size_pot <- length(abund_data$Abundance[abund_data$Gear == "Pots"])
sd(abund_data$Abundance[abund_data$Gear == "Pots"]) / sqrt(sample_size_pot)

mean(abund_data$Abundance[abund_data$Gear == "Lines"])
sample_size_line <- length(abund_data$Abundance[abund_data$Gear == "Lines"])
sd(abund_data$Abundance[abund_data$Gear == "Lines"]) / sqrt(sample_size_line)

## Fixed-gear foraging

# Run GAM model
pots_model <- brm(
  X.CrabPots ~ 
    t2(StartLon, StartLat) +
    s(as.numeric(Date)) +
    s(StartHour, bs = "cc"),
  family = negbinomial(),
  data = boat_data,
  prior = c(
    prior(normal(0, 2), class = "Intercept"),
    prior(exponential(1), class = "sds")
  ),
  chains = 4, 
  cores = 4,
  iter = 4000
)

saveRDS(pots_model, "pots_model.RData")
pots_model <- readRDS("pots_model.RData")
summary(pots_model)

# Predict uncertainty
pp_check(pots_model)

# Create plot for time series
# Create grid
time_grid <- data.frame(
  StartLon = median(boat_data$StartLon),
  StartLat = median(boat_data$StartLat),
  Date = seq(min(boat_data$Date), max(boat_data$Date), length.out = 100),
  StartHour = median(boat_data$StartHour)
)

time_grid$DateNum <- as.numeric(time_grid$Date)

preds_time <- posterior_epred(pots_model, newdata = transform(time_grid, Date = DateNum))

time_grid$mean <- apply(preds_time, 2, mean)
time_grid$sd   <- apply(preds_time, 2, sd)

# Plot
ggplot(time_grid, aes(x = Date)) +
  geom_line(aes(y = mean), color = "#2c7fb8", linewidth = 1.2) +
  geom_ribbon(
    aes(ymin = mean - sd, ymax = mean + sd),
    fill = "#7fcdbb",
    alpha = 0.4
  ) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  theme_minimal() +
  theme(
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

## Scavenging and depredating foraging

# Run GAM model
lines_model <- brm(
  X.Lines ~ 
    t2(StartLon, StartLat) +
    s(as.numeric(Date)) +
    s(StartHour, bs = "cc"),
  family = negbinomial(),
  data = boat_data,
  prior = c(
    prior(normal(0, 2), class = "Intercept"),
    prior(exponential(1), class = "sds")
  ),
  chains = 4, 
  cores = 4,
  iter = 4000
)

saveRDS(lines_model, "lines_model.RData")
lines_model <- readRDS("lines_model.RData")
summary(lines_model)

# Create plot for time series
preds_time <- posterior_epred(lines_model, newdata = transform(time_grid, Date = DateNum))

time_grid$mean <- apply(preds_time, 2, mean)
time_grid$sd   <- apply(preds_time, 2, sd)

# Plot
ggplot(time_grid, aes(x = Date)) +
  geom_line(aes(y = mean), color = "#d81b60", linewidth = 1.2) +
  geom_ribbon(
    aes(ymin = mean - sd, ymax = mean + sd),
    fill = "#f8b4c9",
    alpha = 0.4
  ) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  theme_minimal() +
  theme(
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

