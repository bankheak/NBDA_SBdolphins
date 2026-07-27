## --- Bayesian Multi-network Diffusion Analysis --- ##

####=== Measure HC predictability and abundance ===####

# load packages
if(!require(pacman)){install.packages('pacman'); library(pacman)} # to load all packages
pacman::p_load(stringr, tidyr, abind, ggmap, ggplot2,
               dplyr, tidyverse, sf, brms,
               sp, adehabitatHR, lubridate)

# Set relative path working directory
setwd("../Data") 

# Read in data
boat_data <- read.csv("boat_data.csv")
boat_data$Date <- as.Date(boat_data$Date,
                          format = "%Y-%m-%d")

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

