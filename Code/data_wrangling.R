## --- Bayesian Multi-network Diffusion Analysis --- ##

#### DATA WRANGLING ####===

# load packages
if(!require(pacman)){install.packages('pacman'); library(pacman)} # to load all packages
pacman::p_load(stringr, tidyr, ggmap, STbayes, ggplot2,
               dplyr, lubridate, tidyverse, ggExtra, sf,
               sp, adehabitatHR, hrbrthemes, doParallel, brms)

# Set relative path working directory
setwd("../Data") 
# set working directory

#### PART 1: Wrangling data ####

# create fake data
set.seed(123)

# Demographic data #
n <- 1142

# Create a pool of maternal IDs
mom_pool <- paste0(
  "F",
  sample(100:999, 150, replace = FALSE)
)

ILV_all <- data.frame(
  Alias = as.character(sample(1000:9999, n, replace = FALSE)),
  
  HI_Indiv = sample(
    c(0, 1),
    n,
    replace = TRUE,
    prob = c(0.75, 0.25)
  ),
  
  Mom = sample(
    mom_pool,
    n,
    replace = TRUE
  ),
  
  Sex = sample(
    c("Female", "Male", NA),
    n,
    replace = TRUE,
    prob = c(0.45, 0.40, 0.15)
  ),
  
  BirthYear = as.character(
    sample(
      1980:2020,
      n,
      replace = TRUE
    )
  ),
  
  stringsAsFactors = FALSE
)

write.csv(ILV_all, "ILV_all.csv") # Save data


# Observation data #
n <- 10000  # smaller example dataset

filtered_data <- data.frame(
  X.1 = 1:n,
  X = 1:n,
  
  SightingFID = sample(10000:30000, n, replace = TRUE),
  
  Date = as.character(
    sample(
      seq(as.Date("1995-01-01"),
          as.Date("2023-12-31"),
          by = "day"),
      n,
      replace = TRUE
    )
  ),
  
  StartTime = sprintf(
    "%02d:%02d",
    sample(6:18, n, replace = TRUE),
    sample(c(0,15,30,45), n, replace = TRUE)
  ),
  
  EndTime = sprintf(
    "%02d:%02d",
    sample(6:19, n, replace = TRUE),
    sample(c(0,15,30,45), n, replace = TRUE)
  ),
  
  SurveyNum = sample(1000:6000, n, replace = TRUE),
  
  SurveyType = sample(
    c("Photo-id (Pop Mon)",
      "No Data",
      "Behavior Survey",
      "Health Assessment"),
    n,
    replace = TRUE,
    prob = c(0.7, 0.1, 0.15, 0.05)
  ),
  
  SurveyArea = sample(
    c("SB","NB","CB",""),
    n,
    replace = TRUE,
    prob = c(0.5,0.2,0.2,0.1)
  ),
  
  Subarea = sample(
    c("AMG","RBC","SAR","VEN","LKB","PAS"),
    n,
    replace = TRUE
  ),
  
  StartLat = round(runif(n, 27.1, 27.9), 5),
  StartLon = round(runif(n, -82.9, -82.2), 5),
  
  Sighting = as.character(sample(1:1000, n, replace = TRUE)),
  
  Effort = sample(
    c("Photo-id (Pop Mon)",
      "Focal Follow",
      "No Data"),
    n,
    replace = TRUE
  ),
  
  TDolphinBest = sample(1:30, n, replace = TRUE),
  TCalfBest = sample(0:10, n, replace = TRUE),
  
  Behaviors = sample(
    c("Travel",
      "Forage",
      "Socialize",
      "Rest",
      "Mill"),
    n,
    replace = TRUE
  ),
  
  HumanInteraction = sample(
    c("999","FG","HL","Boating","Feeding"),
    n,
    replace = TRUE,
    prob = c(0.8,0.08,0.08,0.03,0.01)
  ),
  
  Code = paste0(
    sample(c("F","FB","C","D"), n, replace = TRUE),
    sample(100:999, n, replace = TRUE)
  ),
  
  SightingNotes = "",
  
  ConfHI = sample(
    c("0","FG","HL"),
    n,
    replace = TRUE,
    prob = c(0.9,0.05,0.05)
  ),
  
  PossHI = sample(
    c("0","FG","HL"),
    n,
    replace = TRUE,
    prob = c(0.85,0.075,0.075)
  ),
  
  Confirmed_HI = sample(
    c(0,1),
    n,
    replace = TRUE,
    prob = c(0.9,0.1)
  )
)

# Derive time-related variables
filtered_data$Year <- as.integer(format(as.Date(filtered_data$Date), "%Y"))
filtered_data$MonthIndex <- sample(1:348, n, replace = TRUE)

filtered_data$DiffHI <- sample(
  c("None","FG","SD","Both"),
  n,
  replace = TRUE,
  prob = c(0.7,0.15,0.1,0.05)
)

filtered_data$sixmon <- sample(1:58, n, replace = TRUE)

write.csv(filtered_data, "filtered_data.csv") # Save data


#### PART 2: Inspect data on human-centric behavior ####

# Read in data
filtered_data <- read.csv("filtered_data.csv")
group_data_sd <- read.csv("group_data_sd.csv")
group_data_fg <- read.csv("group_data_fg.csv")

# How many total dolphins are there?
length(unique(filtered_data$Code))

# How many dolphins with each HC are there?
length(unique(group_data_sd$Code[group_data_sd$DiffHI == "SD"]))
length(unique(group_data_fg$Code[group_data_fg$DiffHI == "FG"]))

# See which individuals overlap
presence_df <- do.call(
  rbind,
  lapply(
    split(filtered_data$DiffHI, filtered_data$Code),
    function(x) c(
      FG = "FG" %in% x,
      SD = "SD" %in% x
    )
  )
)

presence_df <- as.data.frame(presence_df)
presence_df$Code <- rownames(presence_df)
rownames(presence_df) <- NULL

# Get the number of times either are both true
sum(presence_df$SD == T & presence_df$FG == F) 
# Combine BG and SD but not FG

# What percentage does the behavior take up for foraging
# Subset for only dolphins seen foraging
forage_data <- subset(
  filtered_data,
  grepl("Feed|pFeed", Behaviors)
)

# Separate data before and after HAB
forage_data_Pre <- subset(forage_data, Year < 2005)
forage_data_Post <- subset(forage_data, Year > 2004)

# Catergorize into raw data
## Pre
rawHI_diff_pre <- as.data.frame(table(forage_data_Pre$Code,
                                  forage_data_Pre$DiffHI))
colnames(rawHI_diff_pre) <- c("Code", "DiffHI", "Freq")
## Post
rawHI_diff_post <- as.data.frame(table(forage_data_Post$Code,
                                  forage_data_Post$DiffHI))
colnames(rawHI_diff_post) <- c("Code", "DiffHI", "Freq")

# Categorize ID to Sightings
## Pre
IDbehav_pre <- as.data.frame(table(forage_data_Pre$Code))
colnames(IDbehav_pre) <- c("Code", "Sightings")
## Post
IDbehav_post <- as.data.frame(table(forage_data_Post$Code))
colnames(IDbehav_post) <- c("Code", "Sightings")

# Divide the two
## Pre
HC_props_pre <- rawHI_diff_pre %>%
  left_join(IDbehav_pre %>% select(Code, Sightings), by = "Code") %>%
  mutate(Prop = Freq / Sightings)
# Get rid of individuals with zero
HC_props_pre <- HC_props_pre[HC_props_pre$Prop != 0,]
## Post
HC_props_post <- rawHI_diff_post %>%
  left_join(IDbehav_post %>% select(Code, Sightings), by = "Code") %>%
  mutate(Prop = Freq / Sightings)
# Get rid of individuals with zero
HC_props_post <- HC_props_post[HC_props_post$Prop != 0,]

# Distribution of SD foraging
## Pre
hist(
  HC_props_pre$Prop[HC_props_pre$DiffHI == "SD"],
  breaks = seq(0, 1, by = 0.05)
)
median(HC_props_pre$Prop[HC_props_pre$DiffHI == "SD"])
## Post
hist(
  HC_props_post$Prop[HC_props_post$DiffHI == "SD"],
  breaks = seq(0, 1, by = 0.05)
)
median(HC_props_post$Prop[HC_props_post$DiffHI == "SD"])

# Distribution of FG foraging
## Pre
hist(
  HC_props_pre$Prop[HC_props_pre$DiffHI == "FG"],
  breaks = seq(0, 1, by = 0.05)
)
median(HC_props_pre$Prop[HC_props_pre$DiffHI == "FG"])
## Post
hist(
  HC_props_post$Prop[HC_props_post$DiffHI == "FG"],
  breaks = seq(0, 1, by = 0.05)
)
median(HC_props_post$Prop[HC_props_post$DiffHI == "FG"])

# If it has FG or SD then don't assign None
cleaned <- filtered_data %>%
  group_by(Code, sixmon) %>%
  mutate(has_real_cat = any(DiffHI %in% c("FG", "SD", "BG"))) %>%
  ungroup() %>%
  filter(!(has_real_cat & DiffHI == "None")) %>%  # drop None only where a real category exists
  select(-has_real_cat)

# Categorize DiffHI to IDs
HI_diff <- as.data.frame(table(cleaned$Code, cleaned$DiffHI, cleaned$sixmon))
colnames(HI_diff) <- c("Code", "DiffHI", "sixmon", "Freq")

# Add dates
seq_dates <- seq(from = as.Date("1995-01-01"),
                 to   = as.Date("2014-12-01"),
                 by   = "6 month")

# Filter to HCs only
HI_diff$ID_count <- ifelse(HI_diff$Freq > 0, 1, 0)
fg <- HI_diff %>% filter(DiffHI == "FG")
fg <- aggregate(ID_count ~ sixmon, data = fg, sum)
fg$HC <- "FG"
sd <- HI_diff %>% filter(DiffHI == "SD")
sd <- aggregate(ID_count ~ sixmon, data = sd, sum)
sd$HC <- "SD"
None <- HI_diff %>% filter(DiffHI == "None")
None <- aggregate(ID_count ~ sixmon, data = None, sum)
None$HC <- "None"

# Now combine
dfs <- list(fg, sd, None)
time_series_data <- Reduce(function(x, y) merge(x, y, all = TRUE), dfs)
time_series_data$Date <- seq_dates[time_series_data$sixmon]

# Plot the different HC by year

# FG
ggplot(fg, aes(x = factor(sixmon), y = ID_count, group = 1)) +
  geom_line(color = "#4C78A8", linewidth = 1) +
  geom_point(size = 3, color = "#4C78A8") +
  labs(
    x = "sixmon",
    y = "ID_count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
# Could do 6 months and start 1994

# SD
ggplot(sd, aes(x = factor(sixmon), y = ID_count, group = 1)) +
  geom_line(color = "#4C78A8", linewidth = 1) +
  geom_point(size = 3, color = "#4C78A8") +
  labs(
    x = "sixmon",
    y = "ID_count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
# Could do 6 months and start 1995

# Exclude None
time_series_data <- time_series_data[time_series_data$HC != 'None',]

# Make a line plot of naive dolphins and begging, depredating and fixed-gear 
ggplot(time_series_data, aes(x = Date, y = ID_count, color = HC)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c(
    "springgreen4",
    "lightcoral"
  )) +
  xlab("Date") +
  ylab("Number of Individuals") +
  labs(color = "HC Category") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )


#### PART 3: Look at demographics ####

# Read in data
filtered_data <- read.csv("filtered_data.csv")

# Find the demographics of the population
ILV_all <- read.csv("Individuals_Residency.csv", header=TRUE, sep=",")
ILV_all <- ILV_all[, c("Alias", "HI_Indiv", "Mom", "Sex", "BirthYear")]

Codes <- unique(filtered_data$Code)
ILV_pat <- subset(ILV_all, Alias %in% Codes)

# Assign ages
## Change Alias to Code
colnames(ILV_pat) <- c("Code", "HI", "Mom", "Sex", "BirthYear")

## Add birthyear
filtered_data <- merge(
  filtered_data,
  ILV_pat[, c("Code", "BirthYear")],
  by = "Code",
  all.x = TRUE
)

## Add sex
ILV_pat$Sex <- ifelse(ILV_pat$Sex == "Probable Female", "Female", 
                      ifelse(ILV_pat$Sex == "Probable Male", "Male", ILV_pat$Sex))
filtered_data <- merge(
  filtered_data,
  ILV_pat[, c("Code", "Sex")],
  by = "Code",
  all.x = TRUE
)

## Add ages
filtered_data$age <- as.numeric(filtered_data$Year) - as.numeric(filtered_data$BirthYear)
filtered_data$age_class <- ifelse(
  filtered_data$age >= 10, "adult",
  ifelse(filtered_data$age > 4, "juvenile", "calf")
)

# How many dolphins from the individual analyses
## SD
group_data_sd <- read.csv("group_data_sd.csv")

group_data_sd <- merge(
  group_data_sd,
  ILV_pat[, c("Code", "Sex")],
  by = "Code",
  all.x = TRUE
)

### Sex
length(unique(group_data_sd$Code[group_data_sd$Sex == "Female"]))
length(unique(group_data_sd$Code[group_data_sd$Sex == "Male"]))
length(unique(group_data_sd$Code[is.na(group_data_sd$Sex)]))
### HC
length(unique(group_data_sd$Code[group_data_sd$DiffHI == "SD"]))
### Both
length(unique(group_data_sd$Code[group_data_sd$DiffHI == "SD" & group_data_sd$Sex == "Female"]))
length(unique(group_data_sd$Code[group_data_sd$DiffHI == "SD" & group_data_sd$Sex == "Male"]))
length(unique(group_data_sd$Code[group_data_sd$DiffHI == "SD" & is.na(group_data_sd$Sex)]))

## FG
group_data_fg <- read.csv("group_data_fg.csv")

group_data_fg <- merge(
  group_data_fg,
  ILV_pat[, c("Code", "Sex")],
  by = "Code",
  all.x = TRUE
)

### Sex
length(unique(group_data_fg$Code[group_data_fg$Sex == "Female"]))
length(unique(group_data_fg$Code[group_data_fg$Sex == "Male"]))
length(unique(group_data_fg$Code[is.na(group_data_fg$Sex)]))
### HC
length(unique(group_data_fg$Code[group_data_fg$DiffHI == "FG"]))
### Both
length(unique(group_data_fg$Code[group_data_fg$DiffHI == "FG" & group_data_fg$Sex == "Female"]))
length(unique(group_data_fg$Code[group_data_fg$DiffHI == "FG" & group_data_fg$Sex == "Male"]))
length(unique(group_data_fg$Code[group_data_fg$DiffHI == "FG" & is.na(group_data_fg$Sex)]))

# Count each dolphin behavior
count_data <- as.data.frame(
  table(
    id = filtered_data$Code,
    age_class = filtered_data$age_class,
    behavior  = filtered_data$DiffHI,
    sex = filtered_data$Sex
  )
)

# Remove zero-count combinations
count_data <- subset(count_data, Freq > 0)

# Get rid of frequency
count_data <- count_data[,-5]

# Subset without None
count_data <- subset(count_data, behavior != "None")
count_data <- subset(count_data, behavior != "BG")
count_data <- subset(count_data, !is.na(age_class))

# Count each dolphin behavior per age-class
plot_data <- as.data.frame(
  table(
    age_class = count_data$age_class,
    behavior  = count_data$behavior,
    sex = count_data$sex
  )
)

# Make proportions
plot_data$prop <- ave(
  plot_data$Freq,
  plot_data$behavior,
  FUN = function(x) x / sum(x)
)

plot_data$behavior <- factor(
  plot_data$behavior,
  levels = c("SD", "FG")
)

# Get rid of missing data
plot_data <- na.omit(plot_data)

# Plot the data
ggplot(
  plot_data,
  aes(
    x = age_class,
    y = prop,
    fill = interaction(behavior, sex)
  )
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(~ behavior) +
  labs(
    x = "Age Class",
    y = "Proportion of Behaviors",
    fill = "Sex"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(
    values = c(
      "SD.Female" = "lightcoral",
      "SD.Male"   = "firebrick3",
      "FG.Female" = "springgreen3",
      "FG.Male"   = "springgreen4"
    ),
    labels = c(
      "SD.Female" = "Female",
      "SD.Male"   = "Male",
      "FG.Female" = "Female",
      "FG.Male"   = "Male"
    )
  ) +
  theme_classic()



#### PART 4: Map HAB and HC overlap ####

# Read in data
filtered_data <- read.csv("filtered_data.csv")

# # Subset for only dolphins seen foraging
# forage_data <- subset(
#   filtered_data,
#   grepl("Feed|pFeed", Behaviors)
# )

forage_data <- filtered_data

# Subset for HC behaviors only
forage_data <- subset(
  forage_data,
  DiffHI == "SD" | DiffHI == "FG"
)

# Get coordinates with HC only
HC_coordinates <- forage_data[c("Date", "Year", "Code", "StartLon", "StartLat", "DiffHI")]
# Add HAB count
HC_coordinates$HAB_time <- ifelse(HC_coordinates$Year > 2005, 2, 1)
# Now create a list for each year
HC_list <- split(HC_coordinates, HC_coordinates$HAB_time)
# Save data
write.csv(HC_list[[1]], "HC_map_PreHAB.csv")
write.csv(HC_list[[2]], "HC_map_PostHAB.csv")

# See which individuals overlap
presence_df <- do.call(
  rbind,
  lapply(
    split(forage_data$DiffHI, forage_data$Code),
    function(x) c(
      FG = "FG" %in% x,
      SD = "SD" %in% x
    )
  )
)

presence_df <- as.data.frame(presence_df)
presence_df$Code <- rownames(presence_df)
rownames(presence_df) <- NULL

sum(presence_df$SD == T & presence_df$FG == T) 
# 69 SD, 46 FG, 19 Both

# Add HAB count
forage_data$HAB_time <- ifelse(forage_data$Year > 2006, 2, 1)

# Create HC data only 
HC_data <- subset(forage_data, DiffHI != "None")

# Now create a list for each year
HC_list <- split(HC_data, HC_data$HAB_time)

# save data
write.csv(HC_list[[1]], "HC_data_preHAB.csv")
write.csv(HC_list[[2]], "HC_data_postHAB.csv")

# Now create a list for each year
forage_list <- split(forage_data, forage_data$HAB_time)

# Calculate home ranges of foraging dolphins
kernel_data <- list()

for (i in seq_along(forage_list)) {
  
  clean_data <- forage_list[[i]][
    !is.na(forage_list[[i]]$StartLon) &
      !is.na(forage_list[[i]]$StartLat) &
      forage_list[[i]]$StartLon >= -180 &
      forage_list[[i]]$StartLon <= 180 &
      forage_list[[i]]$StartLat >= -90 &
      forage_list[[i]]$StartLat <= 90,
  ]
  
  # keep dolphins with >=5 sightings
  clean_data <- clean_data[
    clean_data$Code %in% names(which(table(clean_data$Code) >= 5)),
  ]
  
  coords <- clean_data[, c("StartLon", "StartLat")]
  ids <- clean_data$Code   
  coords_sp <- SpatialPointsDataFrame(
    coords = coords,
    data = data.frame(id = ids)
  )
  
  proj4string(coords_sp) <- CRS("+proj=longlat +datum=WGS84")
  
  dolph.sp <- spTransform(
    coords_sp,
    CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
  )
  
  # kernel per individual
  kernel <- kernelUD(
    dolph.sp[, "id"],
    h = "href",
    grid = 100,
    extent = 5
  )
  
  kernel_data[[i]] <- kernel
}

# Save eco dat
saveRDS(kernel_data, "kernel_data.RData")
kernel_data <- readRDS("kernel_data.RData")

# Plot map
poly_list <- list()

for (i in seq_along(kernel_data)) {
  
  # get 40% polygons
  ver95 <- try(getverticeshr(kernel_data[[i]], percent = 40))
  
  # skip if failed
  if (inherits(ver95, "try-error")) next
  
  # convert to sf
  ver95_sf <- st_as_sf(ver95)
  
  # transform to lat/long
  ver95_sf <- st_transform(ver95_sf, crs = 4326)
  
  ver95_sf$dataset <- paste("Dataset", i)
  
  poly_list[[i]] <- ver95_sf
}

# combine all
poly_data <- do.call(rbind, poly_list)

# Create a for loop to store each period's average coordinates
# Extract 50% home range polygons
homerange50 <- lapply(kernel_data, function (kud) getverticeshr(kud, percent = 50))

centroid_list <- list()

for (i in seq_along(homerange50)) {
  
  # Convert entire dataset to sf
  hr_sf <- st_as_sf(homerange50[[i]])
  
  # Ensure it's lat/long
  hr_sf <- st_transform(hr_sf, crs = 4326)
  
  # Compute centroids for ALL polygons at once
  centroids <- st_centroid(hr_sf)
  
  # Extract coordinates
  coords <- st_coordinates(centroids)
  
  # Build dataframe
  centroids_df <- data.frame(
    ID = hr_sf$id,
    Longitude = coords[,1],
    Latitude = coords[,2]
  )
  
  centroid_list[[i]] <- centroids_df
}

# Connect ids to human-centric behavior
# Get unique IDs
ids <- unique(forage_data$Code)

# Build classification per individual
HC_data <- aggregate(DiffHI ~ Code, forage_data, function(x) {
  
  # remove "None"
  x <- x[x != "None"]
  
  # now classify
  if (length(x) == 0) return(NA)   # no behavior at all
  
  if (all(x == "SD")) return("SD")
  if (all(x == "FG")) return("FG")
  if (any(x == "SD") & any(x == "FG")) return("Both")
})

colnames(HC_data) <- c("ids", "HC")

poly_data <- merge(
  poly_data,
  HC_data,
  by.x = "id",
  by.y = "ids",
  all.x = TRUE
)
poly_data$HC <- as.character(poly_data$HC)

for (i in seq_along(centroid_list)) {
  centroid_list[[i]] <- merge(
    centroid_list[[i]],
    HC_data,
    by.x = "ID",
    by.y = "ids",
    all.x = TRUE
  )
}

for (i in seq_along(centroid_list)) {
  centroid_list[[i]]$HC <- as.character(centroid_list[[i]]$HC)
  centroid_list[[i]]$HC[is.na(centroid_list[[i]]$HC)] <- "Unknown"
}

saveRDS(centroid_list, "centroid_list.RData")
centroid_list <- readRDS("centroid_list.RData")

# Seperate saved centroids
write.csv(centroid_list[[1]], "centroids_PreHAB.csv")
write.csv(centroid_list[[2]], "centroids_PostHAB.csv")

# Read in human data
boat_data <- read.csv("human_dolphin_data.csv")

# Subset with only boats present
boat_data <- subset(boat_data, X.Boats > 0 |
                      X.CrabPots > 0 |
                      X.Lines > 0)

# Fix data INF
boat_data$X.Boats <- ifelse(boat_data$X.Boats < 0, 0, boat_data$X.Boats)
boat_data$X.CrabPots <- ifelse(boat_data$X.CrabPots < 0, 0, boat_data$X.CrabPots)
boat_data$X.Lines <- ifelse(boat_data$X.Lines < 0, 0, boat_data$X.Lines)

# Save data
write.csv(boat_data, "boat_data.csv")


