## --- Bayesian Multi-network Diffusion Analysis --- ##

#### DATA WRANGLING ####===

# load packages
## Clean up
if(!require(tidyr)){install.packages('tidyr'); library(tidyr)} 
if(!require(doParallel)){install.packages('doParallel'); library(doParallel)} # Faster computing
if(!require(dplyr)){install.packages('dplyr'); library(dplyr)} 
## Mapping/Graphs
if(!require(ggmap)){install.packages('ggmap'); library(ggmap)} # register API key version = '3.0.0'
if(!require(ggplot2)){install.packages('ggplot2'); library(ggplot2)} 
if(!require(raster)){install.packages('raster'); library(raster)} 
if(!require(sf)){install.packages('sf'); library(sf)} # Convert degrees to meters
if(!require(sp)){install.packages('sp'); library(sp)} # Convert degrees to meters
if(!require(adehabitatHR)){install.packages('adehabitatHR'); library(adehabitatHR)} # Caluculate MCPs and Kernel density 
if(!require(maps)){install.packages('maps'); library(maps)} # To map florida
if(!require(ggOceanMaps)){install.packages('ggOceanMaps'); library(ggOceanMaps)} # To map florida
if(!require(ggspatial)){install.packages('ggspatial'); library(ggspatial)} # To map florida
if(!require(prettymapr)){install.packages('prettymapr'); library(prettymapr)} # To map florida

# Set relative path working directory
setwd("../Data") 
# set working directory

#### PART 1: Wrangling data ####

# Subset original data
data_1 <- read.csv("93_04_data.csv")
data_2 <- read.csv("05_14_data.csv")
orig_data <- merge(data_1, data_2, all = T)

# Clean data
orig_data$Confirmed_HI <- ifelse(orig_data$ConfHI != "0", 1, 0)
orig_data$Date <- as.Date(as.character(orig_data$Date), format="%d-%b-%y")

# Change names to match their earlier names
orig_data$Code <- ifelse(orig_data$Code == "1312", "F222", orig_data$Code)

# Add year
orig_data$Year <- as.numeric(format(orig_data$Date, format = "%Y"))

# Save transformed data
write.csv(orig_data, "orig_data.csv")

# Read orig_data
orig_data <- read.csv("orig_data.csv")

# Filter data to include individuals that were seen at least 10 times 
tab <- table(orig_data$Code)
codes_in_all <- rownames(tab > 10)
filtered_data <- orig_data[orig_data$Code %in% codes_in_all, ]

# If needed subset
#start_year <- 2004
#filtered_data <- subset(filtered_data, Year %in% start_year:2014)

# Add a month rank column
filtered_data$Date <- as.Date(as.character(filtered_data$Date), format="%Y-%m-%d")
filtered_data <- filtered_data[order(filtered_data$Date), ]
ym <- format(filtered_data$Date, "%Y-%m")
filtered_data$MonthIndex <- match(ym, unique(ym))

# Add different HC categories
# Separate HI Behaviors to create weighted HI prop variable
#' BG = Beg: F, G, H
#' SD = Scavenge and Depredation: A, B, C, D, E
#' FG = Fixed Gear Interaction: P
# Change the code using ifelse statements
filtered_data$DiffHI <- ifelse(filtered_data$ConfHI %in% c("F", "G", "H", "A", "B", "C", "D", "E"), "SD",
                               ifelse(filtered_data$ConfHI %in% c("P"), "FG", "None"))

# Make a column for 6 months instead
split <- 6
filtered_data$sixmon <- ceiling(as.numeric(filtered_data$MonthIndex) / split)

write.csv(filtered_data, "filtered_data.csv")


#### PART 2: Inspect data on human-centric behavior ####

# Read in data
filtered_data <- read.csv("filtered_data.csv")

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


# If it has FG or SD then don't assign None
cleaned <- filtered_data %>%
  group_by(Code, sixmon) %>%
  mutate(has_real_cat = any(DiffHI %in% c("FG", "SD"))) %>%
  ungroup() %>%
  filter(!(has_real_cat & DiffHI == "None")) %>%  # drop None only where a real category exists
  select(-has_real_cat)

# Categorize DiffHI to IDs
HI_diff <- as.data.frame(table(cleaned$Code, cleaned$DiffHI, cleaned$sixmon))
colnames(HI_diff) <- c("Code", "DiffHI", "sixmon", "Freq")

# Add dates
seq_dates <- seq(from = as.Date("2004-01-01"),
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
  geom_line() +
  xlab("Date") +
  ylab("No. of Individuals") +
  labs(color = "HC Category") +
  theme_minimal()


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

# Count each dolphin behavior
count_data <- as.data.frame(
  table(
    age_class = filtered_data$age_class,
    behavior  = filtered_data$DiffHI,
    sex = filtered_data$Sex
  )
)

# Remove zero-count combinations
count_data <- subset(count_data, Freq > 0)

# Subset without None
count_data <- subset(count_data, behavior != "None")

# Make proportions
count_data$prop <- ave(
  count_data$Freq,
  count_data$behavior,
  FUN = function(x) x / sum(x)
)

# Plot the data
ggplot(count_data,
       aes(x = age_class, y = prop, fill = behavior)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8)) +
  facet_wrap(~ sex) +
  labs(
    x = "Age Class",
    y = "Proportion of Behaviors",
    fill = "Behavior"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("SD" = "steelblue", "FG" = "darkorange")) +
  theme_classic()



#### PART 4: Check variation in associations ####

# Read in network
nxn <- readRDS("nxn_sd.RData")
gbi <- readRDS("gbi_sd.RData")

# Done in the HPC 

#  Create 1000 random group-by-individual binary matrices
reps<- 1000
n.cores <- detectCores()
registerDoParallel(n.cores)

source("functions.R")
nF <- lapply(gbi, function (group_index) null(group_index, iter=reps))
saveRDS(nF, "nF.RData")

#' Calculate the association and CV for each of the 1000 permuted matrices to
#' create null distribution
nF <- readRDS("nF.RData")

cv_null <- foreach(i = 1:reps, .combine = c) %dopar% {
  SRI_cv(nF[[i]])
}

stopImplicitCluster()

saveRDS(cv_null, "cv_null_sd.RData")

# Next take results from the HPC

# Read in null cv values for one year
cv_null <- readRDS("../data/cv_years.RData")
## Remove NAs, if any
# cv_null = cv_null[!is.na(cv_null)]

# Calculate the CV of the observation association data
# CV = (SD/mean)*100
cv_obs <- lapply(nxn, function (df) {(sd(df) / mean(df)) * 100})  # Very high CV = unexpectedly 
# high or low association indices in the empirical distribution

# Calculate 95% confidence interval, in a two-tailed test
cv_ci = lapply(cv_null, function (df) {quantile(df, probs=c(0.025, 0.975), type=2)})

# Check whether pattern of connections is non-random
par(mfrow=c(2, 1))

# Create a list to store the histograms
hist_cvs <- list()

# Create histograms for each element in cv_null
for (i in seq_along(cv_null)) {
  hist_cvs[[i]] <- hist(cv_null[[i]], 
                        breaks=50,
                        xlim = c(min(cv_null[[i]]), max(cv_obs[[i]] + 10)),
                        col='grey70',
                        main = NULL,
                        xlab="Null CV SRI")
  
  # Add lines for empirical CV, 2.5% CI, and 97.5% CI
  abline(v= cv_obs[[i]], col="red")
  abline(v= cv_ci[[i]], col="blue")
  abline(v= cv_ci[[i]], col="blue")
}

#' This shows whether there are more preferred/avoided 
#' relatioNPhips than we would expect at random



#### PART 5: Map HAB and HC overlap ####

# Read in data
filtered_data <- read.csv("filtered_data.csv")

# Subset for only dolphins seen foraging
forage_data <- subset(
  filtered_data,
  grepl("Feed|pFeed", Behaviors)
)

# Subset for HC behaviors only
forage_data <- subset(
  forage_data,
  Code %in% unique(Code[DiffHI %in% c("SD", "FG")])
)

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

sum(presence_df$SD == T & presence_df$FG == F) 
# 70 SD, 46 FG, 19 Both

# Add HAB count
forage_data$HAB_time <- ifelse(forage_data$Year > 2006, 2, 1)

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

saveRDS(centroid_list, "centroid_list.RData")
centroid_list <- readRDS("centroid_list.RData")

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

# Upload florida map
register_google(key = "AIzaSyCNUfReSv2TSoMrxLnDC0glT8kffvSpLGM")

bbox <- st_bbox(poly_data)

florida_map <- basemap(
  limits = c(
    bbox["xmin"] - 0.2,
    bbox["xmax"] + 0.2,
    bbox["ymin"] - 0.2,
    bbox["ymax"] + 0.2
  ),
  bathymetry = F 
) +
  theme_void()


# Plot data
for (i in unique(poly_data$dataset)) {
  
  idx <- as.numeric(gsub("Dataset ", "", i))
  centroids_df <- centroid_list[[idx]]
  
  p <- florida_map +
    
    # ONLY centroids
    geom_point(
      data = centroids_df,
      aes(x = Longitude, y = Latitude),
      color = "black",
      size = 3.5
    ) +
    geom_point(
      data = centroids_df,
      aes(x = Longitude, y = Latitude, color = HC),
      size = 2.5
    )
  
  print(p)
}

