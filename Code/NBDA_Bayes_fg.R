## --- Bayesian Multi-network Diffusion Analysis --- ##

#### Fixed Gear Foraging ####===

# load packages
#if(!require(devtools)){install.packages('devtools'); library(devtools)} # To load NBDA
#if(!require(remotes)){install.packages('remotes'); library(remotes)} 
#remotes::install_github("stan-dev/cmdstanr") # If STBayes doesn't download
# Install NBDA package
#devtools::install_github("whoppitt/NBDA")
# install devtools if not already
#devtools::install_github("michaelchimento/STbayes")
#cmdstanr::set_cmdstan_path(cmdstanr::install_cmdstan())
## Bayesian
if(!require(tidyr)){install.packages('tidyr'); library(tidyr)} 
if(!require(abind)){install.packages('abind'); library(abind)} # array
if(!require(STbayes)){install.packages('STbayes'); library(STbayes)} 
if(!require(ggplot2)){install.packages('ggplot2'); library(ggplot2)} 
if(!require(dplyr)){install.packages('dplyr'); library(dplyr)} 
if(!require(posterior)){install.packages('posterior'); library(posterior)} 
## Creating networks
if(!require(asnipe)){install.packages('asnipe'); library(asnipe)} # get_group_by_individual
if(!require(sf)){install.packages('sf'); library(sf)} # Convert degrees to meters
if(!require(sp)){install.packages('sp'); library(sp)} # Convert degrees to meters
if(!require(adehabitatHR)){install.packages('adehabitatHR'); library(adehabitatHR)} # Caluculate MCPs and Kernel density 
if(!require(kinship2)){install.packages('kinship2'); library(kinship2)} # genetic relatedness
if(!require(doParallel)){install.packages('doParallel'); library(doParallel)} # for faster computing

# Set relative path working directory
setwd("../Data") 
# set working directory

#### PART 1: Create Networks ####
# Horizontal network -----------------------------------------

# Read in filtered data
filtered_data <- read.csv("filtered_data.csv")

# Add individual data
# Read ILVs
ILV_all <- read.csv("Individuals_Residency.csv", header=TRUE, sep=",")
ILV_all <- ILV_all[, c("Alias", "HI_Indiv", "Mom", "Sex", "BirthYear")]

# Group each individual by date and sighting
group_data <- filtered_data[,c("Date","Sighting","Code","Year","sixmon","DiffHI")]
group_data$Group <- cumsum(!duplicated(group_data[1:2])) # Create sequential group # by date
group_data <- group_data[,c(1, 3:7)] # Subset ID and group #

# Find the first time globally where Confirmed_HI == HC
HC <- "FG"
first_year <- min(filtered_data$Year[filtered_data$DiffHI == HC])

# Get individuals who had Confirmed_HI == HC in that first year
first_year_individuals <- unique(filtered_data$Code[
  filtered_data$DiffHI == HC & filtered_data$Year == first_year
])

# Assign yes/no based on that
ILV_all$Demons_HI_forage <- ifelse(ILV_all$Alias %in% first_year_individuals, "yes", "no")
write.csv(ILV_all, "ILV_FG_subset.csv")

# Subset the data to include observations only from before acquisition
# 1. Identify Codes to exclude
exclude_codes <- ILV_all$Alias[ILV_all$Demons_HI_forage == 'yes']

# 2. Compute first HI index for each Code
first_HI_index <- tapply(seq_len(nrow(group_data)), group_data$Code, function(idx) {
  hi_rows <- idx[group_data$Confirmed_HI[idx] == 1]
  if (length(hi_rows) > 0) hi_rows[1] else Inf
})

# 3. Split rows by Code
split_rows <- split(seq_len(nrow(group_data)), group_data$Code)

# 4. For excluded Codes: keep all rows
# For others: keep up to first HI, or all if no HI (cutoff = Inf)
keep_rows <- mapply(function(idx, cutoff, code) {
  if (code %in% exclude_codes || is.infinite(cutoff)) {
    idx  # keep everything for excluded Codes or codes with no HI
  } else {
    idx[idx <= cutoff]  # keep up to first HI for others
  }
}, split_rows, first_HI_index[names(split_rows)], names(split_rows))

# 5. Flatten and subset
keep_rows <- unlist(keep_rows)

# Subset the data
group_data <- group_data[keep_rows, ]

# Now create a list for each time period
group_data_list <- split(group_data, group_data$Year)

# Calculate Gambit of the group
create_gbi <- function(list_years) {
  gbi <- list()
  for (i in seq_along(list_years)) {
    
    # Gambit of the group index
    gbi[[i]] <- get_group_by_individual(list_years[[i]][,c("Code", "Group")], data_format = "individuals")
  }
  return(gbi)                                      
}

gbi_sd <- create_gbi(group_data_list)
saveRDS(gbi_fg, "gbi_fg.RData")

# Create association matrix
source("../Code/functions.R") # SRI.func
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

nxn_fg <- create_nxn(gbi_fg)

# Save nxn
saveRDS(nxn_fg, "nxn_fg.RData")

# Ecological network -----------------------------------------
# Transform coordinate data into a Spatial Points Dataframe in km

# Read in filtered data
filtered_data <- read.csv("filtered_data.csv")

# Read in vertical network
SRI_vert_all <- readRDS("SRI_vert_all.RData")

# Now create a list for each year
filtered_list <- split(filtered_data, filtered_data$Year)

# Create a list for Home Range
ecol_all <- list()
for (i in seq_along(filtered_list)) {
  
  # Remove NAs and invalid coordinates
  clean_data <- filtered_list[[i]][
    !is.na(filtered_list[[i]]$StartLon) &
      !is.na(filtered_list[[i]]$StartLat) &
      filtered_list[[i]]$StartLon >= -180 & filtered_list[[i]]$StartLon <= 180 &
      filtered_list[[i]]$StartLat >= -90 & filtered_list[[i]]$StartLat <= 90,
  ]
  
  # Keep only individuals with at least 5 observations
  clean_data <- clean_data[clean_data$Code %in% names(which(table(clean_data$Code) >= 5)), ]
  
  ids <- clean_data$Code
  coordinates <- clean_data[, c("StartLon", "StartLat")]
  
  # Create a SpatialPointsDataFrame with coordinates
  coords_sp <- SpatialPointsDataFrame(coords = coordinates, data = data.frame(id = ids))
  
  # Set CRS to WGS84
  proj4string(coords_sp) <- CRS("+proj=longlat +datum=WGS84")
  
  # Transform to a UTM CRS that uses km as the unit
  dolph.sp <- spTransform(coords_sp, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
  
  # Use the calculated extent in kernelUD
  kernel <- kernelUD(dolph.sp, h = 1000)
  
  # Calculate Dyadic HRO Matrix: HRO = (Rij/Ri) * (Rij/Rj)
  kov <- kerneloverlaphr(kernel, method = "HR", lev = 95)
  
  # Order data
  #order_rows <- rownames(SRI_vert_all)
  #order_cols <- colnames(SRI_vert_all)
  
  # Apply the order to each matrix in the list
  #ecol_all[[i]] <- kov[order_rows, order_cols] 
  ecol_all[[i]] <- kov
}

# Save eco dat
saveRDS(ecol_all, "ecol_all_fg.RData")

# Relatedness network -----------------------------------------

# Read in data
ILV_pat <- read.csv("Paternity_data.csv") 
ILV_all <- read.csv("Individuals_Residency.csv", header=TRUE, sep=",")

# Merge to fill in empty data
ILV <- merge(
  ILV_pat,
  ILV_all,
  all.y = TRUE
)

# Subset paternity data
pedigree_df <- data.frame(Alias = ILV$Alias,
                          Mom = ILV$Mom,
                          Dad = ILV$Dad,
                          Sex = ILV$Sex)

# Fix dad data
pedigree_df$Dad <- ifelse(pedigree_df$Dad == "na", NA, pedigree_df$Dad)
pedigree_df$Dad <- ifelse(pedigree_df$Dad == "FB26 or FB66", "FB26", pedigree_df$Dad)
pedigree_df$Dad <- ifelse(pedigree_df$Dad == "FB76 or FB38", "FB76", pedigree_df$Dad)

# Fix sex so that probable is assigned
pedigree_df$Sex <- ifelse(ILV$Sex == "Probable Female", "Female",
                          ifelse(ILV$Sex == "Probable Male", "Male", ILV$Sex))

# Make sex numeric
pedigree_df$Sex <- ifelse(pedigree_df$Sex == "Female", 2, 
                          ifelse(pedigree_df$Sex == "Male", 1, NA))

# Limit data to non-missing paternity IDs
pedigree_subset <- pedigree_df[!is.na(pedigree_df$Mom) | !is.na(pedigree_df$Dad), ]

# Reset row names to be sequential
row.names(pedigree_subset) <- NULL

# Make id numeric
## Moms
pedigree_df$ID <- rownames(pedigree_df)
for (i in 1:nrow(pedigree_df)) {
  pedigree_df$Mom <- ifelse(pedigree_df$Mom %in% pedigree_df$Alias[i], 
                            pedigree_df$ID[i], pedigree_df$Mom)
}

## Dads
for (i in 1:nrow(pedigree_df)) {
  pedigree_df$Dad <- ifelse(pedigree_df$Dad %in% pedigree_df$Alias[i], 
                            pedigree_df$ID[i], pedigree_df$Dad)
}

# Only take the ids that aren't found in the list
missing_moms<- subset(pedigree_df, nchar(Mom) > 3)
missing_dads<- subset(pedigree_df, nchar(Dad) > 3)

# Create the sequence of numbers starting from 118
number_mom <- data.frame(Mom = unique(missing_moms$Mom), 
                         ID = c((nrow(pedigree_df) + 1):(nrow(pedigree_df) + length(unique(missing_moms$Mom)))))

# Fill in numbers
for (i in 1:nrow(missing_moms)) {
  missing_moms$Mom <- ifelse(missing_moms$Mom %in% number_mom$Mom[i], 
                             number_mom$ID[i],
                             missing_moms$Mom)
}

# Make ID numeric
missing_moms$Mom <- as.numeric(missing_moms$Mom)

# Do the same thing with dads
number_dad <- data.frame(Dad = unique(missing_dads$Dad), 
                         ID = c((max(missing_moms$Mom) + 1):(max(missing_moms$Mom) + length(unique(missing_dads$Dad)))))
for (i in 1:nrow(missing_dads)) {
  missing_dads$Dad <- ifelse(missing_dads$Dad %in% number_dad$Dad[i], 
                             number_dad$ID[i],
                             missing_dads$Dad)
}

# Make ID numeric
missing_dads$Dad <- as.numeric(missing_dads$Dad)

# Fill in the rest of the NAs with random numbers
## Moms
missing_moms_match <- subset(pedigree_df, nchar(Mom) > 3)
matching_indices <- match(pedigree_df$Mom, missing_moms_match$Mom)
pedigree_df$Mom <- ifelse(!is.na(matching_indices), missing_moms$Mom[matching_indices], pedigree_df$Mom)

## Dads
missing_dads_match<- subset(pedigree_df, nchar(Dad) > 3)
matching_indices <- match(pedigree_df$Dad, missing_dads_match$Dad)
pedigree_df$Dad <- ifelse(!is.na(matching_indices), missing_dads$Dad[matching_indices], pedigree_df$Dad)

# Now create data for function
pedigree_data <- data.frame(id = as.numeric(pedigree_df$ID),
                            mom = as.numeric(pedigree_df$Mom),
                            dad = as.numeric(pedigree_df$Dad),
                            sex = pedigree_df$Sex)
# Assuming your dataframe is named pedigree_data
pedigree_data$dad[is.na(pedigree_data$dad)] <- 0  # Replace NA with 0 or another appropriate code
pedigree_data$mom[is.na(pedigree_data$mom)] <- 0  # Replace NA with 0 or another appropriate code

# Add Fake Fathers
for (i in which(pedigree_data$mom > 0 & pedigree_data$dad == 0)) {
  pedigree_data$dad[i] <- i + max(pedigree_data$dad)
}

# Create fake individuals
fake_ids <- (nrow(pedigree_df) + 1):(max(pedigree_data$dad) + 1)
fake <- data.frame(id = fake_ids,
                   mom = rep(0, length(fake_ids)),
                   dad = rep(0, length(fake_ids)),
                   sex = rep(3, length(fake_ids)))
pedigree_data <- rbind(pedigree_data, fake)

# Change errors
pedigree_data$sex[pedigree_data$id %in% c(139:270)] <- 1
pedigree_data$sex[pedigree_data$id %in% c(118:138)] <- 2

# For limited data
pedigree_data$sex[pedigree_data$id %in% c(94:112, 117:nrow(pedigree_data))] <- 1
pedigree_data$sex[pedigree_data$id %in% c(58:93)] <- 2

# Create GR matrix
ped <- pedigree(id = pedigree_data$id, 
                dadid = pedigree_data$dad, 
                momid = pedigree_data$mom,
                sex = pedigree_data$sex)

# Calculate kinship matrix
kinship_matrix <- kinship(ped)
relate_all <- kinship_matrix[1:117, 1:117]
saveRDS(kinship_matrix, "kinship_matrix.RData")



#### PART 2: Create acquisition data for model input ####

# Read in all network data
nxn <- readRDS("nxn_fg.RData")
SRI_vert_all <- readRDS("SRI_vert_all.RData")
ecol_all <- as.array(readRDS("ecol_all_fg.RData"))

# Get rid of IDs without data in nxn from vert data
for (i in seq_along(nxn)) {
  target_names <- rownames(SRI_vert_all)
  
  rn <- rownames(nxn[[i]])
  cn <- colnames(nxn[[i]])
  
  keep_r <- target_names[target_names %in% rn]
  keep_c <- target_names[target_names %in% cn]
  
  nxn[[i]] <- nxn[[i]][keep_r, keep_c, drop = FALSE]
}

# Get rid of IDs without data in ecol data from nxn
for (i in seq_along(nxn)) {
  target_names <- rownames(ecol_all[[i]])
  
  rn <- rownames(nxn[[i]])
  cn <- colnames(nxn[[i]])
  
  keep_r <- target_names[target_names %in% rn]
  keep_c <- target_names[target_names %in% cn]
  
  nxn[[i]] <- nxn[[i]][keep_r, keep_c, drop = FALSE]
}

# Get rid of IDs without data in nxn from ecol data
for (i in seq_along(ecol_all)) {
  target_names <- rownames(nxn[[i]])
  
  rn <- rownames(ecol_all[[i]])
  cn <- colnames(ecol_all[[i]])
  
  keep_r <- target_names[target_names %in% rn]
  keep_c <- target_names[target_names %in% cn]
  
  ecol_all[[i]] <- ecol_all[[i]][keep_r, keep_c, drop = FALSE]
}

# Add zeros to rows that don't have all individuals

# Get all unique IDs across all matrices
total_ids <- unique(unlist(lapply(nxn, rownames)))

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


# Do the same for the vert matrix
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
years <- 11
SRI_vert_all <- replicate(years, full_mat, simplify = FALSE)
SRI_vert_all <- as.array(SRI_vert_all)

# Also do this for the ecol matrix
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

# Put matrices into data frame
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

# Read in full data
ILV_all <- read.csv("ILV_FG_subset.csv")

# Subset event data to include only edge_list ids
ILV_all <- subset(ILV_all, Alias %in% unique(edge_list$focal))

# Add order of acquisition data
# Create demonstrator column
filtered_data <- read.csv("filtered_data.csv")

# Create acquisition data
# Step 1: Filter filtered_data for confirmed HI behavior after first date
HC <- "FG"
first_year <- min(filtered_data$Year[filtered_data$DiffHI == HC])

hi_data <- filtered_data[filtered_data$DiffHI == HC & filtered_data$Year != first_year, ]

# Step 2: Get the first year each Alias showed the behavior
first_hi_yr <- aggregate(Year ~ Code, data = hi_data, FUN = min)

# Step 3: Create a new column for order of acquisition
first_hi_yr$HI_Order_acquisition <- as.numeric(first_hi_yr$Year)

# Step 4: Merge this info back into ILV_all
ILV_all <- merge(ILV_all, first_hi_yr[, c("Code", "HI_Order_acquisition")],
                 by.x = "Alias", by.y = "Code", all.x = TRUE)

# Step 5: Replace NA with 0 for individuals who had the behavior in 1995
ILV_all$HI_Order_acquisition[ILV_all$Demons_HI_forage == "yes"] <- 0

# Change individuals who didn't require behavior to t_end +1
ILV_all$time <- ifelse(is.na(ILV_all$HI_Order_acquisition), 
                       max(na.omit(ILV_all$HI_Order_acquisition)) + 1, 
                       ILV_all$HI_Order_acquisition)

# Get unique years excluding 0
years <- sort(unique(ILV_all$time[ILV_all$time != 0]))

# Create a mapping: year → index
yr_map <- setNames(seq_along(years), years)

# Replace years with mapped values, keep 0 as 0
ILV_all$time <- ifelse(ILV_all$time == 0, 0, yr_map[as.character(ILV_all$time)])

# Add end time
ILV_all$t_end <- max(na.omit(ILV_all$time)) - 1

# Edit the other needed columns
ILV_all$id <- ILV_all$Alias
ILV_all$trial <- 1

# Get rid of other columns
event_data <- ILV_all[, c("trial", "id", "time", "t_end")]
write.csv(event_data, "event_data_fg.csv")

# Get rid of unused individuals
edge_list <- subset(edge_list, focal %in% event_data$id)
edge_list <- subset(edge_list, other %in% event_data$id)

saveRDS(edge_list, "edge_list_fg.RData")

# individual level variables
ILV_all$Sex <- ifelse(ILV_all$Sex == "Female", 1, 0)
ILV_all$Sex <- ifelse(is.na(ILV_all$Sex), 1, 0) # Fix NAs

ILV_all$BirthYear <- as.numeric(ILV_all$BirthYear)
ILV_all$BirthYear <- ifelse(is.na(ILV_all$BirthYear), 1985, ILV_all$BirthYear)

# Add HAB
ILV_all$HAB <- ifelse(ILV_all$time > 3, 0, 1)

# Constant ILVs
ILV_c <- data.frame(id = ILV_all$Alias,
                    sex = ILV_all$Sex)

# Time varying ILVs
res <- 11
ILV_tv <- data.frame(
  trial = 1,
  id = rep(ILV_all$Alias, each = res),
  time = rep(1:res, times = length(ILV_all$Alias)),  
  age = (rep(2004:2014, times = length(ILV_all$Alias)) - 
           rep(ILV_all$BirthYear, each = res)) - 1,
  HAB = rep(ILV_all$HAB, each = res)
)

# Shorten it by 1
#ILV_tv <- subset(ILV_tv, time !=res)

# Change age to age groups
ILV_tv$age_group <- ifelse(ILV_tv$age >= 10, "adult", 
                           ifelse(ILV_tv$age > 4, "juvenile",
                                  ifelse(ILV_tv$age > 0, "calf", "unborn")))
ILV_tv <- ILV_tv[, -4]

# Categorize DiffHI to IDs
## Dynamic raw count data
# rawHI_diff <- as.data.frame(table(filtered_data$Code, 
#                                   filtered_data$DiffHI, 
#                                   filtered_data$sixmon))
# colnames(rawHI_diff) <- c("Code", "DiffHI", "Time", "Freq")

## Static proportion data
# Catergorize into raw data
rawHI_diff <- as.data.frame(table(filtered_data$Code,
                                  filtered_data$DiffHI))
colnames(rawHI_diff) <- c("Code", "DiffHI", "Freq")

# Categorize ID to Sightings
IDbehav <- as.data.frame(table(filtered_data$Code))
colnames(IDbehav) <- c("Code", "Sightings")
# Order data
order_rows <- rownames(nxn_full[[1]])

# Now reorder the dataframe
IDbehav <- IDbehav %>%
  arrange(match(Code, order_rows))

# Create a frequency count for each HI behavior
get_IDHI <- function(HI, IDbehav_data, rawHI_diff_data) {
  df <- IDbehav_data
  HI_freq <- rawHI_diff_data$Freq[rawHI_diff_data$DiffHI %in% HI]
  df$Behav <- HI_freq[match(df$Code, rawHI_diff_data$Code)]
  colnames(df) <- c("Code", "Sightings", "Behav")
  return(df)
}

IDbehav_HI <- get_IDHI("FG", IDbehav, rawHI_diff)
# Proportion of Sightings spent in HI
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

prob_HI <- Prop_HI(IDbehav_HI, filtered_data, "FG")

# Convert list of HI prop vectors into a matrix
ids <- unique(prob_HI$Code)
HIProp <- prob_HI$HIprop

# Create trans weights
HI_matrix <- data.frame(
  trial = 1,
  id = ids,
  #time = rawHI_diff$Time,
  t_weight = HIProp)

# Get rid of unused individuals
HI_matrix <- HI_matrix[HI_matrix$id %in% event_data$id, ]

# Read in edge list
edge_list <- readRDS("edge_list_fg.RData")

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

# Rearrange and get rid of the last time
edge_list <- edge_list[, c("focal", "other", "trial", 
                           "assoc", "vert", "ecol", "time")]
edge_list <- subset(edge_list, time != res)

# Input data
data_list <- import_user_STb(
  event_data = event_data,
  networks = edge_list,
  ILV_c = ILV_c,
  ILV_tv = ILV_tv,
  ILVi = c("age_group", "sex", "HAB"),
  ILVs = c("age_group", "sex", "HAB"),
  t_weights = HI_matrix
)

saveRDS(data_list, "data_list_fg.RData")



#### PART 3: Run the model and summary outputs ####

# Input data_list
data_list <- readRDS("data_list_fg.RData")

# Diagnoses notes
#' Detected N_veff = 0, likihood provides no information
#' Parameters that must be positive / bounded
#' Weak or default priors
#' Hierarchical or hazard-based structure

# Generate the model
model_full <- generate_STb_model(data_list, 
                                 data_type = "continuous_time",
                                 gq = T, 
                                 est_acqTime = T)

# Test
# Edit the priors in the model.STAN file
# writeLines(model_full, "model_strong_priors_fg.stan")

# fit_test <- fit_STb(
#   data_list,
#   "model_strong_priors.stan",
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

STb_save(fit, output_dir = "cmdstan_saves", name="my_first_fit")
fit <- readRDS('cmdstan_saves/my_first_fit.rds')

# View output
STb_summary(fit, digits = 3)

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

# Posterior Predictive Checks
event_data <- read.csv("event_data.csv")

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

#' The model reproduces the overall learning dynamics at the population level.

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

#' This is exactly what you expect in a hazard‑based model:
#' Early learners are well constrained
#' Later learners accumulate uncertainty across time
