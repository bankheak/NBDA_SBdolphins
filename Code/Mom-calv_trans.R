## --- Mother-Calve Relationships --- ##

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

# Create a frequency denominator for foraging
filtered_data$Foraging <- "Other"
filtered_data$Foraging[grepl(pattern = 'Feed', x = filtered_data$Behaviors, ignore.case = FALSE)] <- "Feed"

# Fix birthyear
ILV_all$BirthYear <- as.numeric(ILV_all$BirthYear)
ILV_all$BirthYear <- ifelse(is.na(ILV_all$BirthYear), 1985, ILV_all$BirthYear)

# Subset data for those in filtered data
Codes <- unique(filtered_data$Code)
ILV_pat <- subset(ILV_all, Alias %in% Codes)

# Count calves per Mom per BirthYear from ILV_pat
calves_by_mom_year <- ILV_pat[c("Mom", "Alias", "BirthYear")]
calves_by_mom_year$Depart <- ILV_pat$BirthYear + 6

write.csv(calves_by_mom_year, "calves_by_mom_year.csv")

# Records from the period in which mother and calf are associated

## Build a lookup: calf Alias -> Depart year
depart_lookup <- setNames(calves_by_mom_year$Depart, calves_by_mom_year$Alias)

## Keep observations strictly before departure; drop NAs (no known depart)
mom_calf <- subset(filtered_data, Code %in% calves_by_mom_year$Alias)

## Match each row's Code to its depart year
row_depart <- depart_lookup[mom_calf$Code]  # returns NA if Code not in lookup

## Keep observations strictly before departure; drop NAs (no known depart)
mom_calf_data <- subset(mom_calf, !is.na(row_depart) & Year < row_depart)

# Records of the mother after independence of the calf

## Keep observations strictly before departure; drop NAs (no known depart)
mom <- subset(filtered_data, Code %in% calves_by_mom_year$Mom)

## Match each row's Code to its depart year
row_depart_2 <- depart_lookup[mom$Code]  # returns NA if Code not in lookup

## Keep observations strictly before departure; drop NAs (no known depart)
mom_data <- subset(mom, !is.na(row_depart_2) & Year < row_depart_2)

# Records of the calf after independence

## Keep observations strictly before departure; drop NAs (no known depart)
calves_data <- subset(mom_calf, !is.na(row_depart) & Year > row_depart)

# Get rid of data with observations less than 5

## mom-calf
counts <- table(mom_calf_data$Code)
keep <- counts[mom_calf_data$Code] >= 5
mom_calf_data <- mom_calf_data[keep, ]

## mom
counts <- table(mom_data$Code)
keep <- counts[mom_data$Code] >= 5
mom_data <- mom_data[keep, ]

## calf
counts <- table(calves_data$Code)
keep <- counts[calves_data$Code] >= 5
calves_data <- calves_data[keep, ]

# Save the data sets
write.csv(mom_calf_data, "mom_calf_data.csv")
write.csv(mom_data, "mom_data.csv")
write.csv(calves_data, "calves_data.csv")


#### PART 2: Home Range ####

# Load packages
if(!require(sf)){install.packages('sf'); library(sf)} # Convert degrees to meters
if(!require(sp)){install.packages('sp'); library(sp)} # Convert degrees to meters
if(!require(adehabitatHR)){install.packages('adehabitatHR'); library(adehabitatHR)} # KernelUD

# Read in data:
mom_calf_data <- read.csv("mom_calf_data.csv")
mom_data <- read.csv("mom_data.csv")
calves_data <- read.csv("calves_data.csv")

# Create home range function
create_hr_data <- function(df) {
    
    # Get rid of na data
    df <- df[!is.na(df$StartLat) & !is.na(df$StartLon), ]
    df <- df[df$StartLat != 999 & df$StartLon != 999, ]
  
    # Extract IDs and coordinates
    ids <- df$Code
    coordinates <- df[, c("StartLon", "StartLat")]
    
    # Create a SpatialPointsDataFrame with coordinates
    coords_sp <- SpatialPointsDataFrame(coords = coordinates, data = data.frame(id = ids))
    
    # Set CRS to WGS84
    proj4string(coords_sp) <- CRS("+proj=longlat +datum=WGS84")
    
    # Transform to a UTM CRS that uses km as the unit
    coords_sp_utm <- spTransform(coords_sp, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
  
    # Use the calculated extent in kernelUD
    kernel <- kernelUD(coords_sp_utm, h = 1000)
    
  return(kernel)
}

# Calculate home range for each individual before independence
hr_mom_calf <- create_hr_data(mom_calf_data)

# Calculate home range for each mom after independence
hr_mom <- create_hr_data(mom_data)

# Calculate home range for each individual after independence
hr_calf <- create_hr_data(calves_data)

# Save HRO data
saveRDS(hr_mom_calf, "hr_mom_calf.RDS")
saveRDS(hr_mom, "hr_mom.RDS")
saveRDS(hr_calf, "hr_calf.RDS")

#### PART 3: Human-centric data ####

# Read in Data
mom_calf_data <- read.csv("mom_calf_data.csv")
mom_data <- read.csv("mom_data.csv")
calves_data <- read.csv("calves_data.csv")

# Calculate frequency of calf hc before departure
hc_behavs <- table(mom_calf_data$Code, mom_calf_data$DiffHI)
tot_forg <- table(mom_calf_data$Code, mom_calf_data$Foraging)

hc_mom_calf_freq <- data.frame(Code = row.names(hc_behavs),
                             Freq_bg = hc_behavs[, "BG"]/tot_forg[, "Feed"],
                             Freq_fg = hc_behavs[, "FG"]/tot_forg[, "Feed"],
                             Freq_sd = hc_behavs[, "SD"]/tot_forg[, "Feed"])

# Calculate frequency of mom hc after departure
hc_behavs <- table(mom_data$Code, mom_data$DiffHI)
tot_forg <- table(mom_data$Code, mom_data$Foraging)

hc_mom_freq <- data.frame(Code = row.names(hc_behavs),
                               Freq_fg = hc_behavs[,"FG"]/tot_forg[, "Feed"],
                               Freq_sd = hc_behavs[,"SD"]/tot_forg[, "Feed"])

# Calculate frequency of calf hc after departure
hc_behavs <- table(calves_data$Code, calves_data$DiffHI)
tot_forg <- table(calves_data$Code, calves_data$Foraging)

hc_calf_freq <- data.frame(Code = row.names(hc_behavs),
                               Freq_bg = hc_behavs[,"BG"]/tot_forg[, "Feed"],
                               Freq_fg = hc_behavs[,"FG"]/tot_forg[, "Feed"],
                               Freq_sd = hc_behavs[,"SD"]/tot_forg[, "Feed"])

# Save hc_data
write.csv(hc_mom_calf_freq, "hc_mom_calf_freq.csv")
write.csv(hc_mom_freq, "hc_mom_freq.csv")
write.csv(hc_calf_freq, "hc_calf_freq.csv")

#### PART 4: Run Correlation Tests for Home range ####

# Load packages
if(!require(sf)){install.packages('sf'); library(sf)} # Convert degrees to meters
if(!require(sp)){install.packages('sp'); library(sp)} # Convert degrees to meters
if(!require(adehabitatHR)){install.packages('adehabitatHR'); library(adehabitatHR)} # KernelUD

# Read in Data
## HR
hr_mom_calf <- readRDS("hr_mom_calf.RDS")
hr_mom <- readRDS("hr_mom.RDS")
hr_calf <- readRDS("hr_calf.RDS")

## Raw data
calves_by_mom_year <- read.csv("calves_by_mom_year.csv")
mom_calf_data <- read.csv("mom_calf_data.csv")
mom_data <- read.csv("mom_data.csv")
calves_data <- read.csv("calves_data.csv")

# Mom and calf overlaps
calves_by_mom <- na.omit(calves_by_mom_year[c("Alias", "Mom")]) 
names(calves_by_mom) <- c("calf", "mom")

# Ensure we only evaluate pairs that exist in each period
pairs_weaning <- subset(
  calves_by_mom,
  calf %in% names(hr_mom_calf) & mom %in% names(hr_mom_calf)
)
pairs_post <- subset(
  calves_by_mom,
  calf %in% names(hr_calf) & mom %in% names(hr_calf)
)

pairs_mom_after_vs_calf_after <- subset(
  calves_by_mom,
  calf %in% names(hr_calf) & mom %in% names(hr_mom)
)

# Helper to safely extract overlap values from a matrix
get_ovl <- function(M, a, b) {
  if (a %in% rownames(M) && b %in% colnames(M)) {
    M[a, b]
  } else if (b %in% rownames(M) && a %in% colnames(M)) {
    M[b, a]
  } else {
    NA_real_
  }
}

# Overlap matrix within hr_mom_calf (moms & calves during weaning)
ovl_weaning <- kerneloverlaphr(hr_mom_calf, method = "HR")  

pairs_weaning$overlap_weaning <- mapply(
  function(m, c) get_ovl(ovl_weaning, m, c),
  pairs_weaning$mom, pairs_weaning$calf
)

# Overlap matrix within hr_calf (moms & calves after weaning)
ovl_post <- kerneloverlaphr(hr_calf, method = "HR")

pairs_post$overlap_postweaning <- mapply(
  function(m, c) get_ovl(ovl_post, m, c),
  pairs_post$mom, pairs_post$calf
)


# Mom-after vs Calf-after (hr_mom vs hr_calf comparison)
ud_after_combined <- list()
class(ud_after_combined) <- class(hr_mom)  # retain estUDm class

moms_needed   <- unique(pairs_mom_after_vs_calf_after$mom)
calves_needed <- unique(pairs_mom_after_vs_calf_after$calf)

# Add moms (tag as "_momAfter")
for (m in moms_needed) {
  ud_after_combined[[paste0(m, "_momAfter")]] <- hr_mom[[m]]
}

# Add calves (tag as "_calfAfter")
for (c in calves_needed) {
  ud_after_combined[[paste0(c, "_calfAfter")]] <- hr_calf[[c]]
}

# Compute overlap in the combined set
ovl_mom_after_vs_calf_after <- kerneloverlaphr(ud_after_combined, method = "HR")

pairs_mom_after_vs_calf_after$overlap_momAfter_calfAfter <- mapply(
  function(m, c) {
    get_ovl(
      ovl_mom_after_vs_calf_after,
      paste0(m, "_momAfter"),
      paste0(c, "_calfAfter")
    )
  },
  pairs_mom_after_vs_calf_after$mom, pairs_mom_after_vs_calf_after$calf
)


# Combine results
results_overlap <- Reduce(function(x, y) {
  merge(x, y, by = c("calf", "mom"), all = TRUE)
}, list(
  pairs_weaning[, c("calf", "mom", "overlap_weaning")],
  pairs_post[,   c("calf", "mom", "overlap_postweaning")],
  pairs_mom_after_vs_calf_after[, c("calf", "mom", "overlap_momAfter_calfAfter")]
))

results_overlap <- results_overlap[order(results_overlap$mom, results_overlap$calf), ]


# Name overlaps
ids_mom <- names(hr_mom_calf)
ids_calf <- names(hr_calf)
shared_ids <- intersect(ids_mom, ids_calf)

# Tag names so pairs are unambiguous
names(hr_mom_calf) <- paste0(names(hr_mom_calf), "_mom")
names(hr_calf)     <- paste0(names(hr_calf), "_calf")

ud_all <- hr_mom_calf             # class(estUDm)
for (n in names(hr_calf)) {
  ud_all[[n]] <- hr_calf[[n]]     # add each UD; class remains 'estUDm'
}

# Compute overlap matrices
ovl_HR   <- kerneloverlaphr(ud_all, method = "HR")

# Extract overlap for the same individual across periods
shared_ids <- intersect(
  sub("_mom$", "", names(hr_mom_calf)),
  sub("_calf$", "", names(hr_calf))
)

pair_val <- function(M, id) M[paste0(id, "_mom"), paste0(id, "_calf")]

results_overlap <- data.frame(
  id   = shared_ids,
  HR   = sapply(shared_ids, pair_val, M = ovl_HR),
  row.names = NULL
)


#### PART 5: Run Correlation Tests for HC Behavior ####

# Read in Data
hc_mom_calf_freq <- read.csv("hc_mom_calf_freq.csv")
hc_mom_freq <- read.csv("hc_mom_freq.csv")
hc_calf_freq <- read.csv("hc_calf_freq.csv")
calves_by_mom_year <- read.csv("calves_by_mom_year.csv")

# Re-arrange moms and calves
calves_by_mom <- calves_by_mom_year[c("Alias", "Mom")]

# Before Independence

# Join calf freq
hc_mom_calf_freq$X <- NULL
dur_wean <- merge(calves_by_mom, hc_mom_calf_freq, by.x = "Alias", by.y = "Code", all.x = TRUE)
names(dur_wean)[3:5] <- c("freq_calf_bg", "freq_calf_fg", "freq_calf_sd")

# Join mom freq
dur_wean <- merge(dur_wean, hc_mom_calf_freq, by.x = "Mom", by.y = "Code", all.x = TRUE)
names(dur_wean)[6:8] <- c("freq_mom_bg", "freq_mom_fg", "freq_mom_sd")

# Get rid of NAs
dur_wean <- na.omit(dur_wean)

# Euclidean distance (scalar case)
dur_wean$euclid_dist_bg <- sqrt((dur_wean$freq_calf_bg - dur_wean$freq_mom_bg)^2)
dur_wean$euclid_dist_fg <- sqrt((dur_wean$freq_calf_fg - dur_wean$freq_mom_fg)^2)
dur_wean$euclid_dist_sd <- sqrt((dur_wean$freq_calf_sd - dur_wean$freq_mom_sd)^2)

# Transform to similarity
# Normalize Euclidian distances to make it range [0,1] and then simply apply 1-normalized distance
sim_eu <- function(value) {
  m <- max(value, na.rm = TRUE)
  if (is.na(m) || m == 0) return(rep(0, length(value)))
  1 - (value / m)
}

dur_wean$sim_calf_bg <- sim_eu(dur_wean$euclid_dist_bg)
dur_wean$sim_calf_fg <- sim_eu(dur_wean$euclid_dist_fg)
dur_wean$sim_calf_sd <- sim_eu(dur_wean$euclid_dist_sd)


# After Independence

# Merge the data sets for mom and calf after weaning
hc_mom_freq$X <- NULL
hc_calf_freq$X <- NULL
hc_after_wean_freq <- merge(hc_mom_freq, hc_calf_freq, all = T)

# Join calf freq
after_wean <- merge(calves_by_mom, hc_after_wean_freq, by.x = "Alias", by.y = "Code", all.x = TRUE)
names(after_wean)[3:5] <- c("freq_calf_bg", "freq_calf_fg", "freq_calf_sd")

# Join mom freq
after_wean <- merge(after_wean, hc_after_wean_freq, by.x = "Mom", by.y = "Code", all.x = TRUE)
names(after_wean)[6:8] <- c("freq_mom_bg", "freq_mom_fg", "freq_mom_sd")

# Get rid of NAs
after_wean <- na.omit(after_wean)

# Euclidean distance (scalar case)
after_wean$euclid_dist_bg <- sqrt((after_wean$freq_calf_bg - after_wean$freq_mom_bg)^2)
after_wean$euclid_dist_fg <- sqrt((after_wean$freq_calf_fg - after_wean$freq_mom_fg)^2)
after_wean$euclid_dist_sd <- sqrt((after_wean$freq_calf_sd - after_wean$freq_mom_sd)^2)

# Transform to similarity
# Normalize Euclidian distances to make it range [0,1] and then simply apply 1-normalized distance
sim_eu <- function(value) {
  m <- max(value, na.rm = TRUE)
  if (is.na(m) || m == 0) return(rep(0, length(value)))
  1 - (value / m)
}

after_wean$sim_calf_bg <- sim_eu(after_wean$euclid_dist_bg)
after_wean$sim_calf_fg <- sim_eu(after_wean$euclid_dist_fg)
after_wean$sim_calf_sd <- sim_eu(after_wean$euclid_dist_sd)

# Subset for one the overlap of moms and calves
dur_wean <- subset(dur_wean, Alias %in% after_wean$Alias)
after_wean <- subset(after_wean, Alias %in% dur_wean$Alias)

# Calculate the statistical difference
group1 <- after_wean$sim_calf_fg
group2 <- dur_wean$sim_calf_fg

wilcox.test(group1, group2, paired = TRUE)
