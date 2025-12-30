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

# Fix birthyear
ILV_all$BirthYear <- as.numeric(ILV_all$BirthYear)
ILV_all$BirthYear <- ifelse(is.na(ILV_all$BirthYear), 1985, ILV_all$BirthYear)

# Subset data for those in filtered data
Codes <- unique(filtered_data$Code)
ILV_pat <- subset(ILV_all, Alias %in% Codes)

# Count calves per Mom per BirthYear from ILV_pat
calves_by_mom_year <- ILV_pat[c("Mom", "Alias", "BirthYear")]
calves_by_mom_year$Depart <- ILV_pat$BirthYear + 6

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

hc_mom_calf_freq <- data.frame(Code = row.names(hc_behavs),
                             Freq_bg = hc_behavs[,"BG"]/rowSums(hc_behavs),
                             Freq_fg = hc_behavs[,"FG"]/rowSums(hc_behavs),
                             Freq_sd = hc_behavs[,"SD"]/rowSums(hc_behavs))

# Calculate frequency of mom hc after departure
hc_behavs <- table(mom_data$Code, mom_data$DiffHI)

hc_mom_freq <- data.frame(Code = row.names(hc_behavs),
                               Freq_fg = hc_behavs[,"FG"]/rowSums(hc_behavs),
                               Freq_sd = hc_behavs[,"SD"]/rowSums(hc_behavs))

# Calculate frequency of calf hc after departure
hc_behavs <- table(calves_data$Code, calves_data$DiffHI)

hc_calf_freq <- data.frame(Code = row.names(hc_behavs),
                               Freq_bg = hc_behavs[,"BG"]/rowSums(hc_behavs),
                               Freq_fg = hc_behavs[,"FG"]/rowSums(hc_behavs),
                               Freq_sd = hc_behavs[,"SD"]/rowSums(hc_behavs))


#### PART 4: Run Correlation Tests ####

# Read in Data
hr_mom_calf <- readRDS("hr_mom_calf.RDS")
hr_mom <- readRDS("hr_mom.RDS")
hr_calf <- readRDS("hr_calf.RDS")

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

