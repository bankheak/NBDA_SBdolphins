## --- Bayesian Multi-network Diffusion Analysis --- ##

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

# all networks available under https://datadryad.org/review?doi=doi:10.5061/dryad.sc26m6c.
setwd("../Data") # set working directory
source("../Code/functions.R") # nxn

#######################################################################################################################################
####### PART 1: Wrangling data:

# Read ILVs
ILV_all <- read.csv("ILV_dem.csv", header=TRUE, sep=",")
ILV_all <- ILV_all[, c("Alias", "HI_Indiv", "Mom", "Sex", "BirthYear")]

# Subset original data
data_1 <- read.csv("93_04_data.csv")
data_2 <- read.csv("05_14_data.csv")
orig_data <- merge(data_1, data_2, all = T)
write.csv(orig_data, "orig_data.csv")

# Read orig_data
orig_data <- read.csv("orig_data.csv")
orig_data$Confirmed_HI <- ifelse(orig_data$ConfHI != "0", 1, 0)
orig_data$Date <- as.Date(as.character(orig_data$Date), format="%d-%b-%y")

# Add year
orig_data$Year <- as.numeric(format(orig_data$Date, format = "%Y"))

# Subset data to first acquisition event
first_idx <- which(orig_data$Confirmed_HI == 1)[1]
orig_data <- orig_data[first_idx:nrow(orig_data), ]

# Filter data to include individuals that were seen at least 10 times 
tab <- table(orig_data$Code)
codes_in_all <- rownames(tab > 10)
filtered_data <- orig_data[orig_data$Code %in% codes_in_all, ]

write.csv(filtered_data, "filtered_data.csv")

#######################################################################################################################################
####### PART 2: Create Networks:

# Horizontal network -----------------------------------------

# Read in filtered data
filtered_data <- read.csv("filtered_data.csv")

# Add individual data
ILV_all <- read.csv("ILV_all_subset.csv")

# Group each individual by date and sighting
group_data <- filtered_data[,c("Date","Sighting","Code","Year", "ConfHI")]
group_data$Confirmed_HI <- ifelse(group_data$ConfHI != "0", 1, 0)
group_data$Group <- cumsum(!duplicated(group_data[1:2])) # Create sequential group # by date
group_data <- group_data[,c(1, 3, 4, 6, 7)] # Subset ID and group #

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

# Now create a list for each year
group_data_list <- split(group_data, group_data$Year)
saveRDS(group_data_list, "group_data_list.RData")

# Calculate Gambit of the group
create_gbi <- function(list_years) {
  gbi <- list()
  for (i in seq_along(list_years)) {
    
    # Gambit of the group index
    gbi[[i]] <- get_group_by_individual(list_years[[i]][,c("Code", "Group")], data_format = "individuals")
  }
  return(gbi)                                      
}

gbi <- create_gbi(group_data_list)
saveRDS(gbi, "gbi.RData")

# Create association matrix
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

nxn <- create_nxn(gbi)

# Save nxn
saveRDS(nxn, "nxn.RData")

# Vertical network -----------------------------------------

# Read in original data
filtered_data <- read.csv("filtered_data.csv") 

# Read in Horizontal data
nxn <- readRDS("nxn.RData")

# Find the demographics of the population
ILV_pat <- read.csv("Paternity_data.csv")

Codes <- unique(filtered_data$Code)
ILV_pat <- subset(ILV_pat, Alias %in% Codes)

# Subset paternity data
ILV_pat <- data.frame(Code = ILV_pat$Alias,
                      Mom = ILV_pat$Mom)
SRI_vert_all <- vert.func(ILV_pat)

# Save vert
saveRDS(SRI_vert_all, "SRI_vert_all.RData")

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
  
  ids <- filtered_list[[i]]$Code
  coordinates <- filtered_list[[i]][, c("StartLon", "StartLat")]
  
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
  order_rows <- rownames(SRI_vert_all)
  order_cols <- colnames(SRI_vert_all)
  
  # Apply the order to each matrix in the list
  ecol_all[[i]] <- kov[order_rows, order_cols] 
}

# Save eco dat
saveRDS(ecol_all, "ecol_all.RData")

# Read relatedness network -----------------------------------------
ILV_pat <- read.csv("Paternity_data.csv") 

# Order data
order_rows <- rownames(nxn)
order_cols <- colnames(nxn)

# Reorder rows in 'ILV' based on 'order_rows'
ILV <- ILV_pat[ILV_pat$Alias %in% order_rows, ]
ILV <- ILV[match(order_rows, ILV$Alias), ]

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

# Only take the ids that aren't found in the 117 list
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


#######################################################################################################################################
####### PART 3: Check variation in associations:

# Read in network
nxn <- readRDS("nxn.RData")
gbi <- readRDS("gbi.RData")

# Done in the HPC --------------------------------------------------------------

#  Create 1000 random group-by-individual binary matrices
reps<- 1000
n.cores <- detectCores()
registerDoParallel(n.cores)

source("functions.R")
nF <- lapply(gbi, function (group_index) null(group_index, iter=reps))
saveRDS(nF, "nF.RData")

#' Calculate the association and CV for each of the 1000 permuted matrices to
#' create null distribution
cv_null <- rep(NA,reps)

foreach(i = 1:reps, 
        .combine = c) %dopar% { 
          sri_null = as.matrix(SRI.func(nF[[i]]))
          cv_null[i] <- ( sd(sri_null) / mean(sri_null) ) * 100}

stopImplicitCluster()

saveRDS(cv_null, "cv_null.RData")

# Next take results from the HPC ------------------------------------------------

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

#######################################################################################################################################
####### PART 4: Create acquisition data for model input:

# Read in all network data
nxn <- readRDS("nxn.RData")
SRI_vert_all <- readRDS("SRI_vert_all.RData")
ecol_all <- as.array(readRDS("ecol_all.RData"))

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

# Turn vertical matrix into list
SRI_vert_all <- replicate(20, SRI_vert_all, simplify = FALSE)
SRI_vert_all <- as.array(SRI_vert_all)

# Subset nxn
nxn <- nxn[3:22]

# Get rid of IDs without data in nxn
## Vertical network
for (i in seq_along(SRI_vert_all)) {
  # Get the row/column names from the nxn matrix
  target_names <- rownames(nxn[[i]])
  
  # Subset the SRI matrix to match those names
  SRI_vert_all[[i]] <- SRI_vert_all[[i]][target_names, target_names, drop = FALSE]
}

# Ecol network
for (i in seq_along(ecol_all)) {
  # Get the row/column names from the nxn matrix
  target_names <- rownames(nxn[[i]])
  
  # Subset the SRI matrix to match those names
  ecol_all[[i]] <- ecol_all[[i]][target_names, target_names, drop = FALSE]
}

# Put matrices into data frame
edge_list <- do.call(rbind, lapply(seq_along(nxn), function(t) {
  mat <- nxn[[t]]
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
filtered_data <- read.csv("filtered_data.csv")
ILV_all <- read.csv("ILV_dem.csv", header=TRUE, sep=",")
ILV_all <- ILV_all[, c("Alias", "HI_Indiv", "Mom", "Sex", "BirthYear")]

# Subset event data to include only edge_list ids
ILV_all <- subset(ILV_all, Alias %in% unique(edge_list$focal))

# Add order of acquisition data
# Create demonstrator column
ILV_all$Demons_HI_forage <- ifelse(
  ILV_all$Alias %in% unique(filtered_data$Code[filtered_data$Confirmed_HI == 1 & 
                                                 filtered_data$Date == min(filtered_data$Date)]),
  "yes",
  "no"
)

# Create acquisition data
# Step 1: Filter filtered_data for confirmed HI behavior after first date
hi_data <- filtered_data[filtered_data$Confirmed_HI == 1 & filtered_data$Year != min(filtered_data$Year), ]

# Step 2: Get the first year each Alias showed the behavior
first_hi_year <- aggregate(Year ~ Code, data = hi_data, FUN = min)

# Step 3: Create a new column for order of acquisition
first_hi_year$HI_Order_acquisition <- as.numeric(first_hi_year$Year - min(first_hi_year$Year))

# Step 4: Merge this info back into ILV_all
ILV_all <- merge(ILV_all, first_hi_year[, c("Code", "HI_Order_acquisition")],
                 by.x = "Alias", by.y = "Code", all.x = TRUE)

# Step 5: Replace NA with 0 for individuals who had the behavior in 1995
ILV_all$HI_Order_acquisition[ILV_all$Demons_HI_forage == "yes"] <- 0

# Add end time
ILV_all$t_end <- max(na.omit(ILV_all$HI_Order_acquisition))

# Change individuals who didn't require behavior to t_end +1
ILV_all$time <- ifelse(is.na(ILV_all$HI_Order_acquisition), 
                                       max(na.omit(ILV_all$HI_Order_acquisition)) + 1, 
                                       ILV_all$HI_Order_acquisition)

# Edit the other needed columns
ILV_all$id <- ILV_all$Alias
ILV_all$trial <- 1

# Get rid of other columns
event_data <- ILV_all[, c("trial", "id", "time", "t_end")]

ILV_all$Sex <- ifelse(ILV_all$Sex == "Female", 1, 0)
ILV_all$BirthYear <- as.numeric(ILV_all$BirthYear)
ILV_all$BirthYear <- ifelse(is.na(ILV_all$BirthYear), 1985, ILV_all$BirthYear)

# Constant ILVs
ILV_c <- data.frame(id = ILV_all$Alias,
                      sex = ILV_all$Sex)
# Time varying ILVs
ILV_tv <- data.frame(
  trial = 1,
  id = rep(ILV_all$Alias, each = 20),
  time = rep(1:20, times = length(ILV_all$Alias)),       
  year = rep(1995:2014, times = length(ILV_all$Alias)),  
  age = rep(1995:2014, times = length(ILV_all$Alias)) - 
    rep(ILV_all$BirthYear, each = 20)
)

# Separate HI Behaviors to create weighted HI prop variable
#' BG = Beg: F, G, H
#' SD = Scavenge and Depredation: A, B, C, D, E
#' FG = Fixed Gear Interaction: P
# Change the code using ifelse statements
filtered_data <- subset(filtered_data, Year %in% 1995:2014)
filtered_data <- subset(filtered_data, Code %in% unique(edge_list$focal))

filtered_data_list <- split(filtered_data, filtered_data$Year)

subset_HI <- function(aux_data) {
  for (i in seq_along(aux_data)) {
    aux_data[[i]]$DiffHI <- ifelse(aux_data[[i]]$ConfHI %in% c("F", "G", "H"), "BG",
                                   ifelse(aux_data[[i]]$ConfHI %in% c("A", "B", "C", "D", "E"), "SD",
                                          ifelse(aux_data[[i]]$ConfHI %in% c("P"), "FG", "None")))
  }
  return(aux_data)  # Return the modified list of data frames
}

aux <- subset_HI(filtered_data_list)

# Categorize DiffHI to IDs
diff_raw <- function(aux_data) {
  rawHI_diff <- lapply(aux_data, function(df) {
    table_df <- as.data.frame(table(df$Code, df$DiffHI))
    colnames(table_df) <- c("Code", "DiffHI", "Freq")
    
    return(table_df)
  })}

rawHI_diff <- diff_raw(aux)

# Categorize ID to Sightings
ID_sight <- function(aux_data) {
  IDbehav <- lapply(aux_data, function(df) {
    data <- as.data.frame(table(df$Code))
    colnames(data) <- c("Code", "Sightings")
    # Order data
    order_rows <- rownames(nxn[[1]])
    
    # Now reorder the dataframe
    data <- data %>%
      arrange(match(Code, order_rows))
    
  })
  return(IDbehav)
}

IDbehav <- ID_sight(aux)

# Create a frequency count for each HI behavior
get_IDHI <- function(HI, IDbehav_data, rawHI_diff_data) {
  lapply(seq_along(IDbehav_data), function(i) {
    df <- IDbehav_data[[i]]
    HI_freq <- rawHI_diff_data[[i]]$Freq[rawHI_diff_data[[i]]$DiffHI == HI]
    df$Behav <- HI_freq[match(df$Code, rawHI_diff_data[[i]]$Code)]
    colnames(df) <- c("Code", "Sightings", "Behav")
    return(df)
  })
}

IDbehav_HI <- get_IDHI(c("BG", "FG", "SD"), IDbehav, rawHI_diff)

# Proportion of Sightings spent in HI
Prop_HI <- function(IDbehav) {
  lapply(seq_along(IDbehav), function(i) {
    df <- IDbehav[[i]]
    df$HIprop <- as.numeric(df$Behav) / as.numeric(df$Sightings)
    df$HIprop[is.na(df$HIprop)] <- 0
    # Keep only 'Code' and 'HIprop' columns
    df <- df[, c('Code', 'HIprop')]
    df
  })
}

prob_HI <- Prop_HI(IDbehav_HI)

# Convert list of HIprop vectors into a matrix
HI_matrix <- do.call(rbind, lapply(seq_along(prob_HI), function(t) {
  dat <- prob_HI[[t]]
  ids <- dat$Code
  HIProp <- dat$HIprop
  
  data.frame(
    id = ids,
    HIProp = HIProp,
    time = t
  )
}))

# Order edge data
edge_list <- edge_list[order(edge_list$time), ]

# Add vertical network to edge_list
edge_list_vert <- do.call(rbind, lapply(seq_along(SRI_vert_all), function(t) {
  mat <- SRI_vert_all[[t]]
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

edge_list <- merge(edge_list, edge_list_vert, all = T)

# Add transmission weights
ILV_tv$weight <- HI_matrix$HIProp

# Input data
data_list <- import_user_STb(
  event_data = event_data,
  networks = edge_list,
  ILV_c = ILV_c,
  ILV_tv = ILV_tv,
  ILVi = c("age", "weight"),
  ILVs = c("sex")
)

# For test
edge_list_test <- edge_list$assoc

data_list <- import_user_STb(
  event_data = event_data,
  networks = edge_list_test
)

saveRDS(data_list, "data_list.RData")

#######################################################################################################################################
####### PART 5: Run the model and summary outputs:

# Input data_list
data_list <- readRDS("data_list.RData")

# Generate the model
model_full <- generate_STb_model(data_list, gq = T, est_acqTime = T)

# Fit the model
full_fit <- fit_STb(data_list,
                    model_full,
                    parallel_chains = 4,
                    chains = 4,
                    cores = 4,
                    iter = 2000,
                    refresh=1000
)

STb_save(full_fit, output_dir = "cmdstan_saves", name="my_first_fit")

# View output
STb_summary(full_fit, digits = 3)

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

draws_df <- as_draws_df(full_fit$draws(variables = "acquisition_time", inc_warmup = FALSE))

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

# Estimated versus observed 
acqdata = extract_acqTime(full_fit, data_list)

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

# 