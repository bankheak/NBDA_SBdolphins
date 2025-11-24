## --- Bayesian Multi-network Diffusion Analysis --- ##

# load packages
if(!require(devtools)){install.packages('devtools'); library(devtools)} # To load NBDA
# Install NBDA package
devtools::install_github("whoppitt/NBDA")
# install devtools if not already
devtools::install_github("michaelchimento/STbayes")
## Bayesian
if(!require(tidyr)){install.packages('tidyr'); library(tidyr)} # array
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
####### PART 3: Create acquisition data:

# Read in all network data
nxn <- readRDS("nxn.RData")
SRI_vert_all <- readRDS("SRI_vert_all.RData")
ecol_all <- as.array(readRDS("ecol_all.RData"))

# Turn vertical matrix into list
SRI_vert_all <- replicate(20, SRI_vert_all, simplify = FALSE)
SRI_vert_all <- as.array(SRI_vert_all)

# Subset nxn
nxn <- nxn[3:22]
# Get the intersection of row names across all matrices
common_names <- Reduce(intersect, lapply(nxn, rownames))

# Subset each matrix to include only those common names
for (i in seq_along(nxn)) {
  nxn[[i]] <- nxn[[i]][common_names, common_names, drop = FALSE]
}

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

# Prepare random effect for MCMC
num_nodes <- lapply(nxn, function(df) dim(df)[1])
node_names <- lapply(nxn, function(df) colnames(df))

# Separate IDs into i and j
node_ids_i <- lapply(num_nodes, function(df) matrix(rep(1:df, each = df), nrow = df, ncol = df))
node_ids_j <- lapply(node_ids_i, function(df) t(df))

# Abind nxn
upper_tri <- lapply(nxn, function(df) upper.tri(df, diag = TRUE))
edge_nxn <- abind(lapply(nxn, function(mat) mat[upper.tri(mat, diag = TRUE)]), along = 2)

# Create the edge_list
n_mats <- ncol(edge_nxn)
n_edges <- nrow(edge_nxn)

# Flatten all values
values <- as.vector(edge_nxn)

# Create group labels (1 to n_mats repeated for each row)
groups <- rep(1:n_mats, each = n_edges)

# Combine into a data frame
edge_data <- data.frame(value = values, group = groups)

one <- lapply(seq_along(node_ids_i), function(i) factor(as.vector(node_names[[i]][node_ids_i[[i]][upper_tri[[i]]]]), levels = node_names[[i]]))
two <- lapply(seq_along(node_ids_j), function(i) factor(as.vector(node_names[[i]][node_ids_j[[i]][upper_tri[[i]]]]), levels = node_names[[i]]))

# Create the edge_list
edge_list = data.frame(focal = unlist(one),
                       other = unlist(two),
                       trial = 1,
                       assoc = edge_data[, 1],
                       time = edge_data[, 2])
#HRO = unlist(lapply(kov, function (df) df[upper.tri(df, diag = TRUE)])),
#vert = )

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

# Save data
write.csv(ILV_all, "ILV_all_subset.csv")

# Add ILV data
ILV_all <- read.csv("ILV_all_subset.csv")

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

# Prepare individual-level variables
ILV_all <- subset(ILV_all, Alias %in% unique(edge_list$focal))

ILV_all$Sex <- ifelse(ILV_all$Sex == "Female", 1, 0)
ILV_all$BirthYear <- as.numeric(ILV_all$BirthYear)
ILV_all$BirthYear <- ifelse(is.na(ILV_all$BirthYear), 1985, ILV_all$BirthYear)

# Constant ILVs
ILV_c <- data.frame(id = ILV_all$Alias,
                      age = ILV_all$Sex)
# Time varying ILVs
ILV_tv <- data.frame(
  trial = 1,
  id = rep(ILV_all$Alias, each = 20),
  time = rep(1:20, times = length(ILV_all$Alias)),       
  year = rep(1995:2014, times = length(ILV_all$Alias)),  
  age = rep(1995:2014, times = length(ILV_all$Alias)) - 
    rep(ILV_all$BirthYear, each = 20)
)

# Separate HI Behaviors
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
HI_matrix <- do.call(cbind, lapply(prob_HI, function(x) x$HIprop))

# Order edge data
edge_list <- edge_list[order(edge_list$time), ]

# Add transmission weights
ILV_tv$weight <- as.vector(t(HI_matrix))

# Input data
data_list <- import_user_STb(
  event_data = event_data,
  networks = edge_list,
  ILV_c = ILV_c,
  ILV_tv = ILV_tv,
  ILVi = c("age", "weight"),
  ILVs = c("sex")
)
data_list <- import_user_STb(
  event_data = event_data,
  networks = edge_list
)

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
