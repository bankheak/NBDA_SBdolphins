##############  Multi-network Network-Based Diffusion Analysis

#######################################################################################################################################
####### PART 1: Wrangling data:

# load packages
if(!require(devtools)){install.packages('devtools'); library(devtools)} # To load NBDA
if(!require(asnipe)){install.packages('asnipe'); library(asnipe)} # get_group_by_individual
if(!require(sf)){install.packages('sf'); library(sf)} # Convert degrees to meters
if(!require(sp)){install.packages('sp'); library(sp)} # Convert degrees to meters
if(!require(adehabitatHR)){install.packages('adehabitatHR'); library(adehabitatHR)} # Caluculate MCPs and Kernel density 
if(!require(kinship2)){install.packages('kinship2'); library(kinship2)} # genetic relatedness

setwd("../../NBDA")
load_all()

# all networks available under https://datadryad.org/review?doi=doi:10.5061/dryad.sc26m6c.
setwd("../NBDA_SBdolphins/Data") # set working directory
source("../Code/functions.R") # nxn

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

# Subset ILVs to include these filtered individuals
Codes <- unique(filtered_data$Code)
ILV_all <- subset(ILV_all, Alias %in% Codes)

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
hi_data <- filtered_data[filtered_data$Confirmed_HI == 1 & filtered_data$Date != min(filtered_data$Date), ]

# Step 2: Get the first year each Alias showed the behavior
first_hi_year <- aggregate(Date ~ Code, data = hi_data, FUN = min)

# Step 3: Create a new column for order of acquisition
first_hi_year$HI_Order_acquisition <- as.numeric(first_hi_year$Date - min(first_hi_year$Date))

# Step 4: Merge this info back into ILV_all
ILV_all <- merge(ILV_all, first_hi_year[, c("Code", "HI_Order_acquisition")],
                 by.x = "Alias", by.y = "Code", all.x = TRUE)

# Step 5: Replace NA with 0 for individuals who had the behavior in 1995
ILV_all$HI_Order_acquisition[ILV_all$Demons_HI_forage == "yes"] <- 0

# Save data
write.csv(ILV_all, "ILV_all_subset.csv")

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
####### PART 2: Apply NBDA to HI data:

# Add ILV data
ILV_all <- read.csv("ILV_all_subset.csv")

# Extract Confirmed_HI (learners and demonstrators)
Confirmed_HI <- subset(ILV_all, subset = ILV_all$HI_Indiv == 1)
Confirmed_HI <- Confirmed_HI[order(Confirmed_HI$HI_Order_acquisition),]

# get ID codes of all HI
HI_all <- Confirmed_HI$Alias
# extract IDs of all HI that are treated as learners
HI_learners <- as.vector(subset(Confirmed_HI$Alias, subset=Confirmed_HI$Demons_HI_forage=="no"))
# extract IDs of all HI treated as demonstrators
HI_demons <- as.vector(subset(Confirmed_HI$Alias, subset=Confirmed_HI$Demons_HI_forage=="yes"))

# extract ID names from data file
IDs <- sort(as.vector(unique(filtered_data[,"Code"])))

# extract order of acquisition
order <- NULL # create an object to store the vector of acquisition

for (i in 1:length(HI_learners)){ # for each HI, extract the position in the networks and ILV data frame
  order[i] <- which(IDs==HI_learners[i])
}

order <- as.vector(order)
OAc <- order

# extract positions of demonstrators
demons <- NULL # create an object to store the vector of acquistion

for (i in 1:length(HI_demons)){ # for each HI demonstrator, extract the position in the networks and ILV data frame
  demons[i] <- which(IDs==HI_demons[i])
}

# contains positions of all HU demonstrators
demons <- as.vector(demons)

# create vector of length(IDs) with 0 for non-demonstrators and 1 for demonstrators
demons_vector <- c(rep(0,length(IDs)))

for (i in demons){
  demons_vector[i] <- 1
  
}

## prepare individual-level variables
Sex <- ifelse(ILV_all$Sex == "Female", 1, 0)
ILV_all$BirthYear <- as.numeric(ILV_all$BirthYear)
Age <- ifelse(is.na(ILV_all$BirthYear), 1985, ILV_all$BirthYear)

# Read in matrices
SRI_vert_all <- readRDS("SRI_vert_all.RData")
nxn <- readRDS("nxn.RData")
ecol_all <- as.array(readRDS("ecol_all.RData"))

# Get rid of vert data in hor network
SRI_hor_no_vert_all <- lapply(nxn, function(mat) {
  for (r in rownames(mat)) {
    for (c in colnames(mat)) {
      if (!is.na(SRI_vert_all[r, c]) && SRI_vert_all[r, c] == 1) {
        mat[r, c] <- 0
      }
    }
  }
  mat
})

# Turn vertical matrix into list
SRI_vert_all <- replicate(18, SRI_vert_all, simplify = FALSE)
SRI_vert_all <- as.array(SRI_vert_all)

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

# Create array with all matrices
n_indiv <- nrow(SRI_vert_all[[1]])
n_years <- length(SRI_vert_all)
n_networks <- 3

# Initialize the 4D array
assMatrix.B <- array(NA, dim = c(n_indiv, n_indiv, n_networks, n_years))

# Fill the array
for (t in 1:n_years) {
  assMatrix.B[,,1,t] <- SRI_vert_all[[t]]
  assMatrix.B[,,2,t] <- SRI_hor_no_vert_all[[t]]
  assMatrix.B[,,3,t] <- ecol_all[[t]]
}

# ILVs
Sex <- matrix(data = Sex, nrow=length(IDs), byrow=F) # all ILVs need to go into a matrix
Age <- matrix(data = Age, nrow=length(IDs), byrow=F)

ILVs <- c("Sex","Age")

label <- "HIC"

# extract the Confirmed_HI learners with no maternity data available
HI_filter <- subset(Confirmed_HI, subset=Confirmed_HI$Demons_HI_forage=="no")
HI_filter2 <- sort(subset(HI_filter$id_individual, subset=is.na(HI_filter$Mom)))

vec <- NULL
for (i in 1:length(HI_filter2)){
  a <- which(IDs==HI_filter2[i])
  vec[i] <- a
} # get position of HI with no maternity data
# they get set to 0 in the presence matrix (NBDAfilterfunction)

filter <- paste0(label,"_", vec)

# create NBDA Data Object
nbdaDataHI.C <- nbdaData(label=label, assMatrix=assMatrix.B, asoc_ilv=ILVs, 
                         int_ilv=ILVs, multi_ilv=ILVs, orderAcq=OAc,asocialTreatment="constant", 
                         demons = demons_vector) # creates OADA object

# apply filter to exclude individuals without maternity data as learners
nbdaDataHI.C.filter <- filteredNBDAdata(nbdadata=nbdaDataHI.C, filter="id", exclude=filter)

#Then fit the model:
model1_multiNet<-oadaFit(nbdaDataHI.C)
#And get the output:
data.frame(Variable=model1_multiNet@varNames,MLE=model1_multiNet@outputPar,SE=model1_multiNet@se)

# Save filtered data
saveRDS(nbdaDataHI.C.filter, "nbdaDataHI.C.filter.RData")

# the first four positions correspond to the networks (vertical, horizontal, ecology, relatedness),
# the following 4 to int.ILV, then 4 to asoc ILV and then 4 to multi ILV
# all the multi.ILV positions are set to 0


# check ILVs for asoc (ILVs only influence asocial learning), int (ILVs influence asocial and social learning independently)
# and multi (ILVs influence both asocial and social learning to the same extent).
# in this analysis, multiILV will be set to 0 and therefore not estimated (in the constraintsVectMatrix)
nbdaDataHI.C.filter@asoc_ilv
nbdaDataHI.C.filter@int_ilv
nbdaDataHI.C.filter@multi_ilv


# the following part creates a matrix, constraintsVectMatrix, with all possible combinations of networks and ILVs.
# This is then input to the oadaAICtable function below to fit each model
# An explanation of the constraintsVectMatrix is given below

# set number of networks and number of ILVs
num_networks <- 3
num_ILVs <- 2

vector <- seq(1:(num_networks+(2*num_ILVs))) # create a vector for the full model with all networks and ILVs (excluding multiILV slots which will all be set to 0)
count <- 0 # create an object 'count', which starts on 0

constraintsVect <- matrix(nrow = 10000000, ncol=(num_networks+(2*num_ILVs))) # create a matrix to save the combination of parameters in
constraintsVect[1,] <- seq(1:(num_networks+(2*num_ILVs))) # the first row gets filled with a sequence from 1:12 (all parameters will be estimated, none are set to 0)

for (i in 1:(num_networks+(2*num_ILVs)-1)){ # a loop for each number of parameters to be estimated
  array <- combn(vector, i, FUN = NULL, simplify = TRUE) # for each number of paramters to be estiamted (e.g. 2) create all possible combinations of numbers between 1:7 (e.g. 2&7, 1&5 etc)
  
  for (j in 1:length(array[1,])){ # for each of those combinations
    vector2 <- seq(1:((num_networks+(2*num_ILVs))-i)) # create a second vector with 6-i free spaces
    position <- array[,j] # for each created combination
    count <- count+1 # add +1 to the count
    
    for (k in position){ # at each possible position
      vector2 <- append(vector2, 0, after=k-1) # add a 0 (e.g. 1 0 2 3 ...; 1 2 0 3 4 5 ...; 1 2 3 0 4 5 ....)
    }
    constraintsVect[count+1,] <- vector2 # and save the resulting order in a matrix
  }
}


constraintsVect <- na.omit(constraintsVect) # remove all NAs from the matrix
constraintsVect <- rbind(constraintsVect, rep.int(0,(num_networks+2*(num_ILVs)))) # add a last row with all 0

constraintsVect <- cbind(constraintsVect, matrix(0,ncol=num_ILVs, nrow=length(constraintsVect[,1]))) ## add 2 columns at the end with all 0 (multi_ILV)

constraintsVectMatrix<-constraintsVect

# Save vect matrix
saveRDS(constraintsVectMatrix, "constraintsVectMatrix.RData")

# Each line of the resulting object specifies a model
# Each element in the line corresponds to a parameter in the model. When an element is zero, that paramter is constrained
# =0. When two elements have the same value, they are constrained to have the same value (not relevant here).
# For example:
constraintsVectMatrix[1,]
# The first four elements (1-4) are the s parameters for each network- so in this model all networks are included.
# The next four elements (5-8) are the parameters determining the effect each ILV has on asocial learning. So in this model
# all ILVs are assumed to affect asocial learning.
# The next five elements (9-12) are the parameters determining the effect each ILV has on social learning. So in this model
# all ILVs are assumed to affect social learning.
# The final four elements determine the effect each parameter has on both asocial and social learning (the multiplicative
# NBDA model). In our analysis we estimate the effects each ILV has on asocial and social learning independently, so these parameters
# are constrained to be zero for all models fitted.

# run NBDA using the NBDA Data object and the constraitnsVectMatrix
# this fits every model specified by the constrainstsVectMatrix matrix
tableHI.C.filter<-oadaAICtable(nbdadata=nbdaDataHI.C.filter, constraintsVectMatrix=constraintsVectMatrix,writeProgressFile = T)
print(tableHI.C.filter)

save(tableHI.C.filter, file="AIC table HI.Rdata")
load("AIC_table_HI.Rdata")

write.csv(as.data.frame(tableHI.C.filter@printTable), "AIC table HI.csv")

##Create a new object with a printTable that excludes unfitted model
newTableHI<-tableHI.C.filter
newTableHI@printTable<-tableHI.C.filter@printTable[!is.nan(tableHI.C.filter@printTable$aicc)&!is.na(tableHI.C.filter@printTable$aicc),]

tableHI.C.filter@aicc<-tableHI.C.filter@aicc[!is.nan(tableHI.C.filter@aicc)&!is.na(tableHI.C.filter@aicc)]
tableHI.C.filter@MLEs<-tableHI.C.filter@MLEs[!is.nan(tableHI.C.filter@aicc)&!is.na(tableHI.C.filter@aicc),]
newTableHI@MLEilv<-tableHI.C.filter@MLEilv[!is.nan(tableHI.C.filter@aicc)&!is.na(tableHI.C.filter@aicc),]
newTableHI@MLEint<-tableHI.C.filter@MLEint[!is.nan(tableHI.C.filter@aicc)&!is.na(tableHI.C.filter@aicc),]


newTableHI@printTable<-newTableHI@printTable[order(newTableHI@printTable$aicc),]
newTableHI@printTable$deltaAICc<-newTableHI@printTable$aicc-newTableHI@printTable$aicc[1]

# calculate support for fitted models
newTableHI@printTable$RelSupport<- exp(-0.5*newTableHI@printTable$deltaAICc)
newTableHI@printTable$AkaikeWeight<-newTableHI@printTable$RelSupport/sum(newTableHI@printTable$RelSupport)


newTableHI@deltaAIC<-newTableHI@aicc-min(newTableHI@aicc)

newTableHI@RelSupport<- exp(-0.5*newTableHI@deltaAIC)
newTableHI@AkaikeWeight<-newTableHI@RelSupport/sum(newTableHI@RelSupport)


dim(tableHI.C.filter@printTable)[1]-dim(newTableHI@printTable)[1]
dim(tableHI.C.filter@printTable)[1]

##2 models could not be fitted out of 4096- probably too many parameters for the dataset

# save reduced AIC table as csv
write.csv(newTableHI@printTable, file="HI.AIC_table_CORRECTED.csv")


# obtain network support for each network combination
networksSupport_HI<-networksSupport(newTableHI)
networksSupport_HI
write.csv(networksSupport_HI, file="networksSupport_HI.csv")
#83.7% support for vertical social network only (1:0:0:0)

# extract support for each variable
variable_support <- variableSupport(newTableHI, includeAsocial = T)
variable_support
write.csv(variable_support, file="variable_support_HI.csv")

# extract model averaged medians
MLE_med <- modelAverageEstimates(newTableHI,averageType = "median")
MLE_med

write.csv(MLE_med, "MLE_HI.csv")

#######################################################################################################################################
####### PART 3: 95% confidence intervals using profile likelihood:

#This is vital for s parameters since CIs based on SEs will be highly misleading due to frequent assymetry in the profile likelihood
# constraintsVectMatrix[4076,] for best model (vertical social learning + social.sex)

# Read filtered data and matrix
nbdaDataHI.C.filter <- readRDS("nbdaDataHI.C.filter.RData")
constraintsVectMatrix <- readRDS("constraintsVectMatrix.RData")

bestModelData<-constrainedNBDAdata(nbdadata=nbdaDataHI.C.filter,constraintsVect=constraintsVectMatrix[111,]) # ask Sonja how to find this?
model.best.social<-oadaFit(bestModelData)
model.best.social@outputPar
# [1] 1.233004e+10 -4.840486e+00
model.best.social@optimisation
model.best.social@aicc


# extract profile likelihood. which=1 extracts the first parameter (s parameter for vertical social learning)
plotProfLik(which=1,model=model.best.social,range=c(0,1e30), resolution=20)
#Here we can see that we cannot set an upper limit on s (see explanation below)

#Zoom in to locate the lower limit for s
plotProfLik(which=1,model=model.best.social,range=c(0,500), resolution=20)
plotProfLik(which=1,model=model.best.social,range=c(30,40), resolution=20)
profLikCI(which=1,model=model.best.social,lowerRange=c(30,40))

#Lower CI Upper CI
#33.08928       Inf

#To explain why we cannot set an upper limit on s, we extract a table with a row for each acquisition event
#Each row we get the connection of the individual that learned
#The total connections to informed individuals across the population- in this case equal to the number of dolphins with a sponging mother, 
#since connections are binary
#and the maximum connection strength (again this should be 1, but just to confirm)

eventTable<-NULL
for(i in 1:9)
{eventTable<-rbind(eventTable,(c((bestModelData@stMetric[bestModelData@event.id==unique(bestModelData@event.id)[i]])[
  bestModelData@status[bestModelData@event.id==unique(bestModelData@event.id)[i]]==1],
  sum(bestModelData@stMetric[bestModelData@event.id==unique(bestModelData@event.id)[i]]),
  max(bestModelData@stMetric[bestModelData@event.id==unique(bestModelData@event.id)[i]])))
)
}
dimnames(eventTable)[[2]]<-c("Connection of Learner","Sum of Connections","Maximum Connection")
eventTable
#We can see that in every case the individual to learn shelling has an informed mother
#There are a variable number of individuals with sponging mothers across events. But importantly, the dolphin to acquire the behaviour has
#the (joint) maximum connection to informed individuals.

#In any diffusion where, for every event, the individual to acquire the behaviour is always the one with the maximum connection to informed
#individuals (even if it is joint maximum), we cannot set an upper limit on s with OADA. Since the next individual to learn is always the 
#one that the network would predict as being most likely (or joint most likely), a value of s=Inf is plausible.

# which=2 extracts second parameter (gender)
plotProfLik(which=2,model=model.best.social,range=c(-10,-2), resolution=20)
profLikCI(which=2,model=model.best.social,lowerRange=c(-10,-6),upperRange = c(-4,-2))

#Lower CI  Upper CI
#-7.971296 -2.246465
# back-transform
exp(c(4.840486e+00,7.971296, 2.246465))
#[1]  126.530831 2896.608938    9.454256
# Females are an estimated 126x (95% CI= 9.5-2890) faster to learn the behaviour from their mothers than males are


# extract what proportion of spongers are estimated to have learned sponging socially from their mothers
prop.solve.social.byevent <- oadaPropSolveByST.byevent(nbdadata = bestModelData, model=model.best.social) # outputs of 1 for each events means that all sponger offspring have learned socially from their mothers
prop.solve.social <- oadaPropSolveByST(nbdadata = bestModelData, model=model.best.social) # 100% are estimated to have learned socially

prop.solve.social

#To get the estimates for the lower bound we should find the corresponding value of the other parameters to plug in when s1 is constrained to this value
bestModelDataS1LowerBound<-constrainedNBDAdata(nbdadata=nbdaDataHI.C.filter,constraintsVect =constraintsVectMatrix[112,],offset=c(33.08928,rep(0,15)))
bestModelS1LowerBound<-oadaFit(bestModelDataS1LowerBound,type="asocial")
bestModelS1LowerBound@outputPar
#Now plug into the prop solve function in one of these two ways:
prop.solve.social.lower <- oadaPropSolveByST(par= c(33.08928, bestModelS1LowerBound@outputPar),model=NULL, nbdadata = bestModelData)
prop.solve.social.lower <- oadaPropSolveByST(model=bestModelS1LowerBound, nbdadata = bestModelDataS1LowerBound)
prop.solve.social.lower

#AT least 98.9% learned by social transmission from mothers

# In theory repeat for upper limit, but here there is no need since the upper limit is s=inf,
# so the upper limit is 100%

