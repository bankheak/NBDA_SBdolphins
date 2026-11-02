# load packages
if(!require(networkDynamic)){install.packages('networkDynamic'); library(networkDynamic)} # To load NBDA
if(!require(ndtv)){install.packages('ndtv'); library(ndtv)} # get_group_by_individual
if(!require(doParallel)){install.packages('doParallel'); library(doParallel)} # Faster computing
if(!require(network)){install.packages('network'); library(network)} # For assigning coordinates to nodes %v%

setwd("../Data")

# Read in data
nxn <- readRDS("nxn_sd.RData")
ILV_all <- read.csv("Individuals_Residency.csv", header=TRUE, sep=",")
group_data <- read.csv("group_data_sd.csv")

# Choose a HC Behavior
HC <- "SD"

# Get rid of data without individuals
ids_to_keep <- with(group_data, tapply(DiffHI == HC, Code, function(x) any(x, na.rm = TRUE)))
ids_to_keep <- names(ids_to_keep)[ids_to_keep]

# Subset the data to keep only those individuals
HI_ids <- subset(group_data, Code %in% ids_to_keep)

# Only include individuals with IDs in HI_data
ILV_all <- subset(ILV_all, Alias %in% HI_ids$Code)

# Subset data
for (i in seq_along(nxn)) {
  target_names <- as.character(ILV_all$Alias)
  
  rn <- rownames(nxn[[i]])
  cn <- colnames(nxn[[i]])
  
  keep_r <- target_names[target_names %in% rn]
  keep_c <- target_names[target_names %in% cn]
  
  nxn[[i]] <- nxn[[i]][keep_r, keep_c, drop = FALSE]
}

# Get all unique IDs across all matrices
total_ids <- unique(unlist(lapply(nxn, rownames)))

# Update each matrix to include all IDs, filling missing rows/columns with zeros
nxn <- lapply(nxn, function(mat) {
  # Current IDs in this matrix
  current_ids <- rownames(mat)
  
  # Create a full zero matrix with all IDs
  full_mat <- matrix(0, nrow = length(total_ids), ncol = length(total_ids),
                     dimnames = list(total_ids, total_ids))
  
  # Fill in existing values
  full_mat[current_ids, current_ids] <- mat
  
  return(full_mat)
})

# Edgelist: Nodes (i & j) and edge (or link) weight
source("../code/functions.R")

## Create social network
net_list <- lapply(nxn, function (df) {
  as.network(df, matrix.type='adjacency',
             directed = F,
             ignore.eval=FALSE,
             names.eval='weight')
})

# Create a dynamic network object

# Years for each matrix
years <- 1:max(group_data$event_period)

# Create dynamic network with real years
dyn_net <- networkDynamic(
  network.list = net_list,
  onsets = years,
  vertex.pid="vertex.names",
  termini = years + 1  # each network lasts 1 year
)

#saveRDS(dyn_net, "dyn_net_sd.RData")

# Get all edge IDs in the dynamic network
all_edge_ids <- seq_len(network.edgecount(dyn_net))

for (i in seq_along(net_list)) {
  # Get weights for this time slice
  weights <- get.edge.attribute(net_list[[i]], "weight")
  
  # Find which edges are active in this time slice
  active_edges <- which(is.active(dyn_net, e = all_edge_ids, at = years[i]))
  
  # Assign weights to active edges
  activate.edge.attribute(
    dyn_net,
    prefix = "weight",
    value = weights,
    onset = years[i],
    terminus = years[i] + 1,
    e = active_edges
  )
}

# Assume group_data has columns: Alias, Year, Conf_HI
# Get all unique node names
vn <- network.vertex.names(dyn_net)

# Keep track of nodes that have ever been HC
cumulative_hi_nodes <- character(0)

for (year in years) {
  
  # Nodes that become HC this year
  new_hi_nodes <- unique(
    group_data$Code[
      group_data$event_period == year &
        group_data$DiffHI == HC
    ]
  )
  
  # Add to cumulative set
  cumulative_hi_nodes <- union(cumulative_hi_nodes, new_hi_nodes)
  
  # Once red, always red
  colors <- ifelse(vn %in% cumulative_hi_nodes, "red", "gray")
  
  activate.vertex.attribute(
    dyn_net,
    prefix = "color",
    value = colors,
    onset = year,
    terminus = year + 1
  )
}

# Animate the network over time
# Export the dynamic network animation as an HTML file
render.d3movie(dyn_net,
               usearrows = FALSE,
               displaylabels = TRUE,
               vertex.cex = 2,
               edge.lwd = "weight",
               vertex.col = "color",  # Use dynamic color attribute
               main = "Dynamic Network Animation by Year",
               output.mode = "HTML",  # Ensure HTML output
               filename = "dynamic_network_animation.html",  # File name for export
               render.par = list(tween.frames = 10, show.time = TRUE))

