# load packages
if(!require(networkDynamic)){install.packages('networkDynamic'); library(networkDynamic)} # To load NBDA
if(!require(ndtv)){install.packages('ndtv'); library(ndtv)} # get_group_by_individual
if(!require(doParallel)){install.packages('doParallel'); library(doParallel)} # Faster computing
if(!require(network)){install.packages('network'); library(network)} # For assigning coordinates to nodes %v%

setwd("../Data")

# Read in data
nxn <- readRDS("nxn.RData")
ILV_all <- read.csv("ILV_all_subset.csv")

# Subset data
for (i in seq_along(nxn)) {
  # Get the row/column names from the nxn matrix
  target_names <- ILV_all$Alias
  
  # Subset the SRI matrix to match those names
  nxn[[i]] <- nxn[[i]][target_names, target_names, drop = FALSE]
}

# Fake subset
nxn <- nxn[3:22]
# Get the intersection of row names across all matrices
common_names <- Reduce(intersect, lapply(nxn, rownames))

# Subset each matrix to include only those common names
for (i in seq_along(nxn)) {
  nxn[[i]] <- nxn[[i]][common_names, common_names, drop = FALSE]
}

# Edgelist: Nodes (i & j) and edge (or link) weight
n.cores <- detectCores()
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
years <- 1995:2014

# Create dynamic network with real years
dyn_net <- networkDynamic(
  network.list = net_list,
  onsets = years,
  termini = years + 1  # each network lasts 1 year
)

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

# Assume filtered_data has columns: Alias, Year, Conf_HI
# Get all unique node names
all_nodes <- unique(filtered_data$Code)

# Add dynamic color attribute
for (i in seq_along(years)) {
  year <- years[i]
  
  # Nodes with Conf_HI == 1 for this year
  hi_nodes <- filtered_data$Code[filtered_data$Year == year & filtered_data$Confirmed_HI == 1]
  
  # Assign colors: red if in hi_nodes, gray otherwise
  colors <- ifelse(all_nodes %in% hi_nodes, "red", "gray")
  
  # Activate color attribute for this year
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

