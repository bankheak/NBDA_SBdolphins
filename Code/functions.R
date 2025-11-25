# 'Multi-network Network-Based Diffusion Analysis

###########################################################################
# FUNCTIONS
###########################################################################

# load necessary packages
library(asnipe)
library(Matrix)
library(matrixStats)
library(foreach)
library(doParallel)
library(tcltk)
library(survival)
library(raster)

# NULL PERMUTATIONS -------------------------------------------------------
#' @description shuffles binary matrices under different restrictions. 
#' @param mat A quantitative matrix
#' @param iter Number of random matrices to be created 
#' @param model Function to be chosen.
#' @param ... Further arguments from \code{permatswap} or \code{permatfull}
#' @return a list with \code{iter} random matrices
#' @details Totally restricted null model is called. Cell values are permuted restricting all features of the original matrix: column sums, row sums, matrix fill and total sum.
#' @references \code{citation("vegan")}

null <- function (mat, iter, ...){
  require(vegan)
  aux <- permatswap(mat, times=iter, method="quasiswap", fixedmar="both", shuffle="both", mtype="prab")
  return(aux$perm)
}

# VERTICAL LEARNING INDEX ------------------------------------------------------
#' @description This function creates an association matrix using the 1 for mother and calf pairs 0 for all other pairs. 
#' @param dat A data frame depicting IDs and their mothers
#' @return A square matrix in which each cell is an estimate of a dyadic vertical relationship

vert.func <- function(dat) {
  
  # Create an empty matrix with individual codes as row and column names
  matr <- matrix(0, nrow = nrow(dat), ncol = nrow(dat), 
                 dimnames = list(dat$Code, dat$Code))
  
  # Fill the matrix: 1 if individual i's mom is individual j
  for (i in 1:nrow(matr)) {
    for (j in 1:ncol(matr)) {
      matr[i, j] <- ifelse(dat$Mom[i] == dat$Code[j], 1, 0)
    }
  }
  
  #Replace any remaining NA values with 0
  matr[is.na(matr)] <- 0
  
  return(matr)
}


# SIMPLE-RATIO INDEX HORIZONTAL INDEX ------------------------------------------------------
#' @description This function creates an association matrix using the simple-ratio index (SRI). 
#' @param matr A binary matrix depicting individuals in the columns and groups in the rows
#' @return A square matrix in which each cell is an estimate of a dyadic social relationship, from 0 (never seen in the same group) to 1 (always seen in the same group)

SRI.func <-  function (matr) {
  if (any(is.na(matr))) {
    matr <- na.omit(matr)
    cat("The data matrix contains NA, and have been removed.\n")
  }
  matr1 = matr
  N <- nrow(matr1)
  matr1[matr1 > 1] <- 1
  n <- apply(matr1, 2, sum)
  tmatr <- t(matr1)
  df <- as.matrix(t(matr))
  a <- df %*% t(df) # Dyad in same group
  b <- df %*% (1 - t(df)) # A present, B absent
  c <- (1 - df) %*% t(df) # A absent, B present
  d <- ncol(df) - a - b - c # Double absent
  Dice <- data.frame()
  for (i in 1:nrow(a)) {
    for (j in 1:ncol(a)) {
      Dice[i, j] <- a[i, j]/(a[i, j] + b[i, j] + c[i, j])
    }
  }
  rownames(Dice)=colnames(Dice)=colnames(matr)
  Dice
}

# Break apart the numerator and demoninator of SRI
SRI_counts <- function(matr) {
  if (any(is.na(matr))) {
    matr <- na.omit(matr)
    cat("The data matrix contains NA, and have been removed.\n")
  }
  matr1 <- matr
  N <- nrow(matr1)
  matr1[matr1 > 1] <- 1
  n <- apply(matr1, 2, sum)
  tmatr <- t(matr1)
  df <- as.matrix(t(matr))
  
  a <- df %*% t(df)            # both present
  b <- df %*% (1 - t(df))      # A present, B absent
  c <- (1 - df) %*% t(df)      # A absent, B present
  d <- ncol(df) - a - b - c    # both absent (unused)
  
  denom <- a + b + c
  
  list(
    assoc = a,         # numerator
    opps  = denom,     # denominator
    SRI   = a / denom  # same as before, for checking
  )
}

create_counts <- function(gbi) {
  n.cores <- detectCores()
  registerDoParallel(n.cores)
  nxn <- list()
  for (i in seq_along(gbi)) {
    nxn[[i]] <- SRI_counts(gbi[[i]])
  }
  stopImplicitCluster()
  return(nxn)
}

# HRO ------------------------------------------------------
#' @description This function calculateS home range overlap using the formula: 
#' HRO = (Rij/Ri) * (Rij/Rj), where Rij is the overlap in range between individuals i and j, 
#' Ri is the range of individual i, and Rj is the range of individual j. 
#' @param kernel.poly home range polygons for each individual in a SpatialPolygonsDataFrame
#' @return A square matrix in which each cell is an estimate of dyadic homerange overlap

HRO <- function(dolph.kernel.poly) {
  
  # Empty matrix to store HRO values for each pair of individuals
  n_individuals <- length(dolph.kernel.poly@data$id)
  HRO_matrix <- matrix(NA, nrow = n_individuals, ncol = n_individuals, 
                       dimnames = list(dolph.kernel.poly@data$id, dolph.kernel.poly@data$id))
  
  # Loop through each pair of individuals to calculate HRO
  for (i in 1:n_individuals) {
    for (j in (i+1):n_individuals) {
      # Get the home range polygons for individuals i and j
      poly_i <- dolph.kernel.poly@polygons[[i]]
      poly_j <- dolph.kernel.poly@polygons[[j]]
      
      # Calculate the intersection area between the home ranges (Rij)
      intersection_area <- function(poly1, poly2) {
        if (is.null(gIntersects(poly1, poly2))) {
          return(0)  # No intersection, return 0
        } else {
          intersect_poly <- gIntersection(poly1, poly2)
          return(gArea(intersect_poly))
        }
      }
      
      intersect_area <- intersection_area(poly_i, poly_j)
      
      # Calculate the home range area for individuals i and j (Ri and Rj)
      home_range_area_i <- st_area(poly_i)
      home_range_area_j <- st_area(poly_j)
      
      # Calculate HRO using the formula: HRO = (Rij/Ri) * (Rij/Rj)
      HRO_value <- (intersect_area / home_range_area_i) * (intersect_area / home_range_area_j)
      
      # Store the HRO value in the matrix
      HRO_matrix[dolph.kernel.poly@data$id[i], dolph.kernel.poly@data$id[j]] <- HRO_value
      HRO_matrix[dolph.kernel.poly@data$id[j], dolph.kernel.poly@data$id[i]] <- HRO_value
    }
  }
  
}

