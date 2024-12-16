##################################
# RSCR_preamble ##################
# Joshua Cape (jrcape@wisc.edu) ##
##################################

# Load basic packages 
library(Matrix)
library(ggplot2)
library(matrixcalc)
library(RSpectra)
library(reshape2)
library(mclust)
library(ggbreak)
library(plotly)

# For simulations and ARI
library(igraph)
library(matrixStats)
#library(DescTools)

# For connectome data
library(devtools)
#devtools::install_github("youngser/TwoTruth", force=TRUE)
#library(TwoTruth)
library(igraph)
library(tidyverse)
library(sClust)
library(irlba)

# ---------------------------------------------------------------------------------------------------
# Creates a ptr matrix of a symmetric input matrix, w/o considering diagonal, using upper triangular entries.
# NB: Here we scale by (sample size plus one), motivated by literature on theory of rank tests.
# NB: For ties, return the average of rank values.
  ptr_func <- function(input.mtx){
    output.mtx <- matrix(0, nrow=dim(input.mtx)[1], ncol=dim(input.mtx)[2])
    output.mtx[upper.tri(output.mtx,
                         diag=FALSE)] <- rank(input.mtx[upper.tri(input.mtx,
                                                                  diag=FALSE)],
                                              ties.method="average")/(choose(dim(input.mtx)[1],2)+1)
    return(output.mtx + t(output.mtx))
  }
# ---------------------------------------------------------------------------------------------------
  
  
# ---------------------------------------------------------------------------------------------------
# PTR function for symmetric matrices that also ranks the main diagonal elements
  ptr_func_diag <- function(input.mtx){
    output.mtx <- matrix(0, nrow=dim(input.mtx)[1], ncol=dim(input.mtx)[2])
    output.mtx[upper.tri(output.mtx,
                         diag=TRUE)] <- rank(input.mtx[upper.tri(input.mtx,
                                                                  diag=TRUE)],
                                              ties.method="average")/(choose(dim(input.mtx)[1],2) + dim(input.mtx)[1] + 1)
    return(output.mtx + t(output.mtx) - diag(diag(output.mtx)))
  }
# ---------------------------------------------------------------------------------------------------
  
  
# ---------------------------------------------------------------------------------------------------
# PTR function for entire rectangular matrix
  ptr_func_rect <- function(input.mtx){
    input.mtx <- as.matrix(input.mtx)
    return(matrix(rank(input.mtx, ties.method="average")/(dim(input.mtx)[1]*dim(input.mtx)[2] + 1), ncol=dim(input.mtx)[2]))
  }
# ---------------------------------------------------------------------------------------------------

  
# ---------------------------------------------------------------------------------------------------
# PTR applied column-wise to an input matrix
  ptr_func_cols <- function(input.mtx){
    apply(input.mtx, 2, ptr_func_rect)
  }
# ---------------------------------------------------------------------------------------------------
    
  
# ---------------------------------------------------------------------------------------------------
# Function to compute symmetric normalized Laplacian
  norm_Lap <- function(input.mtx){
    deg.vec <- rowSums(input.mtx)
    return(diag(deg.vec^(-1/2)) %*% input.mtx %*% diag(deg.vec^(-1/2)))
  }
# ---------------------------------------------------------------------------------------------------

  
# ---------------------------------------------------------------------------------------------------
# Function converts eigen lists into matrix
  my_eigen_mtx <- function(eigen.list){
    return(eigen.list$vectors %*% diag(eigen.list$values, nrow = length(eigen.list$values)) %*% t(eigen.list$vectors))}
# ---------------------------------------------------------------------------------------------------

  
# ---------------------------------------------------------------------------------------------------  
# Function to align unit norm vectors with equiangular vector sign
quick_sign <- function(vec, equi.vec)
  {return(as.vector(vec * sign(vec * equi.vec)))}
# ---------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------
# Function to compute angle distances
vec_dists <- function(vec1, vec2){
  temp <- c(norm(as.matrix(vec1-vec2), "M"),
            norm(as.matrix(vec1-vec2), "F"),
            sqrt(1-(vec1 %*% vec2)^2))
  names(temp) <- c("max","ell2","sinMax")
  return(temp)}
# ---------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------
# Colors function
ggplotColours <- function(n = 6, h = c(0, 360) + 15) {
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
# ---------------------------------------------------------------------------------------------------