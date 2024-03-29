CPP_smooth.volume.FEM.basis<-function(locations, observations, mesh, lambda, covariates = NULL, incidence_matrix = NULL, ndim, mydim, BC = NULL, GCV,GCVmethod = 2, nrealizations = 100)
{
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  # This is done in C++ now to optimize speed
    
  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }
  
  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = ndim)
  }

  if(is.null(incidence_matrix))
  {
    incidence_matrix<-matrix(nrow = 0, ncol = 1)
  }
  
  if(is.null(BC$BC_indices))
  {
    BC$BC_indices<-vector(length=0)
  }else
  { 
    BC$BC_indices<-as.vector(BC$BC_indices)
  
  }
  
  if(is.null(BC$BC_values))
  {
    BC$BC_values<-vector(length=0)
  }else
  {
    BC$BC_values<-as.vector(BC$BC_values)
  }
  
  ## Set propr type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  data <- as.vector(observations)
  storage.mode(observations) <- "double"
  storage.mode(mesh$order) <- "integer"
  storage.mode(mesh$nnodes) <- "integer"
  storage.mode(mesh$ntetrahedrons) <- "integer"
  storage.mode(mesh$nodes) <- "double"
  storage.mode(mesh$tetrahedrons) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values)  <- "double"
  GCV <- as.integer(GCV)
  storage.mode(GCV) <-"integer"
  
  storage.mode(nrealizations) <- "integer"
  storage.mode(GCVmethod) <- "integer"
  
  ## Call C++ function
  bigsol <- .Call("regression_Laplace", locations, data, mesh, mesh$order, mydim, ndim, lambda, covariates,
                  incidence_matrix, BC$BC_indices, BC$BC_values, GCV, GCVmethod, nrealizations, package = "fdaPDE")

  return(bigsol)
}

CPP_eval.volume.FEM = function(FEM, locations, incidence_matrix, redundancy, ndim, mydim)
{
  FEMbasis = FEM$FEMbasis

  # Imposing types, this is necessary for correct reading from C++
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(FEMbasis$mesh$points) <- "double"
  storage.mode(FEMbasis$mesh$tetrahedrons) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  coeff <- as.matrix(FEM$coeff)
  storage.mode(coeff) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(locations) <- "double"
  storage.mode(redundancy) <- "integer"
  
  #Calling the C++ function "eval_FEM_fd" in RPDE_interface.cpp
  evalmat = matrix(0,nrow(locations),ncol(coeff))
  for (i in 1:ncol(coeff)){
    evalmat[,i] <- .Call("eval_FEM_fd", FEMbasis$mesh, locations, incidence_matrix, coeff[,i],
                          FEMbasis$order, redundancy, mydim, ndim, package = "fdaPDE")
  }

  #Returning the evaluation matrix
  evalmat
}
