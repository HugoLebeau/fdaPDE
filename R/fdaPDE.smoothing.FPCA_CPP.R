CPP_smooth.FEM.FPCA<-function(locations, datamatrix, FEMbasis, incidence_matrix, lambda, ndim, mydim,nPC,validation, NFolds,GCVmethod = 2, nrealizations = 100)
{ 
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  ##TO BE CHANGED SOON: LOW PERFORMANCES, IMPLIES COPY OF PARAMETERS
  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  
  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = ndim)
  }

  if(is.null(incidence_matrix))
  {
    incidence_matrix<-matrix(nrow = 0, ncol = 1)
  }
  
  if(is.null(validation)) 
  {
    validation="NoValidation"
  }

  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(lambda )<- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  nPC <- as.integer(nPC)
  storage.mode(nPC) <- "integer"
  validation <- as.character(validation)
  storage.mode(validation) <- "character"
  NFolds <- as.integer(NFolds)
  storage.mode(NFolds) <- "integer"
  storage.mode(nrealizations) <- "integer"
  storage.mode(GCVmethod) <- "integer"
  
  ## Call C++ function
  bigsol <- .Call("Smooth_FPCA", locations, datamatrix, FEMbasis$mesh, 
                  FEMbasis$order, incidence_matrix, mydim, ndim, 
                  lambda, nPC, validation, NFolds, GCVmethod, nrealizations,
                  package = "fdaPDE")
  
  ## Reset them correctly
  #fdobj$basis$params$mesh$triangles = fdobj$basis$params$mesh$triangles + 1
  #fdobj$basis$params$mesh$edges = fdobj$basis$params$mesh$edges + 1
  #fdobj$basis$params$mesh$neighbors = fdobj$basis$params$mesh$neighbors + 1
  return(bigsol)
}



CPP_smooth.manifold.FEM.FPCA<-function(locations, datamatrix, mesh, incidence_matrix, lambda, ndim, mydim, nPC, validation, NFolds, GCVmethod = 2, nrealizations = 100)
{
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  # This is done in C++ now to optimize speed

  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = ndim)
  }
  
  if(is.null(incidence_matrix))
  {
    incidence_matrix<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(validation))
  {
    validation="NoValidation"
  }
  
  ## Set propr type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  datamatrix <- as.matrix(datamatrix)
  storage.mode(datamatrix) <- "double"
  storage.mode(mesh$order) <- "integer"
  storage.mode(mesh$nnodes) <- "integer"
  storage.mode(mesh$ntriangles) <- "integer"
  storage.mode(mesh$nodes) <- "double"
  storage.mode(mesh$triangles) <- "integer"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  nPC <- as.integer(nPC)
  storage.mode(nPC) <- "integer"
  validation <- as.character(validation)
  storage.mode(validation) <- "character"
  NFolds <- as.integer(NFolds)
  storage.mode(NFolds) <- "integer"
  storage.mode(nrealizations) <- "integer"
  storage.mode(GCVmethod) <- "integer"
  
  ## Call C++ function
  bigsol <- .Call("Smooth_FPCA", locations, datamatrix, mesh, 
                  mesh$order, incidence_matrix, mydim, ndim, lambda,
                  nPC, validation, NFolds, GCVmethod, nrealizations,
                  package = "fdaPDE")

  return(bigsol)
}

CPP_smooth.volume.FEM.FPCA<-function(locations, datamatrix, mesh, incidence_matrix, lambda, ndim, mydim, nPC, validation, NFolds, GCVmethod = 2, nrealizations = 100)
{
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  # This is done in C++ now to optimize speed

  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = ndim)
  }

  if(is.null(incidence_matrix))
  {
    incidence_matrix<-matrix(nrow = 0, ncol = 1)
  }
  
  if(is.null(validation))
  {
    validation="NoValidation"
  }
  
  ## Set propr type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  datamatrix <- as.matrix(datamatrix)
  storage.mode(datamatrix) <- "double"
  storage.mode(mesh$order) <- "integer"
  storage.mode(mesh$nnodes) <- "integer"
  storage.mode(mesh$ntetrahedrons) <- "integer"
  storage.mode(mesh$nodes) <- "double"
  storage.mode(mesh$tetrahedrons) <- "integer"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  nPC <- as.integer(nPC)
  storage.mode(nPC) <- "integer"
  validation <- as.character(validation)
  storage.mode(validation) <- "character"
  NFolds <- as.integer(NFolds)
  storage.mode(NFolds) <- "integer"
  storage.mode(nrealizations) <- "integer"
  storage.mode(GCVmethod) <- "integer"
  
  ## Call C++ function
  bigsol <- .Call("Smooth_FPCA", locations, datamatrix, mesh, 
                  mesh$order, incidence_matrix, mydim, ndim, lambda,
                  nPC, validation, NFolds, GCVmethod, nrealizations,
                  package = "fdaPDE")

  return(bigsol)
}

