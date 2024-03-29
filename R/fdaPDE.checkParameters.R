checkSmoothingParameters<-function(locations = NULL, observations, FEMbasis, lambda, covariates = NULL, incidence_matrix = NULL, BC = NULL, GCV = FALSE, CPP_CODE = TRUE, PDE_parameters_constant = NULL, PDE_parameters_func = NULL,GCVmethod = 2,nrealizations = 100)
{
  #################### Parameter Check #########################
  if(!is.null(locations))
  {
    if(any(is.na(locations)))
      stop("Missing values not admitted in 'locations'.")
    if(any(is.na(observations)))
      stop("Missing values not admitted in 'observations' when 'locations' are specified.")
  }
  
  if (is.null(observations))
    stop("observations required;  is NULL.")
  
  if (is.null(FEMbasis))
    stop("FEMbasis required;  is NULL.")
  if(class(FEMbasis)!= "FEMbasis")
    stop("'FEMbasis' is not class 'FEMbasis'")
  
  if (is.null(lambda))
    stop("lambda required;  is NULL.")
  
  if (!is.null(locations) && !is.null(incidence_matrix))
    stop("Both 'locations' and 'incidence_matrix' are given. In case of pointwise data, set 'incidence_matrix to NULL. In case of areal data, set 'locations' to NULL.")

  if (any(incidence_matrix!=0 & incidence_matrix!=1))
    stop("Value different than 0 or 1 in 'incidence_matrix'.")

  if(!is.null(BC))
  {
    if (is.null(BC$BC_indices)) 
      stop("'BC_indices' required in BC;  is NULL.")
    if (is.null(BC$BC_values)) 
      stop("'BC_indices' required in BC;  is NULL.")
  }
  
  if (is.null(GCV)) 
    stop("GCV required;  is NULL.")
  if(!is.logical(GCV))
    stop("'GCV' is not logical")
  
  if (is.null(CPP_CODE)) 
    stop("CPP_CODE required;  is NULL.")
  if(!is.logical(CPP_CODE))
    stop("'CPP_CODE' is not logical")
  
  if(!is.null(PDE_parameters_constant))
  {
    if (is.null(PDE_parameters_constant$K)) 
      stop("'K' required in PDE_parameters;  is NULL.")
    if (is.null(PDE_parameters_constant$b)) 
      stop("'b' required in PDE_parameters;  is NULL.")
    if (is.null(PDE_parameters_constant$c)) 
      stop("'c' required in PDE_parameters;  is NULL.")
  }
  
  if(!is.null(PDE_parameters_func))
  {
    if (is.null(PDE_parameters_func$K)) 
      stop("'K' required in PDE_parameters;  is NULL.")
    if(!is.function(PDE_parameters_func$K))
      stop("'K' in 'PDE_parameters' is not a function")
    if (is.null(PDE_parameters_func$b)) 
      stop("'b' required in PDE_parameters;  is NULL.")
    if(!is.function(PDE_parameters_func$b))
      stop("'b' in 'PDE_parameters' is not a function")
    if (is.null(PDE_parameters_func$c)) 
      stop("'c' required in PDE_parameters;  is NULL.")
    if(!is.function(PDE_parameters_func$c))
      stop("'c' in 'PDE_parameters' is not a function")
    if (is.null(PDE_parameters_func$u)) 
      stop("'u' required in PDE_parameters;  is NULL.")
    if(!is.function(PDE_parameters_func$u))
      stop("'u' in 'PDE_parameters' is not a function")
  }
  
  if (GCVmethod != 1 && GCVmethod != 2)
    stop("GCVmethod must be either 1(exact calculation) or 2(stochastic estimation)")

  if( !is.numeric(nrealizations) || nrealizations < 1)
    stop("nrealizations must be a positive integer")
}

checkSmoothingParametersSize<-function(locations = NULL, observations, FEMbasis, lambda, covariates = NULL, incidence_matrix = NULL, BC = NULL, GCV = FALSE, CPP_CODE = TRUE, PDE_parameters_constant = NULL, PDE_parameters_func = NULL, ndim, mydim)
{
  #################### Parameter Check #########################
  if(ncol(observations) != 1)
    stop("'observations' must be a column vector")
  if(nrow(observations) < 1)
    stop("'observations' must contain at least one element")
  if(is.null(locations))
  {
    if(class(FEMbasis$mesh) == "MESH2D"){
    	if(nrow(observations) > nrow(FEMbasis$mesh$nodes))
     	 stop("Size of 'observations' is larger then the size of 'nodes' in the mesh")
    }else if(class(FEMbasis$mesh) == "MESH.2.5D" || class(FEMbasis$mesh) == "MESH.3D"){
    	if(nrow(observations) > FEMbasis$mesh$nnodes)
     	 stop("Size of 'observations' is larger then the size of 'nodes' in the mesh")
    }
  }
  if(!is.null(locations))
  {
    if(ncol(locations) != ndim)
      stop("'locations' must be a ndim-columns matrix;")
    if(nrow(locations) != nrow(observations))
      stop("'locations' and 'observations' have incompatible size;")
  }
  if(ncol(lambda) != 1)
    stop("'lambda' must be a column vector")
  if(nrow(lambda) < 1)
    stop("'lambda' must contain at least one element")
  if(!is.null(covariates))
  {
    if(nrow(covariates) != nrow(observations))
      stop("'covariates' and 'observations' have incompatible size;")
  }
  
  if (!is.null(incidence_matrix))
  {
    if (nrow(incidence_matrix) != nrow(observations))
      stop("'incidence_matrix' and 'observations' have incompatible size;")
    if (class(FEMbasis$mesh) == 'MESH2D' && ncol(incidence_matrix) != nrow(FEMbasis$mesh$triangles))
      stop("'incidence_matrix' must be a ntriangles-columns matrix;")
    else if (class(FEMbasis$mesh) == 'MESH.2.5D' && ncol(incidence_matrix) != FEMbasis$mesh$ntriangles)
      stop("'incidence_matrix' must be a ntriangles-columns matrix;")
    else if (class(FEMbasis$mesh) == 'MESH.3D' && ncol(incidence_matrix) != FEMbasis$mesh$ntetrahedrons)
      stop("'incidence_matrix' must be a ntetrahedrons-columns matrix;") 
  }

  if(!is.null(BC))
  {
    if(ncol(BC$BC_indices) != 1)
      stop("'BC_indices' must be a column vector")
    if(ncol(BC$BC_values) != 1)
      stop("'BC_values' must be a column vector")
    if(nrow(BC$BC_indices) != nrow(BC$BC_values))
      stop("'BC_indices' and 'BC_values' have incompatible size;")
     if(class(FEMbasis$mesh) == "MESH2D"){
	    if(sum(BC$BC_indices>nrow(nrow(FEMbasis$mesh$nodes))) > 0)
	      stop("At least one index in 'BC_indices' larger then the number of 'nodes' in the mesh;")
    }else if((class(FEMbasis$mesh) == "MESH.2.5D" || class(FEMbasis$mesh) == "MESH.3D")){
    	if(sum(BC$BC_indices>FEMbasis$mesh$nnodes) > 0)
	      stop("At least one index in 'BC_indices' larger then the number of 'nodes' in the mesh;")
    	}
  }
  
  if(!is.null(PDE_parameters_constant))
  {
    if(!all.equal(dim(PDE_parameters_constant$K), c(2,2)))
      stop("'K' in 'PDE_parameters must be a 2x2 matrix")
    if(!all.equal(dim(PDE_parameters_constant$b), c(2,1)))
      stop("'b' in 'PDE_parameters must be a column vector of size 2")
    if(!all.equal(dim(PDE_parameters_constant$c), c(1,1)))
      stop("'c' in 'PDE_parameters must be a double")
  }
  
  if(!is.null(PDE_parameters_func))
  {
    
    n_test_points = min(nrow(FEMbasis$mesh$nodes), 5)
    test_points = FEMbasis$mesh$nodes[1:n_test_points, ]
    
    try_K_func = PDE_parameters_func$K(test_points)
    try_b_func = PDE_parameters_func$b(test_points)
    try_c_func = PDE_parameters_func$c(test_points)
    try_u_func = PDE_parameters_func$u(test_points)
    
    if(!is.numeric(try_K_func))
      stop("Test on function 'K' in 'PDE_parameters' not passed; output is not numeric")
    if(!all.equal(dim(try_K_func), c(2,2,n_test_points)) )
      stop("Test on function 'K' in 'PDE_parameters' not passed; wrong size of the output")
    
    if(!is.numeric(try_b_func))
      stop("Test on function 'b' in 'PDE_parameters' not passed; output is not numeric")
    if(!all.equal(dim(try_b_func), c(2,n_test_points)))
      stop("Test on function 'b' in 'PDE_parameters' not passed; wrong size of the output")
    
    if(!is.numeric(try_c_func))
      stop("Test on function 'c' in 'PDE_parameters' not passed; output is not numeric")
    if(length(try_c_func) != n_test_points)
      stop("Test on function 'c' in 'PDE_parameters' not passed; wrong size of the output")
    
    if(!is.numeric(try_u_func))
      stop("Test on function 'u' in 'PDE_parameters' not passed; output is not numeric")
    if(length(try_u_func) != n_test_points)
      stop("Test on function 'u' in 'PDE_parameters' not passed; wrong size of the output")
  }
  
}
