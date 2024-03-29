checkSmoothingParametersFPCA<-function(locations = NULL, datamatrix, FEMbasis, incidence_matrix = NULL, lambda, nPC, validation, NFolds, GCVmethod = 2, nrealizations = 100)
{
  #################### Parameter Check #########################
  if(!is.null(locations))
  {
    if(any(is.na(locations)))
      stop("Missing values not admitted in 'locations'.")
    if(any(is.na(datamatrix)))
      stop("Missing values not admitted in 'datamatrix' when 'locations' are specified.")
  }
  if (is.null(datamatrix)) 
    stop("observations required;  is NULL.")
  if (is.null(FEMbasis)) 
    stop("FEMbasis required;  is NULL.")
  if(class(FEMbasis)!= "FEMbasis")
    stop("'FEMbasis' is not class 'FEMbasis'")
  
  if (!is.null(locations) && !is.null(incidence_matrix))
    stop("Both 'locations' and 'incidence_matrix' are given. In case of pointwise data, set 'incidence_matrix to NULL. In case of areal data, set 'locations' to NULL.")

  if (any(incidence_matrix!=0 & incidence_matrix!=1))
    stop("Value different than 0 or 1 in 'incidence_matrix'.")

  if (is.null(lambda)) 
    stop("lambda required;  is NULL.")
  if(is.null(nPC))
    stop("nPC required; is NULL.")
  if(!is.null(validation)){
    if(validation!="GCV" && validation!="KFold")
   	stop("'validation' needs to be 'GCV' or 'KFold'")
    if(validation=="KFold" && is.null(NFolds))
   	stop("NFolds is required if 'validation' is 'KFold'")
   }
   if (GCVmethod != 1 && GCVmethod != 2)
    stop("GCVmethod must be either 1(exact calculation) or 2(stochastic estimation)")

  if( !is.numeric(nrealizations) || nrealizations < 1)
    stop("nrealizations must be a positive integer")

}

checkSmoothingParametersSizeFPCA<-function(locations = NULL, datamatrix, FEMbasis, incidence_matrix, lambda, ndim, mydim, validation, NFolds)
{
  #################### Parameter Check #########################
  if(nrow(datamatrix) < 1)
    stop("'datamatrix' must contain at least one element")
  if(is.null(locations))
  {
    if(class(FEMbasis$mesh) == "MESH2D"){
    	if(ncol(datamatrix) > nrow(FEMbasis$mesh$nodes))
     	 stop("Size of 'datamatrix' is larger then the size of 'nodes' in the mesh")
    }else if(class(FEMbasis$mesh) == "MESH.2.5D" || class(FEMbasis$mesh) == "MESH.3D"){
    	if(ncol(datamatrix) > FEMbasis$mesh$nnodes)
     	 stop("Size of 'datamatrix' is larger then the size of 'nodes' in the mesh")
    }
  }
  if(!is.null(locations))
  {
    if(ncol(locations) != ndim)
      stop("'locations' and the mesh points have incompatible size;")
    if(nrow(locations) != ncol(datamatrix))
      stop("'locations' and 'datamatrix' have incompatible size;")
  }
  if (!is.null(incidence_matrix))
  {
    if (nrow(incidence_matrix) != ncol(datamatrix))
      stop("'incidence_matrix' and 'datamatrix' have incompatible size;")
    if (class(FEMbasis$mesh) == 'MESH2D' && ncol(incidence_matrix) != nrow(FEMbasis$mesh$triangles))
      stop("'incidence_matrix' must be a ntriangles-columns matrix;")
    else if (class(FEMbasis$mesh) == 'MESH.2.5D' && ncol(incidence_matrix) != FEMbasis$mesh$ntriangles)
      stop("'incidence_matrix' must be a ntriangles-columns matrix;")
    else if (class(FEMbasis$mesh) == 'MESH.3D' && ncol(incidence_matrix) != FEMbasis$mesh$ntetrahedrons)
      stop("'incidence_matrix' must be a ntetrahedrons-columns matrix;") 
  }
  if(ncol(lambda) != 1)
    stop("'lambda' must be a column vector")
  if(nrow(lambda) < 1)
    stop("'lambda' must contain at least one element")
  if(nrow(lambda)>1 && is.null(validation))
    stop("If 'lambda' contains more than one element, 'validation' needs to be specified as 'GCV' or 'KFold'")
}
