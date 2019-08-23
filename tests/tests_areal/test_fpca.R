rm(list=objects())
graphics.off()

library(fdaPDE)

path<-"~/Documents/ENSTA Paris 2A/PRe/Projet/fdaPDE/"

## WORKING DIRECTORY
setwd(paste0(path,"tests/tests_areal/Sphere/"))

## READ MESH
node<-"Sfera_vertici.txt"
elem<-"Sfera_tetraedri.txt"
nodes<-read.table(node)[,1:3]
elems<-read.table(elem)[,1:4]
mesh<-create.MESH.3D(nodes,elems)
FEMbasis<-create.FEM.basis(mesh)

## SET PARAMETERS
nIter<-5
lambda<-1
validation<-NULL
nRegions<-500
sizeOfRegions<-c(rep(floor(nrow(elems)/nRegions),nRegions-1),floor(nrow(elems)/nRegions)+nrow(elems)%%nRegions)
seed<-1234
nSamples<-50
sd_score1<-0.1
sd_score2<-0.05
sd_score3<-0.0125
sd_error<-0.075
set.seed(seed)

## READ EIGENFUNCTIONS
eigenf1<-"Sfera_eigenf_1.RData"
eigenf2<-"Sfera_eigenf_2.RData"
eigenf3<-"Sfera_eigenf_3.RData"
load(paste0(eigenf1))
eigenfunc1<-FEM(coeff=sol_esatta,FEMbasis=FEMbasis)
load(paste0(eigenf2))
eigenfunc2<-FEM(coeff=sol_esatta,FEMbasis=FEMbasis)
load(paste0(eigenf3))
eigenfunc3<-FEM(coeff=sol_esatta,FEMbasis=FEMbasis)

truedatarange<-max(c(eigenfunc1$coeff,eigenfunc2$coeff,eigenfunc3$coeff))-min(c(eigenfunc1$coeff,eigenfunc2$coeff,eigenfunc3$coeff))
truecoeff<-cbind(eigenfunc1$coeff,eigenfunc2$coeff,eigenfunc3$coeff)

RMSE.pointwise<-rep(0,nIter)
RMSE.areal<-rep(0,nIter)
for (iter in 1:nIter)
{
  print(paste("Iteration",iter))
  ## REGIONS AND DATAMATRIX
  score1<-rnorm(n=nSamples,sd=sd_score1*truedatarange)
  score2<-rnorm(n=nSamples,sd=sd_score2*truedatarange)
  score3<-rnorm(n=nSamples,sd=sd_score3*truedatarange)
  error<-rnorm(n=nSamples*nrow(nodes),sd=sd_error*truedatarange)
  datamatrix.pointwise<-matrix(score1)%*%t(matrix(eigenfunc1$coeff))+matrix(score2)%*%t(matrix(eigenfunc2$coeff))+matrix(score3)%*%t(matrix(eigenfunc3$coeff))
  datamatrix.pointwise<-datamatrix.pointwise+error
  datamatrix.areal<-matrix(0,nrow=nSamples,ncol=nRegions)
  incidence_matrix<-matrix(0,nrow=nRegions,ncol=nrow(elems))
  regions<-list()
  ind<-sample.int(nrow(elems)) #random regions
  inf<-1
  for (i in 1:nRegions)
  {
    regions[[i]]<-elems[ind[inf:(inf+sizeOfRegions[i]-1)],]
    incidence_matrix[i,ind[inf:(inf+sizeOfRegions[i]-1)]]<-1
    inf<-inf+sizeOfRegions[i]
    # Fill in datamatrix.areal
    regionNodes<-c() #nodes in the region
    for (k in 1:ncol(regions[[i]])) {regionNodes<-c(regionNodes,regions[[i]][,k])}
    regionNodes<-regionNodes[!duplicated(regionNodes)] #remove duplicates
    datamatrix.areal[,i]<-apply(datamatrix.pointwise[,regionNodes],1,sum)/length(regionNodes)
  }
  dm.pointwise.centred<-datamatrix.pointwise-matrix(apply(datamatrix.pointwise,2,mean),ncol=ncol(datamatrix.pointwise),nrow=nrow(datamatrix.pointwise),byrow=TRUE)
  dm.areal.centred<-datamatrix.areal-matrix(apply(datamatrix.areal,2,mean),ncol=ncol(datamatrix.areal),nrow=nrow(datamatrix.areal),byrow=TRUE)
  
  sol.pointwise<-smooth.FEM.FPCA(locations=nodes,datamatrix=dm.pointwise.centred,FEMbasis=FEMbasis,lambda=lambda,nPC=3,validation=validation)
  sol.areal<-smooth.FEM.FPCA(datamatrix=dm.areal.centred,FEMbasis=FEMbasis,incidence_matrix=incidence_matrix,lambda=lambda,nPC=3,validation=validation)
  
  RMSE.pointwise[iter]<-sqrt(mean((sol.pointwise$loadings.FEM$coeff-truecoeff)^2))
  RMSE.areal[iter]<-sqrt(mean((sol.areal$loadings.FEM$coeff-truecoeff)^2))
}

mean(RMSE.pointwise)
mean(RMSE.areal)
