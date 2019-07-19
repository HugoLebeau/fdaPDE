rm(list=objects())
graphics.off()

library(fdaPDE)

dimension<-"2.5"
path<-"~/Documents/ENSTA Paris 2A/PRe/Projet/fdaPDE/"

## WORKING DIRECTORY
setwd(paste(path,"RScripts/",dimension,"D Mesh/",sep=''))

## READ MESH
nodes<-as.matrix(read.csv("node.csv",sep=';',header=FALSE))
elems<-as.matrix(read.csv("elem.csv",sep=';',header=FALSE))
refElems<-as.matrix(read.csv("refElem.csv",sep=';',header=FALSE))

## DEFINE REGIONS
incidence_matrix<-matrix(0,nrow=3,ncol=nrow(elems))
incidence_matrix[1,refElems==2]<-1
incidence_matrix[2,refElems==3]<-1
incidence_matrix[3,refElems==4]<-1

#Z
observations<-c(5,10,15)

## W
#cov<-NULL
cov<-c(9,2,3)
#cov<-cbind(c(9,2,3),c(5,2,4))
#cov<-observations


if (dimension=="2"){
  mesh<-create.MESH.2D(nodes=nodes,triangles=elems)
}else if (dimension=="2.5"){
  mesh<-create.MESH.2.5D(nodes=nodes,triangles=elems)
}else if (dimension=="3"){
  mesh<-create.MESH.3D(nodes=nodes,tetrahedrons=elems)
}

FEMbasis<-create.FEM.basis(mesh)

smooth<-smooth.FEM.basis(observations=observations,FEMbasis=FEMbasis,lambda=1,
                         covariates=cov,incidence_matrix=incidence_matrix,GCV=FALSE)


head(smooth$fit.FEM$coeff)
min(smooth$fit.FEM$coeff)
max(smooth$fit.FEM$coeff)
observations-cov%*%smooth$beta

plot.FEM(smooth$fit.FEM)
