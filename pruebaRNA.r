source('~/.active-rstudio-document')


MLP<-list(R=2,S=c(3,4,1), funAct = c("logsig","pureline","tansig"))

RNA <- createObjRNA(MLP)

a0<-matrix(c(3,7),ncol = 1)
at <- matrix(c(2),ncol = 1)
eps<-0.8

RNA <- propagationAdelante(RNA,a0) 
RNA <- PropagacionAtras(RNA,at,a0,eps)

data<-matrix(c(3,1,5,2,2,4),ncol = 2,byrow =T) 
target <- matrix(c(1,1,0),ncol = 1)