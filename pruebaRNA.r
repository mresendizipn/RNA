MLP<-list(R=2,S=c(3,4,1), funAct = c("tansig","pureline","pureline"))
funcionesActivacion<-c("pureline","logsig","tansig")
funcionesActivacion<-c("purelineDer","logsigDer","tansigDer")

RNA <- createObjRNA(MLP)

a0<-matrix(c(3,7),ncol = 1)
propagationAdelante(RNA,a0) 

data<-matrix(c(3,1,5,2,2,4),ncol = 2,byrow =T) 
target <- matrix(c(1,1,0),ncol = 1)