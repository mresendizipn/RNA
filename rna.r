funcionesActivacion<-c("pureline","logsig","tansig")
funcionesActivacion<-c("purelineDer","logsigDer","tansigDer")


PropagacionAtras <- function(RNA,at,a0,eps){
 
  ##Calculate Matriz S
    namesS <- paste(rep("w",length(MLP$S)),as.character(1:(length(MLP$S))), sep = "")
    S <- vector("list", length(MLP$S))
    names(S) <- namesS
  
    for(i in length(RNA$S):1 ){
      
      if(i == length(RNA$S)){
        S[[i]] <- -2 * (diag(evalFunc(RNA$Nn[[i]],RNA$Sensibilidad[[i]][1])) %*%  (at-RNA$An[[i]]) )
      }else{
        S[[i]] <- diag(as.vector(evalFunc(RNA$Nn[[i]],RNA$Sensibilidad[[i]][1]))) %*% t(RNA$Wn[[i+1]]) %*% S[[i+1]]
      }
    }
    
  ##Update Wn and Bn
    for(i in 1:length(RNA$Wn)){
      if(i == 1){
        RNA$Wn[[i]] = RNA$Wn[[i]] - eps*(S[[i]] %*% t(a0))
        
      }else{
        RNA$Wn[[i]] = RNA$Wn[[i]] - eps*(S[[i]] %*% t(RNA$An[[i-1]])) 
      }
      RNA$Bn[[i]] <- RNA$Bn[[i]] - eps*(S[[i]])
    }
    RNA
}


createObjRNA <- function(MLP){
  ##Create list W
    namesWn <- paste(rep("w",length(MLP$S)),as.character(1:(length(MLP$S))), sep = "")
    Wn <- vector("list", length(MLP$S))
    names(Wn) <- namesWn
    
  ##Create list B
    namesBn <- paste(rep("b",length(MLP$S)),as.character(1:(length(MLP$S))), sep = "")
    Bn <- vector("list", length(MLP$S))
    names(Bn) <- namesBn
 
  ##Create list A
    namesAn <- paste(rep("a",length(MLP$S)),as.character(1:(length(MLP$S))), sep = "")
    An <- vector("list", length(MLP$S))
    names(An) <- namesAn
  
  ##Create list N
  namesNn <- paste(rep("n",length(MLP$S)),as.character(1:(length(MLP$S))), sep = "")
  Nn <- vector("list", length(MLP$S))
  names(Nn) <- namesNn

  #sensibilidad 
    namesSensibilidad <- paste(rep("S",length(MLP$S)),as.character(1:(length(MLP$S))), sep = "")
    Sensibilidad <- vector("list", length(MLP$S))
    names(Sensibilidad) <- namesSensibilidad
  
  ## Load randon values
  rnaConf<-c(MLP$R,MLP$S)
  for(i in 1:length(MLP$S) ){
    Wn[[i]] <- matrix(rnorm(rnaConf[i]*rnaConf[i+1]),nrow = MLP$S[i])
    Bn[[i]] <- matrix((rep(1,rnaConf[i+1])),ncol = 1)
    Sensibilidad[[i]] <- createMatrixSensitive(MLP$S[i],MLP$funAct[i])
  }
  #create list with vector of RNA
  list(Wn=Wn, Bn = Bn,An=An,Nn= Nn, funAct = MLP$funAct,Sensibilidad=Sensibilidad ) 
}

createMatrixSensitive<- function(num,strFun){
  auxTep <- matrix(rep("0",num*num),ncol = num) 
  for(i in 1:num){
    auxTep[i,i]<-paste(strFun,"Der",sep = "")
  }
  auxTep
}

propagationAdelante <- function(RNA,a0){
    RNA$Nn[[1]]<- RNA$Wn[[1]]%*%a0 + RNA$Bn[[1]]
    RNA$An[[1]] <- evalFunc(RNA$Nn[[1]],RNA$funAct[1])
  
  for(i in 2:length(RNA$Wn) ){
    tmp <- RNA$Wn[[i]]%*%RNA$An[[i-1]] + RNA$Bn[[i]]     
    RNA$Nn[[i]]  <- tmp
    RNA$An[[i]]  <- evalFunc(tmp,RNA$funAct[i])
  }  
  RNA
}

evalFunc <- function(x,strFun){
  row<-dim(x)[1]
  col<-dim(x)[2]
  tmp <- matrix(rep(0,col*row),ncol = col, nrow = row)
  for(r in 1:row){
    for(c in 1:col){
      auxText <- paste(strFun,"(",as.character(x[r,c]),")",sep = "")      
      tmp[r,c]<-eval(parse(text = auxText))
    }
  }
  tmp
}

pureline <- function(x){
  x
}

logsig <- function(x){
  1 / (1 + exp(-x))
}

tansig <- function(x){
  2/(1+exp(-2*x))-1
}

purelineDer <- function(x){
  1
}

logsigDer <- function(x){
  logsig(x)*(1-logsig(x))
}

tansigDer <- function(x){
  1 - tansig(x)^2
}