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

  ## Load randon values
  rnaConf<-c(MLP$R,MLP$S)
  for(i in 1:length(MLP$S) ){
    Wn[[i]] <- matrix(rnorm(rnaConf[i]*rnaConf[i+1]),nrow = MLP$S[i])
    Bn[[i]] <- matrix((rep(1,rnaConf[i+1])),ncol = 1)
  }
  #create list with vector of RNA
  list(Wn=Wn, Bn = Bn,An=An, funAct = MLP$funAct) 
}

propagationAdelante <- function(RNA,a0){
  ##Propagacion -->
  RNA$An[[1]] <- evalFunc(RNA$Wn[[1]]%*%a0 + RNA$Bn[[1]],RNA$funAct[length(RNA$Wn)])
  for(i in 2:length(RNA$Wn) ){
    tmp <- RNA$Wn[[i]]%*%RNA$An[[i-1]] + RNA$Bn[[i]]       
    RNA$An[[i]]  <- evalFunc(tmp,RNA$funAct[1-length(RNA$Wn)])
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



