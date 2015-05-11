sigmoide<-function(x,a=1){
  1/(1+exp(-a*x))
}

vec<-c(1,1,1)
########################
m11<-c(1,4,1)
m12<-c(3,2,1)
m21<-c(2,3,1)

#########################
v11<-sum(vec*m11)
v11
v12<-sum(vec*m12)
v12

y11 <- sigmoide(v11)
y11
y12 <- sigmoide(v12)
y12
y1 <- c(y11,y12,1)

v21<-sum(y1*m21)
v21

y21<-sigmoide(v21)
y21
########################


