rm(list = ls(all = TRUE))
library(xtable)
source('/Users/xiaofei/Box Sync/Fuller/Intermidiate design/master/code/sampling_procedure_4.R')

setwd("/Users/xiaofei/Box Sync/Fuller/Intermidiate design/Article/Code")

data <- read.csv("Payrolldata")
data <- data[order(data$EMP_Admin,decreasing = TRUE),]
y <- data$Payroll/10000


N <- nrow(data)
n <- 12
prob <- data$EMP_Admin/sum(data$EMP_Admin)*n


set.seed(1)
sp <- splitp2(prob,bd=1)
pi <- sp[[1]]
Stra_id <- sp[[2]]
Pair_id <- as.integer(ceiling(Stra_id / 2)) 
Set_id <- as.integer(Stra_id)
Group_id <- as.integer(ceiling(Pair_id/3))

TABOUT <- data.frame(ID = 1:nrow(data), Payroll = round(y,digits=2), 
                     EMP = data$EMP_Admin, 
                     Prob = round(prob, digits = 2), 
                     Cond_Prob = round(pi, digits = 2), 
                     Cumu_Prob = round(cumsum(pi), digits = 2), 
                     Set = Set_id, Pair = Pair_id, Group = Group_id)
TABOUT 

print(xtable(TABOUT),include.rownames = FALSE)

# S2
vartwo <- function(prob,y) {
  v <- 0
  N <- length(y)
  for(i in 1:N) 
    for(j in 1:N) {
      p2 <- pij(prob,c(i,j))
      v <- v + (4*prob[i]*prob[j]-p2)*(y[i]/prob[i]-y[j]/prob[j])^2/8
    }
  v
}

S2 <- numeric()
for(i in unique(Stra_id)) {
  S2[i] <- vartwo(prob[Stra_id == i], y[Stra_id == i])
}
S2

# Within mean square
varone <- function(prob,y) {
  v <- 0
  N <- length(y)
  for(i in 1:N) 
    for(j in 1:N) {
      v <- v + prob[i]*prob[j]*(y[i]/prob[i]-y[j]/prob[j])^2
    }
  v/2
}

WMS <- numeric()
for(i in unique(Stra_id)) {
  WMS[i] <- varone(prob[Stra_id == i], y[Stra_id == i])
}
WMS

# Between mean square

Tset <- tapply(y,Set_id,sum)
BMS <- (Tset[(1:12)%%2==1]-Tset[(1:12)%%2==0])^2


# One sample
set.seed(1)
sam <- inter_select(pi,n,2)

sam$Stru_info[,2] <- round(sam$Stru_info[,2] ,2)
for(x in c(1,3:6)) sam$Stru_info[,x] <- as.integer(sam$Stru_info[,x])

print(xtable(sam$Stru_info),include.rownames = FALSE)

print(xtable(rbind(WMS,BMS)),include.rownames = FALSE)

pij <- function(prob,s) {
  sam1 <- s[1]
  sam2 <- s[2]
  4*prob[sam1]*prob[sam2]/(2*(1+sum(prob/(1-2*prob))))*
    (1/(1-2*prob[sam1])+1/(1-2*prob[sam2]))
}

set6 <- sam$Stru_info[sam$Stru_info[,6]==6,1]
prod(pi[set6])/pij(pi[Set_id==6],set6 - (1:50)[Set_id==6][1] + 1) * 
(y[set6[1]]/prob[set6[1]]-y[set6[2]]/prob[set6[2]])^2

set9 <- sam$Stru_info[sam$Stru_info[,6]==9,1]
prod(pi[set9])/pij(pi[Set_id==9],set9 - (1:50)[Set_id==9][1] + 1) * 
  (y[set9[1]]/prob[set9[1]]-y[set9[2]]/prob[set9[2]])^2

var_info <- sam$Var_info
(y[var_info[,1]]/prob[var_info[,1]]-
    y[var_info[,2]]/prob[var_info[,2]])^2

# Theoretical variance
n=12
G=2
4*G*mean(S2) + (n-2*G)*mean(WMS) + 2*G*mean(BMS)

#Estimated variance
var_est1(y,prob,var_info)

# Replicate weights
weight <- c()
w <- 1/prob[sam$Sample]
for(i in 1:n) {
  temp <- w
  pair <- ceiling((1:n) / 2) == ceiling(i/2)
  temp[pair] <- 0
  temp[i] <- 2*w[i]
  weight <- cbind(weight,temp)
}
weight <- rbind(weight,rep(sam$Var_info[,3],each=2))
weight[,(1:12)%%2==1]

x <- pi-prob
z <- (pi-prob)*y/prob
s <- sam$Sample
X <- cbind(rep(1,n),z[s])
wd <- c(N,0) %*% solve(t(X) %*% diag(1/pi[s]) %*% X) %*% t(X) %*% diag(1/pi[s])

weight <- cbind(weight[,(1:12)%%2==1], c(wd+w,1))

row.names(weight) <- c(paste("w_",1:12,'^{(r)}',sep=''),"$c_r$")
colnames(weight) <- c(1:7)

print(xtable(weight))

























