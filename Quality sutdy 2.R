rm(list = ls(all = TRUE))
source('/Users/xiaofei/Box Sync/Fuller/Intermidiate design/master/code/sampling_procedure_4.R')

folder <- '/Users/xiaofei/Box Sync/Fuller/Intermidiate design/Quality Study'
setwd(folder)
iowa_data <- read.csv('iowa_data.csv')
iowa_data <- iowa_data[order(iowa_data$Order),]
y1 <- iowa_data$latitude
y2 <- iowa_data$longitude
y3 <- iowa_data$SAMPLE_CLASS

G <- c(31,20,15,10,5,1)
n <- 124
N <- nrow(iowa_data)
prob <- rep(n/N, N)

m=10000
res <- list()
for(i in 1:length(G)) res[[i]] <- matrix(NA, nrow = m, ncol = 8)

s12 <- s12_c <- matrix(NA, nrow = m, ncol = 3)
ET <- rep(NA,m)
y = y2
set.seed(4321)
for(i in 1:m) {
  # set.seed(i)
  sam <- inter_small_simu(prob,G,y)
  pcd <- sam[[2]]
  for(g in 1:length(G)) {
    var_info <- sam[[1]][[g]]
    s <- unique(c(var_info$Sample,var_info$pair))
    
    y_HT <- sum(y[s]/prob[s])/N
    y_HT2 <- sum(y[s]/pcd[s])/sum(1/pcd[s])
    
    d2 <- y - y_HT2
    
    v_y2_c <- var_est1(d2,pcd,var_info)/N^2
    v_y_c <- var_est1(y,prob,var_info)/N^2
    
    x <- pcd-prob
    z <- (pcd-prob)*y/prob
    fit <- lm(z[s]~x[s], weights = 1/pcd[s])
    bias <- c(N,sum(x)) %*% fit$coef
    e <- z - fit$coef[1] - fit$coef[2]*x
    # Md <- bias^2/N^2 - var_est1(e,pcd,var_info)/N^2

    res[[g]][i,] <- c(y_HT,v_y_c,bias^2/N^2,var_est1(e,pcd,var_info)/N^2,y_HT2,v_y2_c,
                      var_modi(y,prob,var_info))
  }
  s12[i,] <- sam[[3]]
  ET[i] <- sum(y/prob*pcd)
  if(i %% 10 == 0) print(i)
}

temp <- res[[1]][,7]
res[[1]][,7] <- res[[1]][,8]
res[[1]][,8] <- temp

# for(g in 1:length(G)) write.csv(res[[g]],paste('result/Long1_G',G[g],sep=''))
# write.csv(cbind(s12,ET),'result/Long1_TH')
# 
# res <- list()
# for(g in 1:length(G)) res[[g]] <- read.csv(paste('result/Long1_G',G[g],sep=''))[,-1]
# 
# s12 <- read.csv('result/Long1_TH')[,-1]
# ET <- s12[,4]
# s12 <- s12[,1:3]
# 
# df <- sapply(1:length(G), function(x) 2*var(res[[x]][,1])^2/var(res[[x]][,2]+res[[x]][,3]-res[[x]][,4]))
# sapply(1:length(G), function(x) mean(res[[x]][,1]-qt(0.975, df[x])*sqrt(res[[x]][,2]+res[[x]][,3]) < mean(y3) &
#                                        res[[x]][,1]+qt(0.975, df[x])*sqrt(res[[x]][,2]+res[[x]][,3]) > mean(y3)))
# 
# sapply(1:length(G), function(x) mean(res[[x]][,7]))/n
# sapply(1:length(G), function(x) mean(res[[x]][,8]))/2/G
# 
# rf <- qf(0.95,G,61-G)
# var_adj <- function(x) ((res[[x]][,7]/n+rf[x]*res[[x]][,8]/2/G[x]-
#                            abs(res[[x]][,7]/n-rf[x]*res[[x]][,8]/2/G[x]))*n/2+res[[x]][,8])/N^2

res_inter <- rbind((sapply(G, function(G) sum(colMeans(s12)*c(n-2*G,4*G,2*G))) + var(ET))/N^2,
                   sapply(1:length(G), function(x) var(res[[x]][,1])),
                   sapply(1:length(G), function(x) mean(res[[x]][,2])),
                   sapply(1:length(G), function(x) mean(res[[x]][,2]+res[[x]][,3]-res[[x]][,4])),
                   # sapply(1:length(G), function(x) mean(var_adj(x))),
                   sqrt(sapply(1:length(G), function(x) var(res[[x]][,2]+res[[x]][,3]-res[[x]][,4]))),
                   # sapply(1:length(G), function(x) var(res[[x]][,5])),
                   # sapply(1:length(G), function(x) mean(res[[x]][,6])),
                   # sqrt(sapply(1:length(G), function(x) var(res[[x]][,6]))),
                   # sapply(1:length(G), function(x) mean(res[[x]][,1]-mean(y1))/sqrt(var(res[[x]][,1])))*100,
                   # sapply(1:length(G), function(x) sqrt(var(res[[x]][,2]+res[[x]][,3]-res[[x]][,4]))/var(res[[x]][,1])),
                   sapply(1:length(G), function(x) mean(res[[x]][,1]-1.96*sqrt(res[[x]][,2]+res[[x]][,3]) < mean(y) &
                                                          res[[x]][,1]+1.96*sqrt(res[[x]][,2]+res[[x]][,3]) > mean(y))))


res_inter*c(rep(10^2,5),1)

row.names(res_inter) <- c('Theoretical variance of mean',
                          'M.C. variance of mean',
                          'Estimated conditional variance',
                          'Estimated variance',
                          'Standard deviation of est. var.',
                          # 'M.C. variance for $\\bar{y}_{c,q}$',
                          # 'Estimated variance for $\\bar{y}_{c,q}$',
                          # 'S.d. of the estimated varianc for $\\bar{y}_{c,q}$',
                          # 'Bias ratio',
                          # 'CV',
                          'Coverage')
library(xtable)
xtable(res_inter*c(rep(10^2,5),1))

par(mfrow=c(2,3))
for(i in 1:length(G)) {
  hist(res[[i]][,1]-mean(y1),breaks=30, xlab = 'Residual', 
       ylim = c(0,26), freq = FALSE,
       main=paste('G = ', G[i], sep = ''))
  curve(dnorm(x,0,sqrt(var(res[[i]][,1]))),col='red',add=TRUE)
  # curve(dt(x/sqrt(var(res[[i]][,1])),G[i]/2)/sqrt(var(res[[i]][,1])),col='blue',add=TRUE)
}


par(mfrow=c(2,3))
for(i in 1:length(G)) {
  ubv <- res[[i]][,2]+res[[i]][,3]-res[[i]][,4]
  hist(ubv,breaks=15, xlab = 'Estimated variance',
       # ylim = c(0,1200),
       main=paste('G = ', G[i], sep = ''))
}



########## SRS ###########
SRS1 <-SRS2 <- SRS3 <-  matrix(NA, nrow = m, ncol = n)
for(i in 1:m) {
  s <- sample(N,n)
  SRS1[i,] <- y1[s]
  SRS2[i,] <- y2[s]
  SRS3[i,] <- y3[s]
}



