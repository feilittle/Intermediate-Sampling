
### general split procedure to split a sequence at boundary bd
assignset <- function(prob,bd=1) {
  N <- length(prob)
  s2 <- sum(prob)-bd
  set1 <- set2 <- numeric()
  cuprob <- cumsum(prob)
  
  ### find boundary element
  x <- which(cuprob>=bd)[1]
  
  if(cuprob[x]==bd) {
    bdry <- x 
    set1 <- prob[1:x]
    set2 <- prob[(x+1):N]
  } else {
    pbd <- prob[x]
    pbd2 <- cuprob[x]-bd
    pbd1 <- pbd-pbd2
    u <- runif(1,0,1)
    if(u <= pbd1/pbd) {
      ### boundary element is assigned to the left to the boundary
      bdry <- x
      set1 <- c(prob[(1:(bdry-1))]/(bd-pbd1)*(bd-pbd),pbd)
      set2 <- prob[((bdry+1):N)]/(s2-pbd2)*s2
    } else {
      ### boundary element is assigned to the right to the boundary
      bdry <- x-1
      set1 <- prob[(1:bdry)]/(bd-pbd1)*bd
      set2 <- c(pbd,prob[((bdry+2):N)]/(s2-pbd2)*(s2-pbd))
    }
  }
  return(list(set1=set1,set2=set2,bdry=bdry))
}

### Split a sequence of probabilities into sets with 
### sum of probabilities in each group equal to sset
splitp2 <- function(pi,sset=1) {
  
  n <- round(sum(pi)) # sample size 
  pi <- pi/sum(pi)*n
  
  N <- length(pi) # population size
  
  g <- n/sset # number of groups

  bdry <- c(0,rep(NA,g-1),N)  # boundaries of sets 
  
  cumpi <- cumsum(pi)
  
  ### find boundary elements
  B <- rep(NA,g-1)
  for(i in 1:(g-1)) B[i] <- which(cumpi >= i*sset)[1] 
  B <- c(B,N)                     ### boundary elements
  p1B <- (sset*(1:g) - cumpi[B-1])   ### probability to the left
  p2B <- (pi[B] - p1B)            ### probability to the right
  
  prob <- pi
  if(g==2) {
    a <- assignset(prob,bd=sset)
    new_prob <- c(a[[1]],a[[2]])
    bdry <- c(0,a[[3]],N)
  } else {
    for (i in 1:(g-1)) {
      temp <- c(prob[(bdry[i]+1):(B[i+1]-1)],p1B[i+1]) 
      a <- assignset(temp,bd=sset)
      bdry[i+1] <- bdry[i]+a[[3]]
      prob[(bdry[i]+1):B[i+1]] <- c(a[[1]], a[[2]])
      prob[B[i+1]] <- tail(a[[2]], n=1) + p2B[i+1]
    }
    new_prob <- prob
  }
  
  stra <- rep(1:g, diff(bdry))
  return(list(new_prob,stra))
}

split_large <- function(prej) {
  nrej <- round(sum(prej))
  Stra_id <- numeric()
  i <- 1
  a <- 1
  prej_ck <- prej
  nprej <- numeric()
  while(i < round(nrej-1)) {
    if(which(cumsum(prej_ck) >= 1)[1] > 3) {
      temp1 <- assignset(prej_ck,1)
      i = i+1
    } else {
      temp1 <- assignset(prej_ck,2)
      i = i+2
    }
    nprej <- c(nprej,temp1$set1) 
    prej_ck <- temp1$set2
    Stra_id <- c(Stra_id,rep(a,length(temp1$set1)))
    a <- a+1
  }
  nprej <- c(nprej,prej_ck) 
  Stra_id <- c(Stra_id,rep(a,length(prej_ck)))
  return(list(prej=nprej,Stra_id=Stra_id))
}

checksplit <- function(prob,bd=1) {
  N <- length(prob)
  s2 <- sum(prob)-bd
  cuprob <- cumsum(prob)
  x <- which(cuprob>=bd)[1]
  if(x < 3 | (N-x) <= 2) {
    return(check=1)
  } else {
    pbd <- prob[x]
    pbd2 <- cuprob[x]-bd
    pbd1 <- pbd-pbd2
    check <- sum(prob[(1:(x-1))]/(bd-pbd1)*(bd-pbd) > 1) + sum(prob[((x+1):N)]/(s2-pbd2)*s2 > 1) +
      sum(prob[(1:(x-1))]/(bd-pbd1)*bd > 1) + sum(prob[((x+1):N)]/(s2-pbd2)*(s2-pbd) > 1) 
    return(check=check)
  }
}

##### sample selection #####

# One-per-stratum
draw_1 <- function(prob,x) {
  # sum(prob) =1
  u <- runif(1,0,1)
  sam <- which(u < cumsum(prob))[1]
  return(x[sam])
}

# Two-per-stratum
# Brewer's procedure
draw_2 <- function(prob,x) {
  prob <- prob/sum(prob)
  if(length(prob)==2)  return(c(x,1,prob))  else {
    q1 <- prob*(1-prob)/(1-2*prob)/sum(prob*(1-prob)/(1-2*prob))
    sam1 <- which(runif(1,0,1) < cumsum(q1))[1]
    prob2 <- prob
    prob2[sam1] <- 0
    q2 <- prob2/sum(prob2)
    sam2 <- which(runif(1,0,1) < cumsum(q2))[1]
    pij <- 4*prob[sam1]*prob[sam2]/(2*(1+sum(prob/(1-2*prob))))*
      (1/(1-2*prob[sam1])+1/(1-2*prob[sam2]))
    # pij is the joint inclusion probability, 2*prob[sam] is the inclusion prob
    return(c(s1=x[sam1],s2=x[sam2],pij = pij,p1=prob[sam1],p2=prob[sam2]))
  }
}

# Draw large
rejct_1 <- function(prob) {
  N <- length(prob)
  rej <- draw_1(prob,1:N)
}

##### Intermediate design #####

### Sampling procedure for intermediate sampling
intersamp3 <- function(pi,Set_id,G,Idd1,aa) {
  n <- round(sum(pi))
  N <- length(pi)
  Id <- 1:N
  
  Pair_id <- ceiling(Set_id/2)
  
  if(G==1) {
    sg <- list(prob = rep(G/(n/2),n/2), Set_id = rep(1,n/2))
  } else {
    sg <- splitp2(rep(G/(n/2),n/2))
  }
  pig <- sg[[1]]
  K <- tapply(rep(1,n/2),sg[[2]],sum)
  Group_id <- findInterval(Set_id, cumsum(2*K), left.open = TRUE) + 1
  
  two_per <- one_per <- sam_pair <- numeric()
  for(i in 1:G) {
    sam_pair <- c(sam_pair,draw_1(pig[sg[[2]]==i],unique(Pair_id[Group_id==i])))
    two_per <- c(two_per,sample(unique(Set_id[Pair_id==sam_pair[i]]),1))
    one_per <- c(one_per,unique(Set_id[Pair_id!=sam_pair[i] & Group_id==i]))
  }
  
  Sam_two_per <- sapply(two_per,function(x) draw_2(pi[Set_id==x],Id[Set_id==x]))
  Sam_one_per <- sapply(one_per,function(x) draw_1(pi[Set_id==x],Id[Set_id==x]))
  
  d1 <- length(Sam_one_per)
  info_one <- data.frame(Sample=Sam_one_per[seq(d1)%%2==1],pair=Sam_one_per[seq(d1)%%2==0],
                         weight = 2*G/(n-2*G)
  )
  
  info_two <- data.frame(Sample=Sam_two_per[1,],pair=Sam_two_per[2,],
                         weight = (n+ifelse(aa==1,1,0))/G*
                           Sam_two_per[4,]*Sam_two_per[5,]/Sam_two_per[3,] -1 
  )
  
  Sam_info=rbind(info_one,info_two)
  
  Sam_info$Sample <- Idd1[Sam_info$Sample]
  Sam_info$pair <- Idd1[Sam_info$pair]
  
  return(list(Sam_info,Group_id))
}


### Sampling procedure for prob <= 0.5
inter_small <- function(prob2,G) {
  n <- round(sum(prob2))
  N <- length(prob2)
  Id <- 1:N
  
  if (n %% 2 != 0) {
    temp <- assignset(prob2,bd=n-1)
    s1 <- draw_1(temp$set2,(temp$bdry+1):N)
    # s1 <- c(s1,prob[s1])
    prob2 <- temp[[1]]
    pp1 <- temp[[2]]
    Id <- 1:length(prob2)
    
    # temp <- assignset(prob2,bd=1)
    # s1 <- draw_1(temp$set1,1:temp$bdry)
    # # s1 <- c(s1,prob[s1])
    # prob2 <- temp[[2]]
    # pp1 <- temp[[1]]
    # Id <- 1:length(prob2)
  }
  if (round(sum(prob2))==2) {
    pi <- prob2
    Stra_id <- rep(1,length(prob2))
  } else {
    sp <- splitp2(prob2,bd=2)
    pi <- sp[[1]]
    Stra_id <- sp[[2]]
  }
  
  Group <- Pair <- Set <- rep(0,length(pi))
  
  Pair <- Stra_id
  
  keep_two <- numeric()
  for (i in unique(Stra_id)) {
    if(checksplit(pi[Stra_id==i]) == 1) keep_two <- c(keep_two,i)
  }
  
  if(sum(!(unique(Stra_id) %in% keep_two)) <= 2) keep_two <- unique(Stra_id)
  
  Idd1 <- Id[!(Stra_id %in% keep_two)]
  Idd2 <- Id[Stra_id %in% keep_two]
  
  Sam_info <- data.frame()
  prob_cd <- numeric()
  
  if(length(keep_two)!=0) {
    two_per <- sapply(keep_two, function(i) draw_2(pi[Stra_id==i],Id[Stra_id==i]))
    info_twoper <- data.frame(Sample=two_per[1,],pair=two_per[2,],
                              weight = (4*two_per[4,]*two_per[5,] - two_per[3,])/two_per[3,] 
    )
    Sam_info <- rbind(Sam_info,info_twoper)
    prob_cd <- c(prob_cd,pi[Idd2])
  }
  
  if(length(Idd1)!=0) {
    pi_inter <- splitp2(pi[!(Stra_id %in% keep_two)])
    n_int <- round(sum(pi[Idd1]))
    if(G < 1 | G > n_int/4) {
      stop(paste("G should between 1 and", n_int/4))
    }
    samp_inter <- intersamp3(pi_inter[[1]], pi_inter[[2]], G, Idd1, aa = n %% 2)
    Sam_info <- rbind(Sam_info, samp_inter[[1]])
    prob_cd <- c(prob_cd,pi_inter[[1]])[order(c(Idd2,Idd1))]
    Set[Idd1] <- pi_inter[[2]]
    Group[Idd1] <- samp_inter[[2]]
  }
  
  if(n %% 2 != 0) {
    Sam_info <- rbind(Sam_info,c(s1,s1,0))
    prob_cd <- c(prob_cd,pp1)
    Group <- c(Group,rep('one',length(pp1)))
    Pair <- c(Pair,rep('one',length(pp1)))
    Set <- c(Set,rep('one',length(pp1)))
  }
  
  Sam_info <- Sam_info[order(Sam_info[,1]),]
  
  return(list(Sam_info = Sam_info, prob_cd = prob_cd, 
              Stru_info = data.frame(Group,Pair,Set)))
}

### Sampling procedure for prob > 0.5
inter_large <- function(nprej,Stra_id) {
  N <- length(nprej)
  sam <- numeric()
  Samp_info <- data.frame()
  for(i in unique(Stra_id)) {
    pi <- nprej[Stra_id==i] 
    id <- (1:N)[Stra_id==i]
    if(round(sum(pi))==1) {
      a <- draw_1(pi,1:length(id))
      s <- id[-a]
      if(length(s) == 1) {
        info <- data.frame(Sample=s, pair= s, weight = 0)
        Samp_info <- rbind(Samp_info,info)
      } else if(length(s) > 1){
        per <- sample(s)
        per_n <- c(2:(length(pi)-1),1)
        info <- data.frame(Sample=per,pair=per[per_n])
        info <- data.frame(info,
                           weight = (length(id)-2)/2*((1-nprej[info$Sample])*(1-nprej[info$pair])-
                                                        (1-nprej[info$Sample]-nprej[info$pair]))/
                             (1-nprej[info$Sample]-nprej[info$pair]))
        Samp_info <- rbind(Samp_info,info)
      }
    } else if(round(sum(pi))==2) {
      A <- sum(pi/(1-pi)/2)
      N_pi <- length(pi)
      pi_ij <- 2*rep(pi/2,times=N_pi)*2*rep(pi/2,each=N_pi)*(1/(1-2*rep(pi/2,times=N_pi))+
                                                               1/(1-2*rep(pi/2,each=N_pi)))/2/(1+A)
      pi_ij <- matrix(pi_ij,N_pi)
      
      a <- draw_2(pi,1:length(id))
      s <- id[-a[1:2]]
      per <- sample(s)
      per_n <- c(2:(length(pi)-2),1)
      info <- data.frame(Sample=per,pair=per[per_n])
      
      pij <- sapply(1:length(s), function(i) {pi_ij[which(id==info$Sample[i]),
                                                    which(id==info$pair[i])]})
      
      info <- data.frame(info,
                         weight = (length(id)-3)/2*((1-nprej[info$Sample])*(1-nprej[info$pair])-
                                                      (1-nprej[info$Sample]-nprej[info$pair]+pij))/
                           (1-nprej[info$Sample]-nprej[info$pair]+pij))
      Samp_info <- rbind(Samp_info,info)
    }
  }
  Samp_info <- Samp_info[order(Samp_info$Sample),]
  return(Samp_info)
}

inter_select <- function(prob,n,G, method = 'randam') {
  if(method == 'random') assignset <- assignset else assignset <- cumu_round
  
  prob <- prob/sum(prob)*n
  N <- length(prob)
  Id <- 1:N
  
  if(sum(prob>2/3)==0) {
    n1 <- 0
    n2 <- n-n1
    N1 <- sum(prob > 2/3)
    N2 <- N-N1
    
    prej <- (1-prob[prob > 2/3])*(N1-n1)/(N1)
    
    prob2 <- prob[prob <= 2/3]*n2/n
    
    Id1 <- Id[prob > 2/3]
    Id2 <- Id[prob <= 2/3]
  } else {
    s1 <- sum(prob[prob > 0.5])
    n1 <- rbinom(1,1,s1-floor(s1))+floor(s1)
    n2 <- n-n1
    N1 <- sum(prob > 0.5)
    N2 <- N-N1
    
    prej <- (1-prob[prob > 0.5])*(N1-n1)/(N1-s1)
    
    prob2 <- prob[prob <= 0.5]*n2/(n-s1)
    
    Id1 <- Id[prob > 0.5]
    Id2 <- Id[prob <= 0.5]
  }
  
  prob_cd <- numeric()
  Stru_info <- data.frame()
  Part <- rep(0,N)
  
  Part[Id1] <- 1
  Part[Id2] <- 2
  
  if(length(prej)==0) {
    sam1 <- data.frame()
  } else if (N1-n1 == 0) {
    sam1 <- data.frame(Sample = Id1, pair=Id1, weight=0)
    prob_cd <- c(prob_cd,rep(1,n1))
    Stru_id <- c('1-1',rep(1,n1))
  } else {
    temp <- split_large(prej)
    nprej <- temp$prej
    Stra_id <- temp$Stra_id
    
    check1 <- (table(Stra_id)==3) * (tapply(nprej,Stra_id,function(x) round(sum(x)))==2)
    check2 <- (table(Stra_id)==2) * (tapply(nprej,Stra_id,function(x) round(sum(x)))==1)
    check <- check1 + check2
    if(sum(check>0)) {
      a <- which(check==1)
      id <- Id1[Stra_id %in% a]
      Id1 <- Id1[-id]
      
      prob2 <- c(prob2,1-nprej[(Stra_id %in% a)])
      prob2 <- prob2[order(c(Id2,id))]
      Id2 <- sort(c(Id2,id))
      
      nprej <- nprej[!(Stra_id %in% a)]
      Stra_id <- Stra_id[!(Stra_id %in% a)]
    }
    prob_cd <- c(prob_cd,1-nprej)
    Group <- Stra_id
    Pair <- Set <- rep(0,length(Id1))
    sam1 <- inter_large(nprej,Stra_id)
    sam1$Sample <- Id1[sam1$Sample]
    sam1$pair <- Id1[sam1$pair]
    Stru_info <- data.frame(Group,Pair,Set)
  }
  
  
  sam2 <- inter_small(prob2,G)
  sam2[[1]]$Sample <- Id2[sam2[[1]]$Sample]
  sam2[[1]]$pair <- Id2[sam2[[1]]$pair]
  prob_cd <- c(prob_cd,sam2[[2]])[order(c(Id1,Id2))]
  Stru_info <- rbind(Stru_info,sam2[[3]])[order(c(Id1,Id2)),]
  
  Sam_info <- rbind(sam1,sam2[[1]])
  Sam_info <- Sam_info[order(Sam_info[,1]),]
  
  Sample <- sort(unique(c(sam1$Sample,sam2[[1]]$Sample,sam2[[1]]$pair)))
  return(list(Sample = Sample, Var_info = Sam_info, 
              Stru_info = data.frame(Id,prob_cd,Part,Stru_info)[Sample,],
              prob_cd = prob_cd))
}

##### Variance estimation #####

var_est1 <- function(y, pt, var_info) {
  v <- sum((y[var_info[,1]]/pt[var_info[,1]]-
              y[var_info[,2]]/pt[var_info[,2]])^2*var_info[,3])
  return(v)
}

var_est2 <- function(y, pt, var_info, pcd) {
  w <- 1/pt-1/pcd
  
  samp <- unique(c(var_info$Sample,var_info$pair))
  
  yhat <- sum(y[samp]/pt[samp])
  
  v <- sum((y[var_info[,1]]/pt[var_info[,1]]-
              y[var_info[,2]]/pt[var_info[,2]])^2*var_info[,3])
  
  v2 <- v + (yhat-sum(y[samp]/pcd[samp]))^2 -
    sum((y[var_info[,1]]*w[var_info[,1]]-
           y[var_info[,2]]*w[var_info[,2]])^2*var_info[,3])
  
  v2
}

var_modi <- function(y, pt, var_info) {
  within <- which(var_info[,3] > median(var_info[,3]))
  var_w <- var_info[within,]
  var_b <- var_info[-within,]
  s_w <- var_est1(y,prob,var_w)
  s_b <- var_est1(y,prob,var_b)
  
  return(c(s_w, s_b))
}


##### two-per #####

# Data are collected from the survey 
# 'w' are original inclusion probabilities, i.e. weight for element
fvtwo2 <- function(Data,twoper,w_samp) {
  V2 <- (4*twoper[4,]*twoper[5,]/twoper[3,]-1)*(Data[1,]/w_samp[1,]-
                                                  Data[2,]/w_samp[2,])^2
  
  V2[twoper[3,]==1] <- 0
  
  Vhat <- sum(V2)
  
  return(Vhat)
}

fttwo <- function(Data,w_samp) {
  Ttwo <- Data[1,]/w_samp[1,]+Data[2,]/w_samp[2,]
  return(Ttwo=sum(Ttwo))
}

vhat_two <- function(Data,twoper,w_samp,pi_samp) {
  bias2 <- (fttwo(Data,w_samp)-fttwo(Data,pi_samp))^2- fvtwo2(Data,twoper,1/(1/w_samp-1/pi_samp))
  Vtwo <- fvtwo2(Data,twoper,w_samp) + max(bias2,0)  
  return(c(Vtwo=Vtwo,bias=bias2))
}

ftwo <- function(y,prob,n) {
  Id <- 1:length(prob)
  if(n==2) {
    pi <- prob
    stra_id <- rep(1,length(prob))
  } else {
    sp <- splitp2(prob,2)
    pi <- sp[[1]]
    stra_id <- sp[[2]]
  }
  
  twoper <- sapply(unique(stra_id),function(x) draw_2(pi[stra_id==x],Id[stra_id==x]))
  
  Data <- rbind(y[twoper[1,]],y[twoper[2,]])
  w_samp <- rbind(prob[twoper[1,]],prob[twoper[2,]])
  pi_samp <- rbind(pi[twoper[1,]],pi[twoper[2,]])
  
  Ttwo <- fttwo(Data,w_samp)
  Vtwo <- vhat_two(Data,twoper,w_samp,pi_samp)
  return(c(Ttwo,Vtwo,sum(y*prob/pi)))
}

##### one-per #####
fone <- function(y,prob,n) {
  Id <- 1:length(prob)
  sp <- splitp2(prob,n)
  pi <- sp[[1]]
  stra_id <- sp[[2]]
  
  oneper <- sapply(unique(stra_id),function(x) draw_1(pi[stra_id==x],Id[stra_id==x]))
  z_one <- y[oneper]/prob[oneper]*n
  Tone <- sum(z_one)/n
  Vone2 <- var(z_one)/n
  Vone3 <- sum((z_one[(1:floor(n/2))*2]-z_one[(1:floor(n/2))*2-1])^2)/n/n
  Vone4 <- mean(diff(z_one)^2)/2/n
  return(c(Tone,Vone2,Vone3,Vone4))
}



##### systematic #####
fsys <- function(y,prob) {
  a <- runif(1)
  n <- round(sum(prob))
  sys_sam <- sapply(a+(0:(n-1)), function(x) which(cumsum(prob) >= x)[1])
  z <- y[sys_sam]/prob[sys_sam]*n
  Tsys <- sum(z)/n
  c <- sum(prob^2)/n
  Vsys1 <- 0
  for(i in 2:n)
    for(j in 1:i) {
      Vsys1 <- Vsys1 + (1-prob[sys_sam[i]]-prob[sys_sam[j]]+c)*(z[i]/n-z[j]/n)^2/(n-1)
    }
  Vsys2 <- var(z)/n
  Vsys3 <- sum((z[(1:floor(n/2))*2]-z[(1:floor(n/2))*2-1])^2)/n/n
  Vsys4 <- mean(diff(z)^2)/2/n
  return(c(Tsys,Vsys1,Vsys2,Vsys3,Vsys4))
}




##### simulation function #####

inter_small_simu <- function(prob,G,y) {
  n <- round(sum(prob))
  N <- length(prob)
  Id <- 1:N
  
  Sam_info <- list()
  
  sp <- splitp2(prob)
  pi <- sp[[1]]
  g1 <- sp[[2]]
  
  v1 <- mean(sapply(unique(g1), function(x) (V1(y[g1==x]/prob[g1==x],pi[g1==x]))))
  v2 <- mean(sapply(unique(g1), function(x) (V2(y[g1==x]/prob[g1==x],pi[g1==x]))))
  sb <- diff(tapply(y/prob*pi,g1,sum))
  sb <- sb[seq(length(sb)) %% 2 == 1]
  sb <- sum(sb^2)/n
  
  v1_c <- mean(sapply(unique(g1), function(x) (V1(y[g1==x]/pi[g1==x],pi[g1==x]))))
  v2_c <- mean(sapply(unique(g1), function(x) (V2(y[g1==x]/pi[g1==x],pi[g1==x]))))
  sb_c <- diff(tapply(y,g1,sum))
  sb_c <- sb[seq(length(sb)) %% 2 == 1]
  sb_c <- sum(sb^2)/n
  
  for(i in 1:length(G)) {
    samp_inter1 <- intersamp3(pi, g1, G[i], 1:N, aa = 0)
    temp <- samp_inter1[[1]]
    temp <- temp[order(temp[,1]),]
    Sam_info[[i]] <- temp
  }
  
  # samp_inter1 <- intersamp3(pi, g1, G1, 1:N, aa = 0)
  # Sam_info1 <- rbind(Sam_info, samp_inter1[[1]])
  # 
  # samp_inter2 <- intersamp3(pi, g1, G2, 1:N, aa = 0)
  # Sam_info2 <- rbind(Sam_info, samp_inter2[[1]])
  # 
  # samp_inter3 <- intersamp3(pi, g1, G3, 1:N, aa = 0)
  # Sam_info3 <- rbind(Sam_info, samp_inter3[[1]])
  
  # Sam_info1 <- Sam_info1[order(Sam_info1[,1]),]
  # Sam_info2 <- Sam_info2[order(Sam_info2[,1]),]
  # Sam_info3 <- Sam_info3[order(Sam_info3[,1]),]
  
  # return(list(Sam_info1 = Sam_info1, Sam_info2 = Sam_info2, Sam_info3 = Sam_info3,
  #             prob_cd = pi, c(v1,v2,sb)))
  return(list(Sam_info, prob_cd = pi, c(v1,v2,sb), c(v1_c,v2_c,sb_c)))
}

inter_select_simu <- function(prob,n,G,y) {
  prob <- prob/sum(prob)*n
  N <- length(prob)
  Id <- 1:N
  
  s1 <- sum(prob[prob > 0.5])
  n1 <- rbinom(1,1,s1-floor(s1))+floor(s1)
  n2 <- n-n1
  N1 <- sum(prob > 0.5)
  N2 <- N-N1
  
  prej <- (1-prob[prob > 0.5])*(N1-n1)/(N1-s1)
  
  prob2 <- prob[prob <= 0.5]*n2/(n-s1)
  
  Id1 <- Id[prob > 0.5]
  Id2 <- Id[prob <= 0.5]
  
  prob_cd <- numeric()
  
  if(length(prej)==0) {
    sam1 <- data.frame()
  } else if (N1-n1 == 0) {
    sam1 <- data.frame(Sample = Id1, pair=Id1, weight=0)
    prob_cd <- c(prob_cd,rep(1,n1))
  } else {
    temp <- split_large(prej)
    nprej <- temp$prej
    Stra_id <- temp$Stra_id
    
    check1 <- (table(Stra_id)==3) * (tapply(nprej,Stra_id,function(x) round(sum(x)))==2)
    check2 <- (table(Stra_id)==2) * (tapply(nprej,Stra_id,function(x) round(sum(x)))==1)
    check <- check1 + check2
    if(sum(check>0)) {
      a <- which(check==1)
      id <- Id1[Stra_id %in% a]
      Id1 <- Id1[-id]
      
      prob2 <- c(prob2,1-nprej[(Stra_id %in% a)])
      prob2 <- prob2[order(c(Id2,id))]
      Id2 <- sort(c(Id2,id))
      
      nprej <- nprej[!(Stra_id %in% a)]
      Stra_id <- Stra_id[!(Stra_id %in% a)]
    }
    prob_cd <- c(prob_cd,1-nprej)
    sam1 <- inter_large(nprej,Stra_id)
    sam1$Sample <- Id1[sam1$Sample]
    sam1$pair <- Id1[sam1$pair]
  }
  
  sam2 <- inter_small_simu(prob2,G,y[Id2])
  sam2[[1]][[1]]$Sample <- Id2[sam2[[1]][[1]]$Sample]
  sam2[[1]][[1]]$pair <- Id2[sam2[[1]][[1]]$pair]
  sam2[[1]][[2]]$Sample <- Id2[sam2[[1]][[2]]$Sample]
  sam2[[1]][[2]]$pair <- Id2[sam2[[1]][[2]]$pair]
  sam2[[1]][[3]]$Sample <- Id2[sam2[[1]][[3]]$Sample]
  sam2[[1]][[3]]$pair <- Id2[sam2[[1]][[3]]$pair]
  prob_cd <- c(prob_cd,sam2[[2]])[order(c(Id1,Id2))]
  
  Sam_info1 <- rbind(sam1,sam2[[1]][[1]])
  Sam_info2 <- rbind(sam1,sam2[[1]][[2]])
  Sam_info3 <- rbind(sam1,sam2[[1]][[3]])
  
  return(list(Sam_info1 = Sam_info1, Sam_info2 = Sam_info2, Sam_info3 = Sam_info3,
              prob_cd=prob_cd,sam2[[3]]))
}

V1 <- function(y,pi) {
  var1 <- sum(y^2*pi) - sum(rep(y*pi,times=length(y))*rep(y*pi,each=length(y)))
  return(var1)
}

V2 <- function(y,pi) {
  pi <- pi/sum(pi)
  N <- length(y)
  A <- sum(pi/(1-2*pi))
  pi_ij <- 2*rep(pi,times=N)*2*rep(pi,each=N)*(1/(1-2*rep(pi,times=N))+
                                                 1/(1-2*rep(pi,each=N)))/2/(1+A)
  pi_ij[(0:(N-1))*N+1:N] <- 2*pi
  var2 <- sum((pi_ij-4*rep(pi,times=N)*rep(pi,each=N))*rep(y,times=N)*rep(y,each=N)/4)
  return(var2)
} 

