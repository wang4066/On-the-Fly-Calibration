##################################################
#########Bayesian Update + MPWI Item Selection####
##################################################
# Simulation Conditions:
# L: test length
# prop: the proportion of precalibrated items in the item bank
##################################################
# Fixed Parameters for the current study:
# R: the number of replications
# nsess: the number of sessions
# N: sample size
# B: item bank size
#################################################
##### Outputs a list of updated parameters                                                             
##### est_t: ability parameters, a R by N by L+1 array where the first column includes the initial values
##### est_ts: ability parameters' standard deviations, a R by N by L+1 array
##### est_a: the deviance of true and estimated $\alpha$ parameters and their standard deviations,  a R by nsess by B by 2 array  
##### est_b: the deviance of true and estimated $\beta$ parameters and their standard deviations,  a R by nsess by B by 2 array                          
##### tru_t: true ability parameters for each session, a R by N by nsess array
##### exposure: the accumulative times an item is utilized for each session, a R by B by nsess array 
##### adit: the item index administrated for each person, a R by N by L array
##### kit: the precalibrated items' index, a vector of length B*prop
##################################################
set.seed(233)
#varied conditions
#total test length
L <- 120 #80
#precalibrated items
prop <- 0.5 #0.5

#number of session
nsess <- 4
#session length
l <- L/nsess
N <- 3000
#bank size
B <- 500
#starting values and data-generating distributions
mu.b <- 1
mu.a <- 0
mu.t <- 0
sig.b <- 0.60
sig.a <- 1
sig.t <- 1 
#data-generating distributions for beta
sigb <- 0.32

kappa <- 0.0001
gam <- c(-Inf,0,Inf)
#maximal exposure rates
mex <- .4
#number of draws for MPWI
D <- 500
#number of replications
R <- 30

#generate true item parameters
alpha <- rnorm(B,mu.a,sig.a)
repeat{
  beta <- rnorm(B,mu.b,sigb)
  if(min(beta)>0)
    break
}
ip <- cbind(alpha,beta)

resp <- matrix(,N,L)
mu_alpha <- mu_beta <- sig_alpha <- sig_beta <- numeric(B)
mu_theta <- sig_theta <- matrix(,N,L+1)
#index of precalibrated items
kit <- sample.int(B,B*prop)

#load functions to compute parameter update and mpwi
source("par_update.R")
source("mpwi.R")

#objects to store item exposure and item selection indexes of candidate items
exposure <- array(0,c(R,B,4)) 
adit <- array(,c(R,N,L))
tind <- matrix(,B,2)
tind[,1] <- 1:B

est_a <-est_b <- array(,c(R,4,B,2))
est_t <- est_ts <- array(,c(R,N,L+1))
tru_t <- array(,c(R,N,4))

for (r in 1:R){
  #generate thetas for each replication
  theta <- rnorm(N)
  Theta <- cbind(theta,
                 theta+rnorm(N,.2,.1),
                 theta+rnorm(N,.2,.1)+rnorm(N,.2,.1),
                 theta+rnorm(N,.2,.1)+rnorm(N,.2,.1)+rnorm(N,.2,.1))
  #item parameter vectors to compute mpwi
  av <- matrix(rnorm(B*D,mu.a,sig.a),nrow=B,ncol=D)
  bv <-  matrix(rnorm(B*D,mu.b,sig.b),nrow=B,ncol=D)
  #set starting values
  mu_alpha <- rep(mu.a,B)
  mu_beta <- rep(mu.b,B)
  mu_theta[,1] <- rep(mu.t,B)
  sig_alpha <- rep(sig.a,B)
  sig_beta <- rep(sig.b,B)
  sig_theta[,1] <- rep(1,N)
  
  mu_alpha[kit] <- alpha[kit]
  mu_beta[kit] <- beta[kit]
  sig_beta[kit] <- 0 
  sig_alpha[kit] <- 0
  
  av[kit,] <- matrix(alpha[kit],nrow=length(kit),ncol=D)
  bv[kit,] <- matrix(beta[kit],nrow=length(kit),ncol=D)
  
  temp.ad <- matrix(0,N,B)
  for (k in 1:nsess){ 
    #make sure no items can be selected twice in the same session
    temp.select <- matrix(0,N,B)
    temp.select[temp.ad>1] <- 1
    for (i in sample.int(N)) {
      #theta vector to compute mpwi
      tv <- rnorm(D,mu_theta[i,(k-1)*l+1],sig_theta[i,(k-1)*l+1])
      for (jj in 1:l) {
        j <- (k-1)*l+jj
        #compute mpwi for candidate items
        tind[,2] <- mpwi()
        temp.ind <- tind[!temp.select[i,],]
        temp.item <- temp.ind[which.max(temp.ind[,2]),1]
        temp.select[i,temp.item] <- T
        temp.ad[i,temp.item] <- temp.ad[i,temp.item]+1
        #record administered item and # of exposure
        adit[r,i,j] <- temp.item
        exposure[r,temp.item,k] <- exposure[r,temp.item,k]+1
        
        resp[i,j] <- response <- ifelse(pnorm(alpha[temp.item]+
                                                        beta[temp.item]*Theta[i,k])>
                                          runif(1),1,0)
        #compute parameter update
        upd <- par.update()
        mu_alpha[temp.item] <-  mu_alpha[temp.item]+upd[[1]]
        mu_beta[temp.item] <- mu_beta[temp.item]+upd[[2]]
        mu_theta[i,j+1]  <- mu_theta[i,j]+upd[[3]]
        sig_alpha[temp.item] <-  sig_alpha[temp.item]* sqrt(max(1-upd[[4]],kappa))
        sig_beta[temp.item] <-  sig_beta[temp.item]* sqrt(max(1-upd[[5]],kappa))
        sig_theta[i,j+1] <-  sig_theta[i,j]* sqrt(max(1-upd[[6]],kappa))
        #update item parameter vector
        av[temp.item,] <- rnorm(D,mu_alpha[temp.item],sig_alpha[temp.item])
        bv[temp.item,] <- rnorm(D,mu_beta[temp.item],sig_beta[temp.item])
        #add a constant to the posterior variance of theta at the end of a session
        if(j%%l==0&j!=L)
        {sig_theta[i,j+1] <- sig_theta[i,j+1]+.5}
        # update theta vector
        tv <- rnorm(D,mu_theta[i,j+1],sig_theta[i,j+1] )   
      }
    }
    est_a[r,k,,] <- cbind(mu_alpha-alpha,sig_alpha)
    est_b[r,k,,] <- cbind(mu_beta-beta,sig_beta)
    if(k!=nsess)
    {exposure[r,,k+1] <- exposure[r,,k]}
    }
  est_t[r,,] <- mu_theta
  est_ts[r,,] <- sig_theta
  tru_t[r,,] <- Theta
  print(r)
}

#bias and RMSE of thetas at the end of each session
bias_t <- colMeans(cbind(rowMeans(est_t[,,l+1]-tru_t[,,1]),
                         rowMeans(est_t[,,2*l+1]-tru_t[,,2]),
                         rowMeans(est_t[,,3*l+1]-tru_t[,,3]),
                         rowMeans(est_t[,,L+1]-tru_t[,,4])))
rmse_t <- colMeans(cbind(sqrt(rowMeans((est_t[,,l+1]-tru_t[,,1])^2)),
                         sqrt(rowMeans((est_t[,,2*l+1]-tru_t[,,2])^2)),
                         sqrt(rowMeans((est_t[,,3*l+1]-tru_t[,,3])^2)),
                         sqrt(rowMeans((est_t[,,L+1]-tru_t[,,4])^2))))
#bias and RMSE of item parameters
mean(colMeans(est_a[,4,-kit,1]))
mean(colMeans(est_b[,4,-kit,1]))
mean(colMeans(est_a[,4,-kit,2]))
mean(colMeans(est_b[,4,-kit,2]))
mean(sqrt(colMeans((est_a[,4,-kit,1]^2))))
mean(sqrt(colMeans((est_b[,4,-kit,1]^2))))

rs=list(est_t=est_t,est_ts=est_ts,est_a=est_a,est_b=est_b,tru_t=tru_t,exposure=exposure,adit=adit,kit=kit)
saveRDS(rs,"mpwi_120_prop50.rds")