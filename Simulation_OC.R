##################################################
####Online Calibration + Random Item Selection####
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
##### est_t: ability parameters for each session, a R by N by nsess array
##### est_mu: population mean for each session, a R by nsess matrix
##### est_sig: population standard deviation for each session, a R by nsess matrix
##### est_a: the deviance of true and estimated $\alpha$ parameters,  a R by nsess by B array  
##### est_b: the deviance of true and estimated $\beta$ parameters,  a R by nsess by B array                          
##### tru_t: true ability parameters for each session, a R by N by nsess array
##### exposure: the accumulative times an item is utilized for each session, a R by B by nsess array 
##### adit: the item index administrated for each person, a R by N by L array
##### kit: the precalibrated items' index, a vector of length B*prop
##################################################
set.seed(233)
library(mirtCAT)
library(irtoys)
#varied conditions
#total test length
L <- 120
#precalibrated items
prop <- 0.5
source("oc.R")

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
temp.it <- temp.resp <- matrix(,N,l)
mu_alpha <- mu_beta <- numeric(B)
mu_theta <-matrix(,N,nsess)
#index of precalibrated items
kit <- sample.int(B,B*prop)
#number and index of new items
J <- B-length(kit)
nit <- c(1:500)[-kit]


#objects to store item exposure and item selection indexes of candidate items
exposure <- array(0,c(R,B,4))
adit <- array(,c(R,N,L))

est_a <- est_b <- array(,c(R,4,B))
est_t <- array(,c(R,N,nsess))
tru_t <- array(,c(R,N,4))
est.mu <-est.sig <-  matrix(,R,nsess)

for (r in 1:R){
  tryCatch({
  #generate thetas for each replication
  theta <- rnorm(N)
  Theta <- cbind(theta,
                 theta+rnorm(N,.2,.1),
                 theta+rnorm(N,.2,.1)+rnorm(N,.2,.1),
                 theta+rnorm(N,.2,.1)+rnorm(N,.2,.1)+rnorm(N,.2,.1))

  #set starting values
  mu_alpha=rep(mu.a,B)
  mu_beta=rep(mu.b,B)
  
  mu_alpha[kit]=alpha[kit]
  mu_beta[kit]=beta[kit]

  mu <- 0
  sig <- 1
  
  temp.ad <- matrix(0,N,B)
  for (k in 1:nsess){ 
    #make sure no items can be selected twice in the same session
    temp.select <- matrix(0,N,B)
    temp.select[temp.ad>1] <- 1
    for (i in sample.int(N)) {
        temp.item <- sample(c(1:500)[!temp.select[i,]],l)
        
        temp.ad[i,temp.item] <- temp.ad[i,temp.item]+1
        #record administered item and # of exposure
        adit[r,i,((k-1)*l+1):(k*l)] <- temp.it[i,] <- temp.item
        exposure[r,temp.item,k] <- exposure[r,temp.item,k]+1
        resp[i,((k-1)*l+1):(k*l)] <- temp.resp[i,] <- ifelse(pnorm(alpha[temp.item]+
                                                                             beta[temp.item]*Theta[i,k])>
                                                               runif(l),1,0)
       }
    #compute item parameter estimates from mem
        out <- mem(temp.resp,temp.it,mu=mu,sig=sig)
        upd <- out$est
        mu_alpha <- -upd[,1]*upd[,2]
        mu_beta<- upd[,1]
        est.mu[r,k] <-  mu <- out$pop[1]
        est.sig[r,k] <-  sig <- out$pop[2]
        #compute theta estimates 
        for (i in 1:N){
          temp.item <- temp.it[i,]
          pars <- data.frame(a1=mu_beta[temp.item]*1.702,
                             d=mu_alpha[temp.item]*1.702)
          mod <- generate.mirt_object(pars,itemtype = "2PL")
          mu_theta[i,k]<- fscores(mod,response.pattern =temp.resp[i,],
                                  append_response.pattern = F)[1]
        }
        est_a[r,k,] <- mu_alpha-alpha
        est_b[r,k,] <- mu_beta-beta
        if(k!=nsess)
        {exposure[r,,k+1] <- exposure[r,,k]}
  }
    
  
est_t[r,,] <- mu_theta
# est_ts[r,,] <- sig_theta
tru_t[r,,] <- Theta
  }, error=function(e){cat("ERROR:",conditionMessage(e), "\n")})
print(r)
}

#bias and RMSE of thetas at the end of each session
bias_t <- colMeans(cbind(rowMeans(est_t[,,1]-tru_t[,,1]),
                         rowMeans(est_t[,,2]-tru_t[,,2]),
                         rowMeans(est_t[,,3]-tru_t[,,3]),
                         rowMeans(est_t[,,4]-tru_t[,,4])))
rmse_t <- colMeans(cbind(sqrt(rowMeans((est_t[,,1]-tru_t[,,1])^2)),
                         sqrt(rowMeans((est_t[,,2]-tru_t[,,2])^2)),
                         sqrt(rowMeans((est_t[,,3]-tru_t[,,3])^2)),
                         sqrt(rowMeans((est_t[,,4]-tru_t[,,4])^2))))
#bias and RMSE of item parameters
ic <- c(
  mean(colMeans(est_a[,4,-kit])),
  mean(colMeans(est_b[,4,-kit])),
  mean(sqrt(colMeans((est_a[,4,-kit]^2)))),
  mean(sqrt(colMeans((est_b[,4,-kit]^2))))
)

rs=list(est_t=est_t,est.mu=est.mu,est.sig=est.sig,est_a=est_a,est_b=est_b,tru_t=tru_t,exposure=exposure,adit=adit,kit=kit)
saveRDS(rs,"oc_120_prop50.rds")