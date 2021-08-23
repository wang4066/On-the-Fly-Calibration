library(irtoys)

#kit: ID of operational items in the bank
mem <- function(resp,item,Q=21,mu,sig,eps =1e-3,max.iter=500)
{
  
  #One set of quadrature points and weights based on standard normal
  D <- 1.702

  GH <- normal.qu(Q,mu=mu,sigma = sig)
  Xg <- GH$quad.points
  Ag <- GH$quad.weights
  #starting values
  a <- mu_beta
  b <- -mu_alpha/mu_beta
  # a <- rep(1,B)
  # b <- rep(0,B)  
  # a[kit] <- beta[kit]
  # b[kit] <- -alpha[kit]/beta[kit]  
  #intermediate quantities
  ngq <- double(Q)
  rgq <- double(Q)
  # pstar <- matrix(,J,Q)
  
  lik <- function(ip){
    a <- ip[1]
    b <- ip[2]
    p <- 1/(1+exp(-a*D*(Xg-b)))
    p[p==1] <- 1-1e-7
    p[p==0] <- 1e-7
    -sum(rgq*log(p)+(ngq-rgq)*log(1-p))
  }
  
  gr <- function(ip){
    a <- ip[1]
    b <- ip[2]
    p <- 1/(1+exp(-a*D*(Xg-b)))
    c(
      -sum((rgq-ngq*p)*D*(Xg-b)),
      sum((rgq-ngq*p)*a*D)
    )
  }
  
  #### start of the first cycle ####
  for (j in 1:J){
    #current new item ID
    temp.j <- nit[j]
    tempa <- a[temp.j]
    tempb <- b[temp.j]
    #subject ID who gets item j
    jind <- which(item==temp.j,arr.ind = T)
    Nj <- nrow(jind)
    niq <- liq <- matrix(double(Nj*Q),Nj,Q)
    pi <- double(Nj)
    
    for (i in 1:Nj){
      temp.n <- jind[i,1]
      #item IDs administered to subject i
      items <- item[temp.n,]
      #operational item indicator for subject i
      Mi <- items%in%kit
      #number of operational items
      M <- sum(Mi)
      #no operational item was given to this subject
      if(M==0){
        liq[i,] <- NA
      }else{
      #responses to operational items
      mresp <- resp[temp.n,Mi]
      ljq <-  matrix(double(M*Q),M)
      
      for (m in 1:M){
        #current operational item ID
        mi <- items[Mi][m]
        temp.b <- b[mi]
        temp.a <- a[mi]
        tempp <-  1/(1+exp(-temp.a*D*(Xg-temp.b)))
        ljq[m,] <- tempp ^mresp[m]*(1-tempp) ^(1-mresp[m])
      }
      
      liq[i,] <-  apply(ljq,2,prod)*Ag
      }

      # pi[i] <- sum(temp.post)
    } 
    #N
    pi <- apply(liq,1,sum)
    niq <- liq/pi
    ngq <- colSums(niq,na.rm = T)
    #R
    jresp <- resp[item==temp.j]
    riq <-jresp*niq
    rgq <- colSums(riq,na.rm = T) 
    
    #M-step
    temp <- optim(c(tempa,tempb),lik,gr,method = "L-BFGS-B", 
                  lower =c(0.1,-5),upper = c(3,5))$par
    
    a[temp.j] <- temp[1]
    b[temp.j] <- temp[2]
  }
  
#population mean 
  Pjq <-Liq <-  matrix(double(N*Q),N)
  for (i in 1:N) {
    ljq <-  matrix(double(l*Q),l)
    pitem <- item[i,]
    presp <- resp[i,]
    for (j in 1:l){
      mi <- pitem[j]
      temp.b <- b[mi]
      temp.a <- a[mi]
      tempp <-  1/(1+exp(-temp.a*D*(Xg-temp.b)))
      ljq[j,] <- tempp ^presp[j]*(1-tempp) ^(1-presp[j])
    }
    Liq[i,] <-  apply(ljq,2,prod)*Ag
  }
  pi <- rowSums(Liq)
  Pjq <- Liq/pi
  
  mu1 <- sum(Xg*colSums(Pjq))/N
  # sig <- sum((Xg-mu)^2*colSums(Pjq))/N
  #   
  # ts <- c(sig)
  
  #### end of the first cycle (OEM) ####
  
  #start of the second cycle
  df.a <- df.b <-  1
  iter <- 1
  
  while(max(df.a)>eps & max(df.b)>eps  & iter<max.iter)
  {
    iter <- iter+1
    aold <- a
    bold <- b
    
    GH <- normal.qu(Q,mu=mu1,sigma = sig)
    Xg <- GH$quad.points
    Ag <- GH$quad.weights
    
    for (j in 1:J){
      #current new item ID
      temp.j <- nit[j]
      tempa <- a[temp.j]
      tempb <- b[temp.j]
      #subject ID who gets item j
      jind <- which(item==temp.j,arr.ind = T)
      Nj <- nrow(jind)
      niq <- liq <- matrix(double(Nj*Q),Nj,Q)
      pi <- double(Nj)
      
      for (i in 1:Nj){
        temp.n <- jind[i,1]
        #item IDs administered to subject i
        items <- item[temp.n,]

        #responses to all items
        mresp <- resp[temp.n,]
        ljq <-  matrix(double(l*Q),l)
        
        for (m in 1:l){
          mi <- items[m]
          temp.b <- b[mi]
          temp.a <- a[mi]
          tempp <-  1/(1+exp(-temp.a*D*(Xg-temp.b)))
          ljq[m,] <- tempp ^mresp[m]*(1-tempp) ^(1-mresp[m])
        }
        
        liq[i,] <- apply(ljq,2,prod)*Ag
        # pi[i] <- sum(temp.post)
      } 
      
      pi <- apply(liq,1,sum)
      niq <- liq/pi
      ngq <- colSums(niq)
      
      jresp <- resp[item==temp.j]
      riq <-jresp*niq
      rgq <- colSums(riq) 
      
      #M-step
      temp <- optim(c(tempa,tempb),lik,gr,method = "L-BFGS-B", 
                    lower =c(0.1,-5),upper = c(3,5))$par
      
      a[temp.j] <- temp[1]
      b[temp.j] <- temp[2]
    }
    
    #population mean 
    Pjq <-Liq <-  matrix(double(N*Q),N)
    for (i in 1:N) {
      ljq <-  matrix(double(l*Q),l)
      pitem <- item[i,]
      presp <- resp[i,]
      for (j in 1:l){
        mi <- pitem[j]
        temp.b <- b[mi]
        temp.a <- a[mi]
        tempp <-  1/(1+exp(-temp.a*D*(Xg-temp.b)))
        ljq[j,] <- tempp ^presp[j]*(1-tempp) ^(1-presp[j])
      }
      Liq[i,] <-  apply(ljq,2,prod)*Ag
    }
    pi <- rowSums(Liq)
    Pjq <- Liq/pi
    mu1 <- sum(Xg*colSums(Pjq))/N
    # sig <- sum((Xg-mu)^2*colSums(Pjq))/N
    # ts <- c(ts,sig)
    
    df.a <- abs(aold-a)
    df.b <- abs(bold-b)
  }
  
  return(list(est=cbind(a,b),iter=iter,pop=c(mu1,sig),ts=ts))
}

#end of function