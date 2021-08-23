mpwi <- function(){
    z <- bv*tv+av
    #metric
  p <- 1/(1+exp(-1.7*z))
    if(any(p==1))
    {p[p==1] <- 1-1e-5}
    q <- 1-p
    #metric
    tinfo <- (1.7*bv)^2*p*q
    info <-  rowMeans(tinfo) 
  h <- max(info)
  ind <- (1-j/L)*runif(B,max = h)+2*j*info/L
  return(ind*(1-exposure[r,]/(N*mex))) #max exposure 0.4,r is the current replication
}

