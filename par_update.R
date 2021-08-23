Omega <-  function(v1,v2,s){
  (dnorm(v1/s)-dnorm(v2/s))/(pnorm(v1/s)-pnorm(v2/s));
}
Delta <-  function(v1,v2,s,resp){
  if (resp==0)
    return((-(v2/s)*dnorm(v2/s))/(pnorm(v1/s)-pnorm(v2/s)) +
             ((dnorm(v1/s)-dnorm(v2/s))/(pnorm(v1/s)-pnorm(v2/s)))^2)
  if (resp==1)
    return(((v1/s)*dnorm(v1/s))/(pnorm(v1/s)-pnorm(v2/s)) +
             ((dnorm(v1/s)-dnorm(v2/s))/(pnorm(v1/s)-pnorm(v2/s)))^2)
}

par.update <- function(){
  std <- sqrt(1+(sig_alpha[temp.item])^2+
             (sig_beta[temp.item]* mu_theta[i,j])^2 +
             (sig_theta[i,j]*mu_beta[temp.item])^2)
  arg1 <- mu_alpha[temp.item]+mu_beta[temp.item]*mu_theta[i,j]-gam[response+1]
  arg2 <- mu_alpha[temp.item]+mu_beta[temp.item]*mu_theta[i,j]-gam[response+2]
  
  Ometmp <- Omega(arg1,arg2,std);
  Deltmp <- Delta(arg1,arg2,std,response);
  
  Ome.a <- (sig_alpha[temp.item])^2/std *Ometmp;
  Ome.b <-  (sig_beta[temp.item])^2/std * mu_theta[i,j] * Ometmp;
  Ome.t <-  (sig_theta[i,j])^2/std * mu_beta[temp.item] * Ometmp;
  Del.a <-  (sig_alpha[temp.item]/std)^2 *Deltmp;
  Del.b <-  (sig_beta[temp.item] *mu_theta[i,j]/std)^2 *Deltmp;
  Del.t <-  (sig_theta[i,j] *mu_beta[temp.item]/std)^2 *Deltmp;
  
  return(list(Ome.a,Ome.b,Ome.t,Del.a,Del.b,Del.t))
}