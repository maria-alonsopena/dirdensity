
# densities, internal

# cylindrical densities


dms<-function(x, y, mu, kappa, sigma, nu, lambda){


  mu1<-mu[1]
  mu2<-mu[2]

  val<-exp( - ( y - ( mu2 + lambda * cos(x - nu )))^2/( 2*sigma^2 ) + kappa * cos( x - mu1 ))/(sigma * (2*pi)^(3/2) * besselI(kappa,0) )

  return(val)

}


djw<-function(x,y,mu,lambda,kappa){

  val<-(lambda^2-kappa^2)^(1/2)/(2*pi)*exp(-lambda*y+kappa*y*cos(x-mu))
  return(val)
}


dal<-function(x,y,mu,lambda, beta, alpha, kappa){

  val<-alpha*(beta^alpha)*(1+lambda*sin(x-mu))*y^(alpha-1)*exp(-(beta*y)^alpha*(1-tanh(kappa)*cos(x-mu)))/(2*pi*cosh(kappa))
  return(val)
}


#install.packages("Rsolnp")



## 2.1 The general expression of sine-skewed densities.
## Equation (2.1) for the bivariate case
## This function can be employed over any predefined toroidal symmetric density: d"model"=function(x,y,mu1,mu2,kappa1,kappa2,rho)

dsinesk=function(x,y,mu1,mu2,kappa1,kappa2,rho,lambda1,lambda2,model=NULL){

  if (!is.numeric(lambda1)|!is.numeric(lambda2)) {stop("Arguments 'lambda1' and 'lambda2' must be real numbers")}
  if ((length(lambda1)!=1)|(length(lambda2)!=1)) {stop("Arguments 'lambda1' and 'lambda2' must be real numbers")}
  if ((abs(lambda1)+abs(lambda2))>1) {stop("The sum of the absolute values of 'lambda1' and 'lambda2' must be lower or equal than one")}
  if (is.null(model)){
    warning("Argument 'model' is missing. By default, the sine-skewed Sine density is employed, i.e., 'model='sine''.")
    model="sine"
  }


  lamdamod=function(x,y,mu1,mu2,lambda1,lambda2) (1+lambda1*sin(x-mu1)+lambda2*sin(y-mu2))
  devalmodel=eval(parse(text=paste("d",model,"(x,y,mu1,mu2,kappa1,kappa2,rho)",sep="")))


  return(devalmodel*lamdamod(x,y,mu1,mu2,lambda1,lambda2))
}


## Sine distribution (Section 3.2)

dsine=function(x,y,mu1,mu2,kappa1,kappa2,rho,tol=1e-10){

  if (!is.numeric(x)|!is.numeric(y)) {stop("Arguments 'x' and 'y' must be numeric")}
  if (!is.numeric(mu1)|!is.numeric(mu2)) {stop("Arguments 'mu1' and 'mu2' must be real numbers")}
  if ((length(mu1)!=1)|(length(mu2)!=1)) {stop("Arguments 'mu1' and 'mu2' must be real numbers")}
  if (!is.numeric(kappa1)|!is.numeric(kappa2)) {stop("Arguments 'kappa1' and 'kappa2' must be real numbers")}
  if ((length(kappa1)!=1)|(length(kappa2)!=1)) {stop("Arguments 'kappa1' and 'kappa2' must be real numbers")}
  if (!is.numeric(rho)) {stop("Arguments 'rho' must be a real number")}
  if (length(rho)!=1) {stop("Argument 'rho' must be a real number")}

  if (kappa1<0|kappa2<0) {stop("Arguments 'kappa1' and 'kappa2' must be greater or equal than zero")}

  indexes=0:99
  sCvec=0

  Cvec=choose(2*indexes,indexes) *(rho^2/(4*kappa1*kappa2))^indexes*besselI(kappa1,indexes)*besselI(kappa2,indexes)
  sCvec=sCvec+sum( Cvec )
  while(Cvec[length(Cvec)]>tol){
    indexes=indexes+100
    Cvec=choose(2*indexes,indexes) *(rho^2/(4*kappa1*kappa2))^indexes*besselI(kappa1,indexes)*besselI(kappa2,indexes)
    sCvec=sCvec+sum( Cvec )
  }

  C=4*pi^2*sCvec
  val=exp(kappa1*cos(x-mu1)+kappa2*cos(y-mu2)+rho*sin(x-mu1)*sin(y-mu2))/C
  return(val)
}


## Cosine distribution (Section 3.3)

dcosine=function(x,y,mu1,mu2,kappa1,kappa2,rho,tol=1e-10){

  if (!is.numeric(x)|!is.numeric(y)) {stop("Arguments 'x' and 'y' must be numeric")}
  if (!is.numeric(mu1)|!is.numeric(mu2)) {stop("Arguments 'mu1' and 'mu2' must be real numbers")}
  if ((length(mu1)!=1)|(length(mu2)!=1)) {stop("Arguments 'mu1' and 'mu2' must be real numbers")}
  if (!is.numeric(kappa1)|!is.numeric(kappa2)) {stop("Arguments 'kappa1' and 'kappa2' must be real numbers")}
  if ((length(kappa1)!=1)|(length(kappa2)!=1)) {stop("Arguments 'kappa1' and 'kappa2' must be real numbers")}
  if (!is.numeric(rho)) {stop("Arguments 'rho' must be a real number")}
  if (length(rho)!=1) {stop("Argument 'rho' must be a real number")}

  if (kappa1<0|kappa2<0) {stop("Arguments 'kappa1' and 'kappa2' must be greater or equal than zero")}

  indexes=1:100
  sCvec=0
  if(rho<0){
    Cvec= besselI(kappa1,indexes)*besselI(kappa2,indexes)*besselI(abs(rho),indexes)*(-1)^(indexes)
  }else{
    Cvec= besselI(kappa1,indexes)*besselI(kappa2,indexes)*besselI(rho,indexes)
  }
  sCvec=sCvec+sum( Cvec )
  while(abs(Cvec[length(Cvec)])>tol){
    indexes=indexes+100
    if(rho<0){
      Cvec= besselI(kappa1,indexes)*besselI(kappa2,indexes)*besselI(abs(rho),indexes)*(-1)^(indexes)
    }else{
      Cvec= besselI(kappa1,indexes)*besselI(kappa2,indexes)*besselI(rho,indexes)
    }
    sCvec=sCvec+sum( Cvec )
  }



  C=4*pi^2*(besselI(kappa1,0)*besselI(kappa2,0)*besselI(abs(rho),0)+2*sCvec)
  val=exp(kappa1*cos(x-mu1)+kappa2*cos(y-mu2)+rho*cos(x-mu1-y+mu2))/C
  return(val)
}




## Bivariate Wrapped Cauchy (Section 3.4)

dbwc=function(x,y,mu1,mu2,kappa1,kappa2,rho){


  if (!is.numeric(x)|!is.numeric(y)) {stop("Arguments 'x' and 'y' must be numeric")}
  if (!is.numeric(mu1)|!is.numeric(mu2)) {stop("Arguments 'mu1' and 'mu2' must be real numbers")}
  if ((length(mu1)!=1)|(length(mu2)!=1)) {stop("Arguments 'mu1' and 'mu2' must be real numbers")}
  if (!is.numeric(kappa1)|!is.numeric(kappa2)) {stop("Arguments 'kappa1' and 'kappa2' must be real numbers")}
  if ((length(kappa1)!=1)|(length(kappa2)!=1)) {stop("Arguments 'kappa1' and 'kappa2' must be real numbers")}
  if (!is.numeric(rho)) {stop("Arguments 'rho' must be a real number")}
  if (length(rho)!=1) {stop("Argument 'rho' must be a real number")}

  if (kappa1>=1|kappa2>=1|kappa1<0|kappa2<0) {stop("Arguments 'kappa1' and 'kappa2' must be in the interval [0,1)")}
  if ((rho>=1)|(rho<= (-1))) {stop("Argument 'rho' must be in the interval (-1,1)")}


  c0= (1+rho^2)*(1+kappa1^2)*(1+kappa2^2)-8*abs(rho)*kappa1*kappa2
  c1= 2*(1+rho^2)*(1+kappa2^2)*kappa1-4*abs(rho)*(1+kappa1^2)*kappa2
  c2= 2*(1+rho^2)*(1+kappa1^2)*kappa2-4*abs(rho)*(1+kappa2^2)*kappa1
  c3= -4*(1+rho^2)*kappa1*kappa2+2*abs(rho)*(1+kappa1^2)*(1+kappa2^2)
  c4= 2*rho*(1-kappa1^2)*(1-kappa2^2)
  val= (1-rho^2)*(1-kappa1^2)*(1-kappa2^2)/(4*pi^2*(c0-c1*cos(x-mu1)-c2*cos(y-mu2)-c3*cos(x-mu1)*cos(y-mu2)-c4*sin(x-mu1)*sin(y-mu2)))
  return(val)
}

