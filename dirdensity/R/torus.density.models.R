#' dsinsk.torus
#'
#' Function \code{dsinsk.torus} computes the density function of a sine-skewed toroidal distribution and produces an interactive 3D representation of the density on the surface of a cylinder.
#'
#' @param mu Vector with two components representing the mean direction of each marginal model.
#' @param kappa Vector with two components representing the concentration of each marginal model.
#' @param rho Numeric element containing the correlation between the two components.
#' @param lambda Vectorwith two components representing the skewness coefficients.
#' @param model Character value indicating a predefined toroidal symmetric density.  Implemented examples include model="sine", model="cosine", model="bwc".
#' @param plot Logical; if TRUE, the 3D plot is produced.
#' @param axis Logical; if TRUE, axis are plotted.
#' @return An list containing the following components: \item{fx}{ The estimated values of the density.}
#' \item{eval.points}{The points where the estimated density was evaluated.}
#' @references
#' Jose Ameijeiras-Alonso, Christophe Ley, Sine-skewed toroidal distributions and their application in protein bioinformatics, Biostatistics, Volume 23, Issue 3, July 2022, Pages 685–704.
#' @examples
#' ssk<-dsinsk.torus(mu=c(-pi/2,pi/8),kappa=c(5,5),rho=3,lambda=c(-0.9,0.1),model="sine")
#' @export

dsinsk.torus<-function(mu,kappa,rho,lambda,model=NULL,plot=TRUE,axis=TRUE){


  if(!is.numeric(mu) | !is.numeric(kappa) | !is.numeric(rho) | !is.numeric(lambda)){stop("Both mu and k must be numeric")}

  if(length(mu)!=2 | length(kappa)!=2 | length(lambda)!=2){stop("mu, kappa and lambda must have length 2")}

  if (plot!= TRUE & plot!= FALSE){stop("plot must be either TRUE or FALSE")}

  if (axis!= TRUE & axis!= FALSE){stop("axis must be either TRUE or FALSE")}

  mu1<-mu[1]
  mu2<-mu[2]

  kappa1<-kappa[1]
  kappa2<-kappa[2]

  lambda1<-lambda[1]
  lambda2<-lambda[2]

  N<-600

  u = seq(0, 2 * pi, length.out = N)
  v = seq(0, 2 * pi, length.out = N)

  ll<-expand.grid(u,v)
  u2<-ll[,1]
  v2<-ll[,2]


  ftrue<-dsinesk(u2, v2,mu1,mu2,kappa1,kappa2,rho,lambda1,lambda2,model)

  eval.points<-cbind((1 + 0.5 * cos(v2)) * cos(u2), (1 + 0.5 * cos(v2)) * sin(u2), 0.5 * sin(v2))


  if(plot){
    torus.density.plot(eval.points,ftrue,axis)
  }
  
  return(list(fx=ftrue, eval.points=eval.points))
}


#' dsin.torus
#'
#' Function \code{dsin.torus} computes the density function of a bivariate sine model and produces an interactive 3D representation of the density on the surface of a cylinder.
#'
#' @param mu Vector with two components representing the mean direction of each marginal model.
#' @param kappa Vector with two components representing the concentration of each marginal model.
#' @param rho Numeric element containing the correlation between the two components.
#' @param plot Logical; if TRUE, the 3D plot is produced.
#' @return An list containing the following components: \item{fx}{ The estimated values of the density.}
#' \item{eval.points}{The points where the estimated density was evaluated.}
#' @references Singh, H., Hnizdo, V. & Demchuk, E. (2002), ‘Probabilistic model for two dependent circular variables’, Biometrika 89, 719–723.
#' @param axis Logical; if TRUE, axis are plotted.
#' @examples
#' si<-dsin.torus(mu=c(-pi/16,pi/4), kappa=c(4,5), rho=3)
#' @export

dsin.torus<-function(mu,kappa,rho,plot=TRUE,axis=TRUE){


  if(!is.numeric(mu) | !is.numeric(kappa) | !is.numeric(rho) ){stop("Both mu and k must be numeric")}

  if(length(mu)!=2 | length(kappa)!=2 ){stop("mu and kappa must have length 2")}

  if (plot!= TRUE & plot!= FALSE){stop("plot must be either TRUE or FALSE")}

  if (axis!= TRUE & axis!= FALSE){stop("axis must be either TRUE or FALSE")}

  mu1<-mu[1]
  mu2<-mu[2]

  kappa1<-kappa[1]
  kappa2<-kappa[2]

  N<-600

  u = seq(0, 2 * pi, length.out = N)
  v = seq(0, 2 * pi, length.out = N)

  ll<-expand.grid(u,v)
  u2<-ll[,1]
  v2<-ll[,2]


  ftrue<-dsine(u2, v2, mu1, mu2, kappa1, kappa2, rho)

  eval.points<-cbind((1 + 0.5 * cos(v2)) * cos(u2), (1 + 0.5 * cos(v2)) * sin(u2), 0.5 * sin(v2))

  if(plot){
  torus.density.plot(eval.points,ftrue,axis)
  }
  return(list(fx=ftrue, eval.points=eval.points))

}



#' dcos.torus
#'
#' Function \code{dcos.torus} computes the density function of a bivariate cosine model and produces an interactive 3D representation of the density on the surface of a cylinder.
#'
#' @param mu Vector with two components representing the mean direction of each marginal model.
#' @param kappa Vector with two components representing the concentration of each marginal model.
#' @param rho Numeric element containing the correlation between the two components.
#' @param plot Logical; if TRUE, the 3D plot is produced.
#' @param axis Logical; if TRUE, axis are plotted.
#' @return An list containing the following components: \item{fx}{ The estimated values of the density.}
#' \item{eval.points}{The points where the estimated density was evaluated.}
#' @references Mardia, K. V., Taylor, C. C. & Subramaniam, G. K. (2007), ‘Protein bioinformatics and mixtures of bivariate von Mises distributions for angular data’, Biometrics 63, 505–512.
#' @examples
#' co<-dcos.torus(mu=c(-pi/16,pi/4), kappa=c(4,5), rho=3)
#' @export

dcos.torus<-function(mu,kappa,rho,plot=TRUE,axis=TRUE){


  if(!is.numeric(mu) | !is.numeric(kappa) | !is.numeric(rho) ){stop("Both mu and k must be numeric")}

  if(length(mu)!=2 | length(kappa)!=2 ){stop("mu and kappa must have length 2")}

  if (plot!= TRUE & plot!= FALSE){stop("plot must be either TRUE or FALSE")}

  if (axis!= TRUE & axis!= FALSE){stop("axis must be either TRUE or FALSE")}


  mu1<-mu[1]
  mu2<-mu[2]

  kappa1<-kappa[1]
  kappa2<-kappa[2]

  N<-600

  u = seq(0, 2 * pi, length.out = N)
  v = seq(0, 2 * pi, length.out = N)

  ll<-expand.grid(u,v)
  u2<-ll[,1]
  v2<-ll[,2]


  ftrue<-dcosine(u2,v2,mu1,mu2,kappa1,kappa2,rho)

  eval.points<-cbind((1 + 0.5 * cos(v2)) * cos(u2), (1 + 0.5 * cos(v2)) * sin(u2), 0.5 * sin(v2))


  if(plot){
    torus.density.plot(eval.points,ftrue,axis)
  }
  return(list(fx=ftrue, eval.points=eval.points))

}


#' dwc.torus
#'
#' Function \code{dcos.torus} computes the density function of a bivariate wrapped Cauchy and produces an interactive 3D representation of the density on the surface of a cylinder.
#'
#' @param mu Vector with two components representing the mean direction of each marginal model.
#' @param kappa Vector with two components representing the concentration of each marginal model.
#' @param rho Numeric element containing the correlation between the two components.
#' @param plot Logical; if TRUE, the 3D plot is produced.
#' @param axis Logical; if TRUE, axis are plotted.
#' @return An list containing the following components: \item{fx}{ The estimated values of the density.}
#' \item{eval.points}{The points where the estimated density was evaluated.}
#' @references Kato, S. & Pewsey, A. (2015), ‘A Möbius transformation-induced distribution on the torus’, Biometrika 102, 359–370
#' @examples
#' wc<-dwc.torus(mu=c(pi/2,0), kappa=c(0.5,0.15), rho=0.7)
#' @export

dwc.torus<-function(mu,kappa,rho,plot=TRUE,axis=TRUE){

  if(!is.numeric(mu) | !is.numeric(kappa) | !is.numeric(rho) ){stop("Both mu and k must be numeric")}

  if(length(mu)!=2 | length(kappa)!=2 ){stop("mu and kappa must have length 2")}

  if (plot!= TRUE & plot!= FALSE){stop("plot must be either TRUE or FALSE")}

  if (axis!= TRUE & axis!= FALSE){stop("axis must be either TRUE or FALSE")}


  mu1<-mu[1]
  mu2<-mu[2]

  kappa1<-kappa[1]
  kappa2<-kappa[2]

  N<-600

  u = seq(0, 2 * pi, length.out = N)
  v = seq(0, 2 * pi, length.out = N)

  ll<-expand.grid(u,v)
  u2<-ll[,1]
  v2<-ll[,2]


  ftrue<-dbwc(u2, v2, mu1, mu2, kappa1, kappa2, rho)

  eval.points<-cbind((1 + 0.5 * cos(v2)) * cos(u2), (1 + 0.5 * cos(v2)) * sin(u2), 0.5 * sin(v2))

  if(plot){

    torus.density.plot(eval.points,ftrue,axis)

  }

  return(list(fx=ftrue, eval.points=eval.points))

}










