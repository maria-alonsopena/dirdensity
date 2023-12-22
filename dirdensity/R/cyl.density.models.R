#' dMS.cyl
#'
#' Function \code{dMS.cyl} computes the density function of a Mardia-Sutton density and produces an interactive 3D representation of the density on the surface of a cylinder.
#'
#' @param mu Vector with two components: the first one is the circular mean and the second one is the linear mean.
#' @param kappa Concentration parameter (greater or equal than zero).
#' @param sigma Standard deviation (has to be positive).
#' @param nu Parameter between -pi and pi.
#' @param lambda Structure parameter; has to be greater or equal than zero
#' @param rangex Vector with two components giving the interval of the real-line where linear part of the density should be evaluated.
#' @param T If plot=TRUE, parameter controlling the radius of the cylinder.
#' @param plot Logical; if TRUE, the 3D plot is produced.
#' @param axis Logical; if TRUE, axis are plotted.
#' @param orientation If plot=TRUE, character string indicating the orientation of the cylinder. Must be "horiz" or "vert".
#' @return An list containing the following components: \item{fx}{ The estimated values of the density.}
#' \item{eval.points}{The points where the estimated density was evaluated.}
#' @references Mardia, K. V. & Sutton, T. W. (1978), ‘A model for cylindrical variables with applications’, J. Roy. Stat. Soc. Ser. B 40, 229–233.
#' @examples
#' ms<-dMS.cyl(mu=c(pi,0),kappa=5,sigma=0.5,nu=1,lambda=3,rangex=c(-1,6),orientation="vert")
#' @export

dMS.cyl<-function(mu,kappa,sigma,nu,lambda,rangex=NULL,T=NULL,plot=TRUE,axis=TRUE,orientation="horiz"){

  if(!is.numeric(mu) | !is.numeric(kappa) | !is.numeric(sigma) | !is.numeric(nu) | !is.numeric(lambda)){stop("Parameters mu, kappa, sigma, nu and lambda must be numeric")}

  if( length(mu)!=2 ){stop("mu must have length 2")}

  if( length(kappa)!=1 | length(sigma)!=1 | length(nu)!=1 |length(lambda)!=1){stop("kappa, sigma, nu and lambda must have length 1")}

  if(kappa<0){stop("kappa must be greater or equal than 0")}

  if(sigma<=0){stop("sigma must be positive")}

  if(lambda<0){stop("lambda must be greater or equal than 0")}

  if (plot!= TRUE & plot!= FALSE){stop("plot must be either TRUE or FALSE")}

  if (axis!= TRUE & axis!= FALSE){stop("axis must be either TRUE or FALSE")}

  if(orientation!= "horiz" & orientation!= "vert"){stop("orientation must be either 'horiz' or 'vert'")}

  N<-600

  u = seq(0, 2 * pi, length.out = N)
  if(is.null(rangex)){
    v = seq(mu[2]-8*sigma, mu[2]+8*sigma, length.out = N)
  }else{
    if(!is.numeric(rangex)){
      stop("rangex must be numeric")
    }
    if(length(rangex)!=2){
      stop("length of rangex must be 2")
    }
    v = seq(rangex[1],rangex[2], length.out = N)
  }


  # The grid depends on the support of the model...

  ll<-expand.grid(u,v)
  u2<-ll[,1]
  v2<-ll[,2]




  ftrue<-dms(u2, v2, mu, kappa, sigma, nu, lambda)


  if(is.null(T)){T<-diff(range(v))/4}
  if(orientation=="horiz"){
    eval.points<-cbind(v2,  T * cos(u2), T * sin(u2))
  }else{
    eval.points<-cbind(T * cos(u2), T * sin(u2) ,v2 )
  }

  if(plot){
  cyl.density.plot(eval.points,ftrue,axis)
  }

  return(list(fx=ftrue, eval.points=eval.points))

}



#' dJW.cyl
#'
#' Function \code{dJW.cyl} computes the density function of a Jones-Wehrly density and produces an interactive 3D representation of the density on the surface of a cylinder.
#'
#' @param mu Circular location parameter.
#' @param lambda Linear dispersion (positive).
#' @param kappa Dependence parameter (has to be less than lambda).
#' @param rangex Vector with two components giving the interval of the real-line where linear part of the density should be evaluated.
#' @param T If plot=TRUE, parameter controlling the radius of the cylinder.
#' @param plot Logical; if TRUE, the 3D plot is produced.
#' @param axis Logical; if TRUE, axis are plotted.
#' @param orientation If plot=TRUE, character string indicating the orientation of the cylinder. Must be "horiz" or "vert".
#' @return An list containing the following components: \item{fx}{ The estimated values of the density.}
#' \item{eval.points}{The points where the estimated density was evaluated.}
#' @references Johnson, R. A. & Wehrly, T. E. (1978), ‘Some angular-linear distributions and related regression models’, J. Amer. Statist. Assoc. 73, 602–606.
#' @examples
#' jw<-dJW.cyl(mu=pi/2,lambda=2,kappa=1.25,orientation="horiz",rangex=c(0,5))
#' @export

dJW.cyl<-function(mu,lambda,kappa,rangex,T=NULL,plot=TRUE,axis=TRUE,orientation="horiz"){



  if(lambda<=0){stop("lambda must be positive")}

  if(kappa>=lambda){stop("kappa must be less than lambda")}


  if(!is.numeric(mu) | !is.numeric(kappa) | !is.numeric(lambda)){stop("Parameters mu, kappa and lambda must be numeric")}

  if(length(mu)!=1 | length(kappa)!=1 |length(lambda)!=1){stop("mu, kappa and lambda must have length 1")}



  if (plot!= TRUE & plot!= FALSE){stop("plot must be either TRUE or FALSE")}

  if (axis!= TRUE & axis!= FALSE){stop("axis must be either TRUE or FALSE")}

  if(orientation!= "horiz" & orientation!= "vert"){stop("orientation must be either 'horiz' or 'vert'")}

  if(!is.numeric(rangex)){
    stop("rangex must be numeric")
  }
  if(length(rangex)!=2){
    stop("length of rangex must be 2")
  }
  N<-600

  u = seq(0, 2 * pi, length.out = N)
  v = seq(rangex[1],rangex[2], length.out = N)



  ll<-expand.grid(u,v)
  u2<-ll[,1]
  v2<-ll[,2]




  ftrue<-djw(u2, v2, mu,lambda,kappa)


  if(is.null(T)){T<-diff(range(v))/4}
  if(orientation=="horiz"){
    eval.points<-cbind(v2,  T * cos(u2), T * sin(u2))
  }else{
    eval.points<-cbind(T * cos(u2), T * sin(u2) ,v2 )
  }

  if(plot){
  cyl.density.plot(eval.points,ftrue,axis)
  }

  return(list(fx=ftrue, eval.points=eval.points))

}


#' dAL.cyl
#'
#' Function \code{dAL.cyl} computes the density function of a Abe-Ley density and produces an interactive 3D representation of the density on the surface of a cylinder.
#'
#' @param mu Circular location parameter.
#' @param lambda Skewness parameter (must belong to (-1,1)).
#' @param beta Linear dispersion (has to be positive).
#' @param alpha Shape parameter (has to be positive).
#' @param kappa Circular concentration, which acts also as dependence parameter (greater or equal than zero).
#' @param rangex Vector with two components giving the interval of the real-line where linear part of the density should be evaluated.
#' @param T If plot=TRUE, parameter controlling the radius of the cylinder.
#' @param plot Logical; if TRUE, the 3D plot is produced.
#' @param axis Logical; if TRUE, axis are plotted.
#' @param orientation If plot=TRUE, character string indicating the orientation of the cylinder. Must be "horiz" or "vert".
#' @return An list containing the following components: \item{fx}{ The estimated values of the density.}
#' \item{eval.points}{The points where the estimated density was evaluated.}
#' @references Abe, T. and Ley, C. (2017) 'A tractable, parsimonious and flexible model for cylindrical data, with applications', Econometrics and Statistics, 4, 91-104.
#' @examples
#' al<-dAL.cyl(mu=0,lambda=-0.5,kappa=1,beta=2,alpha=3,orientation="vert",rangex=c(0,3),T=1.5)
#' @export


dAL.cyl<-function(mu,lambda,beta,alpha,kappa,rangex,T=NULL,plot=TRUE,axis=TRUE,orientation="horiz"){


  if(!is.numeric(mu) | !is.numeric(kappa) | !is.numeric(lambda) | !is.numeric(beta) | !is.numeric(alpha)){stop("Parameters mu, kappa, beta, alpha and lambda must be numeric")}

  if(length(mu)!=1 | length(kappa)!=1 | length(beta)!=1 | length(alpha)!=1 |length(lambda)!=1){stop("mu, kappa, beta, alpha and lambda must have length 1")}

  if(lambda<=-1 | lambda >=1){stop("lambda must belong to (-1,1)")}

  if(beta<=0){stop("beta must be positive")}

  if(alpha<=0){stop("alpha must be positive")}

  if(kappa<0){stop("kappa must be greater or equal than 0")}

  if (plot!= TRUE & plot!= FALSE){stop("plot must be either TRUE or FALSE")}

  if (axis!= TRUE & axis!= FALSE){stop("axis must be either TRUE or FALSE")}

  if(orientation!= "horiz" & orientation!= "vert"){stop("orientation must be either 'horiz' or 'vert'")}


  if(!is.numeric(rangex)){
    stop("rangex must be numeric")
  }
  if(length(rangex)!=2){
    stop("length of rangex must be 2")
  }

  N<-600

  u = seq(0, 2 * pi, length.out = N)
  v = seq(rangex[1],rangex[2], length.out = N)


  ll<-expand.grid(u,v)
  u2<-ll[,1]
  v2<-ll[,2]


  ftrue<-dal(u2, v2, mu, lambda, beta, alpha, kappa)


  if(is.null(T)){T<-diff(range(v))/4}
  if(orientation=="horiz"){
    eval.points<-cbind(v2,  T * cos(u2), T * sin(u2))
  }else{
    eval.points<-cbind(T * cos(u2), T * sin(u2) ,v2 )
  }
  if(plot){
  cyl.density.plot(eval.points,ftrue,axis)
  }

  return(list(fx=ftrue, eval.points=eval.points))

}




