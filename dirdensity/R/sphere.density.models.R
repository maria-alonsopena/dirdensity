#' dvMF.sph
#'
#' Function \code{dvMF.sph} computes the density function of a von Mises-Fisher density of dimension 2 and produces an interactive 3D representation of the density on the surface of a sphere.
#'
#' @param mu Vector with three components giving the mean of the density in cartesian coordinates.
#' @param kappa Parameter controlling the concentration of the density.
#' @param plot Logical; if TRUE, the 3D plot is produced.
#' @param axis Logical; if TRUE, axis are plotted.
#' @return An list containing the following components: \item{fx}{ The estimated values of the density.}
#' \item{eval.points}{The points where the estimated density was evaluated.}
#' @references  Mardia, K. V. & Jupp, P. E. (2000), Directional Statistics, Wiley, New York.
#' @importFrom DirStats  to_sph  int_sph
#' @importFrom Directional  dmixvmf
#' @examples
#' vm<-dvMF.sph(c(1,0,0),5)
#' @export

dvMF.sph<-function(mu,kappa,plot=TRUE,axis=TRUE){

  if(!is.numeric(mu) | !is.numeric(kappa)){stop("Both mu and kappa must be numeric")}

  if(length(mu)!=3){stop("mu must have length 3")}

  if(length(kappa)!=1){stop("kappa must have length 1")}

  if (plot!= TRUE & plot!= FALSE){stop("plot must be either TRUE or FALSE")}

  if (axis!= TRUE & axis!= FALSE){stop("axis must be either TRUE or FALSE")}

  N<-300
  ele_grid<-seq(0,pi,length=N)
  ori<-list()

  ori[[1]]<-0
  ori[[N]]<-0
  for (i in 2:(N-1)){
    ori[[i]]<- seq(0,2*pi,length=floor(2*N*sin(ele_grid[i])))
  }
  lonxi<-as.numeric(lapply(ori,length))
  both<-cbind(unlist(ori),rep(ele_grid,times=lonxi))
  eval.points<-to_sph(both[,1], both[,2]) # conversion of coordinates



  ftrue<-dmixvmf(eval.points, probs=1, t(mu), kappa, logden = FALSE)
  if(plot){sph.density.plot(eval.points,ftrue,axis)}

  return(list(fx=ftrue, eval.points=eval.points))

}

#' dkent.sph
#'
#' Function \code{dkent.sph} computes the density function of a Kent density of dimension 2 and produces an interactive 3D representation of the density on the surface of a sphere.
#'
#' @param mu Vector with three components giving the mean of the density in cartesian coordinates.
#' @param kappa Parameter controlling the concentration of the density.
#' @param beta Parameter controlling the ovalness of the density.
#' @param G A 3 x 3 matrix with first column equal to the mean direction. The second and third columns are the major and minor axes respectively.
#' @param plot Logical; if TRUE, the 3D plot is produced.
#' @param axis Logical; if TRUE, axis are plotted.
#' @return An list containing the following components: \item{fx}{ The estimated values of the density.}
#' \item{eval.points}{The points where the estimated density was evaluated.}
#' @importFrom DirStats to_sph int_sph
#' @importFrom Directional dkent
#' @references Kent John (1982). The Fisher-Bingham distribution on the sphere. Journal of the Royal Statistical Society, Series B, 44(1): 71–80.
#' Mardia, K. V. & Jupp, P. E. (2000), Directional Statistics, Wiley, New York.
#' @examples
#' ken<-dkent.sph(mu=c(-1,0,0),kappa=5,beta=3,G=cbind(c(-1,0,0),c(0,1,0),c(0,0,1)))
#' @export

dkent.sph<-function(mu,kappa,beta,G,plot=TRUE,axis=TRUE){


  if(!is.numeric(mu) | !is.numeric(kappa) | !is.numeric(beta)){stop("mu, kappa and beta must be numeric")}

  if(!is.matrix(G)){stop("G must be a matrix")}

  if(nrow(G)!=3){stop("G must be a matrix of dimensions 3x3")}
  if(ncol(G)!=3){stop("G must be a matrix of dimensions 3x3")}

  if(length(mu)!=3){stop("mu must have length 3")}

  if(length(kappa)!=1){stop("kappa must have length 1")}

  if(length(beta)!=1){stop("beta must have length 1")}

  if(G[1,1]!=mu[1] & G[2,1]!=mu[2] & G[3,1]!=mu[3]){stop("The first column of G must be mu")}

  if (plot!= TRUE & plot!= FALSE){stop("plot must be either TRUE or FALSE")}

  if (axis!= TRUE & axis!= FALSE){stop("axis must be either TRUE or FALSE")}

  N<-300
  ele_grid<-seq(0,pi,length=N)
  ori<-list()

  ori[[1]]<-0
  ori[[N]]<-0
  for (i in 2:(N-1)){
    ori[[i]]<- seq(0,2*pi,length=floor(2*N*sin(ele_grid[i])))
  }
  lonxi<-as.numeric(lapply(ori,length))
  both<-cbind(unlist(ori),rep(ele_grid,times=lonxi))
  eval.points<-to_sph(both[,1], both[,2]) # conversion of coordinates



  const<-DirStats::int_sph(function(x){Directional::dkent(x,G,c(kappa,beta))})

  ftrue<-Directional::dkent(eval.points,G,c(kappa,beta))/const
  if(plot){sph.density.plot(eval.points,ftrue,axis)}

  return(list(fx=ftrue, eval.points=eval.points))

}


#' dmixvMF.sph
#'
#' Function \code{dmixvMF.sph} computes the density function of a mixture of von Mises-Fisher densities, of dimension 2, and produces an interactive 3D representation of the density on the surface of a sphere.
#'
#' @param mu Matrix with three columns and number of rows equal to the number of components in the mixture, giving the mean vectors for all components.
#' @param kappa Vector of parameters controlling the concentration of each component.
#' @param probs Vector of parameters controlling the mixing probabilities.
#' @param plot Logical; if TRUE, the 3D plot is produced.
#' @param axis Logical; if TRUE, axis are plotted.
#' @return An list containing the following components: \item{fx}{ The estimated values of the density.}
#' \item{eval.points}{The points where the estimated density was evaluated.}
#' @references Kurt Hornik and Bettina Grun (2014). movMF: An R Package for Fitting Mixtures of von Mises-Fisher Distributions http://cran.r-project.org/web/packages/movMF/vignettes/movMF.pdf
#' @importFrom DirStats to_sph int_sph
#' @importFrom Directional dmixvmf
#' @examples
#' mixvm<-dmixvMF.sph(mu=matrix(c(0,0,1,0,0,-1),ncol=3,byrow=TRUE),kappa=c(5,5),probs=c(0.6,0.4))
#' @export

dmixvMF.sph<-function(mu,kappa,probs,plot=TRUE,axis=TRUE){

  if(!is.numeric(kappa) | !is.numeric(probs) ){stop("kappa and probs must be numeric")}

  if(!is.matrix(mu)){stop("mu must be a matrix")}

  if(ncol(mu)!=3){stop("mu must be a matrix with 3 columns")}

  if(length(kappa)!=length(probs)){stop("kappa and probs must have the same length")}

  if(nrow(mu)!=length(probs)){stop("The number of rows of mu must be the same as the length of kappa and probs")}

  if (plot!= TRUE & plot!= FALSE){stop("plot must be either TRUE or FALSE")}

  if (axis!= TRUE & axis!= FALSE){stop("axis must be either TRUE or FALSE")}

  N<-300
  ele_grid<-seq(0,pi,length=N)
  ori<-list()

  ori[[1]]<-0
  ori[[N]]<-0
  for (i in 2:(N-1)){
    ori[[i]]<- seq(0,2*pi,length=floor(2*N*sin(ele_grid[i])))
  }
  lonxi<-as.numeric(lapply(ori,length))
  both<-cbind(unlist(ori),rep(ele_grid,times=lonxi))
  eval.points<-to_sph(both[,1], both[,2]) # conversion of coordinates


  ftrue<-dmixvmf(eval.points, probs, mu, kappa, logden = FALSE)
  if(plot){sph.density.plot(eval.points,ftrue,axis)}

  val<-list(fx=ftrue, eval.points=eval.points)

  return(val)

}






