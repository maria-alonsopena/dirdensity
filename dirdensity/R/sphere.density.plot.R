

#' sph.density.plot
#'
#' Function \code{sph.density.plot} creates an interactive 3D representation of a given spherical density.
#'
#' @param x Matrix with three columns containing the points where the density is evaluated in cartesian coordinates.
#' @param fx Vector containing the values of the density. The length must coincide with the number of rows of x.
#' @param axis Logical; if TRUE, the axis are plotted.
#' @importFrom grDevices colorRampPalette
#' @importFrom misc3d  parametric3d
#' @importFrom rgl open3d points3d  text3d
#' @examples
#' # Constructing a grid on the sphere
#' N<-300
#' ele_grid<-seq(0,pi,length=N)
#' ori<-list()
#' ori[[1]]<-0
#' ori[[N]]<-0
#' for (i in 2:(N-1)){
#'   ori[[i]]<- seq(0,2*pi,length=floor(2*N*sin(ele_grid[i])))
#' }
#' lonxi<-as.numeric(lapply(ori,length))
#' both<-cbind(unlist(ori),rep(ele_grid,times=lonxi))
#' x<-DirStats::to_sph(both[,1], both[,2]) # conversion of coordinates
#' fx<-Directional::iagd(x,mu=c(1,0,0))
#' sph.density.plot(x,fx)
#' @export


sph.density.plot<-function(x, fx, axis=TRUE){

  if(!is.matrix(x)){stop("x must be a matrix")}
  if(!is.numeric(fx)){stop("fx must be numeric")}
  if(ncol(x)!=3){stop("x must have 3 columns")}
  if(nrow(x)!=length(fx)){stop("The number of rows of x but coincide with the length of fx")}
  if (axis!= TRUE & axis!= FALSE){stop("axis must be either TRUE or FALSE")}


  hdr.colors1<- colorRampPalette(c("blue","green" ),alpha=TRUE)(10)
  hdr.colors2<- colorRampPalette(c("green","yellow"),alpha=TRUE)(10)
  hdr.colors3<- colorRampPalette(c("yellow","orange"),alpha=TRUE)(10)
  hdr.colors4<- colorRampPalette(c("orange","red"),alpha=TRUE)(10)
  hdr.colors5<- colorRampPalette(c("red","darkred"),alpha=TRUE)(10)
  mycol<-c(hdr.colors1,hdr.colors2,hdr.colors3,hdr.colors4,hdr.colors5)

  interv<-seq(0,1,length=51)
  vector_int<-findInterval(fx, interv)
  vector_int[vector_int==51]<-50

  colors<-mycol[vector_int]

  zoom<-0.8
  windowRect<-c(500,50,0,0)
  windowRect[3]=windowRect[1]+256*2
  windowRect[4]=windowRect[2]+256*2

  open3d(zoom = zoom,  windowRect=windowRect)

  points3d(x[,1],x[,2],x[,3],col=colors,radius=1,size=3)

  if(axis){
    R<-1;Xaxis1<-seq(-R-0.5,R+0.5,length=100);Xaxis2<-rep(0,100);Xaxis3<-rep(0,100)
    points3d(Xaxis1,Xaxis2,Xaxis3,col="grey")
    text3d(R+0.65,0,0,texts="x",col="grey",cex=.9,font=2)

    R<-1;Yaxis1<-rep(0,100);Yaxis2<-seq(-R-0.5,R+0.5,length=100);Yaxis3<-rep(0,100)
    points3d(Yaxis1,Yaxis2,Yaxis3,col="grey")
    text3d(0,R+0.65,0,texts="y",col="grey",cex=.9,font=2)

    R<-1;Zaxis1<-rep(0,100);Zaxis2<-rep(0,100);Zaxis3<-seq(-R-0.5,R+0.5,length=100)
    points3d(Zaxis1,Zaxis2,Zaxis3,col="grey")
    text3d(0,0,R+0.65,texts="z",col="grey",cex=.9,font=2)
  }

}
