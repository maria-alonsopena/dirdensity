#' torus.density.plot
#'
#' Function \code{torus.density.plot} creates an interactive 3D representation of a given toroidal density.
#'
#' @param x Matrix with three columns containing the points where the density is evaluated in cartesian coordinates.
#' @param fx Vector containing the values of the density. The length must coincide with the number of rows of x.
#' @param axis Logical; if TRUE, the axis are plotted.
#' @importFrom grDevices colorRampPalette
#' @importFrom misc3d  parametric3d
#' @importFrom rgl open3d points3d  text3d
#' @examples
#' N<-600
#' u = seq(0, 2 * pi, length.out = N)
#' v = seq(0, 2 * pi, length.out = N)
#' ll<-expand.grid(u,v)
#' u2<-ll[,1]
#' v2<-ll[,2]
#' x<-cbind((1 + 0.5 * cos(v2)) * cos(u2), (1 + 0.5 * cos(v2)) * sin(u2), 0.5 * sin(v2))
#' fx<-dirdensity:::dbwc(u2,v2,mu1=pi,mu2=0,kappa1=0.5,kappa2=0.15,rho=0.5)
#' torus.density.plot(x, fx)
#' @export

torus.density.plot<-function(x, fx, axis=TRUE){

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

