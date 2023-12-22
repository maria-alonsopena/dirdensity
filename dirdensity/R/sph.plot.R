#' sph.plot
#'
#' Function \code{sph.plot} creates an interactive 3D representation of spherical data on a sphere.
#'
#' @param x Matrix containing the data in cartesian coordinates, where the number of columns must be 3 and the number of rows is the number of observations..
#' @param col.points Vector containing the mean direction of the von Mises-Fisher guide.
#' @param pch Type of points to be plotted.
#' @param size Numerical value containing the concentration of the von Mises-Fisher guide.
#' @param axis Logical; if TRUE, the axis are plotted.
#' @param zoom Parameter controlling the zoom in the plot.
#' @param windowRect aa.
#' @param col.sph Character string indicating the desired color of the sphere.
#' @importFrom misc3d  parametric3d
#' @importFrom rgl open3d points3d  text3d
#' @importFrom movMF rmovMF
#' @examples
#'  library(movMF)
#' n<-200
#' mu<-matrix(c(0,0,1,0,0,-1),ncol=3,byrow=TRUE)
#' k<-c(7,2)
#' probs<-c(0.85,0.15)
#' x<-rmovMF(n,k*mu,alpha=probs)
#' sph.plot(x,2,4)
#' @export



sph.plot<-function(x, col.points = "red", pch=16, size=4, axis = TRUE, zoom=0.8,
                   windowRect=NULL, col.sph="honeydew1"){

  if(!is.matrix(x)){stop("x must be a matrix")}

  d <- ncol(x) - 1

  if(d!=2){stop("Number of columns in x must be 3")}

  if (axis!= TRUE & axis!= FALSE){stop("axis must be either TRUE or FALSE")}


  if (is.null(windowRect)){
    windowRect<-c(500,50,0,0)
    windowRect[3]=windowRect[1]+256*2
    windowRect[4]=windowRect[2]+256*2
  }


  open3d(zoom = zoom,  windowRect=windowRect)

  sphere <- parametric3d(fx = function(u, v) cos(u)*sin(v),
                         fy = function(u, v) sin(u)*sin(v),
                         fz = function(u, v)cos(v),
                         u = seq(0, 2 * pi, length.out = 90),
                         v = seq(0, pi, length.out = 90),
                         color = col.sph,
                         alpha = 0.5)

  points3d(x[,1],x[,2],x[,3],col=col.points,pch=16,radius=1,size=size)

  if(axis==TRUE){
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

