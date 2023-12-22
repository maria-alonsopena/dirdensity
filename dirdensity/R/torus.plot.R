#' torus.plot
#'
#' Function \code{torus.plot} creates an interactive 3D representation of toroidal data on the surface of a torus.
#'
#' @param x Matrix with two columns containing the data, where the number of rows is the number of observations. The first column contains the azimuthal ange, while the second column contains the polar angle (both in radians).
#' @param col.points Vector containing the mean direction of the von Mises-Fisher guide.
#' @param pch Type of points to be plotted.
#' @param size Numerical value containing the concentration of the von Mises-Fisher guide.
#' @param axis Logical; if TRUE, axis are plotted.
#' @param zoom Parameter controlling the zoom in the plot.
#' @param userMatrix A 4 by 4 matrix describing user actions to display the scene. See \code{par3d}.
#' @param windowRect A vector of four values indicating the left, top, right and bottom of the displayed window (in pixels). Applies to the whole device. See \code{par3d}.
#' @param col.torus Character string indicating the desired color of the torus.
#' @param labelsx Numeric vector containing the values to be displayed as labels for the azimuth component.
#' @importFrom misc3d  parametric3d drawScene.rgl
#' @importFrom rgl open3d points3d  text3d
#' @importFrom circular  rvonmises
#' @examples
#' library(circular)
#' a1<-rvonmises(150,0,5)
#' a2<-rvonmises(150,0,15)
#' x<-cbind(a1,a2)
#' torus.plot(x)
#' @export

torus.plot<-function(x,col.points="red",pch=16,size=4,axis = TRUE,zoom=0.7,userMatrix=NULL,
                     windowRect=NULL, col.torus="honeydew1",
                     labelsx=NULL){

  if(!is.matrix(x)){stop("x must be a matrix")}

  d <- ncol(x)

  if(d!=2){stop("Number of columns in x must be 2")}

  if (axis!= TRUE & axis!= FALSE){stop("axis must be either TRUE or FALSE")}




  xx<-x[,2]
  yy<-x[,1]


  if(is.null(userMatrix)){userMatrix<-matrix(c(0.35,-0.55,0.7,0,0.92,0.29,-0.3,0,0.13,0.78,0.57,0,0,0,0,1),ncol=4)}
  if(is.null(windowRect)){
    windowRect<-c(500,50,0,0)
    windowRect[3]=windowRect[1]+256*2
    windowRect[4]=windowRect[2]+256*2
  }

  open3d(zoom = zoom, userMatrix = userMatrix, windowRect=windowRect)
  torus <- parametric3d(fx = function(u, v) (1 + 0.5 * cos(v)) * cos(u),
                        fy = function(u, v) (1 + 0.5 * cos(v)) * sin(u),
                        fz = function(u, v) 0.5 * sin(v),
                        u = seq(0, 2 * pi, length.out = 50),
                        v = seq(0, 2 * pi, length.out = 50),
                        engine = "none",
                        color = "honeydew1",
                        alpha = 0.5)
  drawScene.rgl(torus)


  # Add data points
  xxx<-cos(xx) * (1 + 0.5 * cos(yy))
  yyy<-sin(xx) * (1 + 0.5 * cos(yy))
  zzz<-0.25 * sin(yy)
  points3d(xxx,yyy,zzz,pch=16,size=4,col=col.points)


  # Add labels

  if(!is.null(labelsx)){

    labelsx.pos<-labelsx
    labelsx<-as.character(round(labelsx.pos,2))
    nlx<-length(labelsx)
    text3d(1.55 * cos(labelsx.pos), 1.55 * sin(labelsx.pos), 0, texts=labelsx, cex=1, col=1,font=2)

  }

  if(axis==TRUE){
    R<-1.25;Xaxis1<-seq(-R-0.5,R+0.5,length=100);Xaxis2<-rep(0,100);Xaxis3<-rep(0,100)
    points3d(Xaxis1,Xaxis2,Xaxis3,col="grey")
    text3d(R+0.65,0,0,texts="x",col="grey",cex=.9,font=2)

    R<-1.25;Yaxis1<-rep(0,100);Yaxis2<-seq(-R-0.5,R+0.5,length=100);Yaxis3<-rep(0,100)
    points3d(Yaxis1,Yaxis2,Yaxis3,col="grey")
    text3d(0,R+0.65,0,texts="y",col="grey",cex=.9,font=2)

    R<-1.25;Zaxis1<-rep(0,100);Zaxis2<-rep(0,100);Zaxis3<-seq(-R-0.5,R+0.5,length=100)
    points3d(Zaxis1,Zaxis2,Zaxis3,col="grey")
    text3d(0,0,R+0.65,texts="z",col="grey",cex=.9,font=2)
  }


}



