#' cyl.plot
#'
#' Function \code{cyl.plot} creates an interactive 3D representation of cylindrical data on the surface of a cylinder.
#'
#' @param x Matrix with two columns containing the data in cylindrical coordinates, where the number of rows is the number of observations. The first column contains the azimuthal angle (in radians), while the second column contains the height.
#' @param col.points Vector containing the mean direction of the von Mises-Fisher guide.
#' @param pch Type of points to be plotted.
#' @param size Numerical value containing the concentration of the von Mises-Fisher guide.
#' @param axis Logical; if TRUE, axis are plotted.
#' @param zoom Parameter controlling the zoom in the plot.
#' @param userMatrix A 4 by 4 matrix describing user actions to display the scene. See par3d.
#' @param windowRect A vector of four values indicating the left, top, right and bottom of the displayed window (in pixels). Applies to the whole device. See par3d.
#' @param range Vector with two components giving the interval of the real-line where height of the cylinder should be displayed.
#' @param T If plot=TRUE, parameter controlling the radius of the cylinder.
#' @param col.cyl Character string indicating the desired color of the cylinder.
#' @param labelsx Numeric vector containing the values to be displayed as labels for the circular component.
#' @param labelsy Numeric vector containing the values to be displayed as labels for the linear component.
#' @importFrom misc3d  parametric3d
#' @importFrom rgl open3d points3d  text3d
#' @importFrom circular  rvonmises
#' @examples
#'library(circular)
#' a1<-rvonmises(150,0,5)
#' a2<-runif(150,2,10)
#' x<-cbind(a1,a2)
#' cyl.plot(x)
#' @export




cyl.plot<-function(x,col.points="red",pch=16,size=4,axis=TRUE,zoom=0.9,userMatrix=NULL,
                                windowRect=NULL, range=NULL, T=NULL, col.cyl="honeydew1",
                                labelsx=NULL,labelsy=NULL){




    if(!is.matrix(x)){stop("x must be a matrix")}

    d <- ncol(x)

    if(d!=2){stop("Number of columns in x must be 2")}

    if (axis!= TRUE & axis!= FALSE){stop("axis must be either TRUE or FALSE")}



  xx<-x[,1]
  yy<-x[,2]


  if(is.null(userMatrix)){userMatrix<-matrix(c(0.043,0.88,0.37,0,-0.78,0.23,-0.59,0,-0.74,-0.2,0.6,0,0,0,0,1),ncol=4)}
  if(is.null(windowRect)){
    windowRect<-c(500,50,0,0)
    windowRect[3]=windowRect[1]+256*2
    windowRect[4]=windowRect[2]+256*2
  }
  if(is.null(range)){range<-c(min(yy),max(yy))}
  if(is.null(T)){T<-diff(range)/2}

  open3d(zoom = zoom, userMatrix = userMatrix, windowRect=windowRect)
  cylinder <- parametric3d(fx = function(u, v) v,
                           fy = function(u, v) T * cos(u),
                           fz = function(u, v) T * sin(u),
                           u = seq(0, 2 * pi, length.out = 80),
                           v = seq(range[2], range[1], length.out = 50),
                           color = col.cyl,
                           alpha = 0.5)

  # Add data points
  points3d(yy,T*cos(xx),T*sin(xx),pch=pch,size=size,col=col.points)


  # Add labels
  if(!is.null(labelsy)){

    labelsy.pos<-labelsy
    labelsy<-as.character(round(labelsy.pos,2))
    nly<-length(labelsy)
    text3d(labelsy.pos,rep(T-T/5,nly),rep(T-T/5,nly),texts=labelsy,cex=.9,col=1,font=2)

  }

  if(!is.null(labelsx)){

    labelsx.pos<-labelsx
    labelsx<-as.character(round(labelsx.pos,2))
    nlx<-length(labelsx)
    text3d(rep(range[1]-T/10,nlx),T*cos(labelsx.pos),T*sin(labelsx.pos),texts=labelsx,col=1,cex=.9,font=2)


  }


  if(axis==TRUE){
    Xaxis1<-seq(-T-0.5,T+0.5,length=100);Xaxis2<-rep(0,100);Xaxis3<-rep(0,100)
    points3d(Xaxis1,Xaxis2,Xaxis3,col="grey")
    text3d(T+0.65,0,0,texts="x",col="grey",cex=.9,font=2)

    Yaxis1<-rep(0,100);Yaxis2<-seq(-T-0.5,T+0.5,length=100);Yaxis3<-rep(0,100)
    points3d(Yaxis1,Yaxis2,Yaxis3,col="grey")
    text3d(0,T+0.65,0,texts="y",col="grey",cex=.9,font=2)

    Zaxis1<-rep(0,100);Zaxis2<-rep(0,100);Zaxis3<-seq(-T-0.5,T+0.5,length=100)
    points3d(Zaxis1,Zaxis2,Zaxis3,col="grey")
    text3d(0,0,T+0.65,texts="z",col="grey",cex=.9,font=2)
  }

}


