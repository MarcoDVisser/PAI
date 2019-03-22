##' Single return lidar PAI function
##' 
##' @param(z) 1-D vector of vertical coordinates
##' @param(A) 1-D vector of scan angles
##' @param(zi) 1-D vector of discrete z-axis at which leaf area density will be computed
##' vectors is expected to have a fixed step size (dzi). 
##' @param(lad) leaf inclination model ('planophile','erectophile','spherical', or user defined)
##' @param(ai)  scan angle (optional for re-computing leaf area density with different leaf inclination)
##' @param(verbose) logical: sends signs of life?
##' @import sp
##' @import raster
##' @import raster
##' @import SearchTrees
##' @export
getPAI <- function(z,A,zi,ai,lad="spherical",verbose=FALSE){

    ## build matrix N : matrix MxS number of points per layer and scan angle
    dz <- abs(zi[1]-zi[2])
    da <- abs(ai[1]-ai[2]) 
    M <- length(zi) ## number of layers
    zi <- sort(zi,decreasing=TRUE) ## ensure zi is decreasing
    A <- abs(A) 
    S <- length(ai) ## number of unique absolute angles

    N <- array(NA,dim=c(M,S))
    t0 <- Sys.time()

    ## Height
    for(i in seq_len(M)){
        
        inci <- which(z<zi[i] & z>=zi[i]-dz)

        ## Angle
        for(j in seq_len(S)){
            ## all points at zi[i] % angle S[i]
            incj <- which(A>ai[j]&A<=ai[j]+da)
            ## calculate total passes beyond S[i]
            N[i,j] <- length(intersect(inci,incj))
        }


    t1 <- Sys.time()
    P <- round(i/M,4)
    elap <- round(t1-t0 ,3)
    eta <- round((elap/P)-elap ,3)
    unts <- attr(elap, "units")
        
        if(verbose) {
            cat("\r Code ",P*100, "% Done - ETA:",eta," ",unts,
                     " elapsed: ", elap," ",unts," \t \t \t")
        }
    }
        if(verbose)  cat("\n")

    ## now compute PAI across N's columns

    ## Get I and U matrix
    I <- array(dim=c(M+1,S))
    U <- array(dim=c(M,S))

    ## incoming at each angle
    I0 <- colSums(N)

    ## PAI for each angle
    for(i in seq_len(S)){

        I[1,i]=I0[i]
        I[2:(M+1),i] <- I0[i] - cumsum(N[,i])

        ## PAI over depth (skipping first)
        for(j in seq_len(M)){

            Nom <- -(log(I[j+1,i]/I[j,i])*cos((ai[i]+0.5*da)/180*pi))
            Denom <- Gfunction(lad=lad,ze=ai[i],zi=zi[j])*dz
            Ans <- Nom/Denom
            if(any(c(is.infinite(Ans),is.na(Ans),is.nan(Ans)))) Ans <- 0
            U[j,i] <- Ans

        }
    }

    return(list("N"=N,"incedent"=I,"PAI"=U))
    
}

##'  Get point cloud from a shapefile
##'  lidar subset returned as SpatialPointsDataFrame
##'
##' @param(xyz) matrix point cloud with x,y, and z coordinates
##' @param(qtr) quadtree of xyz
##' @param(shpfile) a SpatialPolygonsDataFrame to subset the point cloud
##' only the first polygon is used if multiple polygons are
##' @param(return.index) logical, return row index for shapefile
##' from xyz
##' present
##' @import sp
##' @import raster
##' @export
shapesXYZ <- function(xyz,qtr,shpfile,return.index=FALSE){

    if(length(shpfile@polygons)>1){
        warning("Only the first polygon selected")
    }

    xy <- apply(shpfile@polygons[[1]]@Polygons[[1]]@coords,2,range)

        tmp <- rectLookup(qtr,xlims=c(xy[1,1],xy[2,1]),
                          ylims=c(xy[1,2],xy[2,2]))
        
        if(length(tmp)>0){

            boxp <- as.data.frame.matrix(xyz[tmp,])
            if(return.index) boxp$ind <- tmp
            coordinates(boxp) <- ~ x + y
            proj4string(boxp) <- proj4string(shpfile)
            final <- boxp[shpfile,]
        } else {


            stop("location contains no points")
        }

    return(final)

}


##'  Extract index of point cloud that are within a shapefile
##'  
##' @param(xyz) matrix point cloud with x,y, and z coordinates
##' @param(qtr) quadtree of xyz
##' @param(shpfile) a SpatialPolygonsDataFrame to subset the point cloud
##' only the first polygon is used if multiple polygons are
##' present
##' @import sp
##' @import raster
##' @export
indexFromShape <- function(xyz,qtr,class="class",shpfile){

    if(length(shpfile@polygons)>1){
        warning("Only the first polygon selected")
    }

    xy <- apply(shpfile@polygons[[1]]@Polygons[[1]]@coords,2,range)

        tmp <- rectLookup(qtr,xlims=c(xy[1,1],xy[2,1]),
                          ylims=c(xy[1,2],xy[2,2]))
        
        if(length(tmp)>0){

            boxp <- as.data.frame.matrix(xyz[tmp,])
            boxp$ind <- tmp
            coordinates(boxp) <- ~ x + y
            proj4string(boxp) <- proj4string(shpfile)
            final <- boxp[shpfile,]@data$tmp
        } else {


            stop("location contains no points")
        }

    return(final)

}



################################################################################
## Translation of the matlab code for calculating plant area index from lidar
## Compute leaf area density profile from multireturn LiDAR cluod point using a
## stochastic radiative transfer model
################################################################################
## Original code by Matteo Detto
##
## Reference:
## Detto, M., G. Asner, Sonnentag, O. and H. C. Muller-Landau. 2015. Using
## stochastic radiative
## transfer models to estimate leaf area density from multireturn LiDAR in
## complex tropical forests. Journal of Geophysical Research.
##
## Princeton, 08.03.2019
## Marco Visser
################################################################################


##'  Compute Ross G function
##' @param(ze) Incoming angle in degrees
##' @param(lad) character: leaf angle distribution (uses radians)
##' @param(zi) vector of vertical coordinates (height)
##' @param(par) = optional parameters to fit the beta
##' @export
Gfunction <- function(lad="spherical",ze,zi, par){ 

    ze <- abs(ze)*pi/180
    th <- seq(0,pi/2,length.out=200)
    dth <- mean(diff(th))
    A <- cos(ze)*cos(th)
    
    J <-   suppressWarnings((1/tan(th))*(1/tan(ze)))
    use <- abs(J)<=1
    phi <- suppressWarnings(acos(J))
    A[use] <- cos(th[use])*cos(ze)*(1+(2/pi)*(tan(phi[use])-phi[use]))

    if(lad=='planophile'){
        f <- 2/pi*(1+cos(2*th))
        G <- sum(A*f)*dth
        
    } else if(lad=='erectophile'){
        
        f <- 2/pi*(1-cos(2*th))
        G <- sum(A*f)*dth
        
    } else if (lad=='spherical'){
        G=0.5
    } else if(lad=='beta') {

        G <- numeric(length(zi))
        
        for(i in seq_along(zi)){
            
            ME <- 60-(60-22)*exp(-(zi[i]/25)^2)
            SD <- 20-(20-13)*exp(-(zi[i]/25)^2)
            tx <- 2*th/pi
            dtx <- diff(tx)[1]

            tbar <- ME/90
            st <- (SD/90)^2
            s0 <- tbar*(1-tbar)
            nu <- tbar*(s0/st-1)
            mu <- (1-tbar)*(s0/st-1)
            f <- 1/beta(mu,nu)*(1-tx)^(mu-1)*tx^(nu-1)
            G[i] <- sum(A*f)*dtx
            
        }} else if(lad=='betafit') {

        G <- numeric(length(zi))

            mu <- par[1]
            nu <- par[2]
            tx <- 2*th/pi
            dtx <- diff(tx)[1]

        for(i in seq_along(zi)){
            
            f <- 1/beta(mu,nu)*(1-tx)^(mu-1)*tx^(nu-1)
            G[i] <- sum(A*f)*dtx
            
        } } else if(lad=='Wirth'){
        G <- numeric(length(zi))
        for(i in seq_len(zi)){
            lad <- (8.0252*log(zi[i]) + 16.949)*pi/180; ## mean angle by zi
            ##Campbell
            chi <- -3+(lad/9.65)^(-0.6061)
            L <- chi+asin(sqrt(1-chi^2))/sqrt(1-chi^2);
            f <- 2*chi^3*sin(th)/(L*(cos(th)^2.+chi^2*sin(th)^2)^2);
            G[i] <- sum(A*f)*dth
        }
        } else if(length(lad)>=1000){
            f <- hist(lad,th,plot=FALSE)
            G <- sum(A*f/sum(f))
        } else if (length(lad)<1000 && length(lad)>10){
            tx <- 2*th/pi
            tbar <-  mean(lad/90)
            st <- var(LIA/90)
            s0 <- tbar*(1-tbar)
            nu <- tbar*(s0/st-1)
            mu <- (1-tbar)*(s0/st-1)
            f <- 1/beta(mu,nu)*(1-t)^(mu-1)*t^(nu-1)
            G <- sum(A*f)*tx
        }

    return(G)
}




