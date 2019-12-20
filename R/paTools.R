##
##################################################################################
## ofm version 1.0
##################################################################################
##
##################################################################################
## Required functions ##
##################################################################################
##
## the run function
run <- function(data = NULL,
                run.model = "FP", 
                grid.input = NULL,
                cov.model = "Exp",
                neighbourhood.points = 50,
                random.seed=1234, # random seed for inside cross-validations
                ...) {
  #				
  set.seed(random.seed)
  start.time <- proc.time()[3]
  #
  if(is.null(data)){ stop("Input a data frame with n rows and 4 columns. \n Col:(1) Easting  (or x axis)\n Col:(2) Northing  (or y axis)\n Col:(3) Observations\n Col:(4) Treatment labels.")  }
  #
  if(!run.model %in% c("AP", "FP")){ stop("\n Define model correctly: 'AP': local adaptive based on points, 'FP': local fixed points. \n")  }
  #
  neighbourhood.size = NULL;
  pts.for.each.treat = 20; # this is for local approach to ensure the number of points in each treatment, default is 20
  subsample.percentage = 10; # for kmeans and bootstrap
  num.of.replications = 100; # for kmeans and bootstrap
  tol.site = 1000;
  #
  if(run.model=="G"){ type <- "global" }
  else if(run.model=="AP"){ type <- "local-processed-p" } # local adaptive based on points cross-val,
  else if(run.model=="FP"){ type <- "local-fixed" }
  else { stop("Correctly define input for 'run.model'\n")}
  #
  if(is.null(grid.input)){  grid.size <- c(10,10); grid.coords <- grid.input; shapefile <- NULL  }
  else if(is.vector(grid.input)){ grid.size <- grid.input; grid.coords <- grid.input; shapefile <- NULL  }
  else if(is.matrix(grid.input) | is.data.frame(grid.input)){ grid.size <- NULL; grid.coords <- grid.input; shapefile <- NULL  }
  else if(class(grid.input)=="SpatialPolygonsDataFrame"){ grid.size <- NULL; grid.coords <- grid.input; shapefile <- grid.input  }
  else{  stop("Correctly define 'grid.input' option. It can take (i) Grid size, e.g., c(10,10); \n (ii) Grid coordinates of easting and northing positions using n x 2 matrix or data frame \n (iii) Shape file loaded to R platform. \n")  }
  #
  out <- coKrig(data = data, type = type, grid.size = grid.size, grid.coords = grid.coords,  shapefile = shapefile,  model = cov.model,
                neighbourhood.size = neighbourhood.size, neighbourhood.points = neighbourhood.points, tol.site = tol.site, 
                sub.replicate = num.of.replications, sub.percentage = subsample.percentage, ...)
  #
  mat <- .z.test.diff(out=out,data=data) 
  #
  out$kriged.output <- mat
  out$run.model <- run.model
  #
  end.time <- proc.time()[3]
  t <- end.time-start.time
  out$comp.time <- .fnc.time_(t)
  #
  class(out) <- "pa"
  out
}
## calculate z stat 
.z.test.diff <- function(out,data){
  ck <- as.matrix(out$kriged.output)
  tr <- sort(unique(data[,4]))
  n <- table(as.numeric(data[,4]))
  if(ncol(ck)==9){
    mat <- matrix(NA,nrow(ck),6+6+3); 
    mat[,1:6] <- ck[,1:6]; mat[,7] <- ck[,1]-ck[,3]; mat[,8] <- abs(ck[,2] + ck[,4] - 2*ck[,7]); 
    mat[,9] <- ck[,1]-ck[,5]; mat[,10] <- abs(ck[,2] + ck[,6] - 2*ck[,8]); 
    mat[,11] <- ck[,3]-ck[,5]; mat[,12] <- abs(ck[,4] + ck[,6] - 2*ck[,9]);
    mat[,13] <- mat[,7]/sqrt(mat[,8]); mat[,14] <- mat[,9]/sqrt(mat[,10]); mat[,15] <- mat[,11]/sqrt(mat[,12]) 	# z test
    colnames(mat) <- c(paste0("tr_",tr[1]),paste0("tr_",tr[1],"_var"),paste0("tr_",tr[2]),paste0("tr_",tr[2],"_var"),paste0("tr_",tr[3]),paste0("tr_",tr[3],"_var"),
                       paste0("tr_diff_",tr[1],"_",tr[2]),paste0("tr_diff_",tr[1],"_",tr[2],"_cov"),paste0("tr_diff_",tr[1],"_",tr[3]),paste0("tr_diff_",tr[1],"_",tr[3],"_cov"),	
                       paste0("tr_diff_",tr[2],"_",tr[3]),paste0("tr_diff_",tr[2],"_",tr[3],"_cov"),paste0("z_",tr[1],"_",tr[2]),paste0("z_",tr[1],"_",tr[3]),paste0("z_",tr[2],"_",tr[3]))
  }
  else if(ncol(ck)==5){
    mat <- matrix(NA,nrow(ck),4+2+1); mat[,1:4] <- ck[,1:4]; mat[,5] <- ck[,1]-ck[,3]; 
    mat[,6] <- abs(ck[,2] + ck[,4] - 2*ck[,5]); 
    mat[,7] <- mat[,5]/sqrt(mat[,6]); # z for 12
    colnames(mat) <- c(paste0("tr_",tr[1]),paste0("tr_",tr[1],"_var"),paste0("tr_",tr[2]),paste0("tr_",tr[2],"_var"),
                       paste0("tr_diff_",tr[1],"_",tr[2]),paste0("tr_diff_",tr[1],"_",tr[2],"_cov"),paste0("z_",tr[1],"_",tr[2]))
  }
  else if(ncol(ck)==2){
    mat <- ck
    colnames(mat) <- c(paste0("tr_",tr[1]),paste0("tr_",tr[1],"_var")) 
  }
  else{
    stop("Error")
  }
  rm(ck); rm(tr);
  mat
}
## calculate the p value
# 1-tail
.p.val.1tailed <- function(x, side="+ve",sig=0.05){
  z.levels <- c(2.3264,1.6449,1.2816)
  p.levels <- c(0.01,0.05,0.1)
  if(side=="+ve"){
    z.levels <- z.levels  
    if (sig==0.01){
      p <- rep(0,length(x))
      p[x>z.levels[1]] <- 1 # reject H0: mu1=m2
    }
    else if (sig==0.05){
      p <- rep(0,length(x))
      p[x>z.levels[2]] <- 1 # reject H0: mu1=m2
    }
    else if (sig==0.1){
      p <- rep(0,length(x))
      p[x>z.levels[3]] <- 1 # reject H0: mu1=m2
    }
    else{
      stop("Takes one of the sig values: c(0.01,0.05,0.1) ")
    }
  }
  else if(side=="-ve"){
    z.levels <- -z.levels  
    if (sig==0.01){
      p <- rep(0,length(x))
      p[x < z.levels[1]] <- 1 # reject H0: mu1=m2
    }
    else if (sig==0.05){
      p <- rep(0,length(x))
      p[x < z.levels[2]] <- 1 # reject H0: mu1=m2
    }
    else if (sig==0.1){
      p <- rep(0,length(x))
      p[x < z.levels[3]] <- 1 # reject H0: mu1=m2
    }
    else{
      stop("Takes one of the sig values: c(0.01,0.05,0.1) ")
    }
  }
  else{
    stop("define argument 'side' properly, can take values '+ve' and '-ve'.")
  }
  #
  p
}
# 2-tailed
.p.val.2tailed <- function(x, sig=0.05){
  x <- abs(x)
  z.levels <- c(3.291,2.576,1.960,1.645,1.281)
  p.levels <- c(0.001,0.01,0.05,0.1,0.2)
  if(sig==0.001){
    p <- rep(0,length(x))
    p[x>3.291] <- 1 # reject H0: mu1=m2
  }
  else if (sig==0.01){
    p <- rep(0,length(x))
    p[x>2.576] <- 1 # reject H0: mu1=m2
  }
  else if (sig==0.05){
    p <- rep(0,length(x))
    p[x>1.960] <- 1 # reject H0: mu1=m2
  }
  else if (sig==0.1){
    p <- rep(0,length(x))
    p[x>1.645] <- 1 # reject H0: mu1=m2
  }
  else if (sig==0.2){
    p <- rep(0,length(x))
    p[x>1.281] <- 1 # reject H0: mu1=m2
  }
  else{
    stop("Takes one of the sig values: c(0.001,0.01,0.05,0.1,0.2) ")
  }
  p
}
## p-val related code




##
set.seed(1234)
# parallel function
coKrig <- 
  function(data = NULL,
           type = "global",
           grid.size = NULL,
           grid.coords = NULL,
           shapefile = NULL,
           model = "Exp",
           neighbourhood.size = NULL,
           neighbourhood.points = NULL,
           pts.for.each.treat = 20, # this is for local approach to ensure the number of points in each treatment, default is 20
           tol.site = 1000,
           sub.replicate = 100,
           num.of.clusters = NULL, 
           sub.percentage = 10,
           parallel = FALSE,
           numCores = NULL,
           ...) {
    #
    #
    # type = "global", "local-adaptive", "local-processed", "local-fixed","kmeans", "boot", "epoch"
    #        "local-processed" => using gam
    #        "local-fixed" => using fixed number of neighbourhood points for each treatment
    #        "epoch" => based on NN stochastic gradient descent, sampling without replacement 
    #
    # local.tuning = only applicable for local-processed to define the number of kmeans sample
    #                varies between 0 to inf; 
    #                default is 2  
    #
    # grid.size = user defined grid.size for kriging, format c(10,10), ie, 10x10 grid 
    # grid.coords = user defined coordinate points
    # shapefile = user defined shapefile for kriging
    # model =  "Exp", "Sph", "Gau"
    # plot.global = if TRUE, plot variograms
    # grid.size = user defined grid.size for kriging, format c(10,10), ie, 10x10 grid 
    # grid.coords = user defined coordinate points
    # shapefile = user defined shapefile for kriging
    #
    # model =  "Exp", "Sph", "Gau"
    #
    # neighbourhood.size => to be defined by the user 
    # neighbourhood.points => to be defined by the user: this is for each treatment
    #
    # For big-n problem: neighbourhood.size = f(maxd,tol.site,n), i.e., neighbourhood.size  <- max(c(d))*(tol.site/nrow(d))
    #
    # plot.local = TRUE will return the local variograms for each kriging points.
    #
    # tol.site = maximum number of observation locations to fit a global model
    #
    # sub.replicate = number of replications: used only for subsampling options
    # num.of.clusters = user defined cluster number: used only for subsampling options
    # sub.percentage = percentage of samples required to identify the cluster size (or boot samples) for each treatment: used only for subsampling options
    #
    # parallel = logical, if TRUE then do parallel computing
    # numCore = only applicable for parallel computing, if NULL then code will find the number of cores of the computer and do the run
    #
    if(is.null(data)){
      stop("\n Define data frame \n")
    }
    if(!(is.data.frame(data) | is.matrix(data))){
      stop("\n Define data frame or data matrix correctly with 4 columns: col(1,2)=(x,y) coordinates, col(3)=yield/observations, col(4)=treatment types/numbers \n")
    }
    if(ncol(data) !=4){
      stop("\n Define data frame correctly with 4 columns: col(1,2)=(x,y) coordinates, col(3)=yield/observations, col(4)=treatment types/numbers \n")
    }
    if (!type %in% c("global", "local-adaptive", "local-processed-r-a", "local-processed-r", "local-processed-p", "local-processed-rp", "local-fixed", "random", "kmeans", "boot", "epoch")) {
      stop("\n Define model correctly: global, local-adaptive, local-processed-radius, local-processed-points, local-processed-radius-points, local-fixed, random, kmeans, boot, epoch \n")
    }
    # check the prediction grid/coordinates 
    if (is.null(grid.size) &
        is.null(grid.coords) & is.null(shapefile)) {
      stop("\n Provide a value for the argument grid.size, grid.coords or shapefile \n")
    }
    else if (!is.null(grid.size)) {
      gr <- spT.grid.coords(c(max(data[,1],na.rm=TRUE),min(data[,1],na.rm=TRUE)), c(max(data[,2],na.rm=TRUE),min(data[,2],na.rm=TRUE)), by = grid.size)
    }
    else if (!is.null(grid.coords)) {
      gr <- as.matrix(grid.coords)
    }
    else if (!is.null(shapefile)) {
      gr <- as.matrix(coordinates(shapefile))
    }
    else {
      stop("\n Provide a value for the argument grid.size, grid.coords or shapefile \n")
    }
    #
    if(isTRUE(parallel)){
      # start the parallel process
      if(is.null(numCores)){
        # Calculate the number of cores
        numCores <- detectCores() - 1
        # Initiate cluster
        cl <- makeCluster(numCores)
      }
      else{
        # Initiate cluster
        cl <- makeCluster(numCores)
      }
      breaks <- round(seq(1,nrow(gr),length.out=(numCores+1)))
      breaks[1] <- 0
      if(breaks[numCores+1] != nrow(gr)){breaks[numCores+1] <- nrow(gr)}
      para_fnc <- function(i){
        gr.ck <- gr[((1)+breaks[i]):breaks[i+1],]
        .coKrig_inner_fnc(data = data,
                          type = type,
                          local.tuning = NULL,
                          grid.size = NULL,
                          grid.coords = gr.ck,
                          shapefile = NULL,
                          model = model,
                          neighbourhood.size = neighbourhood.size,
                          neighbourhood.points = neighbourhood.points,
                          pts = pts.for.each.treat, # this is for local approach to ensure the number of points in each treatment
                          plot.global = FALSE,
                          plot.local = FALSE,
                          tol.site = tol.site,
                          sub.replicate = sub.replicate,
                          num.of.clusters = num.of.clusters,
                          sub.percentage = sub.percentage,
                          ...)
      }
      print(numCores)
      out <- mclapply(1:numCores,para_fnc,mc.cores=numCores)
      n <- unique(unlist(lapply(out, names)))
      names(n) <- n
      #out <- lapply(n, function(i) unlist(lapply(out, `[[`, i)))
      out <- lapply(n, function(i) (lapply(out, `[[`, i)))
      out$kriged.output <- do.call(rbind,out$kriged.output)
      out$grid.coordinates <- do.call(rbind,out$grid.coordinates)
      n <- names(out)[!(names(out)%in%c("kriged.output","grid.coordinates"))]
      out[n] <- lapply(out[n],unlist)
      stopCluster(cl)
      class(out) <- "cokriging"
      out
    }
    else if(!isTRUE(parallel)){
      .coKrig_inner_fnc(data = data,
                        type = type,
                        local.tuning = NULL,
                        grid.size = NULL,
                        grid.coords = gr,
                        shapefile = shapefile,
                        model = model,
                        neighbourhood.size = neighbourhood.size,
                        neighbourhood.points = neighbourhood.points,
                        pts = pts.for.each.treat, # this is for local approach to ensure the number of points in each treatment
                        plot.global = FALSE,
                        plot.local = FALSE,
                        tol.site = tol.site,
                        sub.replicate = sub.replicate,
                        num.of.clusters = num.of.clusters,
                        sub.percentage = sub.percentage,
                        ...)
    }
    else{
      stop("Error: Define the argument 'parallel' correctly.")
    }
  }
#
# main function
.coKrig_inner_fnc <- 
  function(data = NULL,
           type = "global",
           local.tuning = NULL,
           grid.size = NULL,
           grid.coords = NULL,
           shapefile = NULL,
           model = "Exp",
           neighbourhood.size = NULL,
           neighbourhood.points = NULL,
           pts = NULL, 
           plot.global = FALSE,
           plot.local = FALSE,
           tol.site = 1500,
           sub.replicate,
           num.of.clusters,
           sub.percentage,
           ...) {
    #
    #start.time <- proc.time()[3]
    options(warn=-1)
    if(is.null(data)){
      stop("\n Define data frame \n")
    }
    if(!(is.data.frame(data) | is.matrix(data))){
      stop("\n Define data frame or data matrix correctly with 4 columns: col(1,2)=(x,y) coordinates, col(3)=yield/observations, col(4)=treatment types/numbers \n")
    }
    if(ncol(data) !=4){
      stop("\n Define data frame correctly with 4 columns: col(1,2)=(x,y) coordinates, col(3)=yield/observations, col(4)=treatment types/numbers \n")
    }
    if(type=="global"){
      out <- ._Global_fnc(data=data,
                          grid.size=grid.size,
                          grid.coords=grid.coords,
                          shapefile=shapefile,
                          model=model,
                          plot.global=plot.global,
                          tol.site=tol.site,
                          ...)      
    }
    else if(type=="local-adaptive"){
      out <- ._Local_fnc(data=data,
                         model=model,
                         grid.size=grid.size,
                         grid.coords=grid.coords,
                         shapefile=shapefile,
                         pts=pts, # this is for local approach to ensure the number of points in each treatment
                         block.size=neighbourhood.size,
                         plot.local=plot.local,
                         tol.site=tol.site,
                         ...)      
    }
    else if(type=="local-processed-r-a"){ # for radius proposed algorithm 
      coords <- as.matrix(data[, 1:2])
      dimnames(coords) <- NULL
      #x1 <- range(coords[, 1])
      #x2 <- range(coords[, 2])
      if (is.null(grid.size) &
          is.null(grid.coords) & is.null(shapefile)) {
        stop("\n Provide a value for the argument grid.size, grid.coords or shapefile \n")
      }
      else if (!is.null(grid.size)) {
        #gr <- spT.grid.coords(c(max(data[,1],na.rm=TRUE),min(data[,1],na.rm=TRUE)), c(max(data[,2],na.rm=TRUE),min(data[,2],na.rm=TRUE)), by = grid.size)
        #gr <- spT.grid.coords(x1, x2, by = grid.size)
        gr <- spT.grid.coords(range(coords[, 1],na.rm=TRUE), range(coords[, 2],na.rm=TRUE), by = grid.size)
      }
      else if (!is.null(grid.coords)) {
        gr <- as.matrix(grid.coords)
      }
      else if (!is.null(shapefile)) {
        gr <- as.matrix(coordinates(shapefile))
      }
      else {
        stop("\n Provide a value for the argument grid.size, grid.coords or shapefile \n")
      }
      ck <- cbind(gr,1:nrow(gr))
      ch <- chull(gr[,1:2])
      gr.ck <- ck[ch,]
      ck <- ._Local_fnc(data=data,
                        model=model,
                        grid.size=NULL,
                        grid.coords=gr.ck[,1:2],
                        shapefile=NULL,
                        pts = 5, # this is for local approach to ensure the number of points in each treatment
                        block.size=neighbourhood.size,
                        plot.local=FALSE,
                        tol.site=tol.site,
                        ...)$neighbourhood.size      
      #browser()
      names(ck) <- NULL
      gr.ck[,3] <- ck
      dimnames(gr.ck) <- NULL
      print(paste0("Estimating the 'Edge Effects' ... "))
      print(paste0("The initial maximum spatial range parameter is considered as: ", round(max(ck,na.rm=TRUE),4)))
      print(paste0("The initial minimum spatial range parameter is considered as: ", round(min(ck,na.rm=TRUE),4)))
      #browser()
      #
      ck <- cbind(gr,1:nrow(gr))
      ck <- ck[!(ck[,3]%in%ch),]
      ch <- c(ch, ch[1])
      ch <- gr[ch,]
      ch <- Polygon(ch,hole=FALSE) # from sp package
      ch <- ch@area # area of the region 1/2*n*pi*r^2 = A
      if(is.null(local.tuning)){
        local.tuning <- 1 #2
      }
      else if(local.tuning<0){
        local.tuning <- 0
      }
      else{
        local.tuning <- local.tuning
      }
      ch <- round((local.tuning*ch)/(pi*min(gr.ck[,3])^2))
      if(ch <= 25 | ch <= nrow(ck)){
        ch  <- min((nrow(ck)-1),25)
      }
      #browser()
      gr.ck2 <- kmeans(ck[,1:2],ch)$centers
      ck <- ._Local_fnc(data=data,
                        model=model,
                        grid.size=NULL,
                        grid.coords=gr.ck2[,1:2],
                        shapefile=NULL,
                        pts = 5, # this is for local approach to ensure the number of points in each treatment
                        block.size=neighbourhood.size,
                        plot.local=FALSE,
                        tol.site=tol.site,
                        ...)$neighbourhood.size
      names(ck) <- NULL
      gr.ck2 <- cbind(gr.ck2,c(ck))
      dimnames(gr.ck2) <- NULL
      #print(paste0("Number of kmean points: ", ch))
      #print(paste0("The estimated average 'range' from kmean: ", round(mean(ck,na.rm=TRUE),4)))
      print(paste0("Estimating the 'Overall Effects' ... "))
      ## splines 
      gr.ck <- rbind(gr.ck,gr.ck2)
      gr.ck <- data.frame(gr.ck)
      names(gr.ck) <- c("x","y","z")
      #fit <- gam(z~s(x,y,k=nrow(gr.ck)-1),data=gr.ck); 
      fit <- gam(z~s(x,y,k=round(nrow(gr.ck)*0.5)),data=gr.ck); 
      gr <- data.frame(gr)
      names(gr) <- c("x","y")
      pr <- as.vector(predict(fit,gr))
      rm(gr.ck)
      gr <- as.matrix(gr)
      ##
      out <-  ._Local_fnc_processed(data=data,
                                    model=model,
                                    grid.size=grid.size,
                                    grid.coords=grid.coords,
                                    shapefile=shapefile,
                                    pts=pts, # this is for local approach to ensure the number of points in each treatment
                                    #block.size=max(ck,na.rm=TRUE),
                                    block.size=pr,
                                    plot.local=plot.local,
                                    tol.site=tol.site,
                                    ...)       
    }
    else if(type=="local-processed-r"){ # for range - cross val 
      coords <- as.matrix(data[, 1:2])
      dimnames(coords) <- NULL
      if (is.null(grid.size) &
          is.null(grid.coords) & is.null(shapefile)) {
        stop("\n Provide a value for the argument grid.size, grid.coords or shapefile \n")
      }
      else if (!is.null(grid.size)) {
        gr <- spT.grid.coords(range(coords[, 1],na.rm=TRUE), range(coords[, 2],na.rm=TRUE), by = grid.size)
      }
      else if (!is.null(grid.coords)) {
        gr <- as.matrix(grid.coords)
      }
      else if (!is.null(shapefile)) {
        gr <- as.matrix(coordinates(shapefile))
      }
      else {
        stop("\n Provide a value for the argument grid.size, grid.coords or shapefile \n")
      }
      #
      # find train and test data for cross-val
      tr <- unique(data[,4])
      data$id <- 1:nrow(data)
      index <- chull(data[,1:2]); ndat <- data[index,]; dd <- data[-index,]
      index <- chull(dd[,1:2]); ndat <- rbind(ndat,dd[index,]); dd <- dd[-index,]
      index <- chull(dd[,1:2]); ndat <- rbind(ndat,dd[index,]); dd <- dd[-index,]
      index <- chull(dd[,1:2]); ndat <- rbind(ndat,dd[index,]); dd <- dd[-index,]
      index <- chull(dd[,1:2]); ndat <- rbind(ndat,dd[index,]); dd <- dd[-index,]
      ndat <- ndat[order(ndat[,4]),]; rownames(ndat) <- NULL
      if(length(tr)==3){
        size <- 5 #20 # 5
        index <- c(sample(x=1:nrow(ndat[ndat[,4]==tr[1],]),size=size),
                   sample(x=(nrow(ndat[ndat[,4]==tr[1],])+1):(nrow(ndat[ndat[,4]==tr[1],])+nrow(ndat[ndat[,4]==tr[2],])),size=size),
                   sample(x=(nrow(ndat[ndat[,4]==tr[1],])+nrow(ndat[ndat[,4]==tr[2],])+1):(nrow(ndat[ndat[,4]==tr[1],])+nrow(ndat[ndat[,4]==tr[2],])+nrow(ndat[ndat[,4]==tr[3],])),size=size))
        index <- ndat[index,]$id
        ii <- c(sample(x=1:nrow(data[data[,4]==tr[1],]),size=size),
                sample(x=(nrow(data[data[,4]==tr[1],])+1):(nrow(data[data[,4]==tr[1],])+nrow(data[data[,4]==tr[2],])),size=size),
                sample(x=(nrow(data[data[,4]==tr[1],])+nrow(data[data[,4]==tr[2],])+1):(nrow(data[data[,4]==tr[1],])+nrow(data[data[,4]==tr[2],])+nrow(data[data[,4]==tr[3],])),size=size))
        index <- sort(c(index,ii))
      }
      else if(length(tr)==2){
        size <- 7 #25 # 7
        index <- c(sample(x=1:nrow(ndat[ndat[,4]==tr[1],]),size=size),
                   sample(x=(nrow(ndat[ndat[,4]==tr[1],])+1):(nrow(ndat[ndat[,4]==tr[1],])+nrow(ndat[ndat[,4]==tr[2],])),size=size))
        index <- ndat[index,]$id
        ii <- c(sample(x=1:nrow(data[data[,4]==tr[1],]),size=size),
                sample(x=(nrow(data[data[,4]==tr[1],])+1):(nrow(data[data[,4]==tr[1],])+nrow(data[data[,4]==tr[2],])),size=size))
        index <- sort(c(index,ii))	 
      }
      else{
        stop("Check the number of treatments. Currently it takes: 2 or 3 treatments.")
      }
      rm(dd); rm(ndat)
      data$id <- NULL
      test <- data[index,]
      n.pts <- c(table(data[,4])); n.pts <- n.pts[!n.pts%in%0];
      n.pts <- round(seq(20,min(200,c(n.pts-1)),by=5)) # greedy search for the initial sample points 
      test$range <- NA # range 
      test$points <- NA # neighbourhood points 
      print(paste0("Performing cross-validation"))
      pb <- txtProgressBar(min = 0, max = nrow(test), style = 3)   # set progress bar
      for(i in 1:nrow(test)){
        ck <- list()
        j <- 1
        sink(tempfile())
        ck[[j]] <- ._Local_fnc_fixed(data=data[-index[i],], model=model, grid.coords=as.matrix(test[i,1:2]),neighbourhood.points = n.pts[j], tol.site = tol.site)
        sink()
        repeat{
          j <- j+1
          sink(tempfile())
          ck[[j]] <- ._Local_fnc_fixed(data=data[-index[i],], model=model, grid.coords=as.matrix(test[i,1:2]),neighbourhood.points = n.pts[j], tol.site = tol.site)
          sink()
          if(length(tr)==3){
            r1 <- data.frame(t(ck[[j-1]]$kriged.output[,c(1,3,5)]),test[i,3:4])
            bias1 <- abs(r1[,which(tr%in%r1[,ncol(r1)])]-test[i,3])
            r2 <- data.frame(t(ck[[j]]$kriged.output[,c(1,3,5)]),test[i,3:4])
            bias2 <- abs(r2[,which(tr%in%r2[,ncol(r2)])]-test[i,3])
          }
          else if(length(tr)==2){
            r1 <- data.frame(t(ck[[j-1]]$kriged.output[,c(1,3)]),test[i,3:4])
            bias1 <- abs(r1[,which(tr%in%r1[,ncol(r1)])]-test[i,3])
            r2 <- data.frame(t(ck[[j]]$kriged.output[,c(1,3)]),test[i,3:4])
            bias2 <- abs(r2[,which(tr%in%r2[,ncol(r2)])]-test[i,3])
          }
          else{
            stop("Currently can take 2 or 3 treatments.")
          }
          #
          pass <- abs((bias2-bias1)/(bias2+bias1))
          if(!is.na(pass)){
            if (pass <= 0.01 | n.pts[j] == n.pts[length(n.pts)]) {      
              break
            }
          }
        }
        ns1 <- unlist(ck[[j-1]]$est.parameters) # estimated parameters 
        test[i,]$range <- mean(ns1[grep("range",names(ns1))]) # estimated initial range parameters 
        test[i,]$points <- mean(ck[[j-1]]$est.sample) # estimated initial points  
        do.call(file.remove, list(list.files(tempdir(), full.names = TRUE)))
        setTxtProgressBar(pb, i)
      }
      close(pb)  # end progress bar 
      print(paste0("Estimating the 'Overall Effects' ... "))
      ## splines 
      test <- test[,c(1,2,5,6)]
      names(test) <- c("x","y","range","points")
      fit1 <- gam(range~s(x,y,k=round(nrow(test)*0.5)),data=test); 
      gr <- data.frame(gr)
      names(gr) <- c("x","y")
      gr$range <- as.vector(predict(fit1,gr))
      fit2 <- gam(points~s(x,y,range,k=round(nrow(test)*0.5)),data=test); 
      gr$points <- round(as.vector(predict(fit2,gr)))
      gr <- as.matrix(gr)
      # check the thresholds 
      gr[gr[,3] <=0,3] <- c(quantile(gr[,3],prob=0.1))
      gr[gr[,4] <=min(n.pts),4] <- min(n.pts)
      #browser()
      ##
      out <-  ._Local_fnc_processed(data=data,
                                    model=model,
                                    grid.size=NULL,
                                    grid.coords=as.matrix(gr[,1:2]),
                                    shapefile=NULL,
                                    pts=min(n.pts), # this is for local approach to ensure the number of points in each treatment
                                    block.size=gr[,3],
                                    plot.local=FALSE,
                                    tol.site=tol.site,
                                    ...)       
    }
    else if(type=="local-processed-p"){
      coords <- as.matrix(data[, 1:2])
      dimnames(coords) <- NULL
      if (is.null(grid.size) &
          is.null(grid.coords) & is.null(shapefile)) {
        stop("\n Provide a value for the argument grid.size, grid.coords or shapefile \n")
      }
      else if (!is.null(grid.size)) {
        gr <- spT.grid.coords(range(coords[, 1],na.rm=TRUE), range(coords[, 2],na.rm=TRUE), by = grid.size)
      }
      else if (!is.null(grid.coords)) {
        gr <- as.matrix(grid.coords)
      }
      else if (!is.null(shapefile)) {
        gr <- as.matrix(coordinates(shapefile))
      }
      else {
        stop("\n Provide a value for the argument grid.size, grid.coords or shapefile \n")
      }
      #
      # find train and test data for cross-val
      tr <- unique(data[,4])
      data$id <- 1:nrow(data)
      index <- chull(data[,1:2]); ndat <- data[index,]; dd <- data[-index,]
      index <- chull(dd[,1:2]); ndat <- rbind(ndat,dd[index,]); dd <- dd[-index,]
      index <- chull(dd[,1:2]); ndat <- rbind(ndat,dd[index,]); dd <- dd[-index,]
      index <- chull(dd[,1:2]); ndat <- rbind(ndat,dd[index,]); dd <- dd[-index,]
      index <- chull(dd[,1:2]); ndat <- rbind(ndat,dd[index,]); dd <- dd[-index,]
      ndat <- ndat[order(ndat[,4]),]; rownames(ndat) <- NULL
      if(length(tr)==3){
        size <- 5 #20 # 5
        index <- c(sample(x=1:nrow(ndat[ndat[,4]==tr[1],]),size=size),
                   sample(x=(nrow(ndat[ndat[,4]==tr[1],])+1):(nrow(ndat[ndat[,4]==tr[1],])+nrow(ndat[ndat[,4]==tr[2],])),size=size),
                   sample(x=(nrow(ndat[ndat[,4]==tr[1],])+nrow(ndat[ndat[,4]==tr[2],])+1):(nrow(ndat[ndat[,4]==tr[1],])+nrow(ndat[ndat[,4]==tr[2],])+nrow(ndat[ndat[,4]==tr[3],])),size=size))
        index <- ndat[index,]$id
        ii <- c(sample(x=1:nrow(data[data[,4]==tr[1],]),size=size),
                sample(x=(nrow(data[data[,4]==tr[1],])+1):(nrow(data[data[,4]==tr[1],])+nrow(data[data[,4]==tr[2],])),size=size),
                sample(x=(nrow(data[data[,4]==tr[1],])+nrow(data[data[,4]==tr[2],])+1):(nrow(data[data[,4]==tr[1],])+nrow(data[data[,4]==tr[2],])+nrow(data[data[,4]==tr[3],])),size=size))
        index <- sort(c(index,ii))
      }
      else if(length(tr)==2){
        size <- 7 #25 # 7
        index <- c(sample(x=1:nrow(ndat[ndat[,4]==tr[1],]),size=size),
                   sample(x=(nrow(ndat[ndat[,4]==tr[1],])+1):(nrow(ndat[ndat[,4]==tr[1],])+nrow(ndat[ndat[,4]==tr[2],])),size=size))
        index <- ndat[index,]$id
        ii <- c(sample(x=1:nrow(data[data[,4]==tr[1],]),size=size),
                sample(x=(nrow(data[data[,4]==tr[1],])+1):(nrow(data[data[,4]==tr[1],])+nrow(data[data[,4]==tr[2],])),size=size))
        index <- sort(c(index,ii))	 
      }
      else{
        stop("Check the number of treatments. Currently it takes: 2 or 3 treatments.")
      }
      rm(dd); rm(ndat)
      data$id <- NULL
      test <- data[index,]
      n.pts <- c(table(data[,4])); n.pts <- n.pts[!n.pts%in%0];
      n.pts <- round(seq(25,min(200,c(n.pts-1)),by=5)) # greedy search for the initial sample points 
      test$range <- NA # range 
      test$points <- NA # neighbourhood points 
      print(paste0("Performing cross-validation"))
      pb <- txtProgressBar(min = 0, max = nrow(test), style = 3)   # set progress bar
      for(i in 1:nrow(test)){
        ck <- list()
        j <- 1
        sink(tempfile())
        ck[[j]] <- ._Local_fnc_fixed(data=data[-index[i],], model=model, grid.coords=as.matrix(test[i,1:2]),neighbourhood.points = n.pts[j], tol.site = tol.site)
        sink()
        repeat{
          j <- j+1
          sink(tempfile())
          ck[[j]] <- ._Local_fnc_fixed(data=data[-index[i],], model=model, grid.coords=as.matrix(test[i,1:2]),neighbourhood.points = n.pts[j], tol.site = tol.site)
          sink()
          if(length(tr)==3){
            r1 <- data.frame(t(ck[[j-1]]$kriged.output[,c(1,3,5)]),test[i,3:4])
            bias1 <- abs(r1[,which(tr%in%r1[,ncol(r1)])]-test[i,3])
            r2 <- data.frame(t(ck[[j]]$kriged.output[,c(1,3,5)]),test[i,3:4])
            bias2 <- abs(r2[,which(tr%in%r2[,ncol(r2)])]-test[i,3])
          }
          else if(length(tr)==2){
            r1 <- data.frame(t(ck[[j-1]]$kriged.output[,c(1,3)]),test[i,3:4])
            bias1 <- abs(r1[,which(tr%in%r1[,ncol(r1)])]-test[i,3])
            r2 <- data.frame(t(ck[[j]]$kriged.output[,c(1,3)]),test[i,3:4])
            bias2 <- abs(r2[,which(tr%in%r2[,ncol(r2)])]-test[i,3])
          }
          else{
            stop("Currently can take 2 or 3 treatments.")
          }
          #
          pass <- abs((bias2-bias1)/(bias2+bias1))
          if(!is.na(pass)){
            if (pass <= 0.01 | n.pts[j] == n.pts[length(n.pts)]) {      
              break
            }
          }
        }
        ns1 <- unlist(ck[[j-1]]$est.parameters) # estimated range parameters 
        test[i,]$range <- mean(ns1[grep("range",names(ns1))])
        test[i,]$points <- mean(ck[[j-1]]$est.sample)
        do.call(file.remove, list(list.files(tempdir(), full.names = TRUE)))
        setTxtProgressBar(pb, i)
      }
      close(pb)  # end progress bar 
      print(paste0("Estimating the 'Overall Effects' ... "))
      ## splines 
      test <- test[,c(1,2,5,6)]
      names(test) <- c("x","y","range","points")
      #fit1 <- gam(range~s(x,y,k=nrow(test)-1),data=test); 
      fit1 <- gam(range~s(x,y,k=round(nrow(test)*0.5)),data=test); 
      gr <- data.frame(gr)
      names(gr) <- c("x","y")
      gr$range <- as.vector(predict(fit1,gr))
      #fit2 <- gam(points~s(x,y,range,k=nrow(test)-1),data=test); 
      fit2 <- gam(points~s(x,y,range,k=round(nrow(test)*0.5)),data=test); 
      gr$points <- round(as.vector(predict(fit2,gr)))
      gr <- as.matrix(gr)
      # check the thresholds 
      gr[gr[,3] <=0,3] <- c(quantile(gr[,3],prob=0.1))
      gr[gr[,4] <= min(n.pts),4] <- min(n.pts)
      ##
      out <- ._Local_fnc_fixed(data=data, model=model, grid.coords=as.matrix(gr[,1:2]),neighbourhood.points=c(gr[,4]), tol.site=tol.site, ...)       
    }
    else if(type=="local-processed-rp"){
      coords <- as.matrix(data[, 1:2])
      dimnames(coords) <- NULL
      if (is.null(grid.size) &
          is.null(grid.coords) & is.null(shapefile)) {
        stop("\n Provide a value for the argument grid.size, grid.coords or shapefile \n")
      }
      else if (!is.null(grid.size)) {
        gr <- spT.grid.coords(range(coords[, 1],na.rm=TRUE), range(coords[, 2],na.rm=TRUE), by = grid.size)
      }
      else if (!is.null(grid.coords)) {
        gr <- as.matrix(grid.coords)
      }
      else if (!is.null(shapefile)) {
        gr <- as.matrix(coordinates(shapefile))
      }
      else {
        stop("\n Provide a value for the argument grid.size, grid.coords or shapefile \n")
      }
      #
      # find train and test data for cross-val
      tr <- unique(data[,4])
      data$id <- 1:nrow(data)
      index <- chull(data[,1:2]); ndat <- data[index,]; dd <- data[-index,]
      index <- chull(dd[,1:2]); ndat <- rbind(ndat,dd[index,]); dd <- dd[-index,]
      index <- chull(dd[,1:2]); ndat <- rbind(ndat,dd[index,]); dd <- dd[-index,]
      index <- chull(dd[,1:2]); ndat <- rbind(ndat,dd[index,]); dd <- dd[-index,]
      index <- chull(dd[,1:2]); ndat <- rbind(ndat,dd[index,]); dd <- dd[-index,]
      ndat <- ndat[order(ndat[,4]),]; rownames(ndat) <- NULL
      if(length(tr)==3){
        size <- 5 #20 #5
        index <- c(sample(x=1:nrow(ndat[ndat[,4]==tr[1],]),size=size),
                   sample(x=(nrow(ndat[ndat[,4]==tr[1],])+1):(nrow(ndat[ndat[,4]==tr[1],])+nrow(ndat[ndat[,4]==tr[2],])),size=size),
                   sample(x=(nrow(ndat[ndat[,4]==tr[1],])+nrow(ndat[ndat[,4]==tr[2],])+1):(nrow(ndat[ndat[,4]==tr[1],])+nrow(ndat[ndat[,4]==tr[2],])+nrow(ndat[ndat[,4]==tr[3],])),size=size))
        index <- ndat[index,]$id
        ii <- c(sample(x=1:nrow(data[data[,4]==tr[1],]),size=size),
                sample(x=(nrow(data[data[,4]==tr[1],])+1):(nrow(data[data[,4]==tr[1],])+nrow(data[data[,4]==tr[2],])),size=size),
                sample(x=(nrow(data[data[,4]==tr[1],])+nrow(data[data[,4]==tr[2],])+1):(nrow(data[data[,4]==tr[1],])+nrow(data[data[,4]==tr[2],])+nrow(data[data[,4]==tr[3],])),size=size))
        index <- sort(c(index,ii))
      }
      else if(length(tr)==2){
        size <- 7 #25 #7
        index <- c(sample(x=1:nrow(ndat[ndat[,4]==tr[1],]),size=size),
                   sample(x=(nrow(ndat[ndat[,4]==tr[1],])+1):(nrow(ndat[ndat[,4]==tr[1],])+nrow(ndat[ndat[,4]==tr[2],])),size=size))
        index <- ndat[index,]$id
        ii <- c(sample(x=1:nrow(data[data[,4]==tr[1],]),size=size),
                sample(x=(nrow(data[data[,4]==tr[1],])+1):(nrow(data[data[,4]==tr[1],])+nrow(data[data[,4]==tr[2],])),size=size))
        index <- sort(c(index,ii))	 
      }
      else{
        stop("Check the number of treatments. Currently it takes: 2 or 3 treatments.")
      }
      rm(dd); rm(ndat)
      data$id <- NULL
      test <- data[index,]
      n.pts <- c(table(data[,4])); n.pts <- n.pts[!n.pts%in%0];
      n.pts <- round(seq(25,min(200,c(n.pts-1)),by=5)) # greedy search for the initial sample points 
      test$range <- NA # range 
      test$points <- NA # neighbourhood points 
      print(paste0("Performing cross-validation"))
      pb <- txtProgressBar(min = 0, max = nrow(test), style = 3)   # set progress bar
      for(i in 1:nrow(test)){
        ck <- list()
        j <- 1
        sink(tempfile())
        ck[[j]] <- ._Local_fnc_fixed(data=data[-index[i],], model=model, grid.coords=as.matrix(test[i,1:2]),neighbourhood.points = n.pts[j], tol.site = tol.site)
        sink()
        repeat{
          j <- j+1
          sink(tempfile())
          ck[[j]] <- ._Local_fnc_fixed(data=data[-index[i],], model=model, grid.coords=as.matrix(test[i,1:2]),neighbourhood.points = n.pts[j], tol.site = tol.site)
          sink()
          if(length(tr)==3){
            r1 <- data.frame(t(ck[[j-1]]$kriged.output[,c(1,3,5)]),test[i,3:4])
            bias1 <- abs(r1[,which(tr%in%r1[,ncol(r1)])]-test[i,3])
            r2 <- data.frame(t(ck[[j]]$kriged.output[,c(1,3,5)]),test[i,3:4])
            bias2 <- abs(r2[,which(tr%in%r2[,ncol(r2)])]-test[i,3])
          }
          else if(length(tr)==2){
            r1 <- data.frame(t(ck[[j-1]]$kriged.output[,c(1,3)]),test[i,3:4])
            bias1 <- abs(r1[,which(tr%in%r1[,ncol(r1)])]-test[i,3])
            r2 <- data.frame(t(ck[[j]]$kriged.output[,c(1,3)]),test[i,3:4])
            bias2 <- abs(r2[,which(tr%in%r2[,ncol(r2)])]-test[i,3])
          }
          else{
            stop("Currently can take 2 or 3 treatments.")
          }
          #
          pass <- abs((bias2-bias1)/(bias2+bias1))
          if(!is.na(pass)){
            if (pass <= 0.01 | n.pts[j] == n.pts[length(n.pts)]) {      
              break
            }
          }
        }
        ns1 <- unlist(ck[[j-1]]$est.parameters) # estimated range parameters 
        test[i,]$range <- mean(ns1[grep("range",names(ns1))])
        test[i,]$points <- mean(ck[[j-1]]$est.sample)
        #test[i,]$yhat <- r1[,which(tr%in%r1[,ncol(r1)])]
        do.call(file.remove, list(list.files(tempdir(), full.names = TRUE)))
        setTxtProgressBar(pb, i)
      }
      close(pb)  # end progress bar 
      print(paste0("Estimating the 'Overall Effects' ... "))
      ## splines 
      test <- test[,c(1,2,5,6)]
      names(test) <- c("x","y","range","points")
      #fit1 <- gam(range~s(x,y,k=nrow(test)-1),data=test); 
      fit1 <- gam(range~s(x,y,k=round(nrow(test)*0.5)),data=test); 
      gr <- data.frame(gr)
      names(gr) <- c("x","y")
      gr$range <- as.vector(predict(fit1,gr))
      #fit2 <- gam(points~s(x,y,range,k=nrow(test)-1),data=test); 
      fit2 <- gam(points~s(x,y,range,k=round(nrow(test)*0.5)),data=test); 
      gr$points <- round(as.vector(predict(fit2,gr)))
      gr <- as.matrix(gr)
      # check the thresholds 
      gr[gr[,3] <=0,3] <- c(quantile(gr[,3],prob=0.1))
      gr[gr[,4] <= min(n.pts),4] <- min(n.pts)
      ##
      out <-  ._Local_fnc_rp(data=data,
                             model=model,
                             grid.size=NULL,
                             grid.coords=as.matrix(gr[,1:2]),
                             shapefile=NULL,
                             pts=c(gr[,4]), # this is for local approach to ensure the number of points in each treatment
                             block.size=c(gr[,3]),
                             plot.local=plot.local,
                             tol.site=tol.site,
                             ...)       
    }
    else if(type=="local-fixed"){
      out <- ._Local_fnc_fixed(data=data,
                               model=model,
                               grid.size=grid.size,
                               grid.coords=grid.coords,
                               shapefile=shapefile,
                               block.size=NULL,
                               neighbourhood.points = neighbourhood.points,
                               plot.local=plot.local,
                               tol.site=tol.site,
                               ...)      
    }
    else if(type=="kmeans"){
      # sub.replicate = number of replications: used only for subsampling options
      # sub.percentage = percentage of samples required to identify the cluster size for each treatment: used only for subsampling options
      pb <- txtProgressBar(min = 0, max = sub.replicate, style = 3)   # set progress bar
      if(sub.percentage <=0 | sub.percentage >100){stop("\n sub.percentage should be between 0 to 100")}
      out <- NULL
      nms <- unique(data[,4])
      n <- table(data[,4])
      n <- n[!n%in%0]
      if(length(n)==3){m <-9}
      if(length(n)==2){m <-5}
      if(length(n)==1){m <-2}
      #if(length(n)==1){stop("\n Please provide at least 2 treatments \n")}
      out <- NULL
      out$kriged.output <- matrix(0,nrow(grid.coords),m)
      out$est.sample <- NULL
      out$num.of.samples <- NULL
      for(i in 1:sub.replicate){
        data.ck <- NULL
        ss <- NULL
        for(j in 1:length(n)){
          ck <- data[data[,4]==nms[j],]
          if(is.null(num.of.clusters)){
            km.number  <- round(sqrt(nrow(ck)))
          }
          else{
            km.number  <- num.of.clusters
          }
          if(km.number <= 2){ km.number <- 2 }
          sub.no.samples <- round(((sub.percentage)/(100))*as.vector(n)[j]) # getting subsample for each treatments
          km.number <- min(km.number,sub.no.samples)
          if(km.number <= 1){ km.number <- 1 }
          set.seed(round(runif(1,1001,10000001)))
          df <- kmeans(ck[,1:2], km.number) #
          ck <- data.frame(ck,df$cluster)
          sk <- floor(sub.no.samples/km.number)
          sk <- min(sk,min(table(ck$df.cluster)))
          so <- sub.no.samples-sk*km.number
          df <- data.frame(ck %>% group_by(ck$df.cluster) %>% sample_n(sk))[,1:4]
          df <- rbind(df, sample_n(ck[,1:4], so, replace=TRUE))
          df <- df[!duplicated(df),]
          data.ck <- rbind(data.ck,df)
          rm(df); rm(ck); rm(sk); rm(so);
          ss <- c(ss,sub.no.samples)
        }
        row.names(data.ck) <- NULL; data.ck <- data.frame(data.ck)
        names(data.ck) <- names(data)
        #print(dim(data.ck))
        ck <- ._Global_fnc(data=data.ck,
                           grid.size=grid.size,
                           grid.coords=grid.coords,
                           shapefile=shapefile,
                           model=model,
                           plot.global=plot.global,
                           tol.site=tol.site,
                           ...)
        out$kriged.output <- out$kriged.output + ck$kriged.output
        out$est.parameters[[i]] <- ck$est.parameters
        out$est.sample <- c(out$est.sample,ck$est.sample)
        #out$num.of.samples <- rbind(out$num.of.samples,ss)
        out$num.of.samples <- ss
        setTxtProgressBar(pb, i)
      }
      close(pb)  # end progress bar 
      out <- list(kriged.output=out$kriged.output/sub.replicate,
                  est.parameters=out$est.parameters, est.sample=out$est.sample,
                  grid.coordinates=ck$grid.coordinates,
                  model=ck$model,grid.size=ck$grid.size,shapefile=ck$shapefile,type="kmeans",
                  num.of.samples=out$num.of.samples,subsample.percentage=sub.percentage,
                  num.of.replications=sub.replicate)
    }
    else if(type=="boot" | type=="random"){
      if(sub.percentage <=0 | sub.percentage >100){stop("\n sub.percentage should be between 0 to 100")}
      nms <- unique(data[,4])
      n <- table(data[,4])
      n <- n[!n%in%0]
      if(length(n)==3){m <-9}
      if(length(n)==2){m <-5}
      if(length(n)==1){m <-2}
      #if(length(n)==1){stop("\n Please provide at least 2 treatments \n")}
      out <- NULL
      out$kriged.output <- matrix(0,nrow(grid.coords),m)
      out$est.sample <- NULL
      out$num.of.samples <- NULL
      if(type=="random"){ sub.replicate <- 1}
      pb <- txtProgressBar(min = 0, max = sub.replicate, style = 3)   # set progress bar
      for(i in 1:sub.replicate){
        data.ck <- NULL
        ss <- NULL
        for(j in 1:length(n)){
          ck <- data[data[,4]==nms[j],]
          #sub.no.samples <- round(((sub.percentage)/(100*length(n)))*nrow(ck))
          sub.no.samples <- round(((sub.percentage)/(100))*as.vector(n)[j])
          if(sub.no.samples>nrow(ck)){stop("\n Decrease the sub.no.samples \n")}
          if(sub.no.samples<=10){stop("\n Increase the subsample.percentage \n")}
          set.seed(round(runif(1,1001,10000001)))
          data.ck <- rbind(data.ck,ck[sample(nrow(ck),size=sub.no.samples),])
          rm(ck)
          ss <- c(ss,sub.no.samples)
        }
        row.names(data.ck) <- NULL; data.ck <- data.frame(data.ck)
        names(data.ck) <- names(data)
        ck <- ._Global_fnc(data=data.ck,
                           grid.size=grid.size,
                           grid.coords=grid.coords,
                           shapefile=shapefile,
                           model=model,
                           plot.global=plot.global,
                           tol.site=tol.site,
                           ...)
        out$kriged.output <- out$kriged.output + ck$kriged.output
        out$est.parameters[[i]] <- ck$est.parameters
        out$est.sample <- c(out$est.sample,ck$est.sample)
        #out$num.of.samples <- rbind(out$num.of.samples,ss)
        out$num.of.samples <- ss
        setTxtProgressBar(pb, i)
      }
      close(pb)  # end progress bar 
      out <- list(kriged.output=out$kriged.output/sub.replicate,
                  est.parameters=out$est.parameters, est.sample=out$est.sample,
                  grid.coordinates=ck$grid.coordinates,
                  model=ck$model,grid.size=ck$grid.size,shapefile=ck$shapefile,type=type,
                  num.of.samples=out$num.of.samples,subsample.percentage=sub.percentage,
                  num.of.replications=sub.replicate)
    }
    else if(type=="epoch"){
      nms <- unique(data[,4])
      n <- table(data[,4])
      n <- n[!n%in%0]
      if(length(n)==3){m <-9}
      if(length(n)==2){m <-5}
      if(length(n)==1){m <-2}
      out <- NULL
      out$kriged.output <- matrix(0,nrow(grid.coords),m)
      out$est.sample <- NULL
      batch.size <- 10
      sub.replicate <- ceiling(sum(n)/(batch.size*length(n)))
      pb <- txtProgressBar(min = 0, max = sub.replicate, style = 3)   # set progress bar
      dat <- data
      #print(sub.replicate)
      for(i in 1:sub.replicate){
        if(i==sub.replicate){
          data.ck <- dat
          row.names(data.ck) <- NULL; data.ck <- data.frame(data.ck)
          names(data.ck) <- names(data)
          ck <- ._Global_fnc(data=data.ck,
                             grid.size=grid.size,
                             grid.coords=grid.coords,
                             shapefile=shapefile,
                             model=model,
                             plot.global=plot.global,
                             tol.site=tol.site,
                             ...)
        }
        else{
          data.ck <- NULL
          for(j in 1:length(n)){
            ck <- dat[dat[,4]==nms[j],]
            #if(batch.size>nrow(ck)){batch.size<-nrow(ck)}
            data.ck <- rbind(data.ck,ck[sample(nrow(ck),size=batch.size),])
            rm(ck)
          }
          dat <- dat[!(row.names(dat)%in%row.names(data.ck)),]
          #row.names(data.ck) <- NULL; data.ck <- data.frame(data.ck)
          names(data.ck) <- names(data)
          ck <- ._Global_fnc(data=data.ck,
                             grid.size=grid.size,
                             grid.coords=grid.coords,
                             shapefile=shapefile,
                             model=model,
                             plot.global=plot.global,
                             tol.site=tol.site,
                             ...)
          #for(j in 1:m){
          #  out$kriged.output[,j] <- apply(cbind(ck$kriged.output[,j],out$kriged.output[,j]),1,max,na.rm=TRUE)
          #}
        }
        #browser()
        #print(i)
        out$kriged.output <- out$kriged.output + ck$kriged.output
        out$est.parameters[[i]] <- ck$est.parameters
        out$est.sample <- c(out$est.sample,ck$est.sample)
        setTxtProgressBar(pb, i)
      }
      #browser()
      close(pb)  # end progress bar 
      out <- list(kriged.output=out$kriged.output/sub.replicate,
                  est.parameters=out$est.parameters, est.sample=out$est.sample,
                  grid.coordinates=ck$grid.coordinates,
                  model=ck$model,grid.size=ck$grid.size,shapefile=ck$shapefile,type="boot")
    }
    else{
      stop("\n Define model correctly: global, mixed or local \n")
    }
    #
    # check the missing values
    #
    # infill with average smooth for missing values
    # 8 adjacent grids are used in this context
    #
    if(length(unlist(out$kriged.output[complete.cases(out$kriged.output),]))==length(unlist(out$kriged.output))){
      out$kriged.output <- out$kriged.output	
    }
    else{
      check <- as.data.frame(out$kriged.output)
      s <- as.numeric(row.names(check[!complete.cases(check),]))
      d.check <- rdist(out$grid.coordinates)
      #d.check <- d.check[s,]
      for(i in 1:length(s)){
        ss <- which(c(d.check[s[i],]) <= sort(d.check[s[i],])[9]+1)
        check[s[i],] <- apply(check[ss[!ss%in%s[i]],],2,mean,na.rm=TRUE)
      }
      rm(d.check); rm(s); rm(ss);
      out$kriged.output <- as.matrix(check)
    }
    #
    #end.time <- proc.time()[3]
    #t <- end.time-start.time
    #out$comp.time <- .fnc.time_(t)
    #
    out
  }
##
##
._Global_fnc <-
  function(data,
           grid.size,
           grid.coords,
           shapefile,
           model,
           plot.global,
           tol.site,
           ...) {
    #
    # grid.size = user defined grid.size for kriging, format c(10,10), ie, 10x10 grid 
    # grid.coords = user defined coordinate points
    # shapefile = user defined shapefile for kriging
    # model =  "Exp", "Sph", "Gau"
    # plot.global = if TRUE, plot variograms
    #
    if (!model %in% c("Exp", "Sph", "Gau")) {
      stop("\n Define model correctly: Exp, Sph or Gau \n")
    }
    #
    data <- data[order(data[, 4]), ]
    ll <- unique(data[, 4])
    tr <- length(unique(data[, 4]))
    if (tr > 3) {
      stop("\n Cann't take more than 3 treatments \n")
    }
    ck <- NULL
    for (i in 1:tr) {
      ck[i] <- length(data[data[, 4] == ll[i], 4])
    }
    coords <- as.matrix(data[, 1:2])
    dimnames(coords) <- NULL
    #
    x <- range(coords[, 1])
    y <- range(coords[, 2])
    if (is.null(grid.size) &
        is.null(grid.coords) & is.null(shapefile)) {
      stop("\n Provide a value for the argument grid.size, grid.coords or shapefile \n")
    }
    else if (!is.null(grid.size)) {
      gr <- spT.grid.coords(x, y, by = grid.size)
    }
    else if (!is.null(grid.coords)) {
      gr <- as.matrix(grid.coords)
    }
    else if (!is.null(shapefile)) {
      gr <- as.matrix(coordinates(shapefile))
    }
    else {
      stop("\n Provide a value for the argument grid.size, grid.coords or shapefile \n")
    }
    dimnames(gr) <- NULL
    # check the big-n problem
    if(nrow(data) > tol.site) {
      cat("\n For each treatment max tol.site = ",tol.site)
      stop("\n Model stopped due to the big-n problem, use the Local Kriging or check the tol.site argument \n")
    }
    # covariance function and kriging/cokriging
    sink(tempfile())
    if (tr == 1) {
      result <- .TR_1(
        data = data,
        grid = gr,
        model = model
      )
    }
    if (tr == 2) {
      result <- .TR_2(
        data = data,
        grid = gr,
        model = model,
        tr = tr
      )
    }
    if (tr == 3) {
      result <- .TR_3(
        data = data,
        grid = gr,
        model = model,
        tr = tr
      )
    }
    sink()
    do.call(file.remove, list(list.files(tempdir(), full.names = TRUE)))
    dimnames(gr)[[2]] <- names(data)[1:2]
    result <-
      list(
        kriged.output = as.data.frame(result[c(1:(length(result)-2))]),
        est.parameters = result$range_psill,
        est.sample = result$n.sample,
        grid.coordinates = gr,
        model = model,
        grid.size = grid.size,
        shapefile = shapefile,
        type="global"
      )
    if(tr == 1){
      names(result$kriged.output) <-
        c("tr.pred","tr.se")
    }	
    if(tr == 2){
      names(result$kriged.output) <-
        c("tr.01.pred","tr.01.se","tr.02.pred","tr.02.se","tr.12.se")
    }	
    if(tr == 3){
      names(result$kriged.output) <-
        c("tr.01.pred","tr.01.se","tr.02.pred","tr.02.se","tr.03.pred","tr.03.se","tr.12.se","tr.13.se","tr.23.se")
    }
    #
    # check the missing values
    #
    # infill with average smooth for missing values
    # 8 adjacent grids are used in this context
    #
    if(length(unlist(result$kriged.output[complete.cases(result$kriged.output),]))==length(unlist(result$kriged.output))){
      result$kriged.output <- result$kriged.output	
    }
    else{
      check <- as.data.frame(result$kriged.output)
      s <- as.numeric(row.names(check[!complete.cases(check),]))
      d.check <- rdist(result$grid.coordinates)
      #d.check <- d.check[s,]
      for(i in 1:length(s)){
        ss <- which(c(d.check[s[i],]) <= sort(d.check[s[i],])[9]+1)
        check[s[i],] <- apply(check[ss[!ss%in%s[i]],],2,mean,na.rm=TRUE)
      }
      rm(d.check); rm(s); rm(ss);
      result$kriged.output <- as.matrix(check)
    }
    #
    class(result) <- "cokriging"
    result
    #
  }
##
##
._Local_fnc <-
  function(data,
           model,
           grid.size=NULL,
           grid.coords=NULL,
           shapefile=NULL,
           pts, # this is for local approach to ensure the number of points in each treatment
           block.size,
           plot.local=FALSE,
           tol.site,
           ...) {
    #
    # grid.size = user defined grid.size for kriging, format c(10,10), ie, 10x10 grid 
    # grid.coords = user defined coordinate points
    # shapefile = user defined shapefile for kriging
    #
    # model =  "Exp", "Sph", "Gau"
    #
    # block.size => to be defined by user (default is a function of: (1) max-distance (2) range (3) no. of spatial points)
    # i.e., block.size = f(phi,maxd,n); phi = range, maxd = max-distance, n = no. of spatial points.
    # f(phi,maxd,n) = phi * (n/maxd)^(-1), if (n/maxd)^(-1) > 1
    #               = phi, otherwise
    #
    # For big-n problem: block.size = f(maxd,tol.site,n), i.e., block.size  <- max(c(d))*(tol.site/nrow(d))
    #
    # plot.local = TRUE will return the local variograms for each kriging points.
    #
    # tol.site = maximum number of observation locations to fit a global model
    #
    if (!model %in% c("Exp", "Sph", "Gau")) {
      stop("\n Define model correctly: Exp, Sph or Gau \n")
    }
    #
    data <- data[order(data[, 4]), ]
    ll <- unique(data[, 4])
    tr <- length(unique(data[, 4]))
    if (tr > 3) {
      stop("\n Cann't take more than 3 treatments \n")
    }
    ck <- NULL
    for (i in 1:tr) {
      ck[i] <- length(data[data[, 4] == ll[i], 4])
    }
    coords <- as.matrix(data[, 1:2])
    dimnames(coords) <- NULL
    #
    if (is.null(grid.size) &
        is.null(grid.coords) & is.null(shapefile)) {
      stop("\n Provide a value for the argument grid.size, grid.coords or shapefile \n")
    }
    else if (!is.null(grid.size)) {
      x <- range(coords[, 1])
      y <- range(coords[, 2])
      gr <- spT.grid.coords(x, y, by = grid.size)
    }
    else if (!is.null(grid.coords)) {
      gr <- as.matrix(grid.coords)
    }
    else if (!is.null(shapefile)) {
      gr <- as.matrix(coordinates(shapefile))
    }
    else {
      stop("\n Provide a value for the argument grid.size, grid.coords or shapefile \n")
    }
    dimnames(gr) <- NULL
    d <- rdist(coords)
    d.gr <- rdist(coords, gr) # rdist is from fields package
    # check big-n problem
    if(nrow(data) > tol.site) {
      print(paste0("Number of data point exceeds 'tol.site=",tol.site,"'. Initiates the big-n problem."))
    }
    #
    result <- ._local_fnc_bign(data=data,gr=gr,model=model,tr=tr,ck=ck,d=d,d.gr=d.gr,tol.site=tol.site,block.size=block.size,plot.local=plot.local,pts=pts)
    #
    result <-
      list(
        kriged.output = result$result$kriged.output,
        est.parameters = result$result$est.parameters,
        est.sample = result$result$est.sample,
        neighbourhood.size = result$block.size,
        #neighbourhood.points = result$neighbourhood.points,
        grid.coordinates = gr,
        model = model,
        grid.size = grid.size,
        shapefile = shapefile,
        type="local-adaptive"
      )
    class(result) <- "cokriging"
    result
  }
##
##
._Local_fnc_processed <-
  function(data,
           model,
           grid.size,
           grid.coords,
           shapefile,
           pts, # this is for local approach to ensure the number of points in each treatment
           block.size,
           plot.local,
           tol.site,
           ...) {
    #
    # grid.size = user defined grid.size for kriging, format c(10,10), ie, 10x10 grid 
    # grid.coords = user defined coordinate points
    # shapefile = user defined shapefile for kriging
    #
    # model =  "Exp", "Sph", "Gau"
    #
    # block.size => to be defined by user (default is a function of: (1) max-distance (2) range (3) no. of spatial points)
    # i.e., block.size = f(phi,maxd,n); phi = range, maxd = max-distance, n = no. of spatial points.
    # f(phi,maxd,n) = phi * (n/maxd)^(-1), if (n/maxd)^(-1) > 1
    #               = phi, otherwise
    #
    # For big-n problem: block.size = f(maxd,tol.site,n), i.e., block.size  <- max(c(d))*(tol.site/nrow(d))
    #
    # plot.local = TRUE will return the local variograms for each kriging points.
    #
    # tol.site = maximum number of observation locations to fit a global model
    #
    if (!model %in% c("Exp", "Sph", "Gau")) {
      stop("\n Define model correctly: Exp, Sph or Gau \n")
    }
    #
    data <- data[order(data[, 4]), ]
    ll <- unique(data[, 4])
    tr <- length(unique(data[, 4]))
    if (tr > 3) {
      stop("\n Cann't take more than 3 treatments \n")
    }
    ck <- NULL
    for (i in 1:tr) {
      ck[i] <- length(data[data[, 4] == ll[i], 4])
    }
    coords <- as.matrix(data[, 1:2])
    dimnames(coords) <- NULL
    #
    if (is.null(grid.size) &
        is.null(grid.coords) & is.null(shapefile)) {
      stop("\n Provide a value for the argument grid.size, grid.coords or shapefile \n")
    }
    else if (!is.null(grid.size)) {
      x <- range(coords[, 1])
      y <- range(coords[, 2])
      gr <- spT.grid.coords(x, y, by = grid.size)
    }
    else if (!is.null(grid.coords)) {
      gr <- as.matrix(grid.coords)
    }
    else if (!is.null(shapefile)) {
      gr <- as.matrix(coordinates(shapefile))
    }
    else {
      stop("\n Provide a value for the argument grid.size, grid.coords or shapefile \n")
    }
    dimnames(gr) <- NULL
    d <- rdist(coords)
    d.gr <- rdist(coords, gr) # rdist is from fields package
    #
    result <- ._local_fnc_bign_processed(data=data,gr=gr,model=model,tr=tr,ck=ck,d=d,d.gr=d.gr,tol.site=tol.site,block.size=block.size,plot.local=plot.local,pts=pts)
    #
    result <-
      list(
        kriged.output = result$result$kriged.output,
        est.parameters = result$result$est.parameters,
        est.sample = result$result$est.sample,
        neighbourhood.size = result$block.size,
        grid.coordinates = gr,
        model = model,
        grid.size = grid.size,
        shapefile = shapefile,
        type="local-processed"
      )
    class(result) <- "cokriging"
    result
  }
##
##
._Local_fnc_rp <-
  function(data,
           model,
           grid.size,
           grid.coords,
           shapefile,
           pts, # this is for local approach to ensure the number of points in each treatment
           block.size,
           plot.local,
           tol.site,
           ...) {
    #
    if (!model %in% c("Exp", "Sph", "Gau")) {
      stop("\n Define model correctly: Exp, Sph or Gau \n")
    }
    #
    data <- data[order(data[, 4]), ]
    ll <- unique(data[, 4])
    tr <- length(unique(data[, 4]))
    if (tr > 3) {
      stop("\n Cann't take more than 3 treatments \n")
    }
    ck <- NULL
    for (i in 1:tr) {
      ck[i] <- length(data[data[, 4] == ll[i], 4])
    }
    coords <- as.matrix(data[, 1:2])
    dimnames(coords) <- NULL
    #
    if (is.null(grid.size) &
        is.null(grid.coords) & is.null(shapefile)) {
      stop("\n Provide a value for the argument grid.size, grid.coords or shapefile \n")
    }
    else if (!is.null(grid.size)) {
      x <- range(coords[, 1])
      y <- range(coords[, 2])
      gr <- spT.grid.coords(x, y, by = grid.size)
    }
    else if (!is.null(grid.coords)) {
      gr <- as.matrix(grid.coords)
    }
    else if (!is.null(shapefile)) {
      gr <- as.matrix(coordinates(shapefile))
    }
    else {
      stop("\n Provide a value for the argument grid.size, grid.coords or shapefile \n")
    }
    dimnames(gr) <- NULL
    d <- rdist(coords)
    d.gr <- rdist(coords, gr) # rdist is from fields package
    #
    result <- ._local_fnc_bign_processed_rp(data=data,gr=gr,model=model,tr=tr,ck=ck,d=d,d.gr=d.gr,tol.site=tol.site,block.size=block.size,plot.local=plot.local,pts=pts)
    #
    result <-
      list(
        kriged.output = result$result$kriged.output,
        est.parameters = result$result$est.parameters,
        est.sample = result$result$est.sample,
        neighbourhood.size = result$block.size,
        neighbourhood.points = pts,
        grid.coordinates = gr,
        model = model,
        grid.size = grid.size,
        shapefile = shapefile,
        type="local-processed-rp"
      )
    class(result) <- "cokriging"
    result
  }  
##
##  
._Local_fnc_fixed <-
  function(data,
           model,
           grid.size=NULL,
           grid.coords=NULL,
           shapefile=NULL,
           block.size=NULL,
           neighbourhood.points=NULL,
           plot.local=NULL,
           tol.site=NULL,
           ...) {
    #
    if (!model %in% c("Exp", "Sph", "Gau")) {
      stop("\n Define model correctly: Exp, Sph or Gau \n")
    }
    #
    data <- data[order(data[, 4]), ]
    ll <- unique(data[, 4])
    tr <- length(unique(data[, 4]))
    if (tr > 3) {
      stop("\n Cann't take more than 3 treatments \n")
    }
    ck <- NULL
    for (i in 1:tr) {
      ck[i] <- length(data[data[, 4] == ll[i], 4])
    }
    coords <- as.matrix(data[, 1:2])
    dimnames(coords) <- NULL
    #
    if (is.null(grid.size) &
        is.null(grid.coords) & is.null(shapefile)) {
      stop("\n Provide a value for the argument grid.size, grid.coords or shapefile \n")
    }
    else if (!is.null(grid.size)) {
      x <- range(coords[, 1])
      y <- range(coords[, 2])
      gr <- spT.grid.coords(x, y, by = grid.size)
    }
    else if (!is.null(grid.coords)) {
      gr <- as.matrix(grid.coords)
    }
    else if (!is.null(shapefile)) {
      gr <- as.matrix(coordinates(shapefile))
    }
    else {
      stop("\n Provide a value for the argument grid.size, grid.coords or shapefile \n")
    }
    dimnames(gr) <- NULL
    d <- rdist(coords)
    d.gr <- rdist(coords, gr) # rdist is from fields package
    #
    if(length(neighbourhood.points)>1){
      result <- ._local_fnc_bign_p(data=data,gr=gr,model=model,tr=tr,ck=ck,d=d,d.gr=d.gr,tol.site=tol.site,block.size=NULL,neighbourhood.points = neighbourhood.points,plot.local=plot.local)
      result <-
        list(
          kriged.output = result$result$kriged.output,
          est.parameters = result$result$est.parameters,
          est.sample = result$result$est.sample,
          neighbourhood.size = result$block.size,
          neighbourhood.points = neighbourhood.points,
          grid.coordinates = gr,
          model = model,
          grid.size = grid.size,
          shapefile = shapefile,
          type="local-processed-p"
        )
    }
    else{
      result <- ._local_fnc_bign_fixed(data=data,gr=gr,model=model,tr=tr,ck=ck,d=d,d.gr=d.gr,tol.site=tol.site,block.size=NULL,neighbourhood.points = neighbourhood.points,plot.local=plot.local)
      #
      result <-
        list(
          kriged.output = result$result$kriged.output,
          est.parameters = result$result$est.parameters,
          est.sample = result$result$est.sample,
          neighbourhood.size = result$block.size,
          neighbourhood.points = neighbourhood.points,
          grid.coordinates = gr,
          model = model,
          grid.size = grid.size,
          shapefile = shapefile,
          type="local-fixed"
        )
    }
    #
    result
  }  
##
##
._local_fnc_bign <- function(data,gr,model,tr,ck,d,d.gr,tol.site,block.size,plot.local,pts){
  #
  result <- list()
  if (tr == 1) {
    result$kriged.output <- matrix(NA, nrow = ncol(d.gr), ncol = 2)
    dimnames(result$kriged.output)[[2]] <- c("tr.pred","tr.var")
    result$est.parameters <- list()
    result$est.sample <- matrix(NA, nrow = ncol(d.gr), ncol = tr)
    tt <- 5
    #tt <- 1
  }
  if (tr == 2) {
    result$kriged.output <- matrix(NA, nrow = ncol(d.gr), ncol = 5)
    dimnames(result$kriged.output)[[2]] <- c("tr.01.pred","tr.01.var","tr.02.pred","tr.02.var","tr.12.cov")
    result$est.parameters <- list()
    result$est.sample <- matrix(NA, nrow = ncol(d.gr), ncol = tr)
    tt <- 5*tr
    #tt <- 1
  }
  if (tr == 3) {
    result$kriged.output <- matrix(NA, nrow = ncol(d.gr), ncol = 9)
    dimnames(result$kriged.output)[[2]] <- c("tr.01.pred","tr.01.var","tr.02.pred","tr.02.var","tr.03.pred","tr.03.var","tr.12.cov","tr.13.cov","tr.23.cov")
    result$est.parameters <- list()
    result$est.sample <- matrix(NA, nrow = ncol(d.gr), ncol = tr)
    tt <- 5*tr
    #tt <- 1
  }
  #
  # dynamic block size 
  #
  if (is.null(block.size)) {
    if(nrow(d)>tol.site){
      set.seed(1234)
      dmin <- quantile(d[sample(x=1:nrow(d),size=10),],tt/nrow(d))
      increment <- sqrt((tr*max(max(d.gr)^2/ck)+pi*dmin^2)/pi)-dmin
      block.size <- dmin
      print(paste0("To address the big-n problem, we define initial 'neighbourhood.size' as: ",round(block.size,2)))
      print(paste0("The increment factor is estimated as: ",round(increment,2)))
    }
    else{
      dmin <- quantile(d[sample(x=1:nrow(d),size=10),],tt/nrow(d))
      increment <- sqrt((tr*max(max(d.gr)^2/ck)+pi*dmin^2)/pi)-dmin
      block.size <- dmin
      print(paste0("Initial 'neighbourhood.size' is defined as: ",round(block.size,2)))
      print(paste0("The increment factor is estimated as: ",round(increment,2)))
    }
  }
  else{
    dmin <- block.size
    increment <- sqrt((tr*(pi/4)*max(max(d.gr)^2/ck)+pi*dmin^2)/pi)-dmin
    block.size <- dmin
    print(paste0("Initial 'neighbourhood.size' is defined as: ",round(block.size,2)))
    print(paste0("The increment factor is estimated as: ",round(increment,2)))
  }
  #
  pb <-
    txtProgressBar(min = 0,
                   max = ncol(d.gr),
                   style = 3)   # set progress bar
  #
  store.block <- block.size
  store.block.size <- NULL
  if (plot.local == TRUE) {
    print("Error: 'plot.local' option is not avaiable for big-n problem")
  }
  #
  for (i in 1:ncol(d.gr)) { 
    # identify the block size 
    #browser()
    block.size <- store.block[1]	
    tmp <- which(d.gr[, i] <= block.size)
    if(length(tmp) <= (tr*30)){ # ensuring minimum (30*tr) observations in the small data 
      # use repeat function 
      j <- 0
      repeat{
        j <- j+1
        block.size <- block.size + increment
        tmp <- which(d.gr[, i] <= block.size)
        if (length(tmp) >= tr*30) {      
          break
        }
      }
      rm(j)
    }
    block.size <- .wrapper_check_bign(data=data,model=model,tr=tr,tmp=d.gr[,i],block.size=block.size,increment=increment)
    tmp <- which(d.gr[, i] <= block.size)
    # ensure at least 'pts=15' points are coming from each treatment
    if(pts <= 15){ pts <- 15}
    if(tr == 1){
      nb <- which(c(d.gr[1:ck[1],i]) <= sort(d.gr[1:ck[1],i])[pts])	# for treatment one: pts points 
    }
    else if(tr == 2){
      nb1 <- which(c(d.gr[1:ck[1],i]) <= sort(d.gr[1:ck[1],i])[pts])	# for treatment one: pts points 
      nb2 <- ck[1]+which(c(d.gr[(ck[1]+1):(ck[1]+ck[2]),i]) <= sort(d.gr[(ck[1]+1):(ck[1]+ck[2]),i])[pts]) # for treatment two: pts points
      nb <- c(nb1,nb2)
    }
    else{
      nb1 <- which(c(d.gr[1:ck[1],i]) <= sort(d.gr[1:ck[1],i])[pts])	# for treatment one 
      nb2 <- ck[1]+which(c(d.gr[(ck[1]+1):(ck[1]+ck[2]),i]) <= sort(d.gr[(ck[1]+1):(ck[1]+ck[2]),i])[pts])	# for treatment two 
      nb3 <- ck[1]+ck[2]+which(c(d.gr[(ck[1]+ck[2]+1):(ck[1]+ck[2]+ck[3]),i]) <= sort(d.gr[(ck[1]+ck[2]+1):(ck[1]+ck[2]+ck[3]),i])[pts])	# for treatment three 
      nb <- c(nb1,nb2,nb3)
    }
    dat <-  data[c(tmp,nb), ]
    dat <- unique(dat)
    dat <- dat[order(dat[,4]),]
    store.block.size <- c(store.block.size,block.size)
    tr.ck <- length(unique(dat[, 4]))
    ll <- unique(dat[, 4])
    ck.new <- NULL
    for (k in 1:tr.ck) {
      ck.new[k] <- length(dat[dat[, 4] == ll[k], 4])
    }
    ###################
    #tr.ck <- length(unique(dat[, 4]))
    # ensure all treatments are in the new dataset 
    #if(tr.ck < tr){
    #  dmin <- min(block.size)
    #  increment <- sqrt((tr*(pi/4)*max(max(d.gr)^2/ck)+pi*dmin^2)/pi)-dmin
    #  j <- 0
    #  repeat{
    #    j <- j+1
    #    block.size <- block.size + increment
    #    tmp <- which(d.gr[, i] <= block.size[i])
    #    dat <-  data[tmp, ]
    #    tr.ck <- length(unique(dat[, 4]))
    #    if (tr.ck == tr) {      
    #      break
    #    }
    #  }
    #  rm(j)
    #}
    #
    #ll <- unique(dat[, 4])
    #ck.new <- NULL
    #for (k in 1:tr.ck) {
    #  ck.new[k] <- length(dat[dat[, 4] == ll[k], 4])
    #}
    # ensuring at least 100 points for each treatment
    #pts <- pts
    #if(min(ck.new) < pts){
    #  dmin <- min(block.size)
    #  increment <- sqrt((tr*(pi/4)*max(max(d.gr)^2/ck)+pi*dmin^2)/pi)-dmin
    #  j <- 0
    #  repeat{
    #    j <- j+1
    #    block.size <- block.size + increment
    #    tmp <- which(d.gr[, i] <= block.size)
    #    dat <-  data[tmp, ]
    #    tr.ck <- length(unique(dat[, 4]))
    #    ll <- unique(dat[, 4])
    #    ck.new <- NULL
    #    for (k in 1:tr.ck) {
    #      ck.new[k] <- length(dat[dat[, 4] == ll[k], 4])
    #    }
    #    if (min(ck.new) >= pts) {      
    #      break
    #    }
    #  }
    #  rm(j)
    #}
    #
    #print(nrow(dat))
    #neighbourhood.points <- rbind(neighbourhood.points,table(dat[,4]))
    #store.block.size <- c(store.block.size,block.size)
    ###################
    #
    ##
    sink(tempfile())
    if (tr == 1) {
      re <- try(.TR_1(
        data = dat,
        grid = gr[i,],
        model = model
      ), TRUE)
      if (class(re) == "try-error") {
        ree <- .try_processed_TR1(i=i, data=data, grid=gr, model=model, tr=tr, ck=ck, d=d, d.gr=d.gr, block.size)
        re <- ree$out
        store.block.size[i] <- ree$block.size
      }
      result$kriged.output[i, ] <- unlist(re[c(1:(length(re)-2))])
      result$est.parameters[[i]] <- re$range_psill
      result$est.sample[i,] <- re$n.sample
    }
    if (tr == 2) {
      re <- try(.TR_2(
        data = dat,
        grid = gr[i,],
        model = model,
        tr = tr.ck
      ), TRUE)
      if (class(re) == "try-error") {
        ree <- .try_processed_TR2(i=i, data=data, grid=gr, model=model, tr=tr, ck=ck, d=d, d.gr=d.gr, block.size)
        re <- ree$out
        store.block.size[i] <- ree$block.size
      }
      result$kriged.output[i, ] <- unlist(re[c(1:(length(re)-2))])
      result$est.parameters[[i]] <- re$range_psill
      result$est.sample[i,] <- re$n.sample
    }
    if (tr == 3) {
      #browser()
      re <- try(.TR_3(
        data = dat,
        grid = gr[i,],
        model = model,
        tr = tr.ck
      ), TRUE)
      if (class(re) == "try-error") {
        ree <- .try_processed_TR3(i=i, data=data, grid=gr, model=model, tr=tr, ck=ck, d=d, d.gr=d.gr, block.size)
        re <- ree$out
        store.block.size[i] <- ree$block.size
      }
      result$kriged.output[i, ] <- unlist(re[c(1:(length(re)-2))])
      result$est.parameters[[i]] <- re$range_psill
      result$est.sample[i,] <- re$n.sample
    }
    sink()
    do.call(file.remove, list(list.files(tempdir(), full.names = TRUE)))
    #file.remove(tempdir())
    setTxtProgressBar(pb, i)
  }
  close(pb)
  #list(result=result,block.size=store.block.size,neighbourhood.points=neighbourhood.points)
  list(result=result,block.size=store.block.size)
}
##
##
._local_fnc_bign_processed <- function(data,gr,model,tr,ck,d,d.gr,tol.site,block.size,plot.local,pts){
  #
  result <- list()
  if (tr == 1) {
    result$kriged.output <- matrix(NA, nrow = ncol(d.gr), ncol = 2)
    dimnames(result$kriged.output)[[2]] <- c("tr.pred","tr.var")
    result$est.parameters <- list()
    result$est.sample <- matrix(NA, nrow = ncol(d.gr), ncol = tr)
    tt <- 5
    #tt <- 1
  }
  if (tr == 2) {
    result$kriged.output <- matrix(NA, nrow = ncol(d.gr), ncol = 5)
    dimnames(result$kriged.output)[[2]] <- c("tr.01.pred","tr.01.var","tr.02.pred","tr.02.var","tr.12.cov")
    result$est.parameters <- list()
    result$est.sample <- matrix(NA, nrow = ncol(d.gr), ncol = tr)
    tt <- 5*tr
    #tt <- 1
  }
  if (tr == 3) {
    result$kriged.output <- matrix(NA, nrow = ncol(d.gr), ncol = 9)
    dimnames(result$kriged.output)[[2]] <- c("tr.01.pred","tr.01.var","tr.02.pred","tr.02.var","tr.03.pred","tr.03.var","tr.12.cov","tr.13.cov","tr.23.cov")
    result$est.parameters <- list()
    result$est.sample <- matrix(NA, nrow = ncol(d.gr), ncol = tr)
    tt <- 5*tr
    #tt <- 1
  }
  #
  pb <-
    txtProgressBar(min = 0,
                   max = ncol(d.gr),
                   style = 3)   # set progress bar
  #
  if (plot.local == TRUE) {
    print("Warnings: 'plot.local' option is not avaiable for big-n problem")
  }
  #
  store.block.size <- NULL
  #
  for (i in 1:ncol(d.gr)) { 
    tmp <- which(d.gr[, i] <= block.size[i])
    if(length(tmp)==0){
      block.size[i] <- sort(d.gr[, i])[2]
      tmp <- which(d.gr[, i] <= block.size[i])
    }
    # ensure at least 'pts=15' points are coming from each treatment
    if(pts <= 15){ pts <- 15}
    if(tr == 1){
      nb <- which(c(d.gr[1:ck[1],i]) <= sort(d.gr[1:ck[1],i])[pts])	# for treatment one: pts points 
    }
    else if(tr == 2){
      nb1 <- which(c(d.gr[1:ck[1],i]) <= sort(d.gr[1:ck[1],i])[pts])	# for treatment one: pts points 
      nb2 <- ck[1]+which(c(d.gr[(ck[1]+1):(ck[1]+ck[2]),i]) <= sort(d.gr[(ck[1]+1):(ck[1]+ck[2]),i])[pts]) # for treatment two: pts points
      nb <- c(nb1,nb2)
    }
    else{
      nb1 <- which(c(d.gr[1:ck[1],i]) <= sort(d.gr[1:ck[1],i])[pts])	# for treatment one 
      nb2 <- ck[1]+which(c(d.gr[(ck[1]+1):(ck[1]+ck[2]),i]) <= sort(d.gr[(ck[1]+1):(ck[1]+ck[2]),i])[pts])	# for treatment two 
      nb3 <- ck[1]+ck[2]+which(c(d.gr[(ck[1]+ck[2]+1):(ck[1]+ck[2]+ck[3]),i]) <= sort(d.gr[(ck[1]+ck[2]+1):(ck[1]+ck[2]+ck[3]),i])[pts])	# for treatment three 
      nb <- c(nb1,nb2,nb3)
    }
    dat <-  data[c(tmp,nb), ]
    dat <- unique(dat)
    dat <- dat[order(dat[,4]),]
    tr.ck <- length(unique(dat[, 4]))
    store.block.size[i] <- block.size[i]
    ##
    sink(tempfile())
    if (tr == 1) {
      re <- try(.TR_1(
        data = dat,
        grid = gr[i,],
        model = model
      ), TRUE)
      if (class(re) == "try-error") {
        ree <- .try_processed_TR1(i=i, data=data, grid=gr, model=model, tr=tr, ck=ck, d=d, d.gr=d.gr, block.size[i])
        re <- ree$out
        store.block.size[i] <- ree$block.size
      }
      result$kriged.output[i, ] <- unlist(re[c(1:(length(re)-2))])
      result$est.parameters[[i]] <- re$range_psill
      result$est.sample[i,] <- re$n.sample
    }
    if (tr == 2) {
      re <- try(.TR_2(
        data = dat,
        grid = gr[i,],
        model = model,
        tr = tr.ck
      ), TRUE)
      if (class(re) == "try-error") {
        ree <- .try_processed_TR2(i=i, data=data, grid=gr, model=model, tr=tr, ck=ck, d=d, d.gr=d.gr, block.size[i])
        re <- ree$out
        store.block.size[i] <- ree$block.size
      }
      result$kriged.output[i, ] <- unlist(re[c(1:(length(re)-2))])
      result$est.parameters[[i]] <- re$range_psill
      result$est.sample[i,] <- re$n.sample
    }
    if (tr == 3) {
      #browser()
      re <- try(.TR_3(
        data = dat,
        grid = gr[i,],
        model = model,
        tr = tr.ck
      ), TRUE)
      if (class(re) == "try-error") {
        #browser()
        ree <- .try_processed_TR3(i=i, data=data, grid=gr, model=model, tr=tr, ck=ck, d=d, d.gr=d.gr, block.size[i])
        re <- ree$out
        store.block.size[i] <- ree$block.size
      }
      result$kriged.output[i, ] <- unlist(re[c(1:(length(re)-2))])
      result$est.parameters[[i]] <- re$range_psill
      result$est.sample[i,] <- re$n.sample
    }
    sink()
    do.call(file.remove, list(list.files(tempdir(), full.names = TRUE)))
    #file.remove(tempdir())
    setTxtProgressBar(pb, i)
  }
  close(pb)
  #list(result=result,block.size=store.block.size,neighbourhood.points=neighbourhood.points)
  list(result=result,block.size=store.block.size)
}
##
##
._local_fnc_bign_processed_rp <- function(data,gr,model,tr,ck,d,d.gr,tol.site,block.size,plot.local,pts){
  #
  result <- list()
  if (tr == 1) {
    result$kriged.output <- matrix(NA, nrow = ncol(d.gr), ncol = 2)
    dimnames(result$kriged.output)[[2]] <- c("tr.pred","tr.var")
    result$est.parameters <- list()
    result$est.sample <- matrix(NA, nrow = ncol(d.gr), ncol = tr)
    tt <- 5
  }
  if (tr == 2) {
    result$kriged.output <- matrix(NA, nrow = ncol(d.gr), ncol = 5)
    dimnames(result$kriged.output)[[2]] <- c("tr.01.pred","tr.01.var","tr.02.pred","tr.02.var","tr.12.cov")
    result$est.parameters <- list()
    result$est.sample <- matrix(NA, nrow = ncol(d.gr), ncol = tr)
    tt <- 5*tr
  }
  if (tr == 3) {
    result$kriged.output <- matrix(NA, nrow = ncol(d.gr), ncol = 9)
    dimnames(result$kriged.output)[[2]] <- c("tr.01.pred","tr.01.var","tr.02.pred","tr.02.var","tr.03.pred","tr.03.var","tr.12.cov","tr.13.cov","tr.23.cov")
    result$est.parameters <- list()
    result$est.sample <- matrix(NA, nrow = ncol(d.gr), ncol = tr)
    tt <- 5*tr
  }
  #
  pb <-
    txtProgressBar(min = 0,
                   max = ncol(d.gr),
                   style = 3)   # set progress bar
  #
  if (plot.local == TRUE) {
    print("Warnings: 'plot.local' option is not avaiable for big-n problem")
  }
  #
  store.block.size <- NULL
  #
  for (i in 1:ncol(d.gr)) { 
    tmp <- which(d.gr[, i] <= block.size[i])
    if(length(tmp)==0){
      block.size[i] <- sort(d.gr[, i])[2]
      tmp <- which(d.gr[, i] <= block.size[i])
    }
    # ensure at least 'pts=5' points are coming from each treatment
    if(pts[i] <= 5){ pts[i] <- 5}
    if(tr == 1){
      nb <- which(c(d.gr[1:ck[1],i]) <= sort(d.gr[1:ck[1],i])[pts[i]])	# for treatment one: pts[i] points 
    }
    else if(tr == 2){
      nb1 <- which(c(d.gr[1:ck[1],i]) <= sort(d.gr[1:ck[1],i])[pts[i]])	# for treatment one: pts[i] points 
      nb2 <- ck[1]+which(c(d.gr[(ck[1]+1):(ck[1]+ck[2]),i]) <= sort(d.gr[(ck[1]+1):(ck[1]+ck[2]),i])[pts[i]]) # for treatment two: pts[i] points
      nb <- c(nb1,nb2)
    }
    else{
      nb1 <- which(c(d.gr[1:ck[1],i]) <= sort(d.gr[1:ck[1],i])[pts[i]])	# for treatment one 
      nb2 <- ck[1]+which(c(d.gr[(ck[1]+1):(ck[1]+ck[2]),i]) <= sort(d.gr[(ck[1]+1):(ck[1]+ck[2]),i])[pts[i]])	# for treatment two 
      nb3 <- ck[1]+ck[2]+which(c(d.gr[(ck[1]+ck[2]+1):(ck[1]+ck[2]+ck[3]),i]) <= sort(d.gr[(ck[1]+ck[2]+1):(ck[1]+ck[2]+ck[3]),i])[pts[i]])	# for treatment three 
      nb <- c(nb1,nb2,nb3)
    }
    dat <-  data[c(tmp,nb), ]
    dat <- unique(dat)
    dat <- dat[order(dat[,4]),]
    tr.ck <- length(unique(dat[, 4]))
    store.block.size[i] <- block.size[i]
    ##
    sink(tempfile())
    if (tr == 1) {
      re <- try(.TR_1(
        data = dat,
        grid = gr[i,],
        model = model
      ), TRUE)
      if (class(re) == "try-error") {
        ree <- .try_processed_TR1(i=i, data=data, grid=gr, model=model, tr=tr, ck=ck, d=d, d.gr=d.gr, block.size[i])
        re <- ree$out
        store.block.size[i] <- ree$block.size
      }
      result$kriged.output[i, ] <- unlist(re[c(1:(length(re)-2))])
      result$est.parameters[[i]] <- re$range_psill
      result$est.sample[i,] <- re$n.sample
    }
    if (tr == 2) {
      re <- try(.TR_2(
        data = dat,
        grid = gr[i,],
        model = model,
        tr = tr.ck
      ), TRUE)
      if (class(re) == "try-error") {
        ree <- .try_processed_TR2(i=i, data=data, grid=gr, model=model, tr=tr, ck=ck, d=d, d.gr=d.gr, block.size[i])
        re <- ree$out
        store.block.size[i] <- ree$block.size
      }
      result$kriged.output[i, ] <- unlist(re[c(1:(length(re)-2))])
      result$est.parameters[[i]] <- re$range_psill
      result$est.sample[i,] <- re$n.sample
    }
    if (tr == 3) {
      #browser()
      re <- try(.TR_3(
        data = dat,
        grid = gr[i,],
        model = model,
        tr = tr.ck
      ), TRUE)
      if (class(re) == "try-error") {
        #browser()
        ree <- .try_processed_TR3(i=i, data=data, grid=gr, model=model, tr=tr, ck=ck, d=d, d.gr=d.gr, block.size[i])
        re <- ree$out
        store.block.size[i] <- ree$block.size
      }
      result$kriged.output[i, ] <- unlist(re[c(1:(length(re)-2))])
      result$est.parameters[[i]] <- re$range_psill
      result$est.sample[i,] <- re$n.sample
    }
    sink()
    do.call(file.remove, list(list.files(tempdir(), full.names = TRUE)))
    #file.remove(tempdir())
    setTxtProgressBar(pb, i)
  }
  close(pb)
  list(result=result,block.size=store.block.size)
}
##
## for fixed number of neighbourhood points
._local_fnc_bign_fixed <- function(data,gr,model,tr,ck,d,d.gr,tol.site,block.size,neighbourhood.points,plot.local){
  #
  result <- list()
  if (tr == 1) {
    result$kriged.output <- matrix(NA, nrow = ncol(d.gr), ncol = 2)
    dimnames(result$kriged.output)[[2]] <- c("tr.pred","tr.var")
    result$est.parameters <- list()
    result$est.sample <- matrix(NA, nrow = ncol(d.gr), ncol = tr)
    tt <- 5
    #tt <- 1
  }
  if (tr == 2) {
    result$kriged.output <- matrix(NA, nrow = ncol(d.gr), ncol = 5)
    dimnames(result$kriged.output)[[2]] <- c("tr.01.pred","tr.01.var","tr.02.pred","tr.02.var","tr.12.cov")
    result$est.parameters <- list()
    result$est.sample <- matrix(NA, nrow = ncol(d.gr), ncol = tr)
    tt <- 5*tr
    #tt <- 1
  }
  if (tr == 3) {
    result$kriged.output <- matrix(NA, nrow = ncol(d.gr), ncol = 9)
    dimnames(result$kriged.output)[[2]] <- c("tr.01.pred","tr.01.var","tr.02.pred","tr.02.var","tr.03.pred","tr.03.var","tr.12.cov","tr.13.cov","tr.23.cov")
    result$est.parameters <- list()
    result$est.sample <- matrix(NA, nrow = ncol(d.gr), ncol = tr)
    tt <- 5*tr
    #tt <- 1
  }
  #
  # fixed neighbourhood.points
  #
  if (is.null(neighbourhood.points)) {
    neighbourhood.points <- min(ck)-1
    if(length(neighbourhood.points) >= min(ck)){
      stop("neighbourhood.points should be less than the minimum number of points for each treatments.")
    }
    if((length(neighbourhood.points)*tr) >= tol.site){
      stop("(neighbourhood.points x number of treatments) should be less than the 'tol.site'.")
    }
    else{
      print(paste0("The 'neighbourhood.points' input for each treatment is selected: ",neighbourhood.points))
      print(paste0("Total number of points used for ",tr," treatments: ",neighbourhood.points*tr))
      #print(paste0("Please IGNORE any WARNINGS"))
    }
  }
  else{
    if(length(neighbourhood.points) >= min(ck)){
      stop("neighbourhood.points should be less than the minimum number of points for each treatments.")
    }
    if((length(neighbourhood.points)*tr) >= tol.site){
      stop("(neighbourhood.points x number of treatments) should be less than the 'tol.site'.")
    }
    else{
      print(paste0(" 'Neighbourhood.points' input is selected as: ",neighbourhood.points))
      print(paste0("Total number of treatment points used: ",neighbourhood.points*tr))
      #print(paste0("Please IGNORE any WARNINGS"))
    }
  }
  #
  pb <-
    txtProgressBar(min = 0,
                   max = ncol(d.gr),
                   style = 3)   # set progress bar
  #
  for (i in 1:ncol(d.gr)) { 
    #print(paste0("Counter local fixed: ",i))
    #browser()
    ##
    sink(tempfile())
    if (tr == 1) {
      # identify the neighbourhood points 
      nb <- which(c(d.gr[,i]) <= sort(d.gr[,i])[neighbourhood.points])	
      dat <-  data[nb, ]
      re <- try(.TR_1(
        data = dat,
        grid = gr[i,],
        model = model
      ), TRUE)
      if (class(re) == "try-error") {
        stop("Check the neighbourhood.points. Requires possible tuning.")
      }
      result$kriged.output[i, ] <- unlist(re[c(1:(length(re)-2))])
      result$est.parameters[[i]] <- re$range_psill
      result$est.sample[i,] <- re$n.sample
    }
    if (tr == 2) {
      # identify the neighbourhood points by treatments 
      nb1 <- which(c(d.gr[1:ck[1],i]) <= sort(d.gr[1:ck[1],i])[neighbourhood.points])	# for treatment one 
      nb2 <- ck[1]+which(c(d.gr[(ck[1]+1):(ck[1]+ck[2]),i]) <= sort(d.gr[(ck[1]+1):(ck[1]+ck[2]),i])[neighbourhood.points])	# for treatment two 
      nb <- c(nb1,nb2)
      dat <-  data[nb, ]
      dat <- dat[order(dat[,4]),]
      re <- try(.TR_2(
        data = dat,
        grid = gr[i,],
        model = model,
        tr = tr
      ), TRUE)
      if (class(re) == "try-error") {
        stop("Check the neighbourhood.points. Requires possible tuning.")
      }
      result$kriged.output[i, ] <- unlist(re[c(1:(length(re)-2))])
      result$est.parameters[[i]] <- re$range_psill
      result$est.sample[i,] <- re$n.sample
    }
    if (tr == 3) {
      # identify the neighbourhood points by treatments 
      nb1 <- which(c(d.gr[1:ck[1],i]) <= sort(d.gr[1:ck[1],i])[neighbourhood.points])	# for treatment one 
      nb2 <- ck[1]+which(c(d.gr[(ck[1]+1):(ck[1]+ck[2]),i]) <= sort(d.gr[(ck[1]+1):(ck[1]+ck[2]),i])[neighbourhood.points])	# for treatment two 
      nb3 <- ck[1]+ck[2]+which(c(d.gr[(ck[1]+ck[2]+1):(ck[1]+ck[2]+ck[3]),i]) <= sort(d.gr[(ck[1]+ck[2]+1):(ck[1]+ck[2]+ck[3]),i])[neighbourhood.points])	# for treatment three 
      nb <- c(nb1,nb2,nb3)
      dat <-  data[nb, ]
      dat <- dat[order(dat[,4]),]
      #browser()
      re <- try(.TR_3(
        data = dat,
        grid = gr[i,],
        model = model,
        tr = tr
      ), TRUE)
      if (class(re) == "try-error") {
        stop("Check the neighbourhood.points. Requires possible tuning.")
      }
      result$kriged.output[i, ] <- unlist(re[c(1:(length(re)-2))])
      result$est.parameters[[i]] <- re$range_psill
      result$est.sample[i,] <- re$n.sample
    }
    sink()
    do.call(file.remove, list(list.files(tempdir(), full.names = TRUE)))
    #file.remove(tempdir())
    setTxtProgressBar(pb, i)
  }
  close(pb)
  #file.remove(tempdir())
  #list(result=result,neighbourhood.points=neighbourhood.points)
  list(result=result)
}
##
## for fixed number of neighbourhood points from gam output 
._local_fnc_bign_p <- function(data,gr,model,tr,ck,d,d.gr,tol.site,block.size,neighbourhood.points,plot.local){
  #
  result <- list()
  if (tr == 1) {
    result$kriged.output <- matrix(NA, nrow = ncol(d.gr), ncol = 2)
    dimnames(result$kriged.output)[[2]] <- c("tr.pred","tr.var")
    result$est.parameters <- list()
    result$est.sample <- matrix(NA, nrow = ncol(d.gr), ncol = tr)
    tt <- 5
    #tt <- 1
  }
  if (tr == 2) {
    result$kriged.output <- matrix(NA, nrow = ncol(d.gr), ncol = 5)
    dimnames(result$kriged.output)[[2]] <- c("tr.01.pred","tr.01.var","tr.02.pred","tr.02.var","tr.12.cov")
    result$est.parameters <- list()
    result$est.sample <- matrix(NA, nrow = ncol(d.gr), ncol = tr)
    tt <- 5*tr
    #tt <- 1
  }
  if (tr == 3) {
    result$kriged.output <- matrix(NA, nrow = ncol(d.gr), ncol = 9)
    dimnames(result$kriged.output)[[2]] <- c("tr.01.pred","tr.01.var","tr.02.pred","tr.02.var","tr.03.pred","tr.03.var","tr.12.cov","tr.13.cov","tr.23.cov")
    result$est.parameters <- list()
    result$est.sample <- matrix(NA, nrow = ncol(d.gr), ncol = tr)
    tt <- 5*tr
    #tt <- 1
  }
  #
  pb <-
    txtProgressBar(min = 0,
                   max = ncol(d.gr),
                   style = 3)   # set progress bar
  #
  for (i in 1:ncol(d.gr)) { 
    ##
    #print(i)
    #print(neighbourhood.points[i])
    sink(tempfile())
    if (tr == 1) {
      # identify the neighbourhood points 
      nb <- which(c(d.gr[,i]) <= sort(d.gr[,i])[neighbourhood.points[i]])	
      dat <-  data[nb, ]
      re <- try(.TR_1(
        data = dat,
        grid = gr[i,],
        model = model
      ), TRUE)
      if (class(re) == "try-error") {
        stop("Check the neighbourhood.points[i]. Requires possible tuning. TR1")
      }
      result$kriged.output[i, ] <- unlist(re[c(1:(length(re)-2))])
      result$est.parameters[[i]] <- re$range_psill
      result$est.sample[i,] <- re$n.sample
    }
    if (tr == 2) {
      # identify the neighbourhood points by treatments 
      nb1 <- which(c(d.gr[1:ck[1],i]) <= sort(d.gr[1:ck[1],i])[neighbourhood.points[i]])	# for treatment one 
      nb2 <- ck[1]+which(c(d.gr[(ck[1]+1):(ck[1]+ck[2]),i]) <= sort(d.gr[(ck[1]+1):(ck[1]+ck[2]),i])[neighbourhood.points[i]])	# for treatment two 
      nb <- c(nb1,nb2)
      dat <-  data[nb, ]
      dat <- dat[order(dat[,4]),]
      re <- try(.TR_2(
        data = dat,
        grid = gr[i,],
        model = model,
        tr = tr
      ), TRUE)
      if (class(re) == "try-error") {
        stop("Check the neighbourhood.points[i]. Requires possible tuning. TR2")
      }
      result$kriged.output[i, ] <- unlist(re[c(1:(length(re)-2))])
      result$est.parameters[[i]] <- re$range_psill
      result$est.sample[i,] <- re$n.sample
    }
    if (tr == 3) {
      # identify the neighbourhood points by treatments 
      nb1 <- which(c(d.gr[1:ck[1],i]) <= sort(d.gr[1:ck[1],i])[neighbourhood.points[i]])	# for treatment one 
      nb2 <- ck[1]+which(c(d.gr[(ck[1]+1):(ck[1]+ck[2]),i]) <= sort(d.gr[(ck[1]+1):(ck[1]+ck[2]),i])[neighbourhood.points[i]])	# for treatment two 
      nb3 <- ck[1]+ck[2]+which(c(d.gr[(ck[1]+ck[2]+1):(ck[1]+ck[2]+ck[3]),i]) <= sort(d.gr[(ck[1]+ck[2]+1):(ck[1]+ck[2]+ck[3]),i])[neighbourhood.points[i]])	# for treatment three 
      nb <- c(nb1,nb2,nb3)
      dat <-  data[nb, ]
      dat <- dat[order(dat[,4]),]
      #browser()
      re <- try(.TR_3(
        data = dat,
        grid = gr[i,],
        model = model,
        tr = tr
      ), TRUE)
      if (class(re) == "try-error") {
        stop("Check the neighbourhood.points[i]. Requires possible tuning. TR3")
      }
      result$kriged.output[i, ] <- unlist(re[c(1:(length(re)-2))])
      result$est.parameters[[i]] <- re$range_psill
      result$est.sample[i,] <- re$n.sample
    }
    sink()
    do.call(file.remove, list(list.files(tempdir(), full.names = TRUE)))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  list(result=result)
}
##
## function inside try-error for local-processed
.try_processed_TR1 <- function(i, data, grid, model, tr, ck, d, d.gr, block.size){
  tt <- 5
  dmin <- quantile(d[sample(x=1:nrow(d),size=10),],tt/nrow(d))
  increment <- sqrt((tr*max(max(d.gr)^2/ck)+pi*dmin^2)/pi)-dmin
  block.size0 <- .wrapper_check_bign(data=data,model=model,tr=tr,tmp=d.gr[,i],block.size=block.size[i],increment=increment)
  #print(block.size0)
  tmp <- which(d.gr[, i] <= block.size0)
  dat <-  data[tmp, ]
  out <- NULL
  out$block.size <- block.size0
  out$out <- .TR_1(
    data = dat,
    grid = grid[i,],
    model = model
  )
  out
}
##
##
.try_processed_TR2 <- function(i, data, grid, model, tr, ck, d, d.gr, block.size){
  tt <- 5*tr
  #tt <- 1
  #dmin <- quantile(d[sample(x=1:nrow(d),size=10),],0.05)
  dmin <- quantile(d[sample(x=1:nrow(d),size=10),],tt/nrow(d))
  #increment <- sqrt((tr*(pi/4)*max(max(d.gr)^2/ck)+pi*dmin^2)/pi)-dmin
  increment <- sqrt((tr*max(max(d.gr)^2/ck)+pi*dmin^2)/pi)-dmin
  block.size0 <- .wrapper_check_bign(data=data,model=model,tr=tr,tmp=d.gr[,i],block.size=block.size[i],increment=increment)
  #print(block.size0)
  tmp <- which(d.gr[, i] <= block.size0)
  dat <-  data[tmp, ]
  tr.ck <- length(unique(dat[, 4]))
  ll <- unique(dat[, 4])
  ck.new <- NULL
  for (k in 1:tr.ck) {
    ck.new[k] <- length(dat[dat[, 4] == ll[k], 4])
  }
  out <- NULL
  out$block.size <- block.size0
  out$out <- .TR_2(
    data = dat,
    grid = grid[i,],
    model = model,
    tr = tr.ck
  )
  out
}
##
##
.try_processed_TR3 <- function(i, data, grid, model, tr, ck, d, d.gr, block.size){
  tt <- 5*tr
  #tt <- 1
  #dmin <- quantile(d[sample(x=1:nrow(d),size=10),],0.05)
  dmin <- quantile(d[sample(x=1:nrow(d),size=10),],tt/nrow(d))
  increment <- sqrt((tr*max(max(d.gr)^2/ck)+pi*dmin^2)/pi)-dmin
  #browser()
  block.size0 <- .wrapper_check_bign(data=data,model=model,tr=tr,tmp=d.gr[,i],block.size=block.size[i],increment=increment)
  #print(block.size0)
  tmp <- which(d.gr[, i] <= block.size0)
  dat <-  data[tmp, ]
  tr.ck <- length(unique(dat[, 4]))
  ll <- unique(dat[, 4])
  ck.new <- NULL
  for (k in 1:tr.ck) {
    ck.new[k] <- length(dat[dat[, 4] == ll[k], 4])
  }
  out <- NULL
  out$block.size <- block.size0
  out$out <- .TR_3(
    data = dat,
    grid = grid[i,],
    model = model,
    tr = tr.ck
  )
  out
}
################
## sub-routines
################
##
##
.wrapper_cov <- function(data, model) {
  dat <- data
  coordinates(dat) <-
    as.formula(paste0("~", names(data)[1], "+", names(data)[2]))
  v0 <- variogram(as.formula(paste0(names(data)[3], "~1")), data = dat)
  v <- fit.variogram(v0, vgm(model = model))
  #v <- fit.variogram(v0, vgm(psill=1,model = model,nugget=1))
  v
}
##
##
.wrapper_check_bign <- function(data,model,tr,tmp,block.size,increment=NULL){
  if(is.null(increment)){
    increment.factor <- c(seq(1.0,10,0.1))
  }
  else{
    increment.factor <- c(0,rep(increment,99)*1:99)
  }
  #  
  options(warn=-1)
  .wrapper_check_internal <- function(data, model,tr) {
    ll <- unique(data[, 4])
    if(tr ==1){
      dat <- data
      coordinates(dat) <-
        as.formula(paste0("~", names(data)[1], "+", names(data)[2]))
      v0 <- variogram(as.formula(paste0(names(data)[3], "~1")), data = dat)
      v <- fit.variogram(v0, vgm(model = model))
      # range from v 
      # max dist for the data points from v0
      # check range(v) > mad.dist*0.8 => stop 
      v <- unlist(v$model)
      if(max(v[grep("range",names(v))]) > (0.8*max(v0$dist))){
        stop("")
      }
      rm(v0); rm(v); rm(dat)
    }
    else if (tr == 2) {
      dat <- data
      g.dat <- list()
      g.dat[[1]] <- subset(dat, dat[,4]==ll[1])
      names(g.dat[[1]])[3] <- paste0(names(g.dat[[1]])[3],ll[1])
      coordinates(g.dat[[1]]) <-as.formula(paste0("~", names(g.dat[[1]])[1], "+", names(g.dat[[1]])[2]))
      g.dat[[2]] <- subset(dat, dat[,4]==ll[2])
      names(g.dat[[2]])[3] <- paste0(names(g.dat[[2]])[3],ll[2])
      coordinates(g.dat[[2]]) <-as.formula(paste0("~", names(g.dat[[2]])[1], "+", names(g.dat[[2]])[2]))
      g <- gstat(id="tr1",formula=as.formula(paste0(names(g.dat[[1]])[1], "~1")),data=g.dat[[1]],set = list(nocheck = 1))
      g <- gstat(g,id="tr2",formula=as.formula(paste0(names(g.dat[[2]])[1], "~1")),data=g.dat[[2]],set = list(nocheck = 1))
      rm(g.dat); rm(dat)
      v0 <- variogram(g)
      v <- fit.lmc(v=v0, g=g, vgm(model = model))
      # range from v 
      # max dist for the data points from v0
      # check range(v) > mad.dist*0.8 => stop 
      v <- unlist(v$model)
      if(max(v[grep("range",names(v))]) > (0.8*max(v0$dist))){
        stop("")
      }
      rm(g); rm(v0); rm(v);   
    }
    else if (tr == 3) {
      dat <- data
      g.dat <- list()
      g.dat[[1]] <- subset(dat, dat[,4]==ll[1])
      names(g.dat[[1]])[3] <- paste0(names(g.dat[[1]])[3],ll[1])
      coordinates(g.dat[[1]]) <-as.formula(paste0("~", names(g.dat[[1]])[1], "+", names(g.dat[[1]])[2]))
      g.dat[[2]] <- subset(dat, dat[,4]==ll[2])
      names(g.dat[[2]])[3] <- paste0(names(g.dat[[2]])[3],ll[2])
      coordinates(g.dat[[2]]) <-as.formula(paste0("~", names(g.dat[[2]])[1], "+", names(g.dat[[2]])[2]))
      g.dat[[3]] <- subset(dat, dat[,4]==ll[3])
      names(g.dat[[3]])[3] <- paste0(names(g.dat[[3]])[3],ll[3])
      coordinates(g.dat[[3]]) <-as.formula(paste0("~", names(g.dat[[3]])[1], "+", names(g.dat[[3]])[2]))
      g <- gstat(id="tr1",formula=as.formula(paste0(names(g.dat[[1]])[1], "~1")),data=g.dat[[1]],set = list(nocheck = 1))
      g <- gstat(g,id="tr2",formula=as.formula(paste0(names(g.dat[[2]])[1], "~1")),data=g.dat[[2]],set = list(nocheck = 1))
      g <- gstat(g,id="tr3",formula=as.formula(paste0(names(g.dat[[3]])[1], "~1")),data=g.dat[[3]],set = list(nocheck = 1))
      rm(g.dat); rm(dat)
      v0 <- variogram(g) 
      v <- fit.lmc(v=v0, g=g, vgm(model = model))
      # range from v 
      # max dist for the data points from v0
      # check range(v) > mad.dist*0.8 => stop 
      v <- unlist(v$model)
      if(max(v[grep("range",names(v))]) > (0.8*max(v0$dist))){ 
        stop("")
      }
      rm(g); rm(v0); rm(v);   
    }
    else{
      stop("\n Treatments should be <= 3. \n")
    }
  }
  ##
  .internal_fnc <- function(data,model,tr,tmp,block.size0){
    wp <- which(c(tmp) <= block.size0)
    if(length(wp)==0){
      block.size0 <- sort(c(tmp))[2]
      wp <- which(c(tmp) <= block.size0)
    }
    wp <- .wrapper_check_internal(data=data[wp,], model=model, tr=tr)
    wp
  }
  ##
  j <- 0
  repeat{
    j <- j+1
    block.size0 <- block.size + increment.factor[j]
    #print(j)
    #print(block.size0)
    #print(increment.factor[j])
    #browser()
    bign.check <- try(.internal_fnc(data=data,model=model,tr=tr,tmp=tmp,block.size0=block.size0),TRUE)
    if (class(bign.check)=="NULL") {      
      #print(block.size0)
      #print(j)
      break
    }
  }
  #print(block.size0)
  block.size0
}		
##
##
.TR_1 <- function(data, grid, model) {
  g.dat <- as.data.frame(data)
  coordinates(g.dat) <- as.formula(paste0("~", names(g.dat)[1], "+", names(g.dat)[2]))
  v0 <- variogram(as.formula(paste0(names(g.dat)[1], "~1")),g.dat)
  v <- fit.variogram(v0,model=vgm(model=model))
  if(is.matrix(grid) | is.data.frame(grid)){
    grid <- as.data.frame(grid)
    coordinates(grid) <- as.formula(paste0("~",names(grid)[1],"+",names(grid)[2]))
  }
  else{
    grid <- data.frame(matrix(c(grid),1,2))
    coordinates(grid) <- as.formula(paste0("~",names(grid)[1],"+",names(grid)[2]))
  }
  out <- krige(formula=as.formula(paste0(names(g.dat)[1], "~1")),g.dat, grid, model=v)
  list(out.pred = out$var1.pred, out.var = out$var1.var, range_psill=v, n.sample=nrow(data))
}
##
##
.TR_2 <- function(data, grid, model, tr) {
  ## two parameter setting
  nBins    <- 20;      # keep same with the one within vesper
  propMaxD <- 0.333333 # 1/3 is recommended by gstat, and verified by several data sets we have tested
  ##
  g.dat <- list()
  ll <- unique(data[, 4])
  for(j in 1:tr){
    g.dat[[j]] <- subset(data,data[,4]==ll[j])
    names(g.dat[[j]])[3] <- paste0(names(g.dat[[j]])[3],ll[j])
    coordinates(g.dat[[j]]) <-
      as.formula(paste0("~", names(g.dat[[j]])[1], "+", names(g.dat[[j]])[2]))
  }
  g <- gstat(id="tr1",formula=as.formula(paste0(names(g.dat[[1]])[1], "~1")),data=g.dat[[1]],set = list(nocheck = 1))
  g <- gstat(g,id="tr2",formula=as.formula(paste0(names(g.dat[[2]])[1], "~1")),data=g.dat[[2]],set = list(nocheck = 1))
  rm(g.dat)
  ck <- NULL
  for (i in 1:tr) {
    ck[i] <- length(data[data[, 4] == ll[i], 4])
  }
  max.dist <- max(c(rdist(data[,1:2])))*propMaxD;  
  v0 <- variogram(g,cutoff=max.dist,width=max.dist/nBins); 
  ##print(sprintf('print(sum(is propMaxD=%.4f nBins=%d',propMaxD, nBins))#for debugging only
  if(length(unique(v0$id))!=3){
    v0 <- variogram(g,cutoff=max(c(rdist(data[,1:2])))/2)
    if(length(unique(v0$id))!=3){
      v0 <- variogram(g,cutoff=max(c(rdist(data[,1:2]))))
      if(length(unique(v0$id))!=3){
        stop("Error: An expert choice is required to identify the cutoff for the variogram ...")
      }
    }
  }
  v <- fit.lmc(v=v0, g=g, vgm(model = model))
  if(is.matrix(grid) | is.data.frame(grid)){
    grid <- as.data.frame(grid)
    coordinates(grid) <- as.formula(paste0("~",names(grid)[1],"+",names(grid)[2]))
  }
  else{
    grid <- data.frame(matrix(c(grid),1,2))
    coordinates(grid) <- as.formula(paste0("~",names(grid)[1],"+",names(grid)[2]))
  }
  out <- predict(v, grid)
  list(
    out.01.pred = out$tr1.pred,
    out.01.var = out$tr1.var,
    out.02.pred = out$tr2.pred,
    out.02.var = out$tr2.var,
    out.12.cov = out$cov.tr1.tr2,
    range_psill=v$model,
    n.sample=ck
  )
}
##
##
.TR_3 <- function(data, grid, model, tr) {
  ## two parameter setting
  nBins    <- 20;      # keep same with the one within vesper
  propMaxD <- 0.333333 # 1/3 is recommended by gstat, and verified by several data sets we have tested
  ##
  g.dat <- list()
  ll <- unique(data[, 4])
  for(j in 1:tr){
    g.dat[[j]] <- subset(data,data[,4]==ll[j])
    names(g.dat[[j]])[3] <- paste0(names(g.dat[[j]])[3],ll[j])
    coordinates(g.dat[[j]]) <-
      as.formula(paste0("~", names(g.dat[[j]])[1], "+", names(g.dat[[j]])[2]))
  }
  g <- gstat(id="tr1",formula=as.formula(paste0(names(g.dat[[1]])[1], "~1")),data=g.dat[[1]],set = list(nocheck = 1))
  g <- gstat(g,id="tr2",formula=as.formula(paste0(names(g.dat[[2]])[1], "~1")),data=g.dat[[2]],set = list(nocheck = 1))
  g <- gstat(g,id="tr3",formula=as.formula(paste0(names(g.dat[[3]])[1], "~1")),data=g.dat[[3]],set = list(nocheck = 1))
  rm(g.dat)
  ck <- NULL
  for (k in 1:tr) {
    ck[k] <- length(data[data[, 4] == ll[k], 4])
  }
  max.dist <- max(c(rdist(data[,1:2])))*propMaxD;  v0 <- variogram(g,cutoff=max.dist,width=max.dist/nBins); 
  #browser()
  if(length(unique(v0$id))!=6){
    v0 <- variogram(g,cutoff=max(c(rdist(data[,1:2])))/2)
    if(length(unique(v0$id))!=6){
      v0 <- variogram(g,cutoff=max(c(rdist(data[,1:2]))))
      if(length(unique(v0$id))!=6){
        stop("Error: An expert choice is required to identify the cutoff for the variogram ...")
      }
    }
  }
  v <- fit.lmc(v=v0, g=g, vgm(model = model))
  if(is.matrix(grid) | is.data.frame(grid)){
    grid <- as.data.frame(grid)
    coordinates(grid) <- as.formula(paste0("~",names(grid)[1],"+",names(grid)[2]))
  }
  else{
    grid <- data.frame(matrix(c(grid),1,2))
    coordinates(grid) <- as.formula(paste0("~",names(grid)[1],"+",names(grid)[2]))
  }
  out <- predict(v, grid)
  list(
    out.01.pred = out$tr1.pred,
    out.01.var = out$tr1.var,
    out.02.pred = out$tr2.pred,
    out.02.var = out$tr2.var,
    out.03.pred = out$tr3.pred,
    out.03.var = out$tr3.var,
    out.12.cov = out$cov.tr1.tr2,
    out.13.cov = out$cov.tr1.tr3,
    out.23.cov = out$cov.tr2.tr3,
    range_psill=v$model,
    n.sample=ck
  )
}
##
## convert seconds into min. hour. and day
##
.fnc.time_<-function(t)
{
  #
  if(t < 60){
    t <- round(t,2)
    tt <- paste(t," - Sec.")
    cat(paste("##\n# Elapsed time:",t,"Sec.\n##\n"))
  } 
  #
  if(t < (60*60) && t >= 60){
    t1 <- as.integer(t/60)
    t <- round(t-t1*60,2) 
    tt <- paste(t1," - Mins.",t," - Sec.")
    cat(paste("##\n# Elapsed time:",t1,"Min.",t,"Sec.\n##\n"))
  }
  #
  if(t < (60*60*24) && t >= (60*60)){
    t2 <- as.integer(t/(60*60))
    t <- t-t2*60*60
    t1 <- as.integer(t/60)
    t <- round(t-t1*60,2) 
    tt <- paste(t2," - Hour/s.",t1," - Mins.",t," - Sec.")
    cat(paste("##\n# Elapsed time:",t2,"Hour/s.",t1,"Min.",t,"Sec.\n##\n"))
  }
  #
  if(t >= (60*60*24)){
    t3 <- as.integer(t/(60*60*24))
    t <- t-t3*60*60*24
    t2 <- as.integer(t/(60*60))
    t <- t-t2*60*60
    t1 <- as.integer(t/60)
    t <- round(t-t1*60,2)
    tt <- paste(t3," - Day/s.",t2," - Hour/s.",t1," - Mins.",t," - Sec.")
    cat(paste("##\n# Elapsed time:",t3,"Day/s.",t2,"Hour/s.",t1,"Mins.",t,"Sec.\n##\n"))
  }
  #
  tt
}
##
## plot function
##
plot.pa <- function(x, treatment=FALSE, differences=FALSE, zval=FALSE, variance=FALSE, covariance=FALSE, length.out=10, ...){
  .wrapper.spplot <- function(dat,zcol,main,cuts,col.regions, ...){
    spplot(dat[,zcol], col="transparent", scales=list(draw=T),main=main,cuts=cuts,key.space = "right",do.log = TRUE, col.regions=col.regions, ...) 
  }
  if(ncol(x$kriged.output)==15){
    tr <- c(colnames(x$kriged.output)[1],colnames(x$kriged.output)[3],colnames(x$kriged.output)[5])
    tr <- c(substr(tr[1],4,nchar(tr[1])),substr(tr[2],4,nchar(tr[2])),substr(tr[3],4,nchar(tr[3])))
    dat <- as.data.frame(cbind(x$grid.coordinates,x$kriged.output))
    coordinates(dat) <- names(dat[,1:2])
    if(isTRUE(treatment)){
      ct <- round(seq(0,ceiling(max(unlist(data.frame(dat[,c(1,3,5)])[,1:3]))),length.out=length.out+1),2)
      .wrapper.spplot(dat,zcol=c(1,3,5),main="Treatments",cuts=ct,col.regions = viridis(length(ct)))
    }
    else if(isTRUE(differences)){
      ct <- round(seq(floor(min(unlist(data.frame(dat[,c(7,9,11)])[,1:3]))),ceiling(max(unlist(data.frame(dat[,c(7,9,11)])[,1:3]))),length.out=length.out+1),2)
      .wrapper.spplot(dat,zcol=c(7,9,11),main="Treatment Differences",cuts=ct,col.regions = viridis(length(ct)))
    }
    else if(isTRUE(zval)){
      #ct <- round(seq(floor(min(unlist(data.frame(dat[,c(13,14,15)])[,1:3]))),ceiling(max(unlist(data.frame(dat[,c(13,14,15)])[,1:3]))),length.out=length.out+1),2)
	  ct <- c(-5000,-1.96,-1.645,1.6451,1.961,5000)
      .wrapper.spplot(dat,zcol=c(13,14,15),main="Z-value for Treatment Differences",cuts=ct,col.regions = viridis(length(ct)),
                legendEntries=c("< -1.96","[-1.96,-1.645)","[-1.645,1.645]","(1.645,1.96]","> 1.96"))
    }
    else if(isTRUE(variance)){
      ct <- round(seq(floor(min(unlist(data.frame(dat[,c(2,4,6)])[,1:3]))),ceiling(max(unlist(data.frame(dat[,c(2,4,6)])[,1:3]))),length.out=length.out+1),2)
      .wrapper.spplot(dat,zcol=c(2,4,6),main="Variances",cuts=ct,col.regions = viridis(length(ct)))
    }
    else if(isTRUE(covariance)){
      ct <- round(seq(floor(min(unlist(data.frame(dat[,c(8,10,12)])[,1:3]))),ceiling(max(unlist(data.frame(dat[,c(8,10,12)])[,1:3]))),length.out=length.out+1),2)
      .wrapper.spplot(dat,zcol=c(8,10,12),main="Covariances",cuts=ct,col.regions = viridis(length(ct)))
    }
    else{
      stop("Correctly define the input. e.g., treatment=TRUE or differences=TRUE or zval=TRUE or variance=TRUE or covariance=TRUE ...")
    }
  }
  else if(ncol(x$kriged.output)==8){
    tr <- c(colnames(x$kriged.output)[1],colnames(x$kriged.output)[3])
    tr <- c(substr(tr[1],4,nchar(tr[1])),substr(tr[2],4,nchar(tr[2])))
    dat <- as.data.frame(cbind(x$grid.coordinates,x$kriged.output))
    coordinates(dat) <- names(dat[,1:2])
    if(isTRUE(treatment)){
      ct <- round(seq(0,ceiling(max(unlist(data.frame(dat[,c(1,3)])[,1:2]))),length.out=length.out+1),2)
      .wrapper.spplot(dat,zcol=c(1,3),main="Treatments",cuts=ct,col.regions = viridis(length(ct)))
    }
    else if(isTRUE(differences)){
      ct <- round(seq(floor(min(unlist(data.frame(dat[,c(5)])[,1]))),ceiling(max(unlist(data.frame(dat[,c(5)])[,1]))),length.out=length.out+1),2)
      .wrapper.spplot(dat,zcol=c(5),main="Treatment Differences",cuts=ct,col.regions = viridis(length(ct)))
    }
    else if(isTRUE(zval)){
 	  ct <- c(-5000,-1.96,-1.645,1.645,1.96,5000)
      .wrapper.spplot(dat,zcol=c(8),main="Z-value for Treatment Differences",cuts=ct,col.regions = rev(gray(1:length(ct)/length(ct))))
    }
    else if(isTRUE(variance)){
      ct <- round(seq(floor(min(unlist(data.frame(dat[,c(2,4)])[,1:2]))),ceiling(max(unlist(data.frame(dat[,c(2,4)])[,1:2]))),length.out=length.out+1),2)
      .wrapper.spplot(dat,zcol=c(2,4),main="Variances",cuts=ct,col.regions = viridis(length(ct)))
    }
    else if(isTRUE(covariance)){
      ct <- round(seq(floor(min(unlist(data.frame(dat[,c(6)])[,1]))),ceiling(max(unlist(data.frame(dat[,c(6)])[,1]))),length.out=length.out+1),2)
      .wrapper.spplot(dat,zcol=c(6),main="Coariance",cuts=ct,col.regions = viridis(length(ct)))
    }
    else{
      stop("Correctly define the input. e.g., treatment=TRUE or differences=TRUE or zval=TRUE or variance=TRUE or covariance=TRUE ...")
    }
  }
  else{
    stop("Check the treatment numbers. Currently two or three treatments are allowed.")
  }
}
##
## print function
##
print.pa<-function(x, ...) {
  if(x$run.model=="G"){ m <- paste0("Global model.")}
  if(x$run.model=="LA-R-A"){ m <- paste0("Local model with adaptive neighbourhood radius algorithm.")}
  if(x$run.model=="LA-R"){ m <- paste0("Local model with adaptive neighbourhood range using cross-validation.")}
  if(x$run.model=="LA-P"){ m <- paste0("Local model with adaptive neighbourhood points using cross-validation.")}
  if(x$run.model=="LA-RP"){ m <- paste0("Local model with adaptive neighbourhood range and points using cross-validation.")}
  if(x$run.model=="LF"){ m <- paste0("Local model with fixed neighbourhood points. \n",x$neighbourhood.points," points used for each treatment.")}
  if(x$run.model=="SR"){ m <- paste0("Model with random subsampling. Not a random replication method.")}
  if(x$run.model=="SK"){ m <- paste0("Model with subsampling based on k-means algorithm.\nUse no.of.replication argument to increase the replications.")}
  if(x$run.model=="SB"){ m <- paste0("Model with subsampling based on bootstrapping algorithm.\nUse no.of.replication argument to increase the replications.")}
  cat("---------------------------------------------------------"); cat('\n');
  cat("Model: "); cat(x$run.model); cat('\n');
  #cat("---------------------------------------------------------"); cat('\n');
  cat(m);  cat('\n');
  cat("---------------------------------------------------------"); cat('\n');
  cat("Computation time: "); cat(x$comp.time); cat("\n")
}
##
## 
##