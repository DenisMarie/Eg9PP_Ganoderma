#--Function for computing the neighborhood matrix

NeighMatrix <- function(fit =fit, type = "euclidian", nb = 4, level = "PAREXP"){
  if (type == "euclidian"){
    xymap <- createxymap(as.numeric(fit$x), as.numeric(fit$y),
                         districts = fit[[level]],
                         p=2, max.dist =nb)
    ##mapGra is used for Bayesx function
    mapGra = xymap$gra
    mapBnd =xymap$bnd
    colnames(mapGra)  = fit[[level]]
    rownames(mapGra)  = fit[[level]]
    nc.nb = gra2nb(mapGra)
    names(nc.nb) = fit[[level]]
    d <- sapply(nc.nb, length)                   # vector with number of neighbors
    C <- nb2mat(nc.nb, style="B") 
    colnames(C) <- fit[[level]]
    rownames(C) <- fit[[level]]
    return(list(C = C, mapBnd = mapBnd, mapGra = mapGra)  )
  }else{
    coords = coordinates(cbind(fit$x,fit$y))
    k <- knn2nb(knearneigh(coords, k=nb))
    C <- nb2mat(k, style="B")
    colnames(C) <- fit[[level]]
    rownames(C) <- fit[[level]]
    return(C)  # structure matrix
  }
  
}

# Function for computing the CAR matrix
#---------------------------------------
GMRF <- function(W , alpha = 0.80, fit, level){
  Minv <- diag(apply(W,1,sum),ncol=nrow(fit),nrow=nrow(fit))
  M <- solve(Minv)
  Id <- diag(1,nrow(fit))
  ##Ws is the scaled adjacency structure matrix
  Ws <- W/diag(Minv)
  B <- Minv%*%(Id - alpha*Ws)
  Binv <- solve(Id - alpha*Ws)%*%M
  colnames(Binv) <- fit[,colnames(fit)==level]
  rownames(Binv) <- fit[,colnames(fit)==level]
  return(list(Binv= Binv, Ws= Ws, M= M, Id=Id))
}
  

gexchange <- function(varlist,positive=TRUE,pdcheck=TRUE) {
  varlist <- varlist
  pdcheck <- pdcheck
  if (!is.logical(positive)) 
    stop("Invalid value for postive argument")
  positive <- positive
  geinit <- function (vinit, vfixed, intercept, G, X, sparse) 
  {
    vardefault <- list(c(0.02, 0.1, 0.4, 0.8)^2,c(0.1, 0.3, 0.5, 0.9))
    ngroup <- min(length(G), ncol(G))
    nvar <- min(length(X), ncol(X))
    if (ngroup > 0 & nvar > 0) 
      return(list(error = "Mlist cannot have both covariates and grouping"))
    if (!is.list(varlist)) 
      varlist <- list(varlist)
    noname <- all(sapply(varlist, function(x) is.null(dimnames(x)) || 
                           (is.null(dimnames(x)[[1]]) & is.null(dimnames(x)[[2]]))))
    namefun <- function(x, names) {
      if (all(dim(x) == rep(length(names), 2))) 
        dimnames(x) <- list(names, names)
      x
    }
    if (ngroup > 0) {
      n <- nrow(G)
      G <- expand.nested(G)
      groups <- G[[ngroup]]
      bname <- levels(groups)
      if (noname) 
        varlist <- lapply(varlist, namefun, bname)
      if (any(sapply(varlist, function(x) inherits(x, "Matrix")))) 
        varlist <- lapply(varlist, function(x) as(x, "bdsmatrix"))
      tlist <- bdsmatrix.reconcile(varlist, bname)
      imap <- matrix(match(groups, dimnames(tlist[[1]])[[1]]))
      xmap <- NULL
      rname <- names(G)[[ngroup]]
    }
    else {
      n <- nrow(X)
      bname <- dimnames(X)[[2]]
      if (noname) 
        varlist <- lapply(varlist, namefun, bname)
      tlist <- bdsmatrix.reconcile(varlist, bname)
      tlist <- lapply(tlist, as.matrix)
      xmap <- match(dimnames(X)[[2]], bname)
      xmap <- matrix(rep(xmap, n), nrow = n, byrow = T)
      imap <- NULL
      rname <- "(Shrink)"
    }
    ntheta <- length(varlist)*2
    itheta <- vector("list", ntheta)
    #if only one variable 
    for (i in 1:ntheta) itheta[[i]] <- vardefault[[i]]
    if (length(vinit) > 0) {
      if (length(vinit) != ntheta) 
        return(list(error = "Wrong length for initial values"))
      indx = match(names(vinit),c("sigma","rho"))
      #indx <- !is.na(vinit) & vinit != 0
      if (any(indx)) 
        itheta[indx] <- vinit[indx]
    }
    which.fixed <- rep(FALSE, ntheta)
    if (length(vfixed) > 0) {
      indx = match(names(vfixed),c("sigma","rho"))
      if (any(is.na(indx)))
        stop("Unrecognized parameter name in vfixed values")
      if (is.list(vfixed)) {
        if (any(sapply(vfixed, length) !=1))
          stop("Any fixed parameter must have a single value")
        vfixed <- unlist(vfixed)
      }
      which.fixed[indx] <- TRUE
      itheta[which.fixed] <-  vfixed[which(names(vfixed)%in%c("sigma", "rho"))] #vfixed[which.fixed]
    }
    if (any(itheta[[1]]<=0)) stop("Variance must be positive")
    if (any(itheta[[2]]<0) || any(itheta[[2]]>=1)) stop("Correlation must be between 0 and 1")
    
    itheta[[2]] <- itheta[[2]]/(1-itheta[[2]])
    for (j in 1:length(varlist)) {
      kmat <- tlist[[j]]
    }
    
    list(theta = lapply(itheta,log)[!which.fixed], imap = imap, X = X, xmap = xmap, 
         parms = list(varlist = tlist, theta = sapply(itheta,function(x)x[1]), fixed = which.fixed, 
                      bname = bname, rname = rname, positive = positive, 
                      vname = names(varlist)))
  }
  
  
  
  
  gegenerate = function (newtheta, parms) 
  {
    safe.exp <- function(x, emax=20){
      exp(pmax(-emax,pmin(emax,x)))
    }
    theta <- parms$theta
    if ( !all(parms$fixed))
      theta[!parms$fixed] <- safe.exp(newtheta)
    correlation = theta[2]/(1+theta[2])
    variance = min(theta[1],20)
    #try
    varmat = variance*solve(Id-correlation*Ws)%*%M
    varmat
  }
  
  
  
  
  gewrapup <- function (newtheta, b, parms) 
  {
    theta <- parms$theta
    theta[!parms$fixed] <- exp(newtheta)
    correlation = theta[2]/(1+theta[2])
    rtheta = list('Vmat'=c(variance=theta[1],
                           correlation = correlation))
    defaultname <- paste("Vmat", 1:length(theta), sep = ".")
    vname <- parms$vname
    if (length(vname) == 0) 
      vname <- defaultname
    else if (any(vname == "")) {
      indx <- which(vname == "")
      vname[indx] <- defaultname[indx]
    }
    names(theta) <- vname
    theta <- list(theta)
    names(theta) <- parms$rname
    names(rtheta) <- parms$rname
    names(b) <- parms$bname
    b <- list(b)
    names(b) <- parms$rname
    list(theta = rtheta, b = b)
  }
  
  out <- list(initialize = geinit,
              generate = gegenerate,
              wrapup = gewrapup)
  class(out) <- "coxmevar"
  out
}
