#  Modified from File src/library/splines/R/splines.R
#  Part of the R package, https://www.R-project.org

bspline <- function(x, df = NULL, knots = NULL, Boundary.knots = range(x)){
  degree <- 3
  intercept <- TRUE
  ord <- 1L + (degree <- as.integer(degree))
  if(ord <= 1) stop("'degree' must be integer >= 1")
  nx <- names(x)
  x <- as.vector(x)
  nax <- is.na(x)
  if(nas <- any(nax))
    x <- x[!nax]
  outside <- if(!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    (ol <- x < Boundary.knots[1L]) | (or <- x > Boundary.knots[2L])
  } else FALSE
  
  if(!is.null(df) && is.null(knots)) {
    nIknots <- df+2L - ord + (1L - intercept) # ==  #{inner knots}
    if(nIknots < 0L) {
      nIknots <- 0L
      warning(gettextf("'df' was too small; have used %d",
                       ord - (1L - intercept) - 2L), domain = NA)
    }
    knots <-
      if(nIknots > 0L) {
        knots <- seq.int(from = 0, to = 1,
                         length.out = nIknots + 2L)[-c(1L, nIknots + 2L)]
        quantile(x[!outside], knots)
      }
  }
  
  Aknots <- sort(c(rep(Boundary.knots, ord), knots))
  if(any(outside)) {
    warning("some 'x' values beyond boundary knots may cause ill-conditioned bases")
    derivs <- 0:degree
    scalef <- gamma(1L:ord)# factorials
    basis <- array(0, c(length(x), length(Aknots) - degree - 1L))
    e <- 1/4 # in theory anything in (0,1); was (implicitly) 0 in R <= 3.2.2
    if(any(ol)) {
      ## left pivot inside, i.e., a bit to the right of the boundary knot
      k.pivot <- (1-e)*Boundary.knots[1L] + e*Aknots[ord+1]
      xl <- cbind(1, outer(x[ol] - k.pivot, 1L:degree, "^"))
      tt <- splineDesign(Aknots, rep(k.pivot, ord), ord, derivs)
      basis[ol, ] <- xl %*% (tt/scalef)
    }
    if(any(or)) {
      ## right pivot inside, i.e., a bit to the left of the boundary knot:
      k.pivot <- (1-e)*Boundary.knots[2L] + e*Aknots[length(Aknots)-ord]
      xr <- cbind(1, outer(x[or] - k.pivot, 1L:degree, "^"))
      tt <- splineDesign(Aknots, rep(k.pivot, ord), ord, derivs)
      basis[or, ] <- xr %*% (tt/scalef)
    }
    if(any(inside <- !outside))
      basis[inside,  ] <- splineDesign(Aknots, x[inside], ord)
  }
  else basis <- splineDesign(Aknots, x, ord)
  if(!intercept)
    basis <- basis[, -1L , drop = FALSE]
  n.col <- ncol(basis)
  if(nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax,  ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  a <- list(degree = degree, knots = if(is.null(knots)) numeric(0L) else knots,
            Boundary.knots = Boundary.knots, intercept = intercept)
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("bs", "basis", "matrix")
  df0 = df +2L
  x.bs <- basis
  attr(x.bs, "intercept") <- NULL
  x.attr = attributes(x.bs);
  x.attr$dim[2] = x.attr$dim[2]-2
  x.attr$dimnames[[2]] = (x.attr$dimnames[[2]])[-c(df0-1,df0)];
  x.bs2 = x.bs[,-c(1,df0)];
  attributes(x.bs2)=x.attr;
  a <- list(df = df)
  attributes(x.bs2) <- c(attributes(x.bs2), a)
  class(x.bs2) <- c("bspline", "basis", "matrix")
  x.bs2
}

predict.bspline <- function(object, newx, ...){
  if(missing(newx))
    return(object)
  a <- c(list(x = newx), attributes(object)[
    c("df", "knots", "Boundary.knots")])
  do.call("bspline", a)
}

makepredictcall.bspline <- function(var, call){
  if(as.character(call)[1L] != "bspline") return(call)
  at <- attributes(var)[c("df", "knots", "Boundary.knots")]
  xxx <- call[1L:2]
  xxx[names(at)] <- at
  xxx
}
