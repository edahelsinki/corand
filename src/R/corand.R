
## corand
##
## Copyright 2019 Kai Puolamaki <kai.puolamaki@iki.fi>
##
## Permission is hereby granted, free of charge, to any person obtaining
## a copy of this software and associated documentation files (the
## "Software"), to deal in the Software without restriction, including
## without limitation the rights to use, copy, modify, merge, publish,
## distribute, sublicense, and/or sell copies of the Software, and to
## permit persons to whom the Software is furnished to do so, subject to
## the following conditions:
##
## The above copyright notice and this permission notice shall be
## included in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
## EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
## MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
## NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
## LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
## OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
## WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


require(MASS)
require(glmnet)

tiling <- function(n,m,Tm=NULL,Tl=NULL,count=NULL,df=function(x) paste(x,collapse="_")) {
  ## Initial tiling
  ## Tiling matrix as in Alg. 1
  if(is.null(Tm)) Tm <- matrix(as.character(rep(1:m,each=n)),nrow=n,ncol=m)
  ## Additionally, we maintain a hash table of tiles (R,C).
  if(is.null(Tl)) {
    Tl <- new.env(hash=TRUE)
    for(j in 1:m) Tl[[as.character(j)]] <- list(R=1:n,C=j)
  }
  if(is.null(count)) count <- m

  ## Output unique id number
  newid <- function() {
    count <<- count+1
    as.character(count)
  }

  ## Tiling that freezes everything within tile.
  add0tile <- function(R=1:n,C=1:m) for(i in R) addtile(i,C)

  ## Add tile to a structure
  addtile <- function(R=1:n,C=1:m) {
    if(length(R)>0 && length(C)>0) {
      S <- new.env(hash=TRUE)
      ids <- NULL
      for(i in R) {
        raw <- Tm[i,C]
        K <- df(raw) # hash of permutation ids
        if(exists(K,envir=S)) {
          ##if(any(raw!=S[[K]]$raw)) stop("addtile: hash collision.") # This is probably over-cautious, but can however catch some programming errors (e.g., NAs)
          S[[K]]$rows <- c(S[[K]]$rows,i)
        } else {
          id <- unique(raw)
          ids <- union(ids,id)
          S[[K]] <- list(rows=i,id=id,raw=raw)
        }
      }
      for(K in ls(S)) {
        Cp <- NULL
        for(id in S[[K]]$id) Cp <- union(Cp,Tl[[id]]$C)
        id <- newid()
        Tm[S[[K]]$rows,Cp] <<- id
        Tl[[id]] <- list(R=S[[K]]$rows,C=Cp)
      }
      for(id in ids) { # remove overlapping rows from existing tiles
        Rnew <- setdiff(Tl[[id]]$R,R)
        if(length(Rnew)>0)
          Tl[[id]]$R <- Rnew
        else
          rm(list=id,envir=Tl) # if there are no rows left remove empty tile
      }
    }
  }

  copy <- function() {
    newTl <- new.env(hash=TRUE)
    for(x in ls(Tl,all.names=TRUE)) assign(x,get(x,Tl),newTl)
    tiling(n,m,Tm=Tm,Tl=newTl,count=count)
  }

  permute <- function() {
    p <- matrix(1:(n*m),nrow=n,ncol=m)
    for(key in ls(Tl)) {
      R <- Tl[[key]]$R
      C <- Tl[[key]]$C
      if(length(R)>1) p[R,C] <- p[sample(R),C]
    }
    c(p)
  }

  permutedata_matrix <- function(x,p) {
    matrix(x[p],nrow=n,ncol=m,dimnames=dimnames(x))
  }

  permutedata_df <- function(x,p) {
    pp <- p-n*rep(0:(m-1),each=n)
    if(any(pp<1 | n<pp)) stop("permutedata_df: permutations are not within column.")
    as.data.frame(lapply(1:m,function(j) x[pp[((j-1)*n+1):(j*n)],j]),
                  row.names=rownames(x),col.names=colnames(x))
  }

  permutedata_list <- function(x,p) {
    pp <- p-n*rep(0:(m-1),each=n)
    if(any(pp<1 | n<pp)) stop("permutedata_list: permutations are not within column.")
    a <- lapply(1:m,function(j) x[[j]][pp[((j-1)*n+1):(j*n)]])
    names(a) <- names(x)
    a
  }

  permutedata <- function(x,p=permute(),nmin=n) {
    if(is.matrix(x)) {
      if(nmin>n) {
        k <- 1+(nmin-1)%/%dim(x)[1]
        y <- matrix(NA,k*dim(x)[1],dim(x)[2],dimnames=dimnames(x))
        for(i in 1:k) y[((i-1)*dim(x)[1]+1):(i*dim(x)[1]),] <- permutedata_matrix(x,p=permute())
        y
      } else {
        permutedata_matrix(x,p)
      }
    } else if(is.data.frame(x)) {
      permutedata_df(x,p)
    } else {
      permutedata_list(x,p)
    }
  }


  cov_local <- function(x) {
    if(!is.matrix(x)) stop("needs matrix")
    x <- scale(x,center=TRUE,scale=FALSE)
    r <- matrix(NA,m,m,dimnames=list(colnames(x),colnames(x)))
    x0 <- matrix(NA,n,m)
    for(i in 1:m) x0[,i] <- ave(x[,i],Tm[,i])
    for(i in 1:(m-1)) {
      for(j in (i+1):m) {
        a <- x[,i]
        b <- x[,j]
        idx <- which(Tm[,i]!=Tm[,j])
        if(length(idx)>0) {
          a[idx] <- x0[idx,i]
          b[idx] <- x0[idx,j]
        }
        r[i,j] <- r[j,i] <- mean(a*b)
      }
    }
    diag(r) <- apply(x,2,function(y) mean(y^2))
    r
  }

  cor_local <- function(x) {
    cov2cor(cov_local(x))
  }

  list(add0tile=add0tile,
       addtile=addtile,
       copy=copy,
       permute=permute,
       permutedata=permutedata,
       cov=cov_local,
       cor=cor_local,
       status=function() list(n=n,m=m,Tm=Tm,Tl=Tl,count=count))
}

#' sorting step in Astrid
#'
#' @param s set (vector) of indices
#' @param f function that takes in a grouping (list of sets of indices) and outputs
#' a real number.
#'
#' The algorithm goes items in s one by one and always removes item from s and adds it to
#' such that the operation reduces f minimally. The final objective is to find as split
#' partition as possible such that f is small.
#'
#' @return $a sorted indices
#'
#' @export
astrid_sorting <- function(s,f) {
  a <- rep(NA,length(s))
  v <- c(f(list(s)),rep(NA,length(s)))
  B <- c()
  for(i in 1:length(s)) {
    if(length(s)>1) {
      aux <- vapply(s,function(x) f(c(list(setdiff(s,x)),as.list(c(B,x)))),0)
      j <- which.min(aux)
    } else {
      aux <- f(as.list(c(B,s)))
      j <- 1
    }
    a[i] <- s[j]
    v[i+1] <- aux[j]
    s <- setdiff(s,a[i])
    B <- union(B,a[i])
  }
  list(a=a,v=v)
}

#' grouping step in Astrid
#'
#' @param a sorted vector of indices (typically output by astrid sorting)
#' @param f function as in astrid_sorting
#'
#' The function goes through all splits in two and outputs the value of the function for
#' each. The intuition is that split with small cost we can make as it breaks nothing.
#'
#' @return vector of costs per split
#'
#' @export
astrid_grouping <- function(a,f) {
  vapply(1:(length(a)-1),function(i) f(list(a[1:i],a[(i+1):length(a)])),0)
}

#' the splitting step in Astrid
#'
#' @param a sorted vector of indices (typically output by astrid sorting)
#' @param v vector of costs per split (typically output by astrid_grouping)
#' @param f function as in astrid_sorting
#' @param x threshold
#'
#' The splitting uses binary search, resulting to O(log(length(a))) evaluations of f.
#'
#' @return maximal k grouping s such that f(s)<=x
#'
#' @export
astrid_split <- function(a,v,f,x=0.01) {
  p <- order(v)
  makes <- function(k) { # Make grouping of size k
    idx <- if(k>1) c(0,sort(p[1:(k-1)]),length(a)) else c(0,length(a))
    lapply(1:k,function(i) a[(idx[i]+1):idx[i+1]])
  }
  kmin <- 1 # value for which f(s)<=x
  if(f(makes(kmin))<=x) {
    kmax <- length(a) # value for which x<f(s)
    if(x<f(makes(kmax))) {
      while(kmin<kmax-1) {
        k <- floor((kmin+kmax)/2)
        if(f(makes(k))<=x) kmin <- k else kmax <- k
      }
      makes(kmin)
    } else {
      makes(kmax)
    }
  } else {
    makes(kmin)
  }
}

#' The pruning part of the Astrid
#'
#' @param g grouping (typically found by astrid_split)
#' @param f function as in astrid_sorting
#' @param x threshold
#'
#' This function throws away groups from grouping g such that f(g) stays in at most x.
#'
#' @return pruned grouping
#'
#' @export
astrid_prune <- function(g,f,x=0.01) {
  r <- c() # removed tiles
  if(f(g)<=x) {
    idx <- 1:length(g)
    while(length(idx)>0) {
      v <- vapply(idx,function(i) f(g[-union(i,r)]),0)
      j <- which.min(v)
      if(v[j]>x) break
      r <- union(r,idx[j]) # add idx[j] to removed tiles
      idx <- idx[-j] # remove jth element from idx
      idx <- idx[v[-j]<=x] # Don't try to throw away again tiles that are obviously important (they will not be less important later!), remove them from idx as well.
      if(any(v[-j]>x) && f(g[-union(idx,r)])<=x) { # If there are new obviously important tiles then try with only them
        r <- union(r,idx)
        break
      }
    }
  }
  if(length(r)>0) g[-r] else g
}

#' initial pruning
#'
#' @param s set (vector) of items
#' @param f function as in astrid_sorting
#' @param x tolerance for astrid_prune
#'
#' This initial pruning can be used to throw out clearly unnecessary attributes. This finds
#' a minimal grouping consisting of one tile such that f(g) stays in at most x.
#'
#' @return new set of items
#'
#' @export
astrid_prune0 <- function(s,f,x=0.01) {
  r <- c() # removed items
  if(f(list(s))<=x) {
    idx <- 1:length(s)
    while(length(idx)>0) {
      v <- vapply(idx,function(i) f(list(s[-union(i,r)])),0)
      j <- which.min(v)
      if(v[j]>x) break
      r <- union(r,idx[j]) # add idx[j] to removed attributes
      idx <- idx[-j] # remove jth element from idx
      idx <- idx[v[-j]<=x] # Don't try to throw again away attributes that are obviously important (they will not be less important later!), remove them from idx as well.
      if(any(v[-j]>x) && f(list(s[-union(idx,r)]))<=x) { # If there are new obviously important tiles then try with only them
        r <- union(r,idx)
        break
      }
    }
  }
  if(length(r)>0) s[-r] else s
}

#' The astrid algorithm
#'
#' @param s set (vector) of indices
#' @param f function that takes in a grouping (list of sets of indices) and outputs
#' a real number (as in astrid_sorting)
#' @param x threshold for the splitting step.
#' @param y threshold for the pruning step
#'
#' This runs the whole astrid algorithm
#'
#' @return Smallest possible maximum cardinality grouping g such that f(g)<=y.
astrid <- function(s,f,x=0.01,y=x,prune0=TRUE,prune=TRUE) {
  if(prune0) s <- astrid_prune0(s,f,x)
  if(length(s)>0) {
    a <- astrid_sorting(s,f)$a
    g <- astrid_split(a,astrid_grouping(a,f),f,x)
    if(prune) astrid_prune(g,f,y) else g
  } else {
    list()
  }
}



#' Whitening operator
#'
#' @param x1 n1Xm matrix
#' @param x2 n2Xm matrix
#'
#' Matrix x1 is transformed so that any difference in unit spherical distribution
#' signifies difference in distributions of x1 and x2. The operation is defined a a linear
#' operator such that whitening(x,x) gives a n1Xm matrix whose covariance matrix is a
#' unit matrix. Use ZCA-cor as recommended by https://arxiv.org/abs/1512.00809
#'
#' @return n1Xm matrix, transformation of x1.
#'
#' @export
whitening <- function(x1,x2) {
  a <- svd(cor(x2))
  w <- (((1/apply(x2,2,sd)) %o% (1/sqrt(a$d)))*a$u) %*% t(a$u)
  list(y=x1 %*% w,w=w)
}

findproj2 <- function(cov1,cov2,k=2) {
  a <- svd(cov2cor(cov2))
  W <- (((1/sqrt(diag(cov2))) %o% (1/sqrt(a$d)))*a$u) %*% t(a$u) # ZCA-cor
  a <- svd(t(W) %*% cov1 %*% W)
  p <- order(a$d,decreasing=TRUE)[1:k]
  list(w=gs(W %*% a$u[,p,drop=FALSE]),v=a$d[p])
}

findproj <- function(rdata1,rdata2,f=function(s) s,k=2) {
  a <- whitening(rdata1,rdata2)
  s <- svd(cov(a$y))
  u <- s$u
  d <- s$d
  v <- sapply(d,f)
  p <- order(v,decreasing=TRUE)[1:k]
  v <- v[p]
  u <- u[,p,drop=FALSE]
  d <- d[p]
  ## Projection directions are not necessarily orthogonal rotations.
  ## Force orthogonality using gs orthonormalisation.
  list(w=gs(a$w %*% u),v=v)
}




#' L1 norms between two empirical distributions
#'
#' @param x vector of real numbers
#' @param y vector of real numbers
#'
#' The function computes empirical L1 norm between distributions x and y
#'
#' @return Empirical L1 norm between x and y.
#'
#' @export
L1 <- function(x,y) {
  nx <- length(x)
  ny <- length(y)
  if(nx==0 || ny==0) stop("L1: zero length vector.")
  n <- nx+ny
  w <- c(x,y) # All observations
  wt <- c(rep(1/nx,nx),rep(-1/ny,ny)) # weight contributed by observations (sums to zero!)
  p <- order(w) # order observations
  w <- w[p]
  wt <- wt[p]
  sum(abs(cumsum(wt[-n])*(w[-1]-w[-n])))
}

#' Compute distance between two 2D distributions
#'
#' @param x nX2 matrix
#' @param y n'X2 matrix
#'
#' Computes difference of Pearson correlation coefficients
#'
#' @return Difference in correlations.
#'
#' @export
value_cor <- function(x,y) { abs(cor(x[,1],x[,2])-cor(y[,1],y[,2]))/2 }

#' Compute distance between two 2D distributions
#'
#' @param x nX2 matrix
#' @param y n'X2 matrix
#'
#' Computes difference of distributions using distribution estimation.
#'
#' @return Hellinger distance between distributions.
#'
#' @export
value_dens <- function(x,y,n=10) {
  h <- c(bandwidth.nrd((c(x[,1],y[,1]))),bandwidth.nrd((c(x[,2],y[,2]))))
  lims <- c(range(c(x[,1],y[,1])),range(c(x[,2],y[,2])))
  fx <- kde2d(x[,1],x[,2],h=h,n=n,lims=lims)
  fy <- kde2d(y[,1],y[,2],h=h,n=n,lims=lims)
  zx <- fx$z/sum(fx$z)
  zy <- fy$z/sum(fy$z)
  1-sum(sqrt(zx)*sqrt(zy)) # Use Hellinger divergence in [0,1] to compute the distance between distributions
}


#' Make LASSO regression
#'
#' @param x input matrix, nXm matrix
#' @param y response variable, n-vector
#' @param k value of non-zero variables, if null, cross validation is used to choose the best k
#'
#' @return w vector of weights
#'
#' @export
lasso <- function(x,y,k=NULL) {
  if(is.null(k)) {
    fit <- cv.glmnet(x,y)
    s <- fit$lambda.min
  } else {
    fit <- glmnet(x,y)
    s <- min(c(Inf,fit$lambda[fit$df<=k])) # Return lambda = Inf for zero length
  }
  w <- as.vector(coef(fit,s=s))
  names(w) <- c("(Intercept)",colnames(x))
  w
}

#' Computes L2 norm
#'
#' @param x Vector of real numbers
#' @return L2 norm
#'
#' @export
d2 <- function(x) sqrt(sum(x^2))

#' Normalizes vector
#'
#' @param x Vector of real numbers
#' @return Normalized vector
#'
#' @export
norm2 <- function(x) {
  s <- d2(x)
  if(s>0) x/s else x
}

#' Gram-Schmidt orthonormalization process
#'
#' @param x Matrix of size nXm of rank m
#' @return Matrix where columns are orthogonal and spanned by columns of x
#'
#' @export
gs <- function(x,eps=.Machine$double.eps^0.5) {
  m <- dim(x)[2]
  for(i in 1:m) {
    x[,i] <- norm2(x[,i])
    if(i<m) for(j in (i+1):m) x[,j] <- x[,j]-x[,i]*sum(x[,i]*x[,j])
  }
  x
}

#' Find an orthogonal component of x
#'
#' @param x vector of dimension d
#' @param W matrix of size dXn, columns consisting of orthogonal unit vectors
#' @return direction of x orthogonal to columns of x
#'
#' @export
gs1 <- function(x,W) {
  for(i in 1:dim(W)[2]) x <- x-W[,i]*sum(x*W[,i])
  norm2(x)
}

#' Find a random unit vector from the orthogonal subspace of x
#'
#' @param W orthogonal space
#' @return a random projection in orthogonal space
#'
#' @export
ro <- function(W) c(norm2(W %*% rnorm(dim(W)[2])))

#' Makes a random rotation matrix
#'
#' @param m dimensionality of the matrix
#' @param mm columns
#'
#' @return matrix of size
rrot <- function(m,mm=m) gs(matrix(rnorm(m*mm),m,mm))

rots <- function(n,m,U=NULL,fixed=NULL,unit=FALSE) {
  if(is.null(U)) {
    U <- array(NA,c(n,m,m))
    for(i in 1:n) U[i,,] <- if(unit) diag(m) else rrot(m)
  }
  if(is.null(fixed)) fixed <- rep(0,n)

  rot <- function(x) {
    for(i in 1:n) x[i,] <- x[i,] %*% U[i,,]
    x
  }

  irot <- function(x) {
    for(i in 1:n) x[i,] <- x[i,] %*% t(U[i,,])
    x
  }

  fixdir0 <- function(W,R=1:n,C=1:dim(W)[2]) {
    mm <- length(C)
    for(i in R) {
      V <- U[i,,]
      V[,C] <- W
      U[i,,] <<- gs(V)
    }
  }

  fixdir <- function(W,R=1:n) {
    mm <- dim(W)[2]
    if(max(fixed[R])+mm>m) stop("fixdir: too many fixes.")
    idx <- (fixed[R] %o% rep(1,mm))+(rep(1,length(R)) %o% (1:mm))
    for(i in 1:length(R)) {
      V <- U[R[i],,]
      V[,idx[i,]] <- W
      U[R[i],,] <<- gs(V)
    }
    fixed[R] <<- fixed[R]+mm
    idx
  }

  copy <- function() rots(n,m,U,fixed)


  list(
    rot=rot,
    irot=irot,
    fixdir0=fixdir0,
    fixdir=fixdir,
    U=function(R) U[R,,],
    copy=copy,
    status=function() list(U=U,fixed=fixed)
    )
}


#' adds a tile so that the subspace spanned by tile is frozen
#'
#' @param R rows affected
#' @param W mXk matrix whose columns give k directions that span the new subspace
#' @param t1
#' @param t2
#' @param rr
addfrozentile <- function(R,W,t1,t2=NULL,rr,copy=TRUE) {
  if(copy) {
    t1 <- t1$copy()
    if(!is.null(t2)) t2 <- t2$copy()
    rr <- rr$copy()
  }
  idx <- rr$fixdir(W,R)
  for(i in 1:length(R)) {
    t1$addtile(R[i],idx[i,])
    if(!is.null(t2)) t2$addtile(R[i],idx[i,])
  }

  list(t1=t1,t2=t2,rr=rr)
}

canonize <- function(g) {
  g <- lapply(g,sort)
  g[order(sapply(g,function(x) x[1]))]
}

#' Computes gradient of a function
#'
#' @param f function of m parameters
#' @param x0 vector of length m where gradient is computed
#'
#' @return Gradient of the function at x0.
#'
#' @export
grad <- function(f,x0,eps=.Machine$double.eps^(1/3)) {
  g <- vapply(1:length(x0),function(i) {
      x1 <- x2 <- x0
      x1[i] <- x1[i]-eps
      x2[i] <- x2[i]+eps
      (f(x2)-f(x1))/(2*eps)
    },0)
  names(g) <- names(x0)
  g
}

#' Computes gradient matrix
#'
#' @param f function that takes as an input nXm matrix where rows are m-dimensional data items and outputs a real number
#' @param x0 nXm matrix containing points where gradient is computed
#'
#' @return nX(m+1) matrix w where the first column is intercept and the remaining columns are
#' the gradient matrix. The following equality should apply f(x0)=w[,1]+rowSums(w[,-1]*x0).
#'
#' @export
wmatrix <- function(f,x0,eps=.Machine$double.eps^(1/3)) {
  w <- matrix(0,nrow=dim(x0)[1],ncol=dim(x0)[2]+1,
              dimnames=list(rownames(x0),c("(Intercept)",colnames(x0))))
  for(i in 1:dim(x0)[2]) {
    x1 <- x2 <- x0
    x1[,i] <- x1[,i]-eps
    x2[,i] <- x2[,i]+eps
    w[,i+1] <- (f(x2)-f(x1))/(2*eps)
  }
  w[,1] <- f(x0)-rowSums(w[,-1]*x0)
  w
}


###############################################################################################
# Plotting and printing functions follow
###############################################################################################

#' Makes a character string of a numeric vector
#'
#' @param x numeric vector
#' @return string representation of the vector
#'
#' @export
largestdim <- function(x,n=names(x),intercept=TRUE) {
  if(intercept) {
    x0 <- x[1]
    x <- x[-1]
  }
  x <- x[abs(x)>0]
  if(length(x)>0) {
    p <- order(abs(x),decreasing=TRUE)[1:min(5,length(x))]
    s <- paste(mapply(function(y,n) sprintf("%+.2f (%s)",y,n),x[p],n[p]),
               collapse=" ")
  } else {
    s <- ""
  }
  if(intercept) sprintf("%+.2f %s",x0,s) else s
}

#' Makes a plot of subset regression model
#'
#' @param w projection direction, (m+1) vectror
#' @param x input matrix, nXm matrix
#' @param y response variable, n vector
#' @param R subset of points in the selection
#' @param s name of the selection
#'
#' @export
projplot <- function(w,x,y,R,y0=NULL,yn="y",s="selection") {
  v <- (w[1]+x %*% w[-1])[,1]
  plot(v,y,type="n",bty="n",xlab=largestdim(w),ylab=yn,cex.lab=0.5)
  legend("topleft",legend=c(s,"others"),bty="n",pch=c(4,3),col=c("black","gray"))
  points(v[-R],y[-R],pch=3,col="gray")
  abline(a=0,b=1)
  if(!is.null(y0)) for(i in R) lines(c(v[i],v[i]),c(y[i],y0[i]))
  points(v[R],y[R],pch=4,col="black")
}

#' ppairs - make a pairplot
ppairs <- function(x,y,z=NULL,value=value_dens,
                   col=brewer.pal(3,"Set2"),pch=c(".",".","."),high=TRUE) {
  nx <- dim(x)[1]
  ny <- dim(y)[1]
  nz <- if(is.null(z)) 0 else dim(z)[1]
  n <- nx+ny+nz
  mat <- if(is.null(z)) rbind(x,y) else rbind(x,y,z)
  if(is.data.frame(mat)) {
    mat <- matrix(unlist(lapply(mat,as.numeric)),
                  nrow=n,ncol=dim(mat)[2],dimnames=dimnames(mat))
  }
  cols <- c(rep(col[1],nx),rep(col[2],ny),rep(col[3],nz))
  pchs <- c(rep(pch[1],nx),rep(pch[2],ny),rep(pch[3],nz))
  ix <- iy <- rep(FALSE,n)
  ix[1:nx] <- TRUE
  iy[(nx+1):(nx+ny)] <- TRUE
  ## randomize plotting order
  p <- sample.int(n)
  mat <- mat[p,]
  ix <- ix[p]
  iy <- iy[p]
  cols <- cols[p]
  pchs <- pchs[p]

  aux <- function(a,b) {
    value(matrix(c(a[ix],b[ix]),nx,2),
          matrix(c(a[iy],b[iy]),ny,2))
  }

  if(high) {
    highv <- -Inf
    for(i in 1:(dim(x)[2]-1)) {
      for(j in (i+1):dim(x)[2]) {
        highv <- max(highv,aux(mat[,i],mat[,j]))
      }
    }
  } else {
    highv <- Inf
  }

  panel.cor <- function(x,y,digits=3,prefix="",cex.cor,...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr=c(0,1,0,1))
    r <- aux(x,y)
    txt <- format(c(r,0.123456789),digits=digits)[1]
    txt <- paste0(prefix,txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    ##
    text(0.5,0.5,txt,cex=cex.cor,col=if(r>=highv) "red" else "black")
  }

  panel.points <- function(x,y,cex=3,...) {
    if(high) {
      if(aux(x,y)>=highv)
        lines(c(min(x),max(x),max(x),min(x),min(x)),
              c(min(y),min(y),max(y),max(y),min(y)),col="red")
    }
    points(x,y,cex=cex,...)
  }

  panel.hist <- function(x,col=NA,...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2],0,1.5))
    h <- hist(x,plot=FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB],0,breaks[-1],y,col="gray",...)
  }

  pairs(mat,col=cols,pch=pchs,bty="n",xaxt="n",yaxt="n",
        diag.panel=panel.hist,lower.panel=panel.cor,upper.panel=panel.points,
        oma=c(0,0,0,0))
}

xyplot <- function(x,class,y=NULL,
                   col=brewer.pal(length(levels(class)),"Set2"),
                   pch=1:length(levels(class)),
                   p=sample.int(dim(x)[1]),
                   n=10) {
  plot(rbind(x,y),bty="n",type="n")
  ## prepare 2d density estimates for x and y, respectively
  h <- c(bandwidth.nrd((c(x[,1],y[,1]))),bandwidth.nrd((c(x[,2],y[,2]))))
  lims <- c(range(c(x[,1],y[,1])),range(c(x[,2],y[,2])))
  fx <- kde2d(x[,1],x[,2],h=h,n=n,lims=lims)
  if(!is.null(y)) {
    fy <- kde2d(y[,1],y[,2],h=h,n=n,lims=lims)
    ##image(x=fx$x,y=fx$y,z=fx$z-fy$z,col=heat.colors(256))
    contour(x=fx$x,y=fx$y,z=fx$z-fy$z,col="gray",add=TRUE)
  } else {
    ##image(x=fx$x,y=fx$y,z=fx$z,col=heat.colors(256))
    contour(fx,col="gray",add=TRUE)
  }
  legend("topleft",legend=levels(class),box.col="gray",col=col,pch=pch) # bottom
  for(i in 1:2) for(j in p) rug(x[j,i],side=i,col=col[class[j]],lwd=1) # middle, randomized
  points(x[p,],col=col[class[p]],pch=pch[class[p]]) # top, randomized
}


lcat <- function(x,s="") {
  lcat0 <- function(y,s="") {
    s <- sprintf("%s{%d",s,y[1])
    if(length(y)>1) for(z in y[-1]) s <- sprintf("%s,%d",s,z)
    sprintf("%s}",s)
  }
  s <- sprintf("%s{",s)
  if(length(x)>0) {
    s <- sprintf("%s%s",s,lcat0(x[[1]]))
    if(length(x)>1) for(y in x[-1]) { s <- sprintf("%s,%s",s,lcat0(y)) }
  }
  sprintf("%s}",s)
}

#' Makes an ellipse of data
#'
#' @param data nXd data matrix
#' @param sigma the size of ellipse (default 0.95 confidence area)
#' @param ... graphics parameters
#' @returns As a side effect makes an 95% confidence ellipse based on data matrix.
#'
#' @export
plotellipse <- function(data,sigma=qnorm(0.975),...) {
  data <- as.matrix(data[,1:2,drop=FALSE])
  s <- svd(mycov(data))
  x <- seq(from=0,to=2*pi,length.out=100)
  lines((rep(1,length(x)) %o% colMeans(data))
        +sigma*sqrt(s$d[1])*(cos(x) %o% s$u[,1])
        +sigma*sqrt(s$d[2])*(sin(x) %o% s$u[,2]),...)
}


#' Makes a scatterplot of the data
#'
#' @param data nXd data matrix
#' @param rdata nXd matrix of random data (not used)
#' @param currentplot list that contains the current w matrix and name of the axis
#' @param grouped a subset of 1:n that contains the current selection
#' @returns As a side effect makes a pairplot to the directions specified by w.
#'
#' @export
plotdata <- function(data,rdata1=NULL,rdata2=NULL,w,grouped=NULL,s1="",s2="") {
  w <- w[,1:2,drop=FALSE] # Projection vectors
  xy <- data %*% w
  if(!is.null(rdata1)) {
    rxy1 <- rdata1 %*% w
  } else {
    rxy1 <- matrix(0,0,dim(data)[2])
  }
  if(!is.null(rdata2)) {
    rxy2 <- rdata2 %*% w
  } else {
    rxy2 <- matrix(0,0,dim(data)[2])
  }
  if(s1!="") s1 <- sprintf("%s ",s1)
  if(s2!="") s2 <- sprintf("%s ",s2)
  plot(rbind(xy,rxy1,rxy2),
       type="n",
       bty="n",
       cex.lab=0.5,
       xlab=sprintf("%s%s",s1,largestdim(w[,1],n=colnames(data),intercept=FALSE)),
       ylab=sprintf("%s%s",s2,largestdim(w[,2],n=colnames(data),intercept=FALSE)))
  idx <- 1:dim(xy)[1] %in% grouped

  for(i in 1:dim(xy)[1]) lines(rbind(rxy1[i,],xy[i,]),col=if(i %in% grouped) "pink" else "lightgray")
  points(rxy1,pch=1,col=ifelse(idx,"pink","lightpink"))
  points(rxy2,pch=1,col=ifelse(idx,"green","lightgreen"))
  points(xy,pch=20,col=ifelse(idx,"red","black"))
  if(length(grouped)>0) {
    plotellipse(xy[grouped,,drop=FALSE],col="darkred",lwd=3)
    plotellipse(rxy1[grouped,,drop=FALSE],col="deeppink",lty="dotted",lwd=2)
    plotellipse(rxy2[grouped,,drop=FALSE],col="darkgreen",lty="dotted",lwd=2)
  }
}

#' Makes a scatterplot of the data, alternative variant
#'
#' @param data nXd data matrix
#' @param rdata nXd matrix of random data (not used)
#' @param currentplot list that contains the current w matrix and name of the axis
#' @param grouped a subset of 1:n that contains the current selection
#' @returns As a side effect makes a pairplot to the directions specified by w.
#'
#' @export
plotdata2 <- function(data,rdata1=NULL,rdata2=NULL,w,grouped=NULL,s1="",s2="",ca = 0.8, labs = FALSE, 
                      showsamples = TRUE, col_bg="gray", col_sel="red", ...,col_data="black") {
  w <- w[,1:2,drop=FALSE] # Projection vectors
  xy <- data %*% w
  if(!is.null(rdata1)) {
    rxy1 <- rdata1 %*% w
  } else {
    rxy1 <- matrix(0,0,dim(data)[2])
  }
  if(!is.null(rdata2)) {
    rxy2 <- rdata2 %*% w
  } else {
    rxy2 <- matrix(0,0,dim(data)[2])
  }
  if(s1!="") s1 <- sprintf("%s ",s1)
  if(s2!="") s2 <- sprintf("%s ",s2)
  
  xlabel <- ylabel <- ""
  if (labs) {
    xlabel <- sprintf("%s%s",s1,largestdim(w[,1],n=colnames(data),intercept=FALSE))
    ylabel <- sprintf("%s%s",s2,largestdim(w[,2],n=colnames(data),intercept=FALSE))
  }
  
  plot(rbind(xy,rxy1,rxy2),
       type="n",
       bty="n",
       cex.lab = 0.6,
       cex.axis = ca,
       xlab = xlabel,
       ylab = ylabel
  )
  idx <- 1:dim(xy)[1] %in% grouped
  
  ## real data data
  points(xy, pch = 16, col = ifelse(idx, col_data, col_data), cex = 1.25)
  
  if (showsamples) {
    points(rxy2, pch=16, col = col_bg, cex = 1.25)
  }
  
  ## plot selection
  points(xy[idx, ], pch=16, col = col_sel, cex = 1.25)
  
}


#' Makes a scatterplot of the data, alternative variant 2
#'
#' @param data nXd data matrix
#' @param rdata nXd matrix of random data (not used)
#' @param currentplot list that contains the current w matrix and name of the axis
#' @param grouped a subset of 1:n that contains the current selection
#' @returns As a side effect makes a pairplot to the directions specified by w.
#'
#' @export
plotdatafocus <- function(data,rdata1=NULL,rdata2=NULL,w,grouped=NULL,s1="",s2="", 
                          ca = 0.8, labs = FALSE, cluster_ind = NULL, showsamples = TRUE,
                          col_sel="red", col_bg="gray") {
  w <- w[,1:2,drop=FALSE] # Projection vectors
  xy <- data %*% w
  if(!is.null(rdata1)) {
    rxy1 <- rdata1 %*% w
  } else {
    rxy1 <- matrix(0,0,dim(data)[2])
  }
  if(!is.null(rdata2)) {
    rxy2 <- rdata2 %*% w
  } else {
    rxy2 <- matrix(0,0,dim(data)[2])
  }
  if(s1!="") s1 <- sprintf("%s ",s1)
  if(s2!="") s2 <- sprintf("%s ",s2)
  
  xlabel <- ylabel <- ""
  if (labs) {
    xlabel <- sprintf("%s%s",s1,largestdim(w[,1],n=colnames(data),intercept=FALSE))
    ylabel <- sprintf("%s%s",s2,largestdim(w[,2],n=colnames(data),intercept=FALSE))
  }
  
  plot(rbind(xy,rxy1,rxy2),
       type="n",
       bty="n",
       cex.lab = 0.6,
       cex.axis = ca,
       xlab = xlabel,
       ylab = ylabel
  )
  
  idx <- 1:dim(xy)[1] %in% grouped
  
  col_green <- rgb(15 / 255, 120 / 255, 50 / 255)
  
  ## real data
  points(xy[idx,], pch=16, col = "black", cex = 1.25)
  points(xy[!(idx), ], pch=1, col = "black", cex = 1.25)
  
  ## marked cluster
  if (! is.null(cluster_ind)) {
    points(xy[cluster_ind, ], pch=16, col = col_sel, cex = 1.25)
  }
  
  if (showsamples) {
    ## H1
    points(rxy1[idx, ], pch = 0, col = col_green, bg = col_green, cex = 0.7)
    points(rxy1[!(idx), ], pch = 3, col = col_green, cex = 0.7)
    
    ## H2
    points(rxy2[idx, ], pch = 2, col = col_bg, bg = col_blue, cex = 0.7)
    points(rxy2[!(idx), ], pch = 3, col = col_bg, cex = 0.7)
  }
}


## Plot a figure using plotdata2 and save the result to a pdf-file
plotandsave <- function(data, rdata1, rdata2, w, grouped, ca = 2.5, fname = NULL, 
                        showplot = TRUE, labs = FALSE, pfunc = plotdata2, ...) {
  if (showplot) {
    pfunc(data    = data,
          rdata1  = rdata1,
          rdata2  = rdata2,
          w       = w,
          grouped = grouped,
          ca = ca,
          labs = labs, ...)
  }
  
  if (! is.null(fname)) {
    pdf(file = fname, width = 8, height = 6)
    
    pfunc(data    = data,
          rdata1  = rdata1,
          rdata2  = rdata2,
          w       = w,
          grouped = grouped,
          ca = ca,
          labs = labs, ...)
    
    dev.off()
  }
  
}

plotpcpandsave <- function(data, ind, fname, col_sel="red") {
  pcplot(data, grouped = ind, col_sel=col_sel)
  
  if (! is.null(fname)) {
    pdf(file = fname, width = 7, height = 5)
    pcplot(data, grouped = ind, col_sel=col_sel)
    dev.off()
  }
}

mycov <- function(x) {
  x <- x-(rep(1,dim(x)[1]) %o% colMeans(x))
  (t(x) %*% x)/dim(x)[1]
}

clustered <- function(data,grouped) {
  a0 <- apply(data,2,sd)
  a <- apply(data[grouped,],2,sd)/a0
  a[a0==0] <- Inf
  p <- order(a)
  list(p=p,a=a[p])
}

pcplot <- function(data,grouped=c(),selected=c(),fc=clustered, col_sel="red") {
  n <- dim(data)[1]
  m <- dim(data)[2]
  notgrouped <- setdiff(1:n,grouped)
  notselected <- setdiff(1:m,selected)

  ## scale data to [0,1]
  r <- apply(data,2,range)
  nz <- which(r[1,]!=r[2,])
  data[,nz] <- (data[,nz]-(rep(1,n) %o% r[1,nz]))/(rep(1,n) %o% (r[2,nz]-r[1,nz]))
  data[,-nz] <- 0.5
  if(length(grouped)>1) {
    aux <- fc(data,grouped)
    a <- aux$a
    p <- aux$p
    data <- data[,p]
  } else {
    p <- 1:m
  }
  sel <- rep(FALSE,m)
  sel[selected] <- TRUE
  sel <- which(sel[p])


  op <- par()
  par(pin=c(par("pin")[1],4.1))
  plot(c(0,1.3),c(1,m),type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",mar=c(0,0,0,0))
  for(i in 1:m) {
    lines(c(0,1),c(i,i),col=if(i %in% sel) "blue" else "gray")
    text(1,i,
         labels=if(length(grouped)>1) sprintf("%s (%.2f)",colnames(data)[i],a[i]) else colnames(data)[i],
         cex=0.5,col=if(i %in% sel) "blue" else "black",pos=4)
  }
  for(i in notgrouped) lines(data[i,],1:m,col="black")
  for(i in grouped) lines(data[i,],1:m,col=col_sel)
  par(op)
  p
}

#' A version of diag that always returns a matrix with x on diagonal
#'
#' @param x a vector of length n
#' @returns A nXn matrix with x on diagonal.
#'
#' @export
diagm <- function(x) { if(length(x)==0) matrix(0,0,0) else if(length(x)==1) matrix(x,1,1) else diag(x) }


#' Classic whitening: transform data so that the covariance matrix is a unit matrix
#'
#' @param X data matrix
#' @return Whitened data.
#'
#' @export
whitening0 <- function(X) {
  X <- as.matrix(X)
  mu <- colMeans(X)
  s <- svd(cov(X))
  t(s$u %*% diagm(1/sqrt(s$d)) %*% t(s$u) %*% t(X-((rep(1,dim(X)[1]) %o% mu))))
}


#' Performs fastICA on matrix Y
#' @param Y "whitened" data
#' @param g ICA function
#' @param gv Expected value of g(x) if x is a Gaussian variable
#' @return Returns a list of ICA vectors in decreasing order with
#'     respect to the weight and the weights.
#'
#' @export
dofastica <- function(Y,
                      g=function(u) log(cosh(u)),
                      gv=0.3745,
                      iter=1) {
  ##                gv=mean(g(rnorm(1000000)))) {
  if(iter==1) {
    a <- fastICA(Y,dim(Y)[2])
    w <- a$K %*% a$W
    ## normalize w so that the mean of ICA vectors is positive
    w <- w*(rep(1,dim(w)[1]) %o%
              (apply(w,2,function(x) {
                if(mean(x)>=0) {
                  1/sqrt(sum(x^2))
                } else {
                  -1/sqrt(sum(x^2))
                }})))
    ##v <- apply(g(a$S),2,mean)-gv
    v <- apply(g(whitening0(Y) %*% w),2,mean)-gv
    p <- order(abs(v),decreasing=TRUE)
    list(w=w[,p],v=v[p])
  } else {
    res <- lapply(1:iter,function(x) dofastica(Y,g=g,gv=gv,iter=1))
    res[[which.max(sapply(res,
                          function(x) sort(abs(x$v),
                                           decreasing=TRUE)[2]))]]
  }
}

## Function for getting the column set of a tile given the row indices in the data
get_column_indices <- function(data, ind, thr) {
  v_tmp <- clustered(data, ind)
  v_tmp <- names(v_tmp$a)[which(v_tmp$a < thr)]
  which( colnames(data) %in% v_tmp )
}

## --------------------------------------------------
## Functions used to calculate the gain for
## different projection vectors and hypothesis pairs
## --------------------------------------------------

## Functions for finding projection vectors
u_findproj <- function(cov_d_1, cov_d_2) {
  findproj2(cov_d_1, cov_d_2, k = 1)$w
}

u_pca <- function(data) {
  svd(data)$v[, 1, drop = FALSE]
}

u_ica <- function(data) {
  dofastica(data)$w[, 1]
}

get_projection_vectors <- function(data, tile_list_user, tile_focus) {
  t0  <- tiling(nrow(data), ncol(data))
  
  ## (1) explore - empty
  t_ee_h1 <- t0$copy()
  t_ee_h2 <- t0$copy()
  t_ee_h1$addtile( C = seq.int(ncol(data)), R = seq.int(nrow(data)) )
  
  cov_d_ee_h1 <- t_ee_h1$cov(data)
  cov_d_ee_h2 <- t_ee_h2$cov(data)
  
  ## (2) explore - user tiles
  t_et_h1 <- t_ee_h1$copy()
  t_et_h2 <- t_ee_h2$copy()
  
  for (tile in tile_list_user) {
    t_et_h1$addtile( R = tile$rows, C = tile$columns )
    t_et_h2$addtile( R = tile$rows, C = tile$columns )
  }
  
  cov_d_et_h1 <- t_et_h1$cov(data)
  cov_d_et_h2 <- t_et_h2$cov(data)
  
  ## (3) focus - empty
  t_fe_h1 <- t0$copy()
  t_fe_h2 <- t0$copy()
  
  t_fe_h1$addtile( R = tile_focus$rows, C = unlist(tile_focus$column_groups) )
  for (cg in tile_focus$column_groups)
    t_fe_h2$addtile( R = tile_focus$rows, C = cg)
  
  cov_d_fe_h1 <- t_fe_h1$cov(data)
  cov_d_fe_h2 <- t_fe_h2$cov(data)
  
  ## (4) focus - user tiles
  t_ft_h1 <- t_fe_h1$copy()
  t_ft_h2 <- t_fe_h2$copy()
  
  for (tile in tile_list_user) {
    t_ft_h1$addtile( R = tile$rows, C = tile$columns )
    t_ft_h2$addtile( R = tile$rows, C = tile$columns )
  }
  
  cov_d_ft_h1 <- t_ft_h1$cov(data)
  cov_d_ft_h2 <- t_ft_h2$cov(data)
  
  ## --------------------------------------------------
  ## Return projection vectors and the data
  ## --------------------------------------------------
  
  u_list <- list("ee"  = u_findproj(cov_d_ee_h1, cov_d_ee_h2),
                 "et"  = u_findproj(cov_d_et_h1, cov_d_et_h2),
                 "fe"  = u_findproj(cov_d_fe_h1, cov_d_fe_h2),
                 "ft"  = u_findproj(cov_d_ft_h1, cov_d_ft_h2),
                 "pca" = u_pca(data),
                 "ica" = u_ica(data))
  
  cov_d_list <- list("ee" = list("h1" = cov_d_ee_h1, "h2" = cov_d_ee_h2),
                     "et" = list("h1" = cov_d_et_h1, "h2" = cov_d_et_h2),
                     "fe" = list("h1" = cov_d_fe_h1, "h2" = cov_d_fe_h2),
                     "ft" = list("h1" = cov_d_ft_h1, "h2" = cov_d_ft_h2))
  
  list("u" = u_list, "cov_d" = cov_d_list)
}

## Functions for calculating the gain
gain <- function(covd1, covd2, u) {
  (sum(u * (covd1 %*% u))) / (sum(u * (covd2 %*% u)))
}

calculate_gains <- function(data, u_list, cov_d_list) {
  out <- matrix(NA, ncol = 4, nrow = length(u_list), dimnames = list(names(u_list), names(cov_d_list)))
  
  for (i in seq.int(length(u_list))) {
    out[i, 1] <- gain(cov_d_list$ee$h1, cov_d_list$ee$h2, u_list[[i]])
    out[i, 2] <- gain(cov_d_list$et$h1, cov_d_list$et$h2, u_list[[i]])
    out[i, 3] <- gain(cov_d_list$fe$h1, cov_d_list$fe$h2, u_list[[i]])
    out[i, 4] <- gain(cov_d_list$ft$h1, cov_d_list$ft$h2, u_list[[i]])
  }
  
  out
}
