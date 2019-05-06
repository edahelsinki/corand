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


## executed before application launch

set.seed(42)

library("shiny")
library("DT")

source("../../../src/R/corand.R")
source("webdemo.R")

diag0 <- function(x) {
    if(!is.vector(x)) stop("diag0: x must be a vector.")
    if(length(x)==1) matrix(x,1,1) else diag(x)
}

getdata <- function(datafile) {
    origdata <- as.data.frame(readRDS(sprintf("data/%s",datafile)))
    ## Shuffle rows to avoid any spurious structures in plots due to ordering of rows
    origdata <- origdata[sample.int(dim(origdata)[1]),]
    ## r is the real valued part of the data
    ii <- which(!sapply(origdata,is.factor))
    r <- as.matrix(origdata[,ii,drop=FALSE])
    rme <- apply(r,2,mean)
    rsd <- apply(r,2,sd)
    raux <- (rep(1,dim(r)[2]) %o% rme)+diag0(rsd)
    colnames(raux) <- rownames(raux) <- colnames(r)
    ## subsets are the factors
    fn <- colnames(origdata)[-ii]
    if(length(fn)>0) {
        f <- lapply(fn,function(i) origdata[,i])
        names(f) <- fn
        a <- lapply(fn,function(s)
            model.matrix(as.formula(sprintf("~0+%s",s)),origdata))
        s <-  unlist(lapply(a,function(x) lapply(as.data.frame(x),
                                                 function(y) which(y==1))),
                     recursive=FALSE)
        m <- lapply(s,function(x) apply(r[x,,drop=FALSE],2,mean))
        data <- cbind(r,do.call(cbind,a))
        raux <- rbind(raux,do.call(rbind,m))
    } else {
        f <- list()
        s <- list()
        data <- r
    }
    descr <- sprintf("%d rows, %d real attributes, %d factor%s",
                     dim(r)[1],dim(r)[2],length(f),if(length(f)==1) "" else "s")
    cat(sprintf("getdata: file %s, %s.\n",datafile,descr))
    list(r=r,f=f,data=data,s=s,descr=descr,raux=raux)
}

initcurrent <- function(data) {
    a <- dopca(data$r)
    ## tiles <- lapply(data$s,function(R) list(R=R,
    ##                                         C=1:dim(data$data)[2],
    ##                                         H1=FALSE,
    ##                                         H2=FALSE))
    if(length(data$s)>0) {
        tiles <- lapply(1:length(data$s),function(i)
            list(R=data$s[[i]],C=dim(data$r)[2]+i,H1=FALSE,H2=FALSE))
        names(tiles) <- names(data$s)
    } else {
        tiles <- list()
    }
    tiles[["ALL"]] <- list(
      R=1:dim(data$data)[1],
      C=1:dim(data$data)[2],
      H1=TRUE,
      H2=FALSE
    )
    H1 <- tiling(dim(data$data)[1],dim(data$data)[2])
    for(i in which(sapply(tiles,function(x) x$H1))) {
      H1$addtile(R=tiles[[i]]$R,C=tiles[[i]]$C)
    }
    H2 <- tiling(dim(data$data)[1],dim(data$data)[2])
    for(i in which(sapply(tiles,function(x) x$H2))) {
      H2$addtile(R=tiles[[i]]$R,C=tiles[[i]]$C)
    }
    list(R=c(),C=c(),w=a$w,v=a$v,tiles=tiles,
         pc=order(abs(a$w[,1]),decreasing=TRUE)[1:min(10,length(a$w[,1]))],
         pcmode="X1",
         factors=TRUE,
         H1=H1,
         H2=H2,
         label="PC",count=1)
}

makecolumnchoises <- function() {
  a <- as.list(1:dim(data$raux)[1])
  names(a) <- rownames(data$raux)
  a
}


## load data
datafiles <- list.files("data/",pattern=".*\\.rds$")
if(length(datafiles)<1) stop("global.R: no rds file found.")
data <- getdata(datafiles[1]) # load the first file

## initialise status
current <- initcurrent(data)
