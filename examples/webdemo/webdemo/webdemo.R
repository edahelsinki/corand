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

## The actual logic is in this file. Now this contains a PCA and plotting function,
## used by server.R.

plotcolors <- c("#1b9e77","#d95f02","#7570b3")


dopca <- function(data) {
  s <- svd(cov(data))
  w <- s$u
  ## normalize w so that the mean of PCA vectors is positive
  w <- w*(rep(1,dim(w)[1]) %o%
            (apply(w,2,
                   function(x) if(mean(x)>=0) 1/sqrt(sum(x^2)) else -1/sqrt(sum(x^2)))))
  p <- order(s$d,decreasing=TRUE)
  list(w=w[,p],v=s$d[p])
}

plotdata <- function(data,current) {
  w <- current$w[,1:2] # Projection vectors
  v <- current$v[1:2]  # Projection costs
  xy <- data$r %*% w
  idx <- 1:dim(xy)[1] %in% current$R
  plot(xy,main=sprintf("rows (%d of %d)",length(current$R),dim(data$data)[1]),
       col=ifelse(idx,plotcolors[1],"black"),
       pch=ifelse(idx,16,1),
       bty="n",
       xlab=sprintf("%s[%.2g] = %s",current$label,
                    v[1],largestdim(w[,1],names=colnames(data$r))),
       ylab=sprintf("%s[%.2g] = %s",current$label,
                    v[2],largestdim(w[,2],names=colnames(data$r))))
}

plotattr <- function(data,current,factors=TRUE) {
  w <- current$w[,1:2] # Projection vectors
  v <- current$v[1:2]  # Projection costs
  xy <- (if(factors) data$raux else data$raux[1:dim(data$raux)[2],,drop=FALSE]) %*% w
  idx <- 1:dim(xy)[1] %in% current$C
  plot(xy,main=sprintf("columns (%d of %d)",length(current$C),dim(data$data)[2]),
       type="n",
       bty="n",
       xlab=sprintf("%s1[%.2g] = %s",current$label,
                    v[1],largestdim(w[,1],names=colnames(data$r))),
       ylab=sprintf("%s2[%.2g] = %s",current$label,
                    v[2],largestdim(w[,2],names=colnames(data$r))))
  text(xy,
       labels=rownames(xy),
       col=ifelse(idx,plotcolors[2],"black"),
       cex=ifelse(idx,1.4,1))
}

pcorder <- function(data,current,mode,n) {
    if(mode=="X1") {
        rev(order(abs(current$w[,1]),decreasing=TRUE)[1:n])
    } else if(mode=="selection") {
        R <- if(length(current$R)>2) current$R else 1:dim(data$r)[2]
        m <- apply(data$r[R,],2,mean)
        m0 <- apply(data$r,2,mean)
        v <- apply(data$r[R,],2,var)
        v0 <- apply(data$r,2,var)
        ##kl <- log(v0/v)/2+(v+(m-m0)^2)/(2*v0)-0.5
        ##kl <- (v+(m-m0)^2)/(2*v0)-0.5
        kl <- v/v0
        rev(order(kl)[1:n])
    } else {
        stop("pcorder: unknown mode")
    }
}

pcplot <- function(data,current,p=1:dim(data$data)[2]) {
    a <- data.frame(data$data[,p],
                    X2=data$r %*% current$w[,2],
                    X1=data$r %*% current$w[,1])
    a <- as.matrix(as.data.frame(lapply(a,as.numeric)))
    amin <- apply(a,2,min)
    adelta <- apply(a,2,max)-amin
    adelta <- ifelse(adelta>0,adelta,1)
    a <- (a-(rep(1,dim(a)[1]) %o% amin))/(rep(1,dim(a)[1]) %o% adelta)
    plot(c(0,1.3),c(0.5,dim(a)[2]+0.5),
         type="n",main="",xlab="",ylab="",bty="n",xaxt="n",yaxt="n")
    text(x=1.2,y=1:dim(a)[2],colnames(a),pos=4,
         col=c(ifelse(c(p %in% current$C),plotcolors[2],"black"),"black","black"))
    for(i in 1:dim(a)[2]) {
        lines(c(0,1.2),c(i,i),
              col=if(i<=length(p)) { if(p[i] %in% current$C) plotcolors[2] else "gray" } else "black",
              lwd=if(i<=length(p) && p[i] %in% current$C) 2 else 1)
    }
    for(i in 1:dim(a)[1]) {
        lines(a[i,],1:dim(a)[2],
              col=if(i %in% current$R) plotcolors[1] else "gray",
              lty=if(i %in% current$R) "solid" else "dotted")
    }
}

#' Makes a character string of a numeric vector
#'
#' @param x numeric vector
#' @return string representation of the vector
#'
#' @export
largestdim <- function(x,names=as.character(1:length(x)),s="") {
  p <- order(abs(x),decreasing=TRUE)[1:min(5,length(x))]
  paste(mapply(function(y,n) sprintf("%+.2f (%s)",y,n),x[p],names[p]),
        collapse=" ")
}
