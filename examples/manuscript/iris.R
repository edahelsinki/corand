
## This script is used to produce the projections in Sec. 2.5
## Run by Rscript --vanilla iris.R

source("../../src/R/corand.R")

data <- as.matrix(scale(iris[,c("Petal.Length","Petal.Width",
                                "Sepal.Length","Sepal.Width")]))


R1 <- which(iris$Species=="setosa")
R2 <- which(iris$Species=="versicolor")
R3 <- which(iris$Species=="virginica")

print(cor(data[c(R1,R2),1:2]))


## Empty tiling
t0 <- tiling(150,4)

#### Case 1

## H1
H11 <- t0$copy()
H11$addtile()

## H2
H21 <- t0$copy()

a1 <- findproj2(H11$cov(data),H21$cov(data),k=1)
print(a1$w[,1])


#### Case 2

H12 <- H11$copy()
H12$addtile(R=c(R1,R2),C=1:2)

H22 <- H21$copy()
H22$addtile(R=c(R1,R2),C=1:2)

a2 <- findproj2(H12$cov(data),H22$cov(data),k=1)
print(a2$w[,1])


#### Case3


H13 <- t0$copy()
H13$addtile(R=c(R2,R3),C=1:3)
H13$addtile(R=c(R1,R2),C=1:2)

H23 <- t0$copy()
H23$addtile(R=c(R2,R3),C=1)
H23$addtile(R=c(R2,R3),C=2:3)
H23$addtile(R=c(R1,R2),C=1:2)

a3 <- findproj2(H13$cov(data),H23$cov(data),k=1)
print(a3$w[,1])
