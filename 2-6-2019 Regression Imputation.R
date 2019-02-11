
library(mice)
library(foreign)

rm(list = ls())

dataset <- read_sav(file="Backpain50 MI missing.sav")

names(dataset)

dataset <- dataset[, c("Pain", "Tampascale")]

#dataset <- dataset[, c(2, 3)]
#dataset

x <- as.matrix(dataset[, 1])
dimnames(x) <- list(NULL, "Pain")
x


#y	= Incomplete data vector of length n  # vector of variable with missing data
y <- dataset$Tampascale
y

#ry	= Vector of missing data pattern (FALSE=missing, TRUE=observed) # missing data indicator
ry <- !is.na(y)
ry

#x <- na.omit(x)

x

#function (y, ry, x, ridge = 1e-05, ...) 
#{
x <- cbind(1, as.matrix(x))
x
  
xobs <- x[ry, ]
xobs
yobs <- y[ry]
yobs
  
  ridge = 1e-05
  xtx <- t(xobs) %*% xobs
  pen <- ridge * diag(xtx)
  if (length(pen) == 1) 
    pen <- matrix(pen)
  v <- solve(xtx + diag(pen))
  coef <- t(yobs %*% xobs %*% v)
  
  return(x[!ry, ] %*% coef)
#}

imp.regress <- mice(dataset, method="norm.predict", m=1, maxit=3)
imp.regress

imp.regress
mean(imp.regress$imp$Tampascale[, 1])


imp.regress$chainMean

