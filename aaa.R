
###############################
#   Chargement des packages   #
###############################

remove(list= ls())

library(cubt)
library(mclust)
library(doParallel)
library(cluster)
library(dbscan)

detectCores()
if(detectCores() < 32){
  registerDoParallel(3)
} else {
  registerDoParallel(10)
}
ncore_actif <- getDoParWorkers()


data1 <- gendata(mod = 3, N = 4000)
data <- data1[, -1]
yobs <- data1[, 1]

devpart(data, yobs, crit0 = "anova")

a <- kmeans(data, centers = 2, nstart = 20)
devpart(data, a$cluster, crit0 = "anova")

table(yobs)
table(a$cluster)
a$centers

plot(data, col = as.factor(yobs))
plot(data, col = as.factor(a$cluster))

devpart(data[yobs == 1,], yobs[yobs == 1], crit0 = "anova")
devpart(data[yobs == 2,], yobs[yobs == 2], crit0 = "anova")

apply(data, 2, mean)
colMeans(data[yobs ==1,])
colMeans(data[yobs ==2,])




################################
#   Chargement des fonctions   #
################################

pred.NFC <- function(data, centres){
  cc <- pdist(data, centres)
  ccc <- as.matrix(cc)
  cl <- numeric(nrow(ccc))
  for (i in 1:nrow(ccc)){
    cl[i] <- which.min(ccc[i,])
  }
  return(cl)
}

pred.NIC <- function(df, centres, clus.centre, seuil.ncore = 2600, ncore.for.seuil = 10){
  ndata = nrow(df)
  ncentre = nrow(centres)
  if (!is.null(seuil.ncore)){
    if(ncentre < seuil.ncore ){
      ncores <- 0
    } else { ncores <- ncore.for.seuil}
  }
  n = ndata + ncentre
  if (ncores > 1) {
    nearest = NULL
    ff = function(a, b) {which.min(colSums((t(a) - b)^2))}
    nearest = foreach(j = 1:ndata, .combine = "c") %dopar% ff(centres, df[j,])
    gc()
  } else {
    AB = as.matrix(pdist(centres, df))
    nearest = apply(AB, 2, which.min)
    
  }
  return(clus.centre[nearest])
}

cv = function(N,nvc=10) {
  taille = N%/%nvc
  reste = N - nvc * taille
  bloc = sample(c(rep(1:nvc, taille), sample(nvc, reste)))
}

Kmeans.fusion <- function(df, n.bloc, nbcl, recup.method = "centers", alloc.method = "NFC", nstart.bloc = 20, nstart.centre = 20, s.ncore.NIC, ncore.NIC){
  
  n <- nrow(df)
  bloc <- cv(n, n.bloc)
  
  k.all.centres <- foreach(i = 1:n.bloc, .combine = "rbind") %dopar% {
    kmeansD <- kmeans(df[bloc == i,], nstart = nstart.bloc, centers = nbcl)
    if(recup.method == "centers"){return(kmeansD$centers)}
    else if(recup.method == "midpoints"){return(aggregate(df, list(kmeansD$cluster), median)[, -1])}
  }
  
  k.res <- kmeans(k.all.centres, nstart = nstart.centre, centers = nbcl)
  
  if(alloc.method == "NFC"){
    clas <- pred.NFC(df, k.res$centers)
  } else if(alloc.method == "NIC"){
    clas <- pred.NIC(df, k.all.centres, k.res$cluster, seuil.ncore = s.ncore.NIC, ncore.for.seuil = ncore.NIC)
  }
  return(clas)
}

space.tabulate <- "\t"
sptab <- function(x){rep(space.tabulate, x)}
now <- function(){return(substr(Sys.time(), 12,21))}

##################
#   Paramètres   #
##################

data1 <- gendata(1, N = 400, sigma = 0.19)
data <- data1[, -1]
yobs <- data1[, 1]
sil_1 <- silhouette(yobs, dist(data))
sil_2 <- silhouette(sample(1:4, 400, replace = T), dist(data))
plot(sil_1)
plot(sil_2)
mean(sil_1[, 3])
mean(sil_2[, 3])

data1 <- gendata(2, N = 400, sigma = 0.19)
data <- data1[, -1]
yobs <- data1[, 1]
sil_1 <- silhouette(yobs, dist(data))
sil_2 <- silhouette(sample(1:10, 400, replace = T), dist(data))
plot(sil_1)
plot(sil_2)
mean(sil_1[, 3])
mean(sil_2[, 3])

data1 <- gendata(3, N = 400, sigma = 0.19)
data <- data1[, -1]
yobs <- data1[, 1]
k.obj <- kmeans(data, centers = 2)
sil_1 <- silhouette(yobs, dist(data))
sil_2 <- silhouette(k.obj$cluster, dist(data))
plot(sil_1)
plot(sil_2)
mean(sil_1[, 3])
mean(sil_2[, 3])

data1 <- gendata(4, N = 300, sigma = 0.05)
data <- data1[, -1]
yobs <- data1[, 1]
sil_1 <- silhouette(yobs, dist(data))
sil_2 <- silhouette(sample(1:3, 300, replace = T), dist(data))
plot(sil_1)
plot(sil_2)
mean(sil_1[, 3])
mean(sil_2[, 3])

#####

dist_fun <- function(data, yobs){
  vec <- 0
  for(i in 1:(length(unique(yobs)) - 1)){
    for(j in (i+1):length(unique(yobs))){
      vec <- vec + distAB(data, yobs == unique(yobs)[i], yobs == unique(yobs)[j])
    }
  }
  return(vec)
}


data1 <- gendata(1, N = 400, sigma = 0.19)
data <- data1[, -1]
yobs <- data1[, 1]
dist_fun(data, yobs)
dist_fun(data, sample(1:4, 400, replace = T))

data1 <- gendata(2, N = 400, sigma = 0.19)
data <- data1[, -1]
yobs <- data1[, 1]
dist_fun(data, yobs)
dist_fun(data, sample(1:10, 400, replace = T))

data1 <- gendata(3, N = 400, sigma = 0.19)
data <- data1[, -1]
yobs <- data1[, 1]
k.obj <- kmeans(data, centers = 2)
dist_fun(data, yobs)
dist_fun(data, k.obj$cluster)

data1 <- gendata(4, N = 300, sigma = 0.05)
data <- data1[, -1]
yobs <- data1[, 1]
dist_fun(data, yobs)
dist_fun(data, sample(1:3, 300, replace = T))

































