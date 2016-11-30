#
# require(mousetrap)
# require(mta2)
# require(cstab)
# require(Rcpp)
# sourceCpp('~/Dropbox (2.0)/Work/Software/cstab/src/distMat.cpp')
# sourceCpp('~/Dropbox (2.0)/Work/Software/cstab/src/kmeansPredict.cpp')
# sourceCpp('~/Dropbox (2.0)/Work/Software/cstab/src/cStability2.R')
# sourceCpp('~/Dropbox (2.0)/Work/Software/cstab/src/instabCorrection.cpp')
# source('~/Dropbox (2.0)/Work/Software/cstab/R/instabCorrection.R')
#
#
# d = readRDS('~/Dropbox (2.0)/Work/Software/cstab/data/example_Spivey2005_E1_C1_trajectories.RDS')
# data = cbind(d[,2,],d[,3,])
#
# system.time({res_4 = cStability_f(data,2:5,type = 'hierarchical',linkage = 'ward.D')})
# system.time({res_4 = cStability_f(rbind(data,data),2:5,type = 'hierarchical',linkage = 'ward.D')})
# system.time({res_4 = cStability(rbind(data,data),2:5,type = 'hierarchical',linkage = 'ward.D')})
#
# system.time({res_4 = cStability_f(data,2:5,type = 'kmeans',linkage = 'ward.D')})
# system.time({res_4 = cStability_f(rbind(data,data),2:5,type = 'kmeans',linkage = 'ward.D')})
# system.time({res_4 = cStability(rbind(data,data),2:5,type = 'kmeans',linkage = 'ward.D')})
#
#
# require(MASS)
# d = rbind(mvrnorm(100,c(1,1,1),diag(3)),
#           mvrnorm(100,c(5,5,5),diag(3)),
#           mvrnorm(100,c(9,9,9),diag(3)),
#           mvrnorm(100,c(13,13,13),diag(3)),
#           mvrnorm(100,c(17,17,17),diag(3)))
# plot(d)
#
# cDistance(d,2:10,method='hierarchical')
# cStability_f(d,2:10)
#
#
# cDistance(data,2:100,method='hierarchical')
# a = cStability_f(data,2:100,method = 'hierarchical',norm = TRUE, nB = 100)
#
#
# GetMeasures(data,2)
#
# WCD <- Sil <- MSE <- numeric()
# for(k in 1:10) {
#   obj = GetMeasures(data, k, method = 'hierarchical')
#   WCD[k] = obj$WCD
#   Sil[k] = obj$Sil
#   MSE[k] = obj$MSE
#   }
#
#

#require(devtools)
#document('~/Dropbox (2.0)/Work/Software/cstab')
##
#check('~/Dropbox (2.0)/Work/Software/cstab')
#
# require(MASS)
# cluster_example = rbind(mvrnorm(100,c(0,0),diag(2)),mvrnorm(100,c(6,0),diag(2)),
#       mvrnorm(100,c(-2,6),diag(2)),mvrnorm(100,c(8,6),diag(2)),
#       mvrnorm(100,c(3,9.873),diag(2)))
# colnames(cluster_example) = c('x','y')
# save(cluster_example,file='~/Dropbox (2.0)/Work/Software/cstab/data/cluster_example.RData')

#' @useDynLib cstab
#' @importFrom Rcpp sourceCpp
NULL
