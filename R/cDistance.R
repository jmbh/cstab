

cDistance <- function(data, # n x p data matrix
                      kseq, #sequence of ks to be checked
                      method = 'kmeans', # or: hiearchical
                      linkage = 'complete',
                      kmIter = 10, # restarts of k means algorithm
                      RunsGap = 10) # number of simulated datasets in gap statistic
{

  # ----- Input Checks -----

  # On Data
  if(sum(is.na(data))>0) stop('No missing values permitted!')
  # On k-sequence
  if(1 %in% kseq) stop('Please select a k sequence starting with 2: {2,3,...K}!')
  # On type
  if(method %in% c('kmeans', 'hierarchical'))


  # ----- Calculate Aux Variables -----

  n <- nrow(data)
  dims <- ncol(data)


  # ----- run clustering clustering on real data -----

  l_WCD_data <- l_Sil <- l_MSE <-  list()
  for(k in c(1, kseq)) {
    obj <- GetMeasures(data, k, method=method, linkage=linkage, kmIter = kmIter, type = 'data')
    l_WCD_data[[k]] <- obj$WCD
    l_Sil[[k]] <- obj$Sil
    l_MSE[[k]] <- obj$MSE
  }
  l_WCD_data <- unlist(l_WCD_data)
  l_Sil <- unlist(l_Sil)
  l_MSE <- unlist(l_MSE)

  # ----- run clustering clustering on synthetic data (for Gap Statistic) -----

  l_WCD_syndata_runs <- list()
  for(runs in 1:RunsGap) {
    l_WCD_syndata1 <- list()
    data_synt <- UniformData(data)
    for(k in c(1,kseq)) {
      obj <- GetMeasures(data_synt, k, method=method, linkage=linkage, kmIter = kmIter, type = 'simulated')
      l_WCD_syndata1[[k]] <- obj$WCD
    }
    l_WCD_syndata_runs[[runs]] <- unlist(l_WCD_syndata1)
  }
  l_WCD_syntetic <- colMeans(do.call(rbind, l_WCD_syndata_runs))


  # ----- Compute Measures -----

  # Gap Statistic
  WCD_data_log <- log(l_WCD_data)
  WCD_syn_log <- log(l_WCD_syntetic)

  WCD_data_log <- WCD_data_log - WCD_data_log[1]
  WCD_syn_log <- WCD_syn_log - WCD_syn_log[1]

  Gaps <- WCD_syn_log-WCD_data_log

  # Slope Statistic
  sil <- l_Sil
  p <- 1
  slope <- -(sil[-1] - sil[-length(sil)]) * sil[-1]^p
  kopt_slope <- which.max(slope)+1 #add one because we have sequence 2:...

  ## Jump Statistic
  MSE_transf <- l_MSE^(- dims/2)
  jump <- (MSE_transf - c(0, MSE_transf[-length(MSE_transf)]))[-1]
  kopt_Jump <- which.max(jump) + 1 #add one because we have sequence 2:...

  outlist <- list('kOpt_Gap'=which.max(Gaps),
                  'kOpt_Slope'=kopt_slope,
                  'kOpt_Jump'=kopt_Jump,
                  'WCD_data'=l_WCD_data,
                  'WCD_syn'=l_WCD_syntetic,
                  'Gaps'= Gaps,
                  'silhouettes'=sil,
                  'jumps'=jump)

  return(outlist)

} # EoF


