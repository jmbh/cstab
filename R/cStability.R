

cStability <- function(data, # n x p data matrix
                       kseq, # sequence of ks tested
                       B = 10, # number of bootstrap comparisons
                       norm = TRUE, # norm over pw equal assign,FALSE=as in Wang etal
                       prediction = TRUE, # use prediction approach, if FALSE, use brute pair in equal cluster approach
                       type = 'kmeans', # or 'hierarchical'
                       linkage = 'complete', # or average, or ...
                       kmIter = 10,
                       pbar = TRUE) # number of reruns of k-means algorithm
{


  # ---------- Auxilliary Functions ----------

  # Function to calculate the geometric mean
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }

  # ---------- Input Checks ----------

  # On Data
  if(sum(is.na(data))>0) stop('No missing values permitted!')
  # On k-sequence
  if(1 %in% kseq) stop('Please select a k sequence starting with 2: {2,3,...K}!')
  # On B
  if(round(B)!=B) stop('The number of bootstrap comparison has to be a positive integer value.')
  # On type
  if(type %in% c('kmeans', 'hierarchical'))


  # ---------- Calculate Auxilliary Variables ----------

  p <- nrow(data)
  m_instab <- matrix(NA, B, length(kseq)) # storage: no normalization
  m_instab_norm <- matrix(NA, B, length(kseq)) # storage: normalization

  # Compute necesarry bootstrap samples to get B comparisons
  Bsamp <- ceiling(.5 * (sqrt(8*B+1) + 1)) # solution to equation B = [Bsamp(Bsamp-1)]/2


  # ---------- Draw bootstrap samples ----------

  l_ind <- list()
  for(b in 1:Bsamp) {
    l_ind[[b]] <-  sample(1:p, p, replace=T)
    l_ind[[b]] <- l_ind[[b]][order(l_ind[[b]])] # order
  }

  combs <- combn(1:Bsamp ,2) # All possible combinations


  # ---------- compute indices for each pair: which objects are in both samples ----------

  # only necessary for intersection approach
  if(!prediction) {
  l_indices <- list()
  for(bc in 1:B) {
    # overlap
    IS <- intersect(l_ind[[combs[1,bc]]],l_ind[[combs[2,bc]]])
    l_pair <- list()
    count <- 1
    for(r in combs[,bc]) {

      ind_is_r1 <- l_ind[[r]] %in% IS # indicator: indice in both?
      ind2_r1 <- !duplicated(l_ind[[r]])
      ind3_r1 <- ind_is_r1 & ind2_r1

      l_pair[[count]] <- ind3_r1
      count<-count+1
    }
    l_indices[[bc]] <- l_pair
  }
  }


  # ----- Loop over B comparisons -----

  if(pbar)  pb <- txtProgressBar(min=0, max=B, style = 2)

  for(bc in 1:B) {

    for(k in kseq) {

      l_clust <- list()

      if(prediction) {

        # kmeans or spectral clustering?
        if(type=='kmeans') {
          l_km_models <- list()
          count <- 1
          for(r in combs[,bc]) {

            l_km <- list()
            for(km in 1:kmIter) {
              l_km[[km]] <- kcca(data[l_ind[[r]],], k=k, kccaFamily("kmeans")) #save whole model
            }

            WCD <- unlist(lapply(l_km, function(x) mean(x@clusinfo$av_dist)))
            km_model <- l_km[[which.min(WCD)]] # pick k-means clustering with smallest WCD
            l_clust[[count]] <- predict(km_model, newdata=data)  #make predictions
            count <- count+1
          }
        }
        if(type=='hierarchical') {

          count <- 1
          for(r in combs[,bc]) {
            hcobj <- hclust(dist(data[l_ind[[r]],]), method = linkage)
            l_clust[[count]] <- cutree(hcobj, k)
            count <- count+1
          }

        }

        # count pairwise equal assignments
        same_a <- as.numeric(dist(l_clust[[1]])==0)*1
        same_b <- as.numeric(dist(l_clust[[2]])==0)*1

        # compute Instability
        InStab <- mean(abs(same_a - same_b))

        # Normalize = FALSE
        m_instab[bc, which(kseq==k)] <- InStab

        # Normalize = TRUE
        tb1 <- table(l_clust[[1]])
        tb2 <- table(l_clust[[2]])
        norm_val <- instab(tb1, tb2, 100)
        m_instab_norm[bc,which(kseq==k)] <- InStab / norm_val

        # else: no prediction
      } else {

        if(type=='kmeans') {
          l_cl <- list()
          l_pairind <- list() # are two objects in same cluster (1 yes, 0 no)
          count <- 1
          for(r in combs[,bc]) {

            # run k means several times
            l_km <- list()
            for(km in 1:kmIter) {
              l_km[[km]] <- kcca(data[l_ind[[r]],], k=k, kccaFamily("kmeans")) #save whole model
            }
            WCD <- unlist(lapply(l_km, function(x) mean(x@clusinfo$av_dist)))
            km_model <- l_km[[which.min(WCD)]] # pick k-means clustering with smallest WCD

            cl_long <- km_model@cluster

            l_cl[[count]] <- cl_long[l_indices[[bc]][[count]]] # only take the ones in the intersection set
            l_pairind[[count]] <- (as.numeric(dist(l_cl[[count]]))==0)*1
            count <- count+1
          }
        }
        if(type=='hierarchical') {
          l_cl <- list()
          l_pairind <- list() # are two objects in same cluster (1 yes, 0 no)
          count <- 1
          for(r in combs[,bc]) {
            hcobj <- hclust(dist(as.matrix(data[l_ind[[r]],])), method = linkage)
            cl_long <- cutree(hcobj, k)
            l_cl[[count]] <- cl_long[l_indices[[bc]][[count]]] # only take the ones in the intersection set
            l_pairind[[count]] <- (as.numeric(dist(l_cl[[count]]))==0)*1
            count <- count+1
          }
        }

        # compute Instability
        InStab <- mean(abs(l_pairind[[1]] - l_pairind[[2]]))


        # Normalize = FALSE
        m_instab[bc, which(kseq==k)] <- InStab

        # Normalize = TRUE
        tb1 <- table(l_cl[[1]])
        tb2 <- table(l_cl[[2]])
        norm_val <- instab(tb1, tb2, 100)
        m_instab_norm[bc,which(kseq==k)] <- InStab / norm_val

      } # end if: prediction TRUE/FALSE

    } # end for k

    if(pbar) setTxtProgressBar(pb, bc)

  } # end for B

  # taking the mean
  m_instab_M <- apply(m_instab, 2, mean)
  m_instab_norm_M <- apply(m_instab_norm, 2, mean)

  kopt_instab <- which.min(m_instab_M)+(min(kseq)-1)
  kopt_instabN <- which.min(m_instab_norm_M)+(min(kseq)-1)

  # replicate function call
  f_call <- list('kseq'=kseq,
                 'B'=B,
                 'norm'=norm,
                 'prediction'=prediction,
                 'type'=type,
                 'linkage'=linkage,
                 'kmIter'=kmIter)

  outlist <- list("kopt_instab"=kopt_instab,
                  "kopt_instab_norm"=kopt_instabN,
                  "Instab_path"=m_instab_M,
                  "Instab_path_norm"=m_instab_norm_M,
                  "Instab_path_matrix"=m_instab,
                  "Instab_path_nrom_matrix"=m_instab_norm,
                  'call'=f_call)

  class(outlist) <- c('list', 'cstab', 'cStability')

  return(outlist)

} # EoF

