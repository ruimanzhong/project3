library("splines")
library(INLA)
knots = seq(0,1, length.out  = J)

# estimate mlik for kth cluster using INLA 
### initialize log likelihood vector
max(membership)
evalLogLike_each_INLA(7, Y, tPHI, lambda, sigmasq_y, population, cluster)
log_like_ini = evalLogLike_all_parallel2(Y, tPHI, lambda, sigmasq_y, population, cluster, cl)
log_like = log_like_ini$llik_all;log_like_vec=log_like_ini$log_like_vec

for(iter in 1:MCMC) {
  
  if(k == 1) {rb = 0.95; rd = 0; rc = 0
  } else if(k == ns) {rb = 0; rd = 0.85; rc = 0.1
  } else {rb = 0.425; rd = 0.425; rc = 0.1}
  
  move = sample(4, 1, prob = c(rb, rd, rc, rhy))
  
  if(move == 1) { ## Birth move
    # split an existing cluster
    split_res = splitCluster(mstgraph, k, cluster)
    membership_new = split_res$cluster 
    split_res$k=k+1
    
    # compute log-proposal ratio
    if(k == ns-1) {
      rd_new = 0.85
    } else {rd_new = 0.425}
    # if(k_m == n-1) {
    #   rd_new = 0.6
    # } else {rd_new = 0.3}
    log_P = log(rd_new) - log(rb)
    # compute log-likelihood ratio
    
    if(k>1){
      
      lambda_new = rbind(lambda, lambda[split_res$clust_old,])
      
    }else{
      
      lambda_new = rbind(lambda, lambda)
    }
    
    ind_k=which(membership_new==k+1); ns1=length(ind_k)
    Y_clust=Y[ind_k,] 
    vec_lambda_k=Para_re2all(lambda_new[k+1,], J)
    # sigmasq_yk=postSigmasq(Y_clust, X, vec_lambda_k, a0, b0)
    # sigmasq_y_new = c(sigmasq_y, sigmasq_yk)
    
    

    
    # compute log-prior ratio
    log_A = log(1-c)  +
      sum(c0/2*log(d0/2)-log(gamma(c0/2))-(c0/2+1)*log(lambda_new[k+1,])-d0/(2*lambda_new[k+1,])) 
    
    ## calculate loglikelihood ratio by only comparing local likelihood of the two clusters that changed.
    
    log_L_new=evalLogLike.ratio('split',log_like_vec, split_res, Y, X, lambda_new, sigmasq_y_new, membership_new)
    log_L=log_L_new$ratio 
    
    #acceptance probability
    acc_prob = min(0, log_A + log_P + log_L)
    acc_prob = exp(acc_prob)
    if(runif(1) < acc_prob){
      # accept
      cluster = membership_new
      k = k + 1
      
      log_like_vec=log_L_new$log_like_vec
      log_like = sum(log_like_vec);
      
      lambda = lambda_new; 
      sigmasq_y = sigmasq_y_new; 
      ind_gam=ind_gam_new
      edge_status = getEdgeStatus(cluster, mstgraph)
      birth_cnt=birth_cnt+1
      
    }
    
  }
  
  if(move == 2) { ## Death move
    # merge two existing clusters (c1, c2) -> c2
    merge_res = mergeCluster(mstgraph, edge_status, cluster)
    membership_new = merge_res$cluster
    
    ##index of cluster that is removed
    cid_rm = merge_res$cluster_rm
    
    
    # # compute log-proposal ratio
    if(k == 2) {rb_new = 0.85
    }else {rb_new = 0.425}
    
    log_P = log(rb_new) - log(rd) 
    
    sigmasq_y_new=sigmasq_y[-cid_rm]; lambda_new = lambda[-cid_rm,];
    
    
    ##For hetero variances
    ##index of the new merged cluster
    cid_comb = merge_res$cluster_comb
    cid_newid = merge_res$cluster_newid
    
    ind_k=which(membership_new == cid_newid); ns1=length(ind_k)
    Y_clust=Y[ind_k,]
    vec_lambda_k=Para_re2all(lambda[cid_comb,], J)
    
    sigmasq_yk=0.01
    sigmasq_y_new[cid_newid]=sigmasq_yk; 
    
    # Obtain the log prior for ind_gam
    
    # compute log-prior ratio
    
    log_A = -log(1-c)  - sum(c0/2*log(d0/2)-log(gamma(c0/2))-(c0/2+1)*log(lambda[cid_rm,])-d0/(2*lambda[cid_rm,]))
    
    # compute log-likelihood ratio
    
    log_L_new=evalLogLike.ratio('merge',log_like_vec, merge_res, Y, X, lambda_new, sigmasq_y_new, membership_new)
    log_L=log_L_new$ratio 
    
    
    
    #acceptance probability
    acc_prob = min(0, log_A + log_P + log_L)
    acc_prob = exp(acc_prob)
    if(runif(1) < acc_prob){
      # accept
      cluster = membership_new
      k = k - 1
      
      log_like_vec=log_L_new$log_like_vec
      log_like = sum(log_like_vec);
      
      lambda = lambda_new; sigmasq_y = sigmasq_y_new; 
      ind_gam=ind_gam_new;
      
      edge_status = getEdgeStatus(cluster, mstgraph)
      
      death_cnt=death_cnt+1
      
    }
    
  }
  
  
  if(move == 3) { ## change move
    # first perform death move: (c1, c2) -> c2
    merge_res = mergeCluster(mstgraph, edge_status, cluster)
    cid_rm = merge_res$cluster_rm
    membership_new = merge_res$cluster
    
    
    sigmasq_y_new=sigmasq_y[-cid_rm]; lambda_new = lambda[-cid_rm,];
    
    ##For hetero variances
    
    cid_comb = merge_res$cluster_comb
    cid_newid = merge_res$cluster_newid
    
    ind_k=which(membership_new == cid_newid); ns1=length(ind_k)
    
    Y_clust=Y[ind_k,]
    vec_lambda_k=Para_re2all(lambda[cid_comb,], J)
    sigmasq_yk=0.01
    sigmasq_y_new[cid_newid]=sigmasq_yk; 
    
    # Obtain the log prior for ind_gam

    
    log_prior_merge= - sum(c0/2*log(d0/2)-log(gamma(c0/2))-(c0/2+1)*log(lambda[cid_rm,])-d0/(2*lambda[cid_rm,]))
    
    
    k = k-1
    
    log_L_new_merge=evalLogLike.ratio('merge',log_like_vec, merge_res, Y, X, lambda_new, sigmasq_y_new, membership_new)
    
    # then perform birth move
    split_res = splitCluster(mstgraph, k, merge_res$cluster);
    split_res$k=k+1;membership_new=split_res$cluster
    
    k = k+1
    
    
    if(k>2){
      lambda_new = rbind(lambda_new, lambda_new[split_res$clust_old,])
      ind_gam_new = rbind(ind_gam_new, ind_gam_new[split_res$clust_old,])
      
    }else{
      
      lambda_new = rbind(lambda_new, lambda_new)
      ind_gam_new = rbind(ind_gam_new, ind_gam_new)
      
      
    }
    # For hetero variances
    ind_k=which(membership_new==k); ns1=length(ind_k)
    Y_clust=Y[ind_k,]
    vec_lambda_k=Para_re2all(lambda_new[k,], J)
    sigmasq_yk=0.01
    sigmasq_y_new = c(sigmasq_y_new, sigmasq_yk)
    
    
    # Obtain the log prior for ind_gam
    
    log_prior_split=sum(c0/2*log(d0/2)-log(gamma(c0/2))-(c0/2+1)*log(lambda_new[k,])-d0/(2*lambda_new[k,]))
    
    
    # compute log-likelihood ratio
    
    log_L_new=evalLogLike.ratio('split',log_L_new_merge$log_like_vec, split_res, Y, X, lambda_new, sigmasq_y_new, membership_new)
    log_L = log_L_new$ratio + log_L_new_merge$ratio + log_prior_merge + log_prior_split
    
    
    
    
    # acceptance probability
    acc_prob = min(0, log_L)
    acc_prob = exp(acc_prob)
    if(runif(1) < acc_prob){
      # accept
      cluster = membership_new
      
      log_like_vec=log_L_new$log_like_vec
      log_like = sum(log_like_vec);
      
      lambda = lambda_new; sigmasq_y = sigmasq_y_new; 
      ind_gam=ind_gam_new
      
      edge_status = getEdgeStatus(cluster, mstgraph)
      
      change_cnt=change_cnt+1
      
    }
    
  }
  
  if(move == 4) { ## Hyper move
    hyper_cnt = hyper_cnt + 1
    
    #################################################################################      
    
    # update MST
    
    edge_status_G = getEdgeStatus(cluster, graph0) 
    mstgraph = proposeMST(graph0, edge_status_G)
    
    V(mstgraph)$vid=1:ns
    
    edge_status = getEdgeStatus(cluster, mstgraph)
    
    #################################################################################      
    
    # update lambda
    
    #####parallel implementation
    
    # lambda.res=postLambda_par(beta, sigmasq_y, ind_gam, c0, d0, cluster, cl)
    lambda.res = matrix(1, J+1, 1)
    lambda=t(matrix(lambda.res, (J+1), k))
    
    #######################################################################################      
    ##update ind_gam
    # ind_gam_new=ind_gam
    
    # for(i in 1:k){
    #   
    #   temp_ind=which(cluster==i); ns1=length(temp_ind)
    #   
    #   if(k>1){
    #     
    #     lambda_k=lambda[i,]; ind_gam_k=ind_gam[i,];
    #     
    #   }else{
    #     
    #     lambda_k=lambda; ind_gam_k=ind_gam;
    #     
    #   }
    #   
    #   Y_clust=matrix(t(Y[temp_ind,]), ns1*nt, 1)
    #   vec_lam_k=Para_re2all(lambda_k, J);
    #   
    #   res_gam=postGamma_loop(Y_clust, X, vec_lam_k, ind_gam_k, a0, b0)
    #   
    #   if(k>1){
    #     ind_gam_new[i,]=unlist(res_gam)
    #   }else{
    #     
    #     ind_gam_new=unlist(res_gam)
    #     
    #   }
    #   
    # }
    # 
    # ind_gam=ind_gam_new
    
    ##Update log-likelihood value
    log_like_res = evalLogLike_all_parallel2(Y, tPHI, lambda, sigmasq_y, cluster, cl)
    log_like=log_like_res$llik_all;log_like_vec=log_like_res$log_like_vec;
  }
  
  ############################################################################## 
  # update sigma2
  
  ####Consider a parallel implementation      
  sigmasq_y_new=postSigmasq_par(Y, X, lambda, ind_gam, a0, b0, cluster, cl)
  sigmasq_y = sigmasq_y_new
  
  # update beta
  beta_res=postBeta_par(Y, X, lambda, sigmasq_y, ind_gam, cluster, cl)
  beta=t(matrix(beta_res, nt, k))
  
  ##Update log-likelihood value
  log_like_res = evalLogLike_all_parallel2(Y, X, lambda, sigmasq_y, ind_gam, cluster, cl)
  log_like=log_like_res$llik_all;log_like_vec=log_like_res$log_like_vec;
  
  ###store estimates
  if(iter %% 100 == 0) {
    cat('Iteration', iter, 'k', k, 'birth cnt', birth_cnt, 'death cnt', death_cnt, 'change cnt', change_cnt, 'hyper cnt', hyper_cnt, '\n')
    cat('move', move, 'k', k, 'ncluster', length(unique(cluster)), '\n') 
  }
  
  ## save result
  if(iter > BURNIN & (iter - BURNIN) %% THIN == 0) {
    
    beta_out[[(iter-BURNIN)/THIN]] = beta
    lambda_out[[(iter-BURNIN)/THIN]] = lambda
    sigmasq_y_out[[(iter-BURNIN)/THIN]] = sigmasq_y
    
    MST_out[[(iter-BURNIN)/THIN]] = mstgraph
    cluster_out[(iter-BURNIN)/THIN, ] = cluster
    loglike_out[(iter-BURNIN)/THIN] = log_like
  }
  
  if(iter %% 100 == 0) {
    cat('Iteration', iter, 'birth cnt', birth_cnt, 'death cnt', death_cnt, 'change cnt', change_cnt, 'hyper cnt', hyper_cnt, '\n')
    
    save(beta_out, lambda_out, sigmasq_y_out, MST_out, cluster_out, loglike_out, hyperpar, birth_cnt, death_cnt, change_cnt, file=path_save)
    
  }
  
}

