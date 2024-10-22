##The main function for the Fclust-RST method
##Y: the n by m data matrix; 
##X: the m by m inverse discrete wavelet transform matrix;
##graph0: the hyper-parameter graph, G;
##init_val: lists of initial values;
##hyperpar: lists of hyperparameters;
##MCMC: the number of MCMC iterations;
##BURNIN: the burn-in period;
##THIN: the thinning setting;
##path_save: the file path and file name for storing the results;
##seed: the random seed specification;
##n.cores: the number of CPU cores to use.

####Functions for implementing the Fclust-RST method

###Function to obtain m by 1 vector of lambdas from the unique lambda of each resolution
###x: the (J+1) by one vector of parameters;
###J: the number of resolutions
Para_re2all <- function(x, J){
  
  vec_x=array(0, 2^J)
  vec_x[1]=x[1]
  start=2
  for(i in 1:J){
    
    ind_sel=start:(start+2^{J-i}-1)
    vec_x[ind_sel]=x[i+1]
    start=start+2^{J-i}
    
  }
  return(vec_x)
  
}

#########  Return both likelihood and local likelihood vector

evalLogLike_all_parallel2 <- function(Y, X, lambda, sigmasq_y, population, membership, cl,a0 = NULL, b0 = NULL) {
  
  p=max(membership)
  
  llik_res=parLapply(cl, 1:p, evalLogLike_each_INLA, Y, X, lambda, sigmasq_y, population, membership, a0, b0)
  log_like_vec=unlist(llik_res)
  llik_all=sum(log_like_vec)
  
  
  return(list(llik_all=llik_all,log_like_vec=log_like_vec))
  
}

#######  Functions for updating local likelihood vector

evalLogLike.ratio=function(move,log_like_vec, control.move, Y, X, lambda, sigmasq_y, membership,a0 = NULL, b0 = NULL){
  
  
  ### update local likelihoods for split move
  if(move=='split'){ 
    log_like_vec_new=log_like_vec;
    M1 <- evalLogLike_each_INLA (control.move$clust_old, Y, X, lambda, sigmasq_y, population, membership, a0, b0)  ## new clust 1 
    M2 <- evalLogLike_each_INLA (control.move$k,  Y, X, lambda, sigmasq_y, population, membership, a0, b0)  ## new clust 2
    log_like_vec_new[control.move$clust_old]=M1
    log_like_vec_new[control.move$k]=M2
    llratio=M1+M2-log_like_vec[control.move$clust_old]
  }
  
  ### update local likelihoods for merge move
  if(move=='merge'){  
    M1 <- evalLogLike_each_INLA (control.move$cluster_newid, Y, X, lambda, sigmasq_y, population, membership, a0, b0)
    log_like_vec_new=log_like_vec[-control.move$cluster_rm]
    log_like_vec_new[control.move$cluster_newid]=M1
    llratio= M1-sum(log_like_vec[c(control.move$cluster_rm,control.move$cluster_comb)])
  }
  
  return(list(ratio=llratio,log_like_vec=log_like_vec_new))
} 

################################################################################

######Compute the log-likelihood value for a cluster
######k: the cluster ID;
######Y: observation vector;
######X: inverse discrete wavelet transform matrix;
######lambda: the (J+1) by one vector of scale parameters;
######sigmasq_y: the measurement-error variance;
######ind_gam: the m by one vector of indicator variables;
######membership: the n b one vector of cluster memberships

evalLogLike_each_INLA <- function(k, Y, X, lambda, sigmasq_y, population, membership, a0 = NULL, b0 = NULL) {
  
  nt=base::nrow(X); J = ncol(X) ; 
  p=max(membership)
  
  ind=which(membership==k); ns=length(ind); n=ns*nt
  
  Yk= matrix(Y[ind,],ncol = nt)
  vec_Yk = as.integer(as.vector(t(Yk))) 
  
  sigma2_yk=sigmasq_y[k]
  
  if(p>1){
    
    lambda_k=lambda[k,]
    # ind_gam_k=ind_gam[k,]
  }else{
    
    lambda_k=lambda
    # ind_gam_k=ind_gam
    
  }
  COV = cbind(1, do.call(rbind, replicate(nrow(Yk), X, simplify = FALSE)))
  beta.prec = diag(lambda_k)*(1/sigma2_yk)
  data = as.data.frame(cbind(vec_Yk = vec_Yk,COV = COV))
  
  idx = rep(1:nt,nrow(Yk))
    formula= vec_Yk~ 0+ COV + f(idx, model = 'iid',
                              hyper = list(prec = list(prior = "loggamma", param = c(a0,b0))))
    m.bs3 <- INLA::inla(formula, family = "poisson", E = rep(population[ind],each = nt),
                        data = data, 
                        control.fixed = list(prec = beta.prec)
    )
    return(m.bs3[["mlik"]][[1]] + 0.5 * log(det(diag(lambda_k))))
  
}

# ################################################################################


##################################################################################
###Function to get the posterior samples of the wavelet coefficients vector, beta

#### For k-th cluster
postBeta_each<-function(k, Y, X, lambda, sigmasq_y,population,membership, a0 = NULL, b0 = NULL){
  
  nt=base::nrow(X); J = ncol(X) ; 
  p=max(membership)
  
  ind=which(membership==k); ns=length(ind); n=ns*nt
  
  Yk= matrix(Y[ind,],ncol = nt)
  vec_Yk = as.integer(as.vector(t(Yk))) 
  
  sigma2_yk=sigmasq_y[k]
  
  if(p>1){
    
    lambda_k=lambda[k,]
    # ind_gam_k=ind_gam[k,]
  }else{
    
    lambda_k=lambda
    # ind_gam_k=ind_gam
    
  }
  COV = cbind(1, do.call(rbind, replicate(nrow(Yk), X, simplify = FALSE)))
  beta.prec = diag(lambda_k)*(1/sigma2_yk)
  data = as.data.frame(cbind(vec_Yk = vec_Yk,COV = COV))
  idx = rep(1:nt,nrow(Yk))
  if (!is.null(a0) & !is.null(b0)){
    formula= vec_Yk~ 0+ COV + f(idx, model = 'iid',
                                hyper = list(prec = list(prior = "loggamma", param = c(a0,b0))))
    m.bs3 <- INLA::inla(formula, family = "poisson", E = rep(population[ind],each = nt),
                        data = data, control.compute=list(config = TRUE),
                        control.fixed = list(prec = beta.prec)
    )
  } else {
    formula= vec_Yk~ 0+ COV + f(idx, model = 'iid',
                                hyper = list(prec = list(prior = "loggamma", param = c(a0,b0))))
    m.bs3 <- INLA::inla(formula, family = "poisson", E = rep(population[ind],each = nt),
                        data = data, control.compute=list(config = TRUE),
                        control.fixed = list(prec = beta.prec)
    )
  }
  sample = INLA::inla.posterior.sample.eval(c( "COV1", "COV2","COV3","COV4" ), INLA::inla.posterior.sample(1, m.bs3))
  beta_new = matrix(sample, J+1, 1)
  
  return(as.numeric(beta_new))
  
}

#####For all clusters
postBeta_par<-function(k, Y, X, lambda, sigmasq_y,population,membership, a0 = NULL, b0 = NULL){
  
  
  p=max(membership) 
  beta_all=parLapply(cl, 1:p, postBeta_each, Y, X, lambda, sigmasq_y,population,membership, a0, b0)
  
  return(unlist(beta_all))
  
  
  
}

################################################################################
###Function to propose a new sigma2


###For k-th cluster

postSigmasq_each_INLA<-function(k, Y, X, lambda, sigmasq_y,population, a0, b0, membership){
  
  nt=base::nrow(X); J = ncol(X) ; 
  p=max(membership)
  
  ind=which(membership==k); ns=length(ind); n=ns*nt
  
  Yk= matrix(Y[ind,],ncol = nt)
  vec_Yk = as.integer(as.vector(t(Yk))) 
  
  sigma2_yk=sigmasq_y[k]
  
  if(p>1){
    
    lambda_k=lambda[k,]
    # ind_gam_k=ind_gam[k,]
  }else{
    
    lambda_k=lambda
    # ind_gam_k=ind_gam
    
  }
  COV = cbind(1, do.call(rbind, replicate(nrow(Yk), X, simplify = FALSE)))
  beta.prec = diag(lambda_k)*(1/sigma2_yk)
  data = as.data.frame(cbind(vec_Yk = vec_Yk,COV = COV))
  
  idx = rep(1,length(vec_Yk))
  formula= vec_Yk~ 0+ COV + f(idx, model = 'iid',
                            hyper = list(prec = list(prior = "loggamma", param = c(a0,b0))))
  m.bs3 <- INLA::inla(formula, family = "poisson", E = rep(population[ind],each = nt),
                      data = data, 
                      control.fixed = list(prec = beta.prec)
  )
  sigmasq_new <- INLA::inla.hyperpar.sample(1,m.bs3)
  return(as.numeric(sigmasq_new))
  
}

######For all clusters

postSigmasq_par<-function(k, Y, X, lambda, sigmasq_y,population, a0, b0, membership, cl){
  
  p=max(membership)  
  sigmasq_all=parLapply(cl, 1:p, postSigmasq_each_INLA, Y, X, lambda, sigmasq_y, population, membership, a0 = a0, b0 = b0)
  
  return(as.numeric(unlist(sigmasq_all)))
  
}

################################################################################
##Function to propose a new lambda

##propose a single lambda for a cluster

###################################################################################
####Obtain the posterior log-likelihood value
################################################################################
################################################################################
######################For spanning-tree operations

####Generate clusters for the simulation data

GenerateClust=function(coord,clustsize,seed=seed){
  set.seed(seed);
  ns=nrow(coord);
  graph0=ConstructGraph0(coord,method='knn', 10)
  E(graph0)$weights=runif(length(E(graph0)));
  mstgraph=mst(graph0);
  
  ###For balanced clusters
  bsize=floor(ns/clustsize/2)
  membership=1:ns;
  while(min(table(membership))< bsize){
    delid=sample(1:(vcount(mstgraph)-1),clustsize-1)
    mstgraph2=delete.edges(mstgraph,delid) 
    membership=components(mstgraph2)$membership 
  }
  
  return(list(graph0=graph0,truemst=mstgraph,membership=components(mstgraph2)$membership)) 
}

###Unbalanced clusters

GenerateClust_unbalanced=function(coord,clustsize,seed=seed){
  set.seed(seed);
  ns=nrow(coord);
  graph0=ConstructGraph0(coord,method='knn', 10)
  E(graph0)$weights=runif(length(E(graph0)));
  mstgraph=mst(graph0);
  
  ###For balanced clusters
  ###bsize=floor(ns/clustsize/2)
  bsize=6
  membership=1:ns;
  while(min(table(membership))> bsize| min(table(membership))<2 ){
    delid=sample(1:(vcount(mstgraph)-1),clustsize-1)
    mstgraph2=delete.edges(mstgraph,delid) 
    membership=components(mstgraph2)$membership 
  }
  
  return(list(graph0=graph0,truemst=mstgraph,membership=components(mstgraph2)$membership)) 
}

# function to get whether an edge is within a cluster or bewteen two clusters
getEdgeStatus <- function(membership, graph) {
  inc_mat = get.edgelist(graph, names = F)
  membership_head = membership[inc_mat[, 1]]
  membership_tail = membership[inc_mat[, 2]]
  edge_status = rep('w', ecount(graph))
  edge_status[membership_head != membership_tail] = 'b'
  return(edge_status)
}

# function to split an existing cluster given MST

splitCluster <- function(mstgraph,k,membership) { 
  
  
  tcluster=table(membership) 
  clust.split=sample.int(k,1,prob=tcluster-1,replace=TRUE)
  edge_cutted=sample.int(tcluster[clust.split]-1,1)
  
  mst_subgraph=igraph::induced_subgraph(mstgraph, membership==clust.split)
  mst_subgraph=delete.edges(mst_subgraph, edge_cutted)
  connect_comp=components(mst_subgraph) 
  cluster_new = connect_comp$membership
  vid_new = (V(mst_subgraph)$vid)[cluster_new == 2]  # vid for vertices belonging to new cluster
  cluster=membership
  cluster[vid_new] = k + 1 
  return(list(cluster = cluster, vid_new = vid_new,
              clust_old = clust.split))
}



# function to merge two existing clusters

mergeCluster <- function(mstgraph, edge_status, membership) {
  # candidate edges for merging
  ecand = E(mstgraph)[edge_status == 'b']
  edge_merge = ecand[sample.int(length(ecand), 1)]
  # update cluster information
  # endpoints of edge_merge, note v1$vid > v2$vid
  v1 = head_of(mstgraph, edge_merge); v2 = tail_of(mstgraph, edge_merge)
  # clusters that v1, v2 belonging to
  cluster = membership
  
  c1 = cluster[v1]; c2 = cluster[v2]
  
  ###merge the cluster with a larger label into the one with a smaller label
  if(c1<c2){
    c_rm=c2
    c_mer=c1
  }
  else{
    
    c_rm=c1
    c_mer=c2
  }
  
  idx_rm = (cluster == c_rm)
  
  # vid of vertices in c_rm
  vid_old = (V(mstgraph))[idx_rm]
  
  # now drop c_rm
  cluster[idx_rm] = c_mer
  cluster[cluster > c_rm] = cluster[cluster > c_rm] - 1
  
  cluster_newid=c_mer
  
  # return the membership, the indices of the merged vertices
  return(list(cluster = cluster, vid_old = vid_old, cluster_rm = c_rm, cluster_comb = c_mer, cluster_newid = cluster_newid))
  
}


# function to propose a new MST
proposeMST <- function(graph0, edge_status_G) {
  nedge = length(edge_status_G)
  nb = sum(edge_status_G == 'b')
  nw = nedge - nb
  weight = numeric(nedge)
  weight[edge_status_G == 'w'] = runif(nw)
  weight[edge_status_G == 'b'] = runif(nb, 100, 200)
  mstgraph = mst(graph0, weights = weight)
  return(mstgraph)
}

################################################################################
#####Plot the MST

plotGraphData=function(coords,graph,Data,title, index=NULL){
  
  p=length(unique(Data))
  
  ##Delete the between-cluster edges
  edge_status=getEdgeStatus(Data, graph)
  edge_cutted=which(edge_status=='b')
  mst_subgraph=delete.edges(graph, edge_cutted)
  
  edgelist<-get.edgelist(mst_subgraph) 
  edgedata <- data.frame(coords[edgelist[,1],],coords[edgelist[,2],])
  colnames(edgedata) <- c("X1","Y1","X2","Y2")
  
  
  if(missing(index)){
    ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edgedata, size = 0.5, colour="grey") + 
      geom_point(data=data.frame(coords), aes(lon, lat, colour=as.factor(Data)))+scale_colour_hue()+
      ggtitle(title)+
      theme(legend.title =element_blank(),plot.title = element_text(hjust = 0.5))+ labs(fill = "Cluster ID") + xlab("X") + ylab("Y")
  }
  else{
    mtsub=data.frame(coords[index,]);colnames(mtsub)=c('X1','Y1');
    ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edgedata, size = 0.5, colour="grey") + 
      geom_point(data=data.frame(coords), aes(lon, lat, colour = as.factor(Data)))+scale_colour_hue()+ggtitle(title)+
      theme(legend.title =element_blank(),plot.title = element_text(hjust = 0.5))+
      geom_text(data = mtsub, aes(x = X1,y=Y1, label =
                                    rownames(mtsub)), hjust = 0, size = 4) + labs(fill = "Cluster ID") + xlab("X") + ylab("Y")
    
  }
  
}

#####Add circles to the mis-classified locations based on "plotGraphData" function
plotGraphData_comp=function(coords,graph, Data, Data_true, title, index=NULL){
  
  p=length(unique(Data))
  
  ##Delete the between-cluster edges
  edge_status=getEdgeStatus(Data, graph)
  edge_cutted=which(edge_status=='b')
  mst_subgraph=delete.edges(graph, edge_cutted)
  
  edgelist<-get.edgelist(mst_subgraph) 
  edgedata <- data.frame(coords[edgelist[,1],],coords[edgelist[,2],])
  colnames(edgedata) <- c("X1","Y1","X2","Y2")
  
  ###Find the location indices associated with a mis-classification
  IDX_miss=which(Data!=Data_true)
  
  
  if(missing(index)){
    ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edgedata, size = 0.5, colour="grey") + 
      geom_point(data=data.frame(coords), aes(lon, lat, color = as.factor(Data)))+scale_colour_hue() +ggtitle(title)+
      theme(legend.title =element_blank(),plot.title = element_text(hjust = 0.5)) + annotate("path",
                                                                                             x=coords[IDX_miss,1]+0.03*cos(seq(0,2*pi,length.out=100)),
                                                                                             y=coords[IDX_miss,2]+0.03*sin(seq(0,2*pi,length.out=100)))+ xlab("X") + ylab("Y")
    
    
  }
  else{
    mtsub=data.frame(coords[index,]);colnames(mtsub)=c('X1','Y1');
    ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edgedata, size = 0.5, colour="grey") + 
      geom_point(data=data.frame(coords), aes(lon, lat, color = Data))+scale_color_gradientn(colors 
                                                                                             = color_plot) +ggtitle(title)+
      theme(legend.title =element_blank(),plot.title = element_text(hjust = 0.5))+
      geom_text(data = mtsub, aes(x = X1,y=Y1, label =
                                    rownames(mtsub)), hjust = 0, size = 4) + annotate("path",
                                                                                      x=coords[IDX_miss,1]+0.03*cos(seq(0,2*pi,length.out=100)),
                                                                                      y=coords[IDX_miss,2]+0.03*sin(seq(0,2*pi,length.out=100)))+ xlab("X") + ylab("Y")
    
  }
  
}

##########Construct the hyper-parameter graph G

ConstructGraph0=function(coords,method='knn',para){
  if(method=='knn')
  {
    k=para
    graph0=nng(coords,k=k);
    edgelist=ends(graph0,E(graph0))
    weight= sqrt(apply((coords[edgelist[,1],]-coords[edgelist[,2],])^2,1,sum));
    E(graph0)$weight=weight;}
  
  
  if(method=='knn_geo')
  {
    k=para
    graph0=nng(coords,k=k);
    edgelist=ends(graph0,E(graph0))
    weight= rdist.earth(coords[edgelist[,1],], coords[edgelist[,2],], miles=FALSE);
    E(graph0)$weight=weight;}
  
  if(method=='rnn')
  {     
    eps=para
    adj=rdist(coords) 
    adj[adj>eps]=0; 
    graph0=graph_from_adjacency_matrix(adj,mode='upper',weighted=TRUE)
    
  }
  
  if(method=='rnn_geo')
  {     
    eps=para
    adj=rdist.earth(coords, miles=FALSE) 
    adj[adj>eps]=0; 
    graph0=graph_from_adjacency_matrix(adj,mode='upper',weighted=TRUE)
    
  }
  
  if(method=='deltri'){
    threshold=para;
    graph0=dentrigraph(coords,threshold=threshold)
    
  }
  return(graph0)
}


ConstructGraph0_revised=function(coords,method='knn',para){
  if(method=='knn')
  {
    k=para
    graph0=nng(coords,k=k);
    edgelist=ends(graph0,E(graph0))
    weight= sqrt(apply((coords[edgelist[,1],]-coords[edgelist[,2],])^2,1,sum));
    E(graph0)$weight=weight;}
  
  
  if(method=='knn_geo')
  {
    k=para
    graph0=nng(coords,k=k);
    edgelist=ends(graph0,E(graph0))
    weight= rdist.earth(coords[edgelist[,1],], coords[edgelist[,2],], miles=FALSE);
    E(graph0)$weight=weight;}
  
  if(method=='rnn')
  {
    eps=para
    adj=rdist(coords)
    adj[adj>eps]=0;
    graph0=graph_from_adjacency_matrix(adj,mode='upper',weighted=TRUE)
    
  }
  
  if(method=='rnn_geo')
  {
    eps=para
    adj=rdist.earth(coords, miles=FALSE)
    adj[adj>eps]=0;
    graph0=graph_from_adjacency_matrix(adj,mode='upper',weighted=TRUE)
    
  }
  
  if(method=='deltri'){
    threshold=para;
    graph0=dentrigraph(coords,threshold=threshold)
    
  }
  return(graph0)
}

ord.mat = function(M, decr = F, cols = NULL){
  if(is.null(cols))
    cols = 1: ncol(M)
  out = do.call( "order", as.data.frame(M[,cols]))
  if (decr)
    out = rev(out)
  return(M[out,])
}


dentrigraph=function(coords,threshold=1000){
  triangulation <- deldir(coords[,1], coords[,2])
  distances <- abs(triangulation$delsgs$x1 - triangulation$delsgs$x2) +
    abs(triangulation$delsgs$y1 - triangulation$delsgs$y2)
  #remove edges that are greater than .001
  edge.list <- as.matrix(triangulation$delsgs[distances < threshold,5:6])
  
  sorted <- sort(c(edge.list), index.return = TRUE)
  run.length <- rle(sorted$x)
  indices <- rep(1:length(run.length$lengths),times=run.length$lengths)
  edge.list.reformatted <- edge.list
  edge.list.reformatted[sorted$ix] <- indices
  edge.list.reformatted=ord.mat(edge.list.reformatted);
  #create graph from list of edges
  distances <- sqrt((coords[edge.list.reformatted[,1],1] - coords[edge.list.reformatted[,2],1])^2 +
                      (coords[edge.list.reformatted[,1],2] - coords[edge.list.reformatted[,2],2])^2)
  graph0 <- graph_from_edgelist(edge.list.reformatted, directed = TRUE)
  E(graph0)$weight=distances
  return(graph0)
}

#####################################################################################
#####################################################################################
#####################################################################################

##The main function for the Fclust-RST method
##Y: the n by m data matrix; 
##X: the m by m inverse discrete wavelet transform matrix;
##graph0: the hyper-parameter graph, G;
##init_val: lists of initial values;
##hyperpar: lists of hyperparameters;
##MCMC: the number of MCMC iterations;
##BURNIN: the burn-in period;
##THIN: the thinning setting;
##path_save: the file path and file name for storing the results;
##seed: the random seed specification;
##n.cores: the number of CPU cores to use.

Fun_Wave_Clust <- function(Y, X, graph0, init_val, hyperpar, MCMC, BURNIN, THIN, path_save, seed = 1234, n.cores) {
  
  cl <- snow:::makeCluster(n.cores)
  clusterExport(cl, "Para_re2all")
  
  set.seed(seed)
  
  nt = nrow(X); ns = nrow(Y); 
  
  # initial values
  
  mstgraph = init_val[['trees']]
  beta = init_val[['beta']]
  cluster = init_val[['cluster']]  
  sigmasq_y = init_val[['sigmasq_y']]
  lambda = init_val[['lambda']]
  # ind_gam=init_val[['ind_gam']]
  
  k = max(cluster) # number of clusters
  
  #For recording the acceptance ratio
  hyper_cnt = 0; birth_cnt=0; death_cnt=0; change_cnt=0;
  
  ################# RJMCMC ####################
  
  # hyper-parameter
  rhy = 0.05
  c = hyperpar$c
  a0 = hyperpar$a0; b0 = hyperpar$b0
  c0 = hyperpar$c0; d0 = hyperpar$d0
  hyper = list(a0, b0, c0 ,d0, c)
  
  
  ### initialize log likelihood vector
  log_like_ini = evalLogLike_all_parallel2(Y, X, lambda, sigmasq_y, population, cluster, cl, a0 = a0, b0 = b0)
  log_like = log_like_ini$llik_all;log_like_vec=log_like_ini$log_like_vec

  
  # whether an edge in graph0 is within a cluster or bewteen two clusters
  # n*p matrix
  edge_status = getEdgeStatus(cluster, mstgraph)
  
  ## MCMC results
  beta_out = list(); sigmasq_y_out=list();
  lambda_out = list(); 
  mstgraph_out=list(); cluster_out = array(0, dim = c((MCMC-BURNIN)/THIN, ns))
  MST_out = list()
  loglike_out = numeric((MCMC-BURNIN)/THIN)
  
  #time_start = Sys.time()
  ## MCMC iteration
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
      # vec_lambda_k=Para_re2all(lambda_new[k+1,], J) for ind_gamma
      sigmasq_yk= postSigmasq_each_INLA(k +1, Y, X, lambda_new, sigmasq_y,population, a0, b0, membership_new)
      sigmasq_y_new = c(sigmasq_y, sigmasq_yk)
      
      
      
      
      # compute log-prior ratio
     log_A = log(1-c) 
      # + (a0/2)*log(b0/2) - log(gamma(a0/2)) -(a0/2+1)*log(sigmasq_yk) - b0/(2*sigmasq_yk) +
      #   sum(c0/2*log(d0/2)-log(gamma(c0/2))-(c0/2+1)*log(lambda_new[k+1,])-d0/(2*lambda_new[k+1,]))
      ## calculate loglikelihood ratio by only comparing local likelihood of the two clusters that changed.
      
      log_L_new=evalLogLike.ratio('split',log_like_vec, split_res, Y, X, lambda_new, sigmasq_y_new, membership_new, a0, b0)
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
        # ind_gam=ind_gam_new
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
      # vec_lambda_k=Para_re2all(lambda[cid_comb,], J)
      
      sigmasq_yk=postSigmasq_each_INLA(k, Y, X, lambda, sigmasq_y,population, a0, b0, membership)
      sigmasq_y_new[cid_newid]=sigmasq_yk; 
      
      # compute log-prior ratio
      
      log_A = -log(1-c)
      
      # compute log-likelihood ratio
      
      log_L_new=evalLogLike.ratio('merge',log_like_vec, merge_res, Y, tPHI, lambda_new, sigmasq_y_new, membership_new)
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
        # ind_gam=ind_gam_new;
        
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
      # ind_gam_new=ind_gam[-cid_rm,]
      
      ##For hetero variances
      
      cid_comb = merge_res$cluster_comb
      cid_newid = merge_res$cluster_newid
      
      ind_k=which(membership_new == cid_newid); ns1=length(ind_k)
      
      Y_clust=Y[ind_k,]
      # vec_lambda_k=Para_re2all(lambda[cid_comb,], J)
      # sigmasq_yk=postSigmasq(Y_clust, X, vec_lambda_k, ind_gam[cid_comb,], a0, b0)
      sigmasq_yk = 0.01
      sigmasq_y_new[cid_newid]=sigmasq_yk; 
      
      # Obtain the log prior for ind_gam
      # log_prior_gam=0; curr=2
      # 
      # for(j in 1:J){
      #   
      #   ind_sel=curr:(curr+2^{J-j}-1)
      #   
      #   temp_gam=sum(ind_gam[cid_rm, ind_sel])
      #   log_prior_gam=log_prior_gam + log(gamma(temp_gam+1)) + log(gamma(2^{J-j}+1-temp_gam)) - log(gamma(2^{J-j}+2))
      #   
      #   curr=curr+2^{J-j}
      #   
      # }
      
      # log_prior_merge=- (a0/2)*log(b0/2) + log(gamma(a0/2)) + (a0/2+1)*log(sigmasq_y[cid_rm]) +  b0/(2*sigmasq_y[cid_rm]) - 
      #   sum(c0/2*log(d0/2)-log(gamma(c0/2))-(c0/2+1)*log(lambda[cid_rm,])-d0/(2*lambda[cid_rm,])) - log_prior_gam
      log_prior_merge = 0
      
      k = k-1
      
      log_L_new_merge=evalLogLike.ratio('merge',log_like_vec, merge_res, Y, tPHI, lambda_new, sigmasq_y_new, membership_new)
      
      # then perform birth move
      split_res = splitCluster(mstgraph, k, merge_res$cluster);
      split_res$k=k+1;membership_new=split_res$cluster
      
      k = k+1
      
      
      if(k>2){
        lambda_new = rbind(lambda_new, lambda_new[split_res$clust_old,])
        # ind_gam_new = rbind(ind_gam_new, ind_gam_new[split_res$clust_old,])
        
      }else{
        
        lambda_new = rbind(lambda_new, lambda_new)
        # ind_gam_new = rbind(ind_gam_new, ind_gam_new)
        
        
      }
      # For hetero variances
      ind_k=which(membership_new==k); ns1=length(ind_k)
      Y_clust=Y[ind_k,]
      # vec_lambda_k=Para_re2all(lambda_new[k,], J)
      sigmasq_yk= 0.01
      sigmasq_y_new = c(sigmasq_y_new, sigmasq_yk)
      
      
      # Obtain the log prior for ind_gam
      # log_prior_gam=0; curr=2
      # 
      # for(j in 1:J){
      #   
      #   ind_sel=curr:(curr+2^{J-j}-1)
      #   
      #   temp_gam=sum(ind_gam_new[k, ind_sel])
      #   log_prior_gam=log_prior_gam + log(gamma(temp_gam+1)) + log(gamma(2^{J-j}+1-temp_gam)) - log(gamma(2^{J-j}+2))
      #   
      #   curr=curr+2^{J-j}
      # }
      
      log_prior_split=0
      
      
      # compute log-likelihood ratio
      
      log_L_new=evalLogLike.ratio('split',log_L_new_merge$log_like_vec, split_res, Y, tPHI, lambda_new, sigmasq_y_new, membership_new)
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
        # ind_gam=ind_gam_new
        
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
      lambda.res = matrix(1, J+1, 1)
      # lambda.res=postLambda_par(beta, sigmasq_y, c0, d0, cluster, cl)
      
      lambda=t(matrix(lambda.res, (J+1), k))
      
      #######################################################################################      
  
      ##Update log-likelihood value
      log_like_res = evalLogLike_all_parallel2(Y, tPHI, lambda, sigmasq_y, population, cluster, cl)
      log_like=log_like_res$llik_all;log_like_vec=log_like_res$log_like_vec;
    }
    
    ############################################################################## 
    # update sigma2
    sigmasq_y_new=postSigmasq_par(k, Y, X, lambda, sigmasq_y,population, a0, b0, cluster, cl)
    # sigmasq_y_new = rep(0.01, k)
    sigmasq_y = sigmasq_y_new
    
    ###store estimates
    if(iter %% 10 == 0) {
      cat('Iteration', iter, 'k', k, 'birth cnt', birth_cnt, 'death cnt', death_cnt, 'change cnt', change_cnt, 'hyper cnt', hyper_cnt, '\n')
      cat('move', move, 'k', k, 'ncluster', length(unique(cluster)), '\n') 
    }
    
    ## save result
    if(iter > BURNIN & (iter - BURNIN) %% THIN == 0) {
      
      beta_out[[(iter-BURNIN)/THIN]] = beta
      # lambda_out[[(iter-BURNIN)/THIN]] = lambda
      # sigmasq_y_out[[(iter-BURNIN)/THIN]] = sigmasq_y
      
      MST_out[[(iter-BURNIN)/THIN]] = mstgraph
      cluster_out[(iter-BURNIN)/THIN, ] = cluster
      loglike_out[(iter-BURNIN)/THIN] = log_like
    }
    
    if(iter %% 100 == 0) {
      cat('Iteration', iter, 'birth cnt', birth_cnt, 'death cnt', death_cnt, 'change cnt', change_cnt, 'hyper cnt', hyper_cnt, '\n')
      
      save(beta_out, lambda_out, sigmasq_y_out, MST_out, cluster_out, loglike_out, hyperpar, birth_cnt, death_cnt, change_cnt, file=path_save)
      
    }
    
  }
  
  ####Consider a parallel implementation      
  # # sigmasq_y_new=postSigmasq_par(Y, X, lambda, ind_gam, a0, b0, cluster, cl)
  # sigmasq_y_new = rep(0.01, k)
  # sigmasq_y = sigmasq_y_new
  
  # update beta
  beta_res=postBeta_par(k, Y, X, lambda, sigmasq_y,population,membership, a0, b0)
  beta=t(matrix(beta_res, J+1, k))
  
  ##Update log-likelihood value
  log_like_res = evalLogLike_all_parallel2(Y, tPHI, lambda, sigmasq_y, population, cluster, cl)
  log_like=log_like_res$llik_all;log_like_vec=log_like_res$log_like_vec;
  
  stopCluster(cl)
  
}





