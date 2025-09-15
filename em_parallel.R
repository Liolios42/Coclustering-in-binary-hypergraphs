EM_parallel <- function(constants = constants2, T_max = 30, epsilon = 0.001, cores, initialization = 'random'){
  
  # create cluster with 6 cores:
  clust = parallel::makeCluster(cores, outfile = "")
  doParallel::registerDoParallel(clust)
  
  hypergraph_matrix <- constants[1]
  M <- as.numeric(constants[3])
  N <- as.numeric(constants[2])
  
  K <- as.numeric(constants[4])
  G <- as.numeric(constants[5])
  
  
  indices.m <- seq(M)
  indices.n <- seq(N)
  
  index_sets.m <- split(indices.m, cut(indices.m, cores, labels = FALSE))
  index_sets.n <- split(indices.n, cut(indices.n, cores, labels = F))
  
  
  
  # Initialization of parameters 
  if (initialization == 'random'){
  gamma <- rdirichlet(1, c(rep(1, K)))
  delta <- rdirichlet(1, c(rep(1, G)))
  prob <- as.matrix(runif(K*G, 0, 1), nrow = K, ncol = G)
  }
  
  initial_d <- rdirichlet(M, rep(1,G))
  initial_c <- rdirichlet(N, rep(1,K))
  
  c_t <- initial_c
  s = 1
  s_max = 2
  fixed_converged = F
  while(s < s_max & !fixed_converged){
  	  if (s > 1){prev_dt <- d_t}
	  d_t = foreach(i = seq(length(index_sets.m)), .combine = 'rbind', .packages = c('gtools', 'Rmpfr', 'HyperG'), .errorhandling="pass") %dopar% {
	    source('d_posterior.R')
	    new.M = index_sets.m[[i]]
	    
	    d_results_temp <- d_parallel(prob, delta, gamma, c_t,constants, subset = new.M)$result[new.M,]
	    d_results_temp
	    
	    # c_results = rbind(c_results, c_results_temp)
	    # c_results
	  }
	  d_t
	  denom <- rowSums(d_t)
	  d_t <- d_t/denom
	  
	 
	  
	  prev_ct <- c_t
	  
	  c_t = foreach(i = seq(length(index_sets.n)), .combine = 'rbind', .packages = c('gtools', 'Rmpfr', 'HyperG'),.errorhandling="pass") %dopar% {
	    source('c_posterior.R')
	    new.N = index_sets.n[[i]]
	    
	    c_results_temp <- c_parallel(prob, delta, gamma, asNumeric(d_t),constants, subset = new.N)$result[new.N,]
	    c_results_temp
	    # c_results = rbind(c_results, c_results_temp)
	    # c_results
	  }
	  
	  denom <- rowSums(c_t)
	  c_t <- c_t/denom
	   
	  if (s > 1){
	  fixed_converged <- all(abs(prev_dt - d_t) < epsilon) & 
      		all(abs(prev_ct - c_t) < epsilon) 
      	  if (is.na(fixed_converged)){fixed_converged <- T}}
      	  s = s + 1 
	  }
  delta_est <- colSums(d_t) / M
  gamma_est <- colSums(c_t)/N
  et <- eta(prob, delta, gamma, d_t, c_t)
  
  delta_update <- delta_est
  gamma_update <- gamma_est
  
  
  t = 1
  converged = F
  while(t < T_max & !converged){

    
    params <- c(prob)
    
    temp_delta_update <- delta_est
    temp_gamma_update <- gamma_est
    temp_prob <- et
    
    s = 1
    fixed_converged <- F
    while(s < s_max & !fixed_converged){
  	  prev_dt <- d_t
	  d_t = foreach(i = seq(length(index_sets.m)), .combine = 'rbind', .packages = c('gtools', 'Rmpfr', 'HyperG'), .errorhandling="pass") %dopar% {
	    source('d_posterior.R')
	    new.M = index_sets.m[[i]]
	    
	    d_results_temp <- d_parallel(prob, delta, gamma, c_t,constants, subset = new.M)$result[new.M,]
	    d_results_temp
	    
	    # c_results = rbind(c_results, c_results_temp)
	    # c_results
	  }
	  d_t
	  denom <- rowSums(d_t)
	  
	  d_t <- d_t/denom
	  
	 
	  
	  
	  prev_ct <- c_t
	  c_t = foreach(i = seq(length(index_sets.n)), .combine = 'rbind', .packages = c('gtools', 'Rmpfr', 'HyperG'),.errorhandling="pass") %dopar% {
	    source('c_posterior.R')
	    new.N = index_sets.n[[i]]
	    
	    c_results_temp <- c_parallel(prob, delta, gamma, asNumeric(d_t),constants, subset = new.N)$result[new.N,]
	    c_results_temp
	    # c_results = rbind(c_results, c_results_temp)
	    # c_results
	  }
	  
	  denom <- rowSums(c_t)
	  c_t <- c_t/denom
	  s = s + 1  
	  fixed_converged <- all(abs(prev_dt - d_t) < epsilon) & 
      		all(abs(prev_ct - c_t) < epsilon) 
      	  if (is.na(fixed_converged)){fixed_converged <- T}
	  }
	  
    
    delta_est <- colSums(d_t) / M
    gamma_est <- colSums(c_t)/N
    
    et <- eta(params, delta_update, gamma_update, d_t, c_t)
    
    
    
    delta_update <- delta_est
    gamma_update <- gamma_est
    prob <- et
    
    
    converged <- all(abs(temp_delta_update - delta_update) < epsilon) & 
      all(abs(temp_gamma_update - gamma_update) < epsilon) & 
      all(abs(temp_prob - prob) < epsilon)
    if(is.na(converged)) {converged = T}
    
    t <- t + 1
    print(paste('t =', t))
    
    
    # print(paste('p_kg =', prob))
    # print(paste('delta =', delta_update))
    # print(paste('gamma =',gamma_update))
  }
  parallel::stopCluster(clust)
  return(list('p_kg' = prob, 'delta' = asNumeric(delta_update), 'gamma' = asNumeric(gamma_update), 'iterations' = t, 'd_t' = d_t, 'c_t'= c_t))
}
