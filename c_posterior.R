

c_parallel<- function(p_kg, pi, tau,  d, constants = constants2, cores = 1, subset){
  
  
  hypergraph_matrix <- constants[1]
  N <- as.numeric(constants[2])
  M <- as.numeric(constants[3])
  K <- as.numeric(constants[4])
  G <- as.numeric(constants[5])
  
  
  X <- as.matrix(hypergraph_matrix$M)
  
  
  
  p_kg <- matrix(p_kg, nrow = K, ncol = G)
  
  u_ig <- t(t(mpfr(d, 1024))%*%X)
  n_g <- colSums(d)
  
  gamma_est <- mpfrArray(rep(0, N*K), dim = c(N,K), 1024)
  numerator <- matrix(rep(0,G), nrow = G, ncol = 1)
  # u_ig = mpfrArray(rep(0, N*G), dim = c(N,G), 156)
  # u_ig = mpfrArray(rep(0, N*G), dim = c(N,G), 512)
  
  
  tau <- mpfrArray(tau, dim = c(K, 1), 1300)
  
  
  for (i in subset){
    
    for (k in seq(K)){
      
      
      for (g in seq(G)){
        
        to_replace <- sum(t(d[,g])*X[,i])
        if (typeof(to_replace) == 'list') {to_replace <- asNumeric(to_replace)}
        # u_ig[i,g] <- to_replace
        # n_k <- sum(d[,g])
        p_choose <-  p_kg[k,g]
        numerator[g] <- p_choose ** (u_ig[i,g]) *(1- p_choose) ** (n_g[g] - u_ig[i,g])
        
        
        
      }
      
      numerator <- mpfr(numerator, 1024)
      gamma_est[i,k] <- mpfr(tau[k] * prod(numerator), 1024)
      
    }
    
  }
  result <- gamma_est
  # 
  # denom <- rowSums(gamma_est)
  # 
  # result <- gamma_est/denom
  # 
  # 
  # 
  # print(paste('Tau 1 > tau 2', sum(result[,1] >= result[,2])))
  # tau_est <- colSums(result) / N
  # print(tau_est)
  # print(delta_est/denom)
  # print(sum(delta_est[,1] >= delta_est[,2]))
  
  
  return(list('result' = result, 'gamma_estimate' = asNumeric(gamma_est)))
  
}




