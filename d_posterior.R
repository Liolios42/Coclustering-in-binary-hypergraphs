


d_parallel <- function(p_kg, pi, tau, c, constants = constants2, subset){
  
  start.time <- Sys.time()
  
  hypergraph_matrix <- constants[1]
  N <- as.numeric(constants[2])
  M <- as.numeric(constants[3])
  K <- as.numeric(constants[4])
  G <- as.numeric(constants[5])
  
  
  X <- as.matrix(hypergraph_matrix$M)
  
  
  tau <- c(tau)
  
  p_kg <- matrix(p_kg, nrow = K, ncol = G)
  pi <- mpfrArray(pi, dim = c(G,1), 1024)
  
  # u_jk = matrix(0, nrow = M, ncol = K)
  
  u_jk <- X%*%mpfr(c, 1024)
  n_k <- colSums(c)
  
  d_est <- mpfrArray(rep(0, M*G), dim = c(M,G), 1024)
  
  numerator <- mpfrArray(rep(0, K), dim = c(K,1), 1024)
  
  
  for (j in subset){
    for (g in seq(G)){
      
      
      for (k in seq(K)){
        
        #to_replace <- sum(t(c[,k])*X[j,])
        #if (typeof(to_replace) == 'list') {to_replace <- asNumeric(to_replace)}
        # u_jk[j,k] = to_replace
        # n_k <- sum(c[,k])
        p_choose <-  p_kg[k,g]
        numerator[k] <- p_choose ** (u_jk[j,k]) *(1- p_choose) ** (n_k[k] - u_jk[j,k])
        
        
        
        
      }
      d_est[j,g] <- mpfr(pi[g] * prod(numerator), 1100)
      
    }
  }
  result <- d_est
  
  # denom <- rowSums(d_est)
  # result <- d_est/denom
  # 
  # 
  # d_est <- colSums(result) / M
  # print(d_est)
  # print(d_est/denom)
  # print(paste('Pi 1 > Pi 2', sum(result[,1] >= result[,2])))
  
  return(list('result' = result, 'delta_estimate' = asNumeric(d_est)))
  
}
