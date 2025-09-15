eta <- function(p_kg, pi, tau, delta, gamma, constants = constants2){
  
  
  hypergraph_matrix <- constants[1]
  N <- as.numeric(constants[2])
  M <- as.numeric(constants[3])
  K <- as.numeric(constants[4])
  G <- as.numeric(constants[5])
  
  
  X <- as.matrix(hypergraph_matrix$M)
  
  
  p_kg <- matrix(p_kg, nrow = K, ncol = G)
  
  gamma <- matrix(asNumeric(gamma), nrow = N, ncol = K)
  
  delta <- matrix(asNumeric(delta), nrow = M, ncol = G)
  
  
  sum <- array(0, dim = c(K,G))
  total_sum <- array(0, dim = c(K,G))
  for (k in seq(K)){
    
    for (g in seq(G)){
      
      for (i in seq(N)){
        
        for (j in seq(M)){
          
          
          sum[k,g] <-  sum[k,g] + delta[j, g] * gamma[i, k] * X[j,i]
          total_sum[k,g] <- total_sum[k,g] + delta[j, g] * gamma[i, k]
          
        }
      }
    }
    
  }
  sum <- mpfrArray(c(as.numeric(sum)), 512)
  total_sum <- mpfrArray(c(as.numeric(total_sum)), 512)
  result <- sum/total_sum
  
  result <- matrix(asNumeric(result), nrow = K, ncol = G)
  
  
  return(result)
}
