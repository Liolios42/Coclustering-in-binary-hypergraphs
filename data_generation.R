data_generation <- function(N,M,K,G,gamma,delta,p_kg){
  
  
  node_group <- sample(1:K, N, prob = gamma, replace = T)
  hyper_group <- sample(1:G, M,  prob = delta, replace = T)
  
  prob_gen <- matrix(p_kg, nrow = K, ncol = G)
  
  x_ij <- matrix(nrow = M, ncol = N)
  
  for (i in seq(M)){
    for (j in seq(N)){
      x_ij[i, j] <- rbinom(prob = prob_gen[node_group[j], hyper_group[i]], n = 1, size = 1)
    }
  }
  
  
  return(list('incidence_matrix' = x_ij, 'node_clusters' = node_group, 'hyper_clusters' = hyper_group))
}