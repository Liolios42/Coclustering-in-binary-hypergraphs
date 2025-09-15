setwd('/data1/s3645177/tidy')

library(HyperG)
library(boot)
library(car)
library(gtools)
library(ggplot2)
library(foreach)
library(Rmpfr)
library(mclust)
library(bikm1)
library(dplyr)
library(Metrics)




source('d_posterior.R')
source('c_posterior.R')
source('eta.R')
source('data_generation.R')
source('em_parallel.R')
source('optimal_perm.R')



start.time <- Sys.time()
mc.replicates <- 500
n.cores <- 32
em.runs <- 3
N <- c(200)# Number of nodes
M <- c(200) # Number of hyperedges
K <- c(2) # Number of node-groups
G <- c(2)  # Number of hyperedge-groups



gamma <- c(0.8, 0.2)
delta <- c(0.7, 0.3)
p_kg <- c(0.9, 0.6, 0.4, 0.2)








results <- data.frame(N = integer(),
                      M = integer(),
                      run = integer(),
                      theta = numeric(),
                      gamma = numeric(),
                      delta = numeric(),
                      iterations = numeric(),
                      stringsAsFactors = FALSE)


em_list <- list()
condition_gamma <- c()
condition_delta <- c()


for(n in N) {
  for(m in M){
    
    for (k in K) {
      for (g in G) {
        
        
        
        for(replic in seq(mc.replicates)){
          
          
          # gamma <- rdirichlet(1, c(rep(1, k))) # Node prior assignment probability
          # delta <- rdirichlet(1, c(rep(1, g))) # Hyperedge prior assignment probability
          # p_kg <- runif(k*g) # Prior connection probability
          
          synthetic_data <- data_generation(n, m, k, g, gamma, delta, p_kg)
          
          x_ij <- synthetic_data$incidence_matrix
          random_hypergraph <- hypergraph_from_incidence_matrix(x_ij)
          
          
          node.clusters.true <- synthetic_data$node_clusters
          hyper.clusters.true <- synthetic_data$hyper_clusters
          
          
          for (run in seq(em.runs)) {
            
            
            #x_ij <- data_generation(n, m, 3, 3, gamma, delta, p_kg)
            #random_hypergraph<- hypergraph_from_incidence_matrix(x_ij)
            
            constants2 <- c(random_hypergraph, n, m, k, g)
            
            
            
            em_results <- data.frame(
              theta = numeric(),
              gamma = numeric(),
              delta = numeric(),
              stringsAsFactors = FALSE)
            
            
            
            e <- EM_parallel(constants = constants2, T_max = 20, cores = n.cores)
            mapping <- optimal_perm(node.clusters.true, hyper.clusters.true, e$c_t, e$d_t, constants = constants2)
	    e$gamma <- e$gamma[mapping$node_labels]
	    e$delta <- e$delta[mapping$hyper_labels]
	    e$p_kg <- e$p_kg[mapping$node_labels, mapping$hyper_labels]
	    conf_hyper <- mapping$confusion_mat_hyper
	    conf_node <- mapping$confusion_mat_node
            em_results <- rbind(em_results, data.frame(theta_1 = e$p_kg[1],
                                                       theta_2 = e$p_kg[2],
                                                       theta_3 = e$p_kg[3],
                                                       theta_4 = e$p_kg[4],
                                                       gamma_1 = e$gamma[1],
                                                       gamma_2 = e$gamma[2],
                                                       delta_1 = e$delta[1],
                                                       delta_2 = e$delta[2],
                                                       iterations = e$iterations))
            
            
            

	    
	    node.clusters.est <- as.numeric(apply(as.matrix(e$c_t), 1, function(row) which.max(row)))
	    node.clusters.est <- mapping$node_labels[node.clusters.est]
	    hyper.clusters.est <- as.numeric(apply(as.matrix(e$d_t), 1, function(row) which.max(row)))
	    hyper.clusters.est <- mapping$hyper_labels[hyper.clusters.est]
            
            if(!any(is.na(node.clusters.est)) & !any(is.na(hyper.clusters.est))){
              ICL <- BinBlocICL_LBM(a = 4, b = 1, x_ij, hyper.clusters.est, node.clusters.est)
              rand_nodes <- adjustedRandIndex(node.clusters.est,node.clusters.true)
              rand_edges <- adjustedRandIndex(hyper.clusters.est,hyper.clusters.true)
              acc_nodes <- mean(node.clusters.est == node.clusters.true)
              acc_hyperedges <- mean(hyper.clusters.est == hyper.clusters.true)
            }
            else {
              ICL <- 1
              rand_nodes <- -1
              rand_edges < -1
            }
            list2 <- list('ICL' = ICL,
                          'rand_nodes' = rand_nodes,
                          'rand_edges' = rand_edges,
                          'N' = n,
                          'M' = m,
                          'replicate'= replic,
                          'em_run' = run,
                          'iterations' = e$iterations,
                          'K' = k,
                          'G' = g,
                          'theta' = e$p_kg,
                          'gamma' = e$gamma,
                          'delta' = e$delta,
                          confusion_matrix_node = conf_node,
                          confusion_matrix_hyper = conf_hyper)
            
            em_list[[em.runs*(replic - 1) + run]] <- list2
            results <- rbind(results, data.frame(N = n, M = m, replicate = replic, em_run = run, ICL = ICL,
                                                 theta_1 = em_results$theta_1,
                                                 theta_2 = em_results$theta_2,
                                                 theta_3 = em_results$theta_3,
                                                 theta_4 = em_results$theta_4,
                                                 gamma_1 = em_results$gamma_1,
                                                 gamma_2 = em_results$gamma_2,
                                                 delta_1 = em_results$delta_1,
                                                 delta_2 = em_results$delta_2,
                                                 iterations = em_results$iterations,
                                                 acc.nodes <- acc_nodes,
                                                 acc.hyperedges <- acc_hyperedges)) 
            print(results)
            
            
            
          }
        }
        
      }
      
      
      print(paste('replic =',replic))
      
      
    }
    print('run finished')
  }
}

end.time <- Sys.time()

df <- results



time_taken <- end.time - start.time
print(time_taken)
# write.csv(results, file = 'test_output/em_estimates.csv')
# write.csv(df, file = 'test_output/em_sorted_estimates.csv')
# save(results, file = 'test_output/em_hypergraph.RData')



write.csv(results, file = 'auc_outputs/em_estimates_200_200.csv')
saveRDS(em_list, file = 'auc_outputs/em_hypergraph_200.RDS')
