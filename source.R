
######### R libraries #########

library(igraph)
library(HelpersMG)
library(graphkernels)
library(Matrix)

################################################################################
############### Simulating network from an IRG model ###########################
################################################################################

#input
# p : edge probability matrix
sample_IRG = function(p){
  n = nrow(p)
  
  # Creating the n x n adjacency matrix  
  adj <- matrix(0, n, n)
  
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      adj[i, j] = rbinom(1, 1, p[i,j]) # We include the edge with probability p
    }
  }
  
  adjsymm = symmetricize(adj, method = "ld")
  
  # graph from the adjacency matrix
  G = igraph::graph_from_adjacency_matrix(adjsymm, mode = "undirected", weighted = NULL)
  
  return(G)
}

######### Planting a clique in a network simulated from an ER model ############

network_with_clique_ER = function(G, K){
  G_pc = G
  smpl = sample(1:gorder(G_pc),K)
  new_edges = combinations(smpl, k = 2, layout = "column")
  ap = as_ids(E(G_pc)[new_edges[1,] %--% new_edges[2,]])
  to_del = c(1:gsize(G_pc))
  to_del = to_del[!to_del %in% ap]
  del_eid = sample(to_del,ncol(new_edges)-length(ap))
  G_pc = delete_edges(G_pc, del_eid)
  G_pc <- G_pc + edges(as.vector(new_edges))
  G_pc = simplify(G_pc, remove.multiple = TRUE, remove.loops = TRUE)
  return(G_pc)
}

######## Planting a clique in a network simulated from an ERMM model ###########

network_with_clique = function(G, C, K, max.rep = 20){
  G_pc = G
  E(G_pc)$label = edge.labels(G, C)
  counter = 0
  repeat{
    counter = counter + 1
    smpl = sample(1:gorder(G_pc),K)
    new_edges = combinations(smpl, k = 2, layout = "column")
    
    L= length(unique(C))
    edge.types = g_labels(L)
    new.edge.labels = c()
    for (i in 1:ncol(new_edges)) {
      new.edge.labels[i] = edge.types[C[new_edges[1,i]],C[new_edges[2,i]]]
    }
    new_type = as.data.frame(table(new.edge.labels, dnn = list("type")), responseName = "type_n")
  
    ap = as_ids(E(G_pc)[new_edges[1,] %--% new_edges[2,]])
    
    if(length(ap) == 0){ 
      to_del_type = new_type
      to_del_type$to_del = to_del_type$type_n
    }
    else{
      ap_type = as.data.frame(table(E(G_pc)[ap]$label, dnn = list("type")), responseName = "type_n")
      to_del_type = merge(new_type, ap_type,by=c("type"), all=TRUE)
      to_del_type[is.na(to_del_type)] <- 0
      to_del_type$to_del = to_del_type$type_n.x - to_del_type$type_n.y
    }
    
    to_del = c(1:gsize(G_pc))
    to_del = to_del[!to_del %in% ap]
    
    ###### Checking if we have enough edges to delete ##########
    cond = c()
    for (j in 1:length(to_del_type$type)) {
      if(length(to_del[E(G_pc)[to_del]$label == to_del_type$type[j]]) < to_del_type$to_del[j]) cond[j] = 1
      else cond[j] = 0
    }
    if(sum(cond) == 0) break
    if(counter > max.rep){ #print("Error: Network too sparse. Not enought edges to delete")
      break}
  }
  
  if (counter > max.rep) {
    print("Error: Network too sparse. Not enought edges to delete.")
    list("network" = make_empty_graph(n = 0, directed = FALSE),"rep" = counter)
  }
  else{
    ####### Sampling the edges to delete from the reduced edge list #######
    
    del_eid = c()
    for (j in 1:length(to_del_type$type)) {
      set = to_del[E(G_pc)[to_del]$label == to_del_type$type[j]]
      if(to_del_type$to_del[j] == 0) del_eid = del_eid
      else{
        if(length(set) == 1 && to_del_type$to_del[j] == 1){del_eid = c(del_eid, set)}
        else del_eid = c(del_eid,sample(set,to_del_type$to_del[j],,replace = FALSE))
      }
    }

    G_pc = delete_edges(G_pc, c(ap,del_eid))
    G_pc <- G_pc + edges(as.vector(new_edges))
    
    list("network" = G_pc, "rep" = counter)
  }
}

########## Planting hubs in a network simulated from an ERMM model #############

planted_hub <- function(G, d_m, k, C) {
  d <- degree(G)
  s_d <- sd(d)
  
  u_max <- which(d == d_m)[1]
  if (is.na(u_max)) return(G)
  
  d_star <- max(1, d_m + ceiling(k * s_d))
  n_new_nei <- d_star - degree(G, u_max)
  
  target_v <- setdiff(V(G), neighbors(G, u_max))
  target_v <- setdiff(target_v, u_max)
  
  if (length(target_v) == 0 || n_new_nei == 0) return(G)
  
  # Filter target_v by type
  filtered_targets <- c()
  for (v in target_v) {
    nei_v <- neighbors(G, v)
    if (any(C[nei_v] == C[u_max])) {
      filtered_targets <- c(filtered_targets, v)
    }
  }
  
  if (length(filtered_targets) == 0) return(G)
  
  if (length(filtered_targets) <= n_new_nei) {
    pot_nei <- filtered_targets
  } else {
    pot_nei <- sample(filtered_targets, n_new_nei)
  }
  
  # Loop to build to_del while avoiding duplicates
  
  to_del <- c()
  new_nei <- c()
  for (v in pot_nei) {
    to_del_temp = to_del
    nei_v <- neighbors(G, v)
    same_type <- nei_v[C[nei_v] == C[u_max]]
    
    if (length(same_type) == 0) next
    
    # Try deleting a unique edge without creating duplicates
    found <- FALSE
    to_sample <- same_type
    while (length(to_sample) > 0 && !found) {
      del_node <- if (length(to_sample) == 1) to_sample else sample(to_sample, 1)
      del_edge <- sort(c(v, del_node))
      to_del_temp <- c(to_del, del_edge)
      
      # Split into edge list and check for duplicates
      edge_list <- split(to_del_temp, ceiling(seq_along(to_del_temp) / 2))
      if (!any(duplicated(edge_list))) {
        to_del <- to_del_temp
        new_nei <- c(new_nei, v)
        found <- TRUE
      } else {
        to_sample <- setdiff(to_sample, del_node)
      }
    }
  }
  
  if (length(to_del) %% 2 != 0) to_del <- to_del[-length(to_del)]
  
  if (length(to_del) == 0) return(G)
  if (length(to_del) != length(new_nei) * 2) {
    message("Mismatch: to_del has ", length(to_del), 
            ", expected ", length(new_nei) * 2)
    return(G)
  }
  
  if (length(to_del) > 0) {
    eid <- get.edge.ids(G, to_del)
    G_hub <- delete_edges(G, eid[eid > 0])
  }
  
  new_edges <- as.vector(rbind(rep(u_max, length(new_nei)), new_nei))
  G_hub <- G_hub + edges(new_edges)
  
  return(G_hub)
}

planted_hubs <- function(G, Rep = 5, k = 3, C) {
  G_h <- G
  d_m = max(degree(G))
  i=1
  while (!is.na(d_m) && i <= Rep) {
    G_h <- planted_hub(G_h, d_m, k, C)
    i = i+1
    d_unique <- sort(unique(degree(G_h)), decreasing = TRUE)
    d_m <- d_unique[i]
  } 
  
  list("Graph" = G_h, "Iter" = i-1)
}

################################################################################
#################### IRG-gKSS and the GOF test setup ###########################
################################################################################

# Different kernel choices
# CalculateWLKernel(P,level = 3)
# CalculateGeometricRandomWalkKernel(P,level = 3)
# compute.sp.kernel(P,level = 1)
# CalculateVertexEdgeHistGaussKernel(P,level = 0.1)
# CalculateGraphletKernel(P,level = 3)
# CalculateConnectedGraphletKernel(P,level = 3)

# test_G: An igraph object containing the observed network
# C: A vector of group labels for all the vertices
# p: The edge probability matrix from the null model
# g.kernel: choice of graph kernel for IRG-gKSS 
# alpha: level of significance for the test 
# M: Number of iterations in Monte Carlo setting to compute the null set

GOF_IRG = function(test_G, C , p_0, M = 200, g.kernel= CalculateWLKernel, level = 3, alpha = 0.05){
  
  ######## Test Statistics from observed network ###########
  
  test_stat = generate.one.gKSS(test_G, C, p_0, g.kernel,level)$stats.value
  
  ############## Simulations from null model ERMM-gKSS ##############

  statistic = c()
  for (i in 1:M){
    g = sample_IRG(p_0)
    statistic[i] = generate.one.gKSS(g, C, p_0, g.kernel,level)$stats.value
  }
  
  critical_value_lower = quantile(statistic,alpha/2)
  critical_value_upper = quantile(statistic,(1-alpha/2))
  
  ######## Decision and P-value #############
  
  pval <- two_tailed_pval(test_stat, statistic)
  
  ################### Output ERMM-gKSS#########################
  list("obs.test.stat" = test_stat, "lower_quantile" = critical_value_lower, "upper_quantile" = critical_value_upper, "P_value" = pval, "null_set" = statistic)
}

generate.one.gKSS=function(G, C, p, g.kernel=CalculateWLKernel,level=3){
  # We use group labels as vertex labels for WL kernel computation
  V(G)$name = C
  
  X = as_adjacency_matrix(G, type = "both")
  S = X[upper.tri(X)]
  n = length(S)
  S.t =  p[upper.tri(p)]
  S.t.vec = abs(S - S.t)
  S.mat = S.t.vec %*% t(S.t.vec)
  
  P = compute.transition.list(X)
  K = compute.kernel(P, g.kernel, level)
  
  J.kernel = S.mat * K
  stats.value = mean(J.kernel) 
  
  #Return:
  #stats.value: IRG-gKSS^2
  #J.kernel: Matrix whose enteries correspond to h(s,s')
  list("stats.value" = stats.value, "h(s,s')" = J.kernel)
}

######################## Test with edge resampling #############################

GOF_IRG_resamp = function(test_G, C , p_0, s_size, M = 200, g.kernel= CalculateWLKernel, level = 3, alpha = 0.05){
  
  ######## Test Statistics from observed network ###########
  
  sample.index = sample_index(n, s.size = s_size, replace = TRUE)
  test_stat = generate.one.gKSS.sampled(test_G, C, p_0, sample.index, g.kernel,level)$stats.value
  
  ############## Simulations from null model ERMM-gKSS ##############
  
  statistic = c()
  for (i in 1:M){
    g = sample_IRG(p_0)
    sample.index = sample_index(n, s.size = s_size, replace = TRUE)
    statistic[i] = generate.one.gKSS.sampled(g, C, p_0, sample.index, g.kernel,level)$stats.value
  }
  
  critical_value_lower = quantile(statistic,alpha/2)
  critical_value_upper = quantile(statistic,(1-alpha/2))
  
  ######## Decision and P-value #############
  
  pval <- two_tailed_pval(test_stat, statistic)
  
  ################### Output ERMM-gKSS#########################
  list("obs.test.stat" = test_stat, "lower_quantile" = critical_value_lower, "upper_quantile" = critical_value_upper, "P_value" = pval, "null_set" = statistic)
}

generate.one.gKSS.sampled=function(G, C, p, sample.index, g.kernel=CalculateWLKernel,level=3){
  V(G)$name = C
  
  X = as_adjacency_matrix(G, type = "both")
  S = matrix(X,byrow=TRUE)[sample.index]
  S.t =  p[sample.index]
  S.t.vec = abs(S - S.t)
  S.mat = S.t.vec %*% t(S.t.vec)
  
  P = compute.sampled.list(X, sample.index)
  K = compute.kernel(P, g.kernel, level)
  
  J.kernel = S.mat * K
  stats.value = mean(J.kernel)
  
  #Return:
  #stats.value: IRG-gKSS^2
  #J.kernel: Matrix whose enteries correspond to h(s_b,(s_b)')
  list("stats.value" = stats.value, "IRG-gKSS^2" = J.kernel)
}

################################################################################
###################### Other Helper functions ##################################
################################################################################

####### Computes Kernel and inner product matrix #######

compute.kernel=function(P, g.kernel = CalculateWLKernel, level=3){
  n = length(P) - 1 
  kernel.matrix = g.kernel(P, level)
  K = kernel.matrix[1:n,1:n] + kernel.matrix[n+1,n+1]
  K.vec = kernel.matrix[1:n,n+1]
  K = K - outer(K.vec, K.vec, '+')
  return(K)
}

compute.transition.list = function(X){
  n <- nrow(X)
  P <- vector("list", n * (n - 1) / 2 + 1)  # Preallocate list
  counter <- 1
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      Y <- X
      Y[i, j] <- abs(1 - X[i, j])
      Y[j, i] <- Y[i, j]
      P[[counter]] <- graph_from_adjacency_matrix(Y, mode = "undirected")
      counter <- counter + 1
    }
  }
  # Original graph as the last element
  P[[counter]] <- graph_from_adjacency_matrix(X, mode = "undirected")
  
  return(P)
}

compute.sampled.list = function(X, sample.index){
  P=list()
  l = length(sample.index)
  for (w in 1:l){
    x = X
    x[sample.index[w]] = abs(1 - X[sample.index[w]])
    x = symmetricize(x ,method = c("ud"), adjacencyList = FALSE)
    G = graph_from_adjacency_matrix(x, mode = "undirected")
    P[[w]] = G
  }
  P[[l+1]] = graph_from_adjacency_matrix(X, mode = "undirected")
  P
}

two_tailed_pval <- function(test_stat, statistic) {
  tstat_set = c(test_stat, statistic)
  p_lower <- mean(tstat_set <= test_stat)
  p_upper <- mean(tstat_set >= test_stat)
  pval <- 2 * min(p_lower, p_upper)
  return(min(pval, 1))  # Cap at 1
}

##### Maximum likelihood estimates of ER and ERMM parameters #######
MLE.est = function(G, C)
{
  L = length(unique(C)) # no of block
  bs = as.numeric(table(C))
  X = as_adjacency_matrix(G, type = "both")
  N_ij = matrix(0,nrow = L,ncol = L)
  for (i in 1:L) {
    for (j in 1:L) {
      if(i==j) N_ij[i,i] = (bs[i]*(bs[i]-1)/2)
      else {N_ij[i,j] = bs[i]*bs[j]
      }
    }
  }
  N_ij
  n_e = matrix(0,nrow = L,ncol = L)
  for (i in 1:L) {
    for (j in i:L) {
      n_e[i,j] = sum(X[C == i, C == j])
    }
  }
  diag(n_e) = diag(n_e)/2
  n_e = symmetricize(n_e, method = "ud") 
  
  p_hat = matrix(0,nrow = L, ncol = L)
  for (i in 1:L) {
    for (j in i:L) {
      if(i==j) p_hat[i,j] = n_e[i,i]/(bs[i]*(bs[i] - 1)/2)
      else p_hat[i,j] = n_e[i,j]/(bs[i]*bs[j]) 
    }
  }
  p_hat = symmetricize(p_hat, method = "ud")
  
  list("estimates" = p_hat, "E" = n_e, "N_ij" = N_ij )
}

##### Sample upper triangular indices #####
sample_index = function(n, s.size, replace = TRUE){
  Non = matrix(seq(1,n*n,1), nrow = n, ncol = n, byrow = FALSE)
  N_set = Non[upper.tri(Non)]
  sample.index = sample(N_set, size = s.size, replace)
  sample.index
}
