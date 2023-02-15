############################################################
#                                                          #
#                    Graph Simulations                     #
#                                                          #
############################################################

source("manuscript/global.R")

############################################################
# Create A Simulation Grid

BernDF = expand.grid(name = "Bernoulli", p = 3:15, prob = seq(0, 1, 0.1))
TreeDF = expand.grid(name = c("star","line","binary","rev_binary"), p = 3:15, prob = NA)
CyclDF = expand.grid(name = "dag_cycle", p = 3:15, prob = NA)
sim_grid = 
  rbind(BernDF, TreeDF, CyclDF) %>% 
  filter(!(name %in% c("binary","rev_binary") & p %in% setdiff(3:15, c(3,5,15))))



for(i in 1:nrow(sim_grid)){
  
  p = sim_grid$p[i]
  prob = sim_grid$prob[i]
  type = sim_grid$name[i]
  
  ############################################################
  # 0. Generate Plot Title
  if(type == "Bernoulli"){
    title = sprintf("%s n:%s p:%s", type, p, prob)  
  } else{
    title = sprintf("%s n:%s", type, p)  
  }
  
  
  ############################################################
  # 1. Generate Lowe Triangular Adjacency Matrix
  
  if(type == "Bernoulli"){
    w = bernouli_dag(p, prob)  
  } else if(type == "dag_cycle"){
    w = dag_cycle(p)
  } else{
    w = tree_dag(p, type = type)
  }
  
  ############################################################
  # 2. Plot Dag Related To Adjacency Matrix
  dag <- adjacency_to_dag(w)
  p_graph <- ggdag(dag, layout = "circle") + 
    theme_dag()

  
  ############################################################
  # 3. Generate SCM corresponding to Noise and Adjacency Matrix
  scm <- lin_scm(w = w, noise = normal_noise, n = 100000)
  
  
  ############################################################
  # 4. Estimate Tree by Chatterjee Algorithms
  xi <- cor_mat(normal_mat(scm))
  
  
  ############################################################
  # 5. plot the correlation matrix VS Adjacency
  p_true <- cor_adj_plot(xi, w)
  
  ############################################################
  # 6. Estimate Tree
  e <- estimate_tree(xi)
  
  e_dag <- adjacency_to_dag(e)
  e_graph <- ggdag(e_dag, layout = "circle") + theme_dag()
  
  
  ############################################################
  # 7. find the difference between true and estimated graph
  p_fault <- cor_fault_plot(xi, w + t(w), e)
  
  
  ############################################################
  # 8. find the minimum spanning Tree
  g  <- graph.adjacency(round(1-xi,2), weighted = TRUE, 
                        mode = c("undirected"))
  g %>% mst() %>% as_adjacency_matrix() %>% as.matrix() -> mst_w
  mst_dag <- adjacency_to_dag(mst_w)
  mst_graph <- ggdag(mst_dag, layout = "circle") + theme_dag()
  p_mst_fault <- cor_fault_plot(xi, mst_w, e)
  
  
  ############################################################
  # 
  cp <- ggarrange(p_graph, e_graph, mst_graph, p_true, p_fault,p_mst_fault,
                  labels = c(title, "Estimated Tree", 
                             "MSP",
                             "Xi Correlation Matrix", 
                             "True VS Estimated","Estimated vs MST"),
                  ncol = 3, nrow = 2,
                  font.label = list(size = 12)
                  )+ bgcolor("white")   
  
  ############################################################
  #
  w = 40 + 3   * (p - 10)*(p > 10)
  h = 30 + 2.25 * (p - 10)*(p > 10)
  ggsave(sprintf("plots/pdf/%s.pdf",gsub("\\s+","_",title)), cp, 
         device = cairo_pdf, width = w, height = h, unit ="cm")
  ggsave(sprintf("plots/%s.png",gsub("\\s+","_",title)), cp, 
         device = "png", width = w, height = h, unit ="cm")
  print(i)
}




