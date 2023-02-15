############################################################
#                                                          #
#                        SCM Models                        #
#                                                          #
############################################################

lin_scm = function(w, noise, n){
  p = ncol(w)
  exogenous = noise(n, p)
  endogenous <- matrix(0, nrow = n, ncol = p)
  
  for (i in 1:p) {
    if(i ==1){
      endogenous[,i] = endogenous[,1:i] * w[1, 1] + exogenous[, i]
    }else{
      endogenous[,i] = endogenous[,1:i] %*% w[1:i,i] + exogenous[, i]  
    }
  }
  colnames(endogenous) = colnames(w)
  return(endogenous)
}

############################################################
#                                                          #
#                       Noise Model                        #
#                                                          #
############################################################

# All Noise have mu = 0 and sigma = 1
# Totally Normal Noise
normal_noise = function(n,p) matrix(rnorm(n*p),ncol = p)

# One Categorical and Another is Normal
one_cat_noise = function(n,p,prob = 0.5) matrix(
  c(2*rbinom(n,1,prob)-1,rnorm(n*(p-1)))
  ,ncol = p, byrow = F)

############################################################
#                                                          #
#                Generate Causal Structure                 #
#                                                          #
############################################################

label_it <- function(w){
  p <- ncol(w)
  nodes <- paste0("X",stringr::str_pad(1:p,2,"left","0"))
  rownames(w) <- nodes
  colnames(w) <- nodes
  return(w)
}
# Generate binary lower triangular matrix as adjacent matrix
bernouli_dag <- function(p, prob = 0.5){
  binaries <- rbinom(p*(p-1)/2, 1, prob)
  w = matrix(0,p,p)
  k = 1
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      w[i,j] = binaries[k]
      k = k+1
    }
  }
  label_it(w)
}

# tree
tree_dag <- function(p, type = "line"){
  w = matrix(0, p, p)
  
  if(type == "line"){
    for (i in 1:(p-1)) {
      w[i,i+1] = 1
    }    
  }
  
  if(type == "star"){
    for (i in 1:(p-1)) {
      w[1,i+1] = 1
    }    
  }
  
  if(type == "binary"){
    k = log(p + 1, base = 2)
    k = ceiling(k)
    p = 2^k - 1
    w = matrix(0, p, p)
    for (i in 1:(k-1)) {
      for (j in 1:(2^i)) {
        w[2^(i-1) - 1 + ceiling(j/2), 2^i - 1 + j] = 1
      }
    }    
  }

  if(type == "rev_binary"){
    k = log(p + 1, base = 2)
    k = ceiling(k)
    p = 2^k - 1
    w = matrix(0, p, p)
    for (i in rev(0:(k-2))) {
      for (j in 1:(2^i)) {
        w[2*(2^i - 1 + j),     2^i - 1 + j] = 1
        w[2*(2^i - 1 + j) + 1, 2^i - 1 + j] = 1
      }
    }
  }
  
  label_it(w)
}


dag_cycle <- function(p){
  w = matrix(0, p, p)
  for (i in 1:(p-1)) {
    w[i,i+1] = 1
  }
  w[1,p] = 1
  label_it(w)
}

############################################################
#                                                          #
#                 Convert Adjacency To Dag                 #
#                                                          #
############################################################

adjacency_to_dag <- function(w){
  edge = list()
  p <- nrow(w)
  rnames = rownames(w)
  cnames = colnames(w)
  for(i in 1:(p)){
    ind <- which(w[i,]!=0)
    if(length(ind)>0){
      edge[[i]] = paste0( rnames[i]," -> ", rnames[ind])
    } else{
      edge[[i]] = paste0( rnames[i]," -> ", rnames[i])
      }
  }
  dagitty::dagitty(sprintf("dag {%s}",paste(unlist(edge), collapse = ";")))
}



cor_mat = function(x) {
  p = ncol(x)
  n = nrow(x)
  d = matrix(0, nrow = p, ncol = p)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      d[i,j] = xicorln(x[,i], x[,j])
      d[j,i] = xicorln(x[,j], x[,i])
    }
  }
  rownames(d) <- colnames(x)
  colnames(d) <- colnames(x)
  return(d)
}




normal_mat <- function(m){
  for(i in 1:ncol(m)){
    m[,i] = (m[,i] - mean(m[,i],na.rm = T))/sd(m[,i],na.rm = T)
  }
  return(m)
}

cor_adj_plot <- function(xi, w){
  p <- ggcorrplot(xi , show.legend = F, p.mat = (w + t(w)), 
                  lab = TRUE, show.diag = FALSE,  pch = 1,
                  colors = c("white","#ffff00", "#FF0000"),
                  outline.color = "transparent", pch.col = "black",
                  pch.cex = 15) +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title=element_blank()
          )
  return(p)
}

estimate_tree <- function(d){
  m = d
  p = ncol(m)
  maxval = max(d)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      u = d[,i]
      v = d[,j]
      w1 = setdiff(which(u >= d[j,i]), c(i,j))
      w2 = setdiff(which(v >= d[i,j]), c(i,j))
      if (length(intersect(w1, w2)) > 0) {
        m[i,j] = 0
        m[j,i] = 0
      }
      if (m[i,j] != 0 & m[j,i] != 0) {
        u1 = maxval + 1 - m[i,j]
        u2 = maxval + 1 - m[j,i]
        u = max(u1,u2)
        m[i,j] = u
        m[j,i] = u
      } 
    }
  }
  m[m!=0] = 1
return(m) 
}

cor_fault_plot <- function(xi, w, e){

  p <- ggcorrplot(xi , show.legend = F, p.mat = abs(w - e), 
                  lab = TRUE, show.diag = FALSE,  pch = 4,
                  colors = c("white","#ffff00","#FF0000"),
                  outline.color = "transparent", pch.col = "#005AFF",pch.cex = 10)+
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title=element_blank()
    )
  return(p)
}
